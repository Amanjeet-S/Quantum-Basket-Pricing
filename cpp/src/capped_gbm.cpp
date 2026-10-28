#include "qbp/capped_gbm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>

namespace qbp {
namespace {

constexpr Real pi = 3.141592653589793238462643383279502884L;
constexpr Real sqrt_two = 1.414213562373095048801688724209698079L;
constexpr Real log_sqrt_two_pi = 0.91893853320467274178032973640561764L;

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

Real finite_result(Real x, const char* message) {
    if (!std::isfinite(x)) throw std::overflow_error(message);
    return x;
}

void finite_greeks(const Greeks& values) {
    finite_result(values.price, "Price exceeds the supported floating-point range.");
    for (const auto& value : {values.delta, values.gamma, values.vega, values.rho, values.theta})
        if (value) finite_result(*value, "A Greek exceeds the supported floating-point range.");
}

Real discount(Real rate, Real maturity) {
    const Real value = finite_result(std::exp(-rate * maturity), "Discount factor is outside the supported floating-point range.");
    if (value == 0) throw std::range_error("A positive discount factor underflowed; scaled evaluation is required outside this implementation's supported range.");
    return value;
}

Real log_ratio(Real a, Real b) {
    // The split form also remains valid where a/b would overflow or underflow.
    // log1p preserves a small relative difference when the two inputs are close.
    const Real relative = (a - b) / b;
    if (std::isfinite(relative) && std::abs(relative) < 0.5L) return std::log1p(relative);
    return std::log(a) - std::log(b);
}

Real kink_width(const Contract& c) {
    const Real ratio = c.cap / c.strike;
    if (std::isfinite(ratio)) return std::log1p(ratio);
    return std::log(c.cap) - std::log(c.strike) + std::log1p(c.strike / c.cap);
}

// Integrate phi from upper-width to upper. Tail-aware erfc subtraction avoids
// subtracting two CDFs close to one. A narrow interval uses positive quadrature.
Real normal_interval(Real upper, Real width, Real lower) {
    require(width >= 0 && std::isfinite(width), "Normal interval width is invalid.");
    if (width == 0) return 0;
    const Real scale = std::max({1.0L, std::abs(lower), std::abs(upper)});
    if (width * scale < 0.05L) {
        constexpr std::array<Real, 4> nodes = {
            0.183434642495649804939476142360184L, 0.525532409916328985817739049189246L,
            0.796666477413626739591553936475831L, 0.960289856497536231683560868569473L};
        constexpr std::array<Real, 4> weights = {
            0.362683783378361982965150449277195L, 0.313706645877887287337962201986601L,
            0.222381034453374470544355994426240L, 0.101228536290376259152531354309962L};
        const Real middle = upper - width / 2;
        Real sum = 0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            sum += weights[i] * (normal_pdf(middle - width * nodes[i] / 2)
                               + normal_pdf(middle + width * nodes[i] / 2));
        }
        return width * sum / 2;
    }
    if (lower >= 0) return (std::erfc(lower / sqrt_two) - std::erfc(upper / sqrt_two)) / 2;
    if (upper <= 0) return (std::erfc(-upper / sqrt_two) - std::erfc(-lower / sqrt_two)) / 2;
    return (std::erf(upper / sqrt_two) - std::erf(lower / sqrt_two)) / 2;
}

// phi(upper)-phi(upper-width), preserving a small difference with expm1.
Real density_difference(Real upper, Real width, Real lower) {
    const Real exponent_difference = width * (upper + lower) / 2;
    if (exponent_difference <= 0) return normal_pdf(upper) * -std::expm1(exponent_difference);
    return normal_pdf(lower) * std::expm1(-exponent_difference);
}

Real log_cdf_upper_bound(Real coordinate) {
    // Mills' upper bound Phi(-x)<=phi(x)/x, x>0. This is used solely
    // as a supported-range diagnostic, never to return a tail price.
    if (coordinate >= 0) return 0;
    return std::min(0.0L, -coordinate * coordinate / 2
                         - std::log(-coordinate) - log_sqrt_two_pi);
}

Real log_interval_upper_bound(Real upper, Real lower, Real width) {
    const Real nearest = lower <= 0 && upper >= 0 ? 0
                       : std::min(std::abs(lower), std::abs(upper));
    return std::min(0.0L, std::log(width) - nearest * nearest / 2 - log_sqrt_two_pi);
}

void reject_recovered_underflow(Real observed, Real log_probability_bound,
                                std::initializer_list<Real> log_prefactors) {
    if (observed >= std::numeric_limits<Real>::min()) return;
    const Real smallest = std::numeric_limits<Real>::denorm_min() > 0
                        ? std::numeric_limits<Real>::denorm_min() : std::numeric_limits<Real>::min();
    for (Real log_prefactor : log_prefactors) {
        if (log_prefactor > 0 && log_prefactor + log_probability_bound >= std::log(smallest) - 4)
            throw std::range_error("An underflowed or subnormal normal tail/density may regain representable scale in a price/Greek; logarithmically scaled evaluation is required outside this implementation's supported range.");
    }
}

struct Quadrature {
    Real integral = 0;
    Real error = 0;
    std::size_t evaluations = 0;
};

// Paired 8-/16-point Gauss-Legendre rules provide a numerical error indicator,
// not an outward-rounded or rigorous quadrature error certificate.
std::pair<Real, Real> gauss_pair(const std::function<Real(Real)>& f, Real left, Real right) {
    constexpr std::array<Real, 8> nodes = {
        0.095012509837637440185319335424958L, 0.281603550779258913230460501460496L,
        0.458016777657227386342419442983578L, 0.617876244402643748446671764048791L,
        0.755404408355003033895101194847442L, 0.865631202387831743880467897712393L,
        0.944575023073232576077988415534608L, 0.989400934991649932596154173450333L};
    constexpr std::array<Real, 8> weights = {
        0.189450610455068496285396723208283L, 0.182603415044923588866763667969220L,
        0.169156519395002538189312079030359L, 0.149595988816576732081501730547479L,
        0.124628971255533872052476282192016L, 0.095158511682492784809925107602246L,
        0.062253523938647892862843836994378L, 0.027152459411754094851780572456018L};
    constexpr std::array<Real, 4> nodes8 = {
        0.183434642495649804939476142360184L, 0.525532409916328985817739049189246L,
        0.796666477413626739591553936475831L, 0.960289856497536231683560868569473L};
    constexpr std::array<Real, 4> weights8 = {
        0.362683783378361982965150449277195L, 0.313706645877887287337962201986601L,
        0.222381034453374470544355994426240L, 0.101228536290376259152531354309962L};
    const Real middle = (left + right) / 2;
    const Real half = (right - left) / 2;
    Real high = 0, low = 0;
    for (std::size_t i = 0; i < nodes.size(); ++i)
        high += weights[i] * (f(middle - half * nodes[i]) + f(middle + half * nodes[i]));
    for (std::size_t i = 0; i < nodes8.size(); ++i)
        low += weights8[i] * (f(middle - half * nodes8[i]) + f(middle + half * nodes8[i]));
    return {high * half, std::abs(high - low) * half};
}

Quadrature adaptive_integral(const std::function<Real(Real)>& f, Real left,
                             Real right, Real tolerance, unsigned depth = 0) {
    const auto [value, error] = gauss_pair(f, left, right);
    const Real roundoff = 64 * std::numeric_limits<Real>::epsilon() * std::abs(value);
    if (error <= tolerance + roundoff) return {value, error, 24};
    if (depth >= 24) throw std::runtime_error("Positive digital-price quadrature did not meet its numerical tolerance.");
    const Real middle = (left + right) / 2;
    if (!(left < middle && middle < right))
        throw std::runtime_error("Positive digital-price quadrature exhausted floating-point interval resolution.");
    const auto first = adaptive_integral(f, left, middle, tolerance / 2, depth + 1);
    const auto second = adaptive_integral(f, middle, right, tolerance / 2, depth + 1);
    return {first.integral + second.integral, first.error + second.error,
            24 + first.evaluations + second.evaluations};
}

Quadrature digital_price(const Contract& c, Real d2, Real v, Real bound) {
    // V/(H D)=integral_0^1 Phi(d2-log(1+u H/K)/v) du. Every integrand is
    // non-negative. Breakpoints follow normal-score levels so a narrow digital
    // transition cannot be entirely missed by the adaptive quadrature nodes.
    const Real ratio = c.cap / c.strike;
    require(std::isfinite(ratio) && ratio > 0,
            "The positive digital fallback cannot resolve the normalised strike scale; H/K must be finite and positive.");
    const Real transition_scale = v / ratio;
    require(std::isfinite(transition_scale)
            && transition_scale >= std::numeric_limits<Real>::min(),
            "The positive digital fallback cannot resolve its normalised transition in the supported floating-point range.");
    const auto log_strike_increment = [&](Real u) {
        return std::log1p(u * ratio);
    };
    const auto f = [&](Real u) { return normal_cdf(d2 - log_strike_increment(u) / v); };
    std::vector<Real> points{0, 1};
    for (int score = -16; score <= 16; score += 2) {
        const Real log_relative = v * (d2 - static_cast<Real>(score));
        if (log_relative <= 0 || log_relative >= kink_width(c)) continue;
        const Real u = std::expm1(log_relative) / ratio;
        if (std::isfinite(u) && u > 0 && u < 1) points.push_back(u);
    }
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    Quadrature result;
    // Tolerance is relative to the positive integral's coarse estimate, with a
    // tiny absolute floor. It is a numerical target, not a rigorous bound.
    Real estimate = 0;
    for (std::size_t i = 1; i < points.size(); ++i) estimate += gauss_pair(f, points[i-1], points[i]).first;
    const Real tolerance = std::max(1e-30L, 2e-15L * std::abs(estimate));
    for (std::size_t i = 1; i < points.size(); ++i) {
        const auto part = adaptive_integral(f, points[i-1], points[i],
                                             tolerance * (points[i] - points[i-1]));
        result.integral += part.integral;
        result.error += part.error;
        result.evaluations += part.evaluations + 24;
    }
    result.integral *= bound;
    result.error *= bound;
    return result;
}

Greeks central_greeks(const Contract& c, const Bumps& bumps,
                      const std::function<Real(const Contract&)>& price,
                      std::optional<Real> base_price = std::nullopt) {
    Greeks out;
    out.price = base_price ? *base_price : price(c);
    for (unsigned i = 0; i < 4; ++i) {
        Contract lower = c, upper = c;
        const Real step = i == 0 ? bumps.spot : i == 1 ? bumps.volatility : i == 2 ? bumps.rate : bumps.maturity;
        Real Contract::* field = i == 0 ? &Contract::spot : i == 1 ? &Contract::volatility
                                    : i == 2 ? &Contract::rate : &Contract::maturity;
        lower.*field -= step;
        upper.*field += step;
        const Real low_price = price(lower), high_price = price(upper);
        const Real difference = (high_price - low_price) / (2 * step);
        if (i == 0) {
            out.delta = difference;
            out.gamma = (high_price - 2 * out.price + low_price) / (step * step);
        } else if (i == 1) out.vega = difference;
        else if (i == 2) out.rho = difference;
        else out.theta = -difference;
    }
    finite_greeks(out);
    return out;
}

Real interpolated_price(const Contract& c, const HeatGrid& grid, const HeatSolution& solution) {
    const Real y = log_ratio(c.spot, c.strike)
                 + (c.rate - c.dividend - c.volatility * c.volatility / 2) * c.maturity;
    // validate() rejects exterior queries. Endpoint equality is handled exactly.
    if (y == grid.upper) return discount(c.rate, c.maturity) * solution.values.back();
    if (y == grid.lower) return discount(c.rate, c.maturity) * solution.values.front();
    const Real coordinate = (y - grid.lower) / grid.spacing();
    const auto left = static_cast<std::size_t>(std::floor(coordinate));
    if (left >= grid.intervals) throw std::runtime_error("Interior interpolation coordinate is not representable on the requested grid.");
    const Real fraction = coordinate - static_cast<Real>(left);
    return finite_result(discount(c.rate, c.maturity)
                         * ((1 - fraction) * solution.values[left] + fraction * solution.values[left + 1]),
                         "PDE price is outside the supported floating-point range.");
}

} // namespace

void Contract::validate() const {
    for (Real x : {spot, strike, cap, maturity, rate, dividend, volatility})
        require(std::isfinite(x), "All contract parameters must be finite.");
    require(spot > 0, "Spot must be strictly positive.");
    require(strike >= 0 && cap >= 0 && maturity >= 0 && volatility >= 0,
            "Strike, cap, maturity and volatility must be non-negative.");
    require(std::isfinite(strike + cap), "Upper strike exceeds the supported floating-point range.");
}

Real normal_cdf(Real x) { return std::erfc(-x / sqrt_two) / 2; }
Real normal_pdf(Real x) { return std::exp(-x * x / 2) / std::sqrt(2 * pi); }

Real payoff(Real terminal_spot, const Contract& c) {
    return std::min(std::max(terminal_spot - c.strike, 0.0L), c.cap);
}

AnalyticResult analytic(const Contract& c) {
    c.validate();
    AnalyticResult out;
    out.derivative_convention = "two-sided derivatives; calendar theta; unit vega/rho";
    auto& g = out.values;
    if (c.cap == 0) {
        g = {0, 0, 0, 0, 0, 0};
        out.price_method = "zero_cap";
        return out;
    }
    if (c.maturity == 0 || c.volatility == 0) {
        // At sigma=0 the analytic extension is the deterministic discounted
        // payoff. At tau=0 theta is its right-maturity limit, if it exists.
        const Real forward = finite_result(c.spot * std::exp((c.rate - c.dividend) * c.maturity),
                                            "Deterministic terminal spot exceeds the supported floating-point range.");
        const Real dr = discount(c.rate, c.maturity), dq = discount(c.dividend, c.maturity);
        g.price = dr * payoff(forward, c);
        const bool lower_kink = forward == c.strike;
        const bool upper_kink = forward - c.strike == c.cap;
        if (!lower_kink && !upper_kink) {
            const bool active = forward > c.strike && forward - c.strike < c.cap;
            g.delta = active ? dq : 0;
            g.gamma = 0;
            g.vega = 0;
            g.rho = c.maturity * (active ? c.strike * dr : -g.price);
            g.theta = active ? c.dividend * c.spot * dq - c.rate * c.strike * dr : c.rate * g.price;
            out.derivative_convention = "deterministic limit away from kinks; one-sided tau=0/sigma=0 derivatives";
        } else {
            // The deterministic spot derivative is undefined at either kink.
            // At expiry the risk/volatility derivatives are identically zero;
            // calendar theta is not finite for positive sigma at a payoff kink.
            if (c.maturity == 0) {
                g.vega = 0;
                g.rho = 0;
                if (c.volatility == 0) {
                    const Real movement = (c.rate - c.dividend) * c.spot;
                    g.theta = lower_kink ? -std::max(movement, 0.0L)
                                        : c.rate * g.price - std::min(movement, 0.0L);
                }
            }
            else {
                // Price as a function of non-negative sigma has a finite
                // right derivative at sigma=0, even though spot derivatives do not.
                g.vega = (lower_kink ? 1 : -1) * c.spot * dq * std::sqrt(c.maturity) / std::sqrt(2 * pi);
                if (c.rate == c.dividend) g.theta = c.rate * g.price;
            }
            out.derivative_convention = "undefined kink derivatives are null; vega is right-sided at sigma=0; rho=vega=0 at tau=0";
        }
        out.price_method = c.maturity == 0 ? "expiry_payoff" : "deterministic_zero_volatility";
        finite_greeks(g);
        return out;
    }

    const Real t = c.maturity, root_t = std::sqrt(t), v = c.volatility * root_t;
    const Real dr = discount(c.rate, t), dq = discount(c.dividend, t);
    const Real sdq = finite_result(c.spot * dq, "Dividend-adjusted spot exceeds the supported floating-point range.");
    const Real bound = finite_result(c.cap * dr, "Discounted cap exceeds the supported floating-point range.");
    require(sdq > 0 && bound > 0,
            "A positive discounted spot or cap underflowed; scaled evaluation is required outside this implementation's supported range.");
    const Real log_stock = std::log(c.spot) - c.dividend * t;
    const Real log_cap = std::log(c.cap) - c.rate * t;
    const Real log_q = c.dividend == 0 ? -std::numeric_limits<Real>::infinity() : std::log(std::abs(c.dividend));
    const Real log_r = c.rate == 0 ? -std::numeric_limits<Real>::infinity() : std::log(std::abs(c.rate));
    const Real log_gamma_factor = -c.dividend * t - std::log(c.spot) - std::log(v);
    const Real log_vega_factor = log_stock + std::log(root_t);
    const Real log_theta_density_factor = log_stock + std::log(c.volatility) - std::log(2.0L) - std::log(root_t);
    if (c.strike == 0) {
        // min(S_T,H) = S_T - call(H). This equivalent put representation
        // avoids subtracting two prices of order S when S is large.
        const Real d1 = (log_ratio(c.spot, c.cap) + (c.rate - c.dividend + c.volatility * c.volatility / 2) * t) / v;
        const Real d2 = d1 - v;
        const Real below1 = normal_cdf(-d1), above2 = normal_cdf(d2), density = normal_pdf(d1);
        reject_recovered_underflow(below1, log_cdf_upper_bound(-d1),
                                  {log_stock, -c.dividend * t, log_q + log_stock});
        reject_recovered_underflow(above2, log_cdf_upper_bound(d2),
                                  {log_cap, std::log(t) + log_cap, log_r + log_cap});
        reject_recovered_underflow(density, -d1 * d1 / 2 - log_sqrt_two_pi,
                                  {log_gamma_factor, log_vega_factor, log_theta_density_factor});
        g.price = sdq * below1 + bound * above2;
        g.delta = dq * below1;
        g.gamma = -dq * density / (c.spot * v);
        g.vega = -sdq * density * root_t;
        g.rho = -t * bound * above2;
        g.theta = sdq * density * c.volatility / (2 * root_t)
                  + c.dividend * sdq * below1 + c.rate * bound * above2;
        out.price_method = "zero_strike_positive_terms";
        finite_greeks(g);
        return out;
    }
    const Real d1 = (log_ratio(c.spot, c.strike)
                  + (c.rate - c.dividend + c.volatility * c.volatility / 2) * t) / v;
    const Real d2 = d1 - v;
    const Real width = kink_width(c) / v;
    require(width > 0, "The positive normalised payoff-kink separation underflowed to zero; a differently scaled implementation is required.");
    const Real log_upper = c.cap > c.strike
        ? log_ratio(c.spot, c.cap) - std::log1p(c.strike / c.cap)
        : log_ratio(c.spot, c.strike) - kink_width(c);
    const Real upper_d1 = (log_upper + (c.rate - c.dividend + c.volatility * c.volatility / 2) * t) / v;
    const Real upper_d2 = upper_d1 - v;
    require(std::isfinite(d1) && std::isfinite(d2) && std::isfinite(width),
            "Normalised strike coordinates exceed the supported floating-point range.");
    // For a broad strike interval, computing the upper-strike coordinate from
    // its own stable log ratio avoids subtracting two enormous d-coordinates.
    // For a narrow interval the separately retained width carries its signal.
    const Real lower1 = width * std::max(1.0L, std::abs(d1)) < 0.05L ? d1 - width : upper_d1;
    const Real lower2 = width * std::max(1.0L, std::abs(d2)) < 0.05L ? d2 - width : upper_d2;
    const Real interval1 = normal_interval(d1, width, lower1);
    const Real interval2 = normal_interval(d2, width, lower2);
    const Real high_tail = normal_cdf(lower2);
    const Real density = density_difference(d1, width, lower1);
    const Real log_strike = std::log(c.strike) - c.rate * t;
    const Real discounted_strike = finite_result(c.strike * dr, "Discounted strike exceeds the supported floating-point range.");
    require(discounted_strike > 0,
            "A positive discounted strike underflowed; scaled evaluation is required outside this implementation's supported range.");
    reject_recovered_underflow(interval1, log_interval_upper_bound(d1, lower1, width),
                              {log_stock, -c.dividend * t, log_q + log_stock});
    reject_recovered_underflow(interval2, log_interval_upper_bound(d2, lower2, width),
                              {log_strike, std::log(t) + log_strike, log_r + log_strike});
    reject_recovered_underflow(high_tail, log_cdf_upper_bound(lower2),
                              {log_cap, std::log(t) + log_cap, log_r + log_cap});
    const Real density_scale = std::max(normal_pdf(d1), normal_pdf(lower1));
    const Real log_density_bound = std::log(2.0L) - std::min(d1 * d1, lower1 * lower1) / 2 - log_sqrt_two_pi;
    reject_recovered_underflow(density_scale, log_density_bound,
                              {log_gamma_factor, log_vega_factor, log_theta_density_factor});
    const Real stock_term = sdq * interval1;
    const Real strike_term = discounted_strike * interval2;
    const Real cap_term = bound * high_tail;
    g.price = std::fma(-discounted_strike, interval2, stock_term) + cap_term;
    const Real total = std::abs(stock_term) + std::abs(strike_term) + std::abs(cap_term);
    out.cancellation_ratio = total > 0 ? std::abs(g.price) / total : 1;
    out.price_method = "tail_aware_closed_form";
    if (g.price < 0 || g.price > bound || out.cancellation_ratio < 1e-7L) {
        const auto quadrature = digital_price(c, d2, v, bound);
        g.price = quadrature.integral;
        out.quadrature_estimated_error = quadrature.error;
        out.quadrature_evaluations = quadrature.evaluations;
        out.price_method = "positive_digital_quadrature";
    }
    if (!(g.price >= 0 && g.price <= bound))
        throw std::runtime_error("Analytic price violates the discounted payoff bound; no clipping was applied.");
    const Real rho_term = std::fma(c.strike, interval2, -c.cap * high_tail);
    g.delta = dq * interval1;
    g.gamma = dq * density / (c.spot * v);
    g.vega = sdq * root_t * density;
    g.rho = t * dr * rho_term;
    g.theta = -sdq * c.volatility * density / (2 * root_t)
              + c.dividend * sdq * interval1 - c.rate * dr * rho_term;
    finite_greeks(g);
    return out;
}

Real HeatGrid::spacing() const { return (upper - lower) / static_cast<Real>(intervals); }

void HeatGrid::validate(const Contract& c) const {
    c.validate();
    require(c.strike > 0 && c.cap > 0 && c.maturity > 0 && c.volatility > 0,
            "The heat PDE requires strictly positive strike, cap, maturity and volatility; use analytic mode for degenerate cases.");
    require(std::isfinite(lower) && std::isfinite(upper) && upper > lower,
            "Heat boundaries must be finite and strictly ordered.");
    require(intervals >= 2 && time_steps >= 1, "At least two space intervals and one time step are required.");
    require(intervals < std::vector<Real>().max_size() - 1,
            "The requested grid exceeds the vector size limit.");
    require(lower < 0 && kink_width(c) < upper, "The heat domain must contain both payoff kinks strictly.");
    const Real y = log_ratio(c.spot, c.strike)
                 + (c.rate - c.dividend - c.volatility * c.volatility / 2) * c.maturity;
    require(std::isfinite(y) && lower <= y && y <= upper,
            "The evaluation coordinate lies outside the fixed heat domain; no clipping is permitted.");
    require(std::isfinite(spacing()) && spacing() > 0,
            "Heat spacing is not representable for the requested domain/grid.");
}

void Bumps::validate(const Contract& c) const {
    for (Real bump : {spot, volatility, rate, maturity})
        require(std::isfinite(bump) && bump > 0, "Every finite-difference bump must be positive and finite.");
    require(spot < c.spot && volatility < c.volatility && maturity < c.maturity,
            "Spot, volatility and maturity bumps must leave strictly positive lower parameters.");
    for (const auto [value, bump] : {std::pair<Real, Real>{c.spot, spot}, {c.volatility, volatility},
                                    {c.rate, rate}, {c.maturity, maturity}})
        require(value - bump < value && value < value + bump,
                "A requested parameter bump is below floating-point resolution.");
}

Bumps default_bumps(const Contract& c, const HeatGrid& grid) {
    grid.validate(c);
    const Real h = grid.spacing(), factor = std::pow(h, 0.25L);
    return {std::min(0.10L * c.spot, c.spot * (0.006L * std::pow(h, 0.125L) + h)),
            std::min(0.25L * c.volatility, 0.10L * c.volatility * factor),
            0.01L * factor, std::min(0.25L * c.maturity, 0.05L * c.maturity * factor)};
}

HeatSolution heat_backward_euler(const Contract& c, const HeatGrid& grid, bool capture_trace) {
    grid.validate(c);
    require(!capture_trace || grid.time_steps < 1000000 / (grid.intervals + 1),
            "Trace capture is limited to one million grid entries.");
    // Prevent size_t overflow and expose the allocation cost before allocating.
    require(grid.intervals < std::vector<Real>().max_size() - 1,
            "The requested grid exceeds the vector size limit.");
    const std::size_t n = grid.intervals - 1;
    HeatSolution out;
    out.nodes.resize(grid.intervals + 1);
    out.values.resize(grid.intervals + 1);
    const Real h = grid.spacing(), kink = kink_width(c);
    for (std::size_t j = 0; j <= grid.intervals; ++j) {
        const Real y = j == grid.intervals ? grid.upper : grid.lower + static_cast<Real>(j) * h;
        out.nodes[j] = y;
        out.values[j] = y <= 0 ? 0 : y >= kink ? c.cap : c.strike * std::expm1(y);
    }
    out.values.front() = 0;
    out.values.back() = c.cap;
    if (capture_trace) out.trace.push_back(out.values);
    const Real time_step = c.maturity / static_cast<Real>(grid.time_steps);
    require(std::isfinite(time_step) && time_step > 0,
            "The positive mathematical time step is not representable.");
    out.lambda = c.volatility * c.volatility / 2 * time_step / (h * h);
    require(std::isfinite(out.lambda) && out.lambda > 0 && std::isfinite(1 + 2 * out.lambda),
            "Backward-Euler coefficients exceed the supported floating-point range.");
    const Real lam = out.lambda;
    const Real boundary_forcing = lam * c.cap;
    require(std::isfinite(boundary_forcing), "The boundary forcing exceeds the supported floating-point range.");
    std::vector<Real> cp(n, 0), inverse(n, 0), rhs(n, 0);
    inverse[0] = 1 / (1 + 2 * lam);
    if (n > 1) cp[0] = -lam * inverse[0];
    for (std::size_t j = 1; j < n; ++j) {
        const Real pivot = 1 + 2 * lam + lam * cp[j-1];
        if (!(pivot > 0 && std::isfinite(pivot)))
            throw std::runtime_error("Thomas factorisation produced an invalid pivot.");
        inverse[j] = 1 / pivot;
        if (j < n - 1) cp[j] = -lam * inverse[j];
    }
    for (std::size_t step = 0; step < grid.time_steps; ++step) {
        for (std::size_t j = 0; j < n; ++j) rhs[j] = out.values[j+1];
        rhs.back() += boundary_forcing;
        finite_result(rhs.back(), "The boundary right-hand side exceeds the supported floating-point range.");
        rhs.front() *= inverse.front();
        for (std::size_t j = 1; j < n; ++j) {
            rhs[j] = (rhs[j] + lam * rhs[j-1]) * inverse[j];
            finite_result(rhs[j], "The Thomas forward recurrence exceeds the supported floating-point range.");
        }
        for (std::size_t j = n - 1; j-- > 0;) rhs[j] -= cp[j] * rhs[j+1];
        for (std::size_t j = 0; j < n; ++j) out.values[j+1] = rhs[j];
        if (capture_trace) out.trace.push_back(out.values);
    }
    for (Real x : out.values) finite_result(x, "Heat solver produced a non-finite value.");
    return out;
}

Real pde_price(const Contract& c, const HeatGrid& grid) {
    return interpolated_price(c, grid, heat_backward_euler(c, grid));
}

Real domain_price_bound(const Contract& c, const HeatGrid& grid) {
    grid.validate(c);
    const Real y = log_ratio(c.spot, c.strike)
                 + (c.rate - c.dividend - c.volatility * c.volatility / 2) * c.maturity;
    const Real sd = c.volatility * std::sqrt(c.maturity), h = grid.spacing();
    const Real coordinate = (y - grid.lower) / h;
    const auto left = coordinate >= static_cast<Real>(grid.intervals) ? grid.intervals - 1
                     : static_cast<std::size_t>(std::floor(coordinate));
    Real exit_bound = 0;
    for (std::size_t j : {left, left + 1}) {
        const Real node = grid.lower + static_cast<Real>(j) * h;
        const Real probability_bound = std::min(1.0L, 2 * normal_cdf(-(node - grid.lower) / sd)
                                                + 2 * normal_cdf(-(grid.upper - node) / sd));
        exit_bound = std::max(exit_bound, probability_bound);
    }
    return discount(c.rate, c.maturity) * c.cap * exit_bound;
}

PDEResult pde(const Contract& c, const HeatGrid& grid, const std::optional<Bumps>& requested) {
    grid.validate(c);
    PDEResult out;
    out.bumps = requested.value_or(default_bumps(c, grid));
    out.bumps.validate(c);
    const auto heat = heat_backward_euler(c, grid);
    out.lambda = heat.lambda;
    out.heat_min = *std::min_element(heat.values.begin(), heat.values.end());
    out.heat_max = *std::max_element(heat.values.begin(), heat.values.end());
    out.base_domain_bound = domain_price_bound(c, grid);
    out.values = central_greeks(c, out.bumps, [&](const Contract& bumped) { return pde_price(bumped, grid); },
                                interpolated_price(c, grid, heat));
    out.analytic_same_bumps = central_greeks(c, out.bumps, [](const Contract& bumped) { return analytic(bumped).values.price; });
    return out;
}

} // namespace qbp
