#include "qbp/capped_gbm.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
unsigned assertions = 0;

void check(bool condition, const std::string& message) {
    ++assertions;
    if (!condition) throw std::runtime_error(message);
}

void close(qbp::Real actual, qbp::Real expected, qbp::Real absolute,
           const std::string& message) {
    check(std::abs(actual - expected) <= absolute, message);
}

void rejects(const std::function<void()>& action, const std::string& message) {
    try { action(); } catch (const std::exception&) { ++assertions; return; }
    throw std::runtime_error(message);
}

void analytic_checks() {
    qbp::Contract c;
    auto a = qbp::analytic(c).values;
    close(a.price, 6.3057373074164L, 1e-12L, "Default price reference");
    const qbp::Real residual = *a.theta + c.volatility * c.volatility * c.spot * c.spot * *a.gamma / 2
                            + (c.rate - c.dividend) * c.spot * *a.delta - c.rate * a.price;
    close(residual, 0, 1e-12L, "Dividend-paying Black-Scholes PDE identity");
    for (qbp::Real spot : {0.0001L, 100.0L, 150.0L, 1e20L}) {
        c.spot = spot;
        a = qbp::analytic(c).values;
        check(a.price >= 0 && a.price <= c.cap * std::exp(-c.rate * c.maturity), "Capped payoff price bounds");
        check(*a.delta >= 0, "Delta non-negativity");
    }
    c.spot = 1e20L;
    close(qbp::analytic(c).values.price, c.cap * std::exp(-c.rate * c.maturity), 1e-12L,
          "Deep in-the-money capped price avoids call-price cancellation");
    c.spot = 150;
    check(*qbp::analytic(c).values.gamma < 0 && *qbp::analytic(c).values.vega < 0,
          "Capped gamma and vega may be negative");
    c = qbp::Contract{};
    c.cap = 1e-14L;
    a = qbp::analytic(c).values;
    check(a.price > 0 && *a.delta > 0, "Narrow cap survives rounded upper-strike cancellation");
    c.cap = 20;
    c.spot = 99.999995L;
    c.rate = c.dividend = 0;
    c.volatility = 1e-8L;
    const auto fallback = qbp::analytic(c);
    check(fallback.price_method == "positive_digital_quadrature",
          "Cancellation stress uses the disclosed positive-integral fallback");
    check(fallback.values.price > 0 && fallback.quadrature_evaluations > 0,
          "Fallback retains a small positive price and evaluation count");
    check(fallback.quadrature_estimated_error >= 0 && fallback.quadrature_estimated_error < 1e-25L,
          "Fallback reports its numerical error indicator");
    c = qbp::Contract{};
    c.spot = c.cap = 100;
    c.strike = 1e-300L;
    c.rate = c.dividend = 0;
    c.volatility = 1e-8L;
    const auto broad = qbp::analytic(c).values;
    close(*broad.delta, qbp::normal_cdf(-0.5e-8L), 1e-14L,
          "Tiny lower strike retains the independently computed upper-strike coordinate");
    c = qbp::Contract{};
    c.spot = c.strike = 1e-100L;
    c.cap = 1e300L;
    c.rate = c.dividend = 0;
    c.volatility = 1e-8L;
    rejects([&] { qbp::analytic(c); },
            "Fallback rejects an unrepresentable cap/strike normalisation");
    c.spot = c.strike = 1;
    c.volatility = 1e-100L;
    rejects([&] { qbp::analytic(c); },
            "Fallback rejects transition mass below normalised floating-point resolution");
    c = qbp::Contract{};
    c.spot = 1e300L;
    c.strike = 1e308L;
    c.cap = 1e307L;
    c.rate = c.dividend = 0;
    c.volatility = 0.46L;
    rejects([&] { qbp::analytic(c); },
            "Underflowed normal probabilities must not silently erase representable scaled prices or Greeks");
    c = qbp::Contract{};
    c.spot = c.strike = 1e300L;
    c.cap = 1e300L;
    c.rate = 1000;
    c.dividend = 0;
    rejects([&] { qbp::analytic(c); },
            "Underflowed discount factors must not silently erase a representable discounted large cap");
    c = qbp::Contract{};
    c.spot = 109;
    c.rate = -0.01L;
    c.dividend = 0.045L;
    c.volatility = 0.31L;
    a = qbp::analytic(c).values;
    for (unsigned parameter = 0; parameter < 4; ++parameter) {
        qbp::Contract low = c, high = c;
        qbp::Real qbp::Contract::* field = parameter == 0 ? &qbp::Contract::spot : parameter == 1 ? &qbp::Contract::volatility
                                       : parameter == 2 ? &qbp::Contract::rate : &qbp::Contract::maturity;
        const qbp::Real step = parameter == 0 ? 0.001L : 0.000001L;
        low.*field -= step;
        high.*field += step;
        const auto low_price = qbp::analytic(low).values.price, high_price = qbp::analytic(high).values.price;
        const qbp::Real derivative = (high_price - low_price) / (2 * step);
        const qbp::Real expected = parameter == 0 ? *a.delta : parameter == 1 ? *a.vega : parameter == 2 ? *a.rho : -*a.theta;
        close(derivative, expected, 2e-7L, "Analytic Greek against independent price bump");
        if (parameter == 0) close((high_price - 2 * a.price + low_price) / (step * step), *a.gamma,
                                  5e-7L, "Analytic gamma against independent price bump");
    }
}

void degenerate_checks() {
    qbp::Contract c;
    c.cap = 0;
    auto a = qbp::analytic(c).values;
    check(a.price == 0 && *a.delta == 0 && *a.gamma == 0 && *a.vega == 0,
          "Zero-cap claim is identically zero");
    c = qbp::Contract{};
    c.maturity = 0;
    c.spot = 110;
    a = qbp::analytic(c).values;
    close(a.price, 10, 0, "Expiry payoff");
    close(*a.delta, 1, 0, "Expiry active-region delta");
    close(*a.theta, -1.9L, 1e-15L, "Expiry active-region right-maturity theta");
    c.spot = c.strike;
    a = qbp::analytic(c).values;
    check(!a.delta && !a.gamma && !a.theta && a.vega == 0 && a.rho == 0,
          "Expiry kink has undefined spot and calendar-time derivatives");
    c = qbp::Contract{};
    c.volatility = 0;
    c.rate = c.dividend = 0;
    a = qbp::analytic(c).values;
    check(!a.delta && !a.gamma && !a.rho && a.theta == 0,
          "Deterministic kink has undefined spot/rate derivatives and an exact constant-drift theta");
    close(*a.vega, 100 / std::sqrt(2 * std::acos(-1.0L)), 1e-12L,
          "Deterministic lower-kink right vega is positive");
    c.spot = 120;
    close(*qbp::analytic(c).values.vega, -120 / std::sqrt(2 * std::acos(-1.0L)), 1e-12L,
          "Deterministic upper-kink right vega is negative");
    c = qbp::Contract{};
    c.strike = 0;
    a = qbp::analytic(c).values;
    check(a.price > 0 && a.price <= c.cap * std::exp(-c.rate * c.maturity), "Zero-strike positive price");
    const qbp::Real residual = *a.theta + c.volatility * c.volatility * c.spot * c.spot * *a.gamma / 2
                            + (c.rate - c.dividend) * c.spot * *a.delta - c.rate * a.price;
    close(residual, 0, 1e-12L, "Zero-strike PDE identity");
    rejects([&] { qbp::pde(c, {}); }, "PDE rejects zero strike outside its transformation assumptions");
}

void pde_checks() {
    qbp::Contract c;
    c.spot = 110;
    qbp::HeatGrid coarse{-2, 2, 64, 512}, fine{-2, 2, 256, 2048};
    const auto reference = qbp::analytic(c).values;
    check(std::abs(qbp::pde_price(c, fine) - reference.price)
          < std::abs(qbp::pde_price(c, coarse) - reference.price) / 3,
          "Independent Thomas implementation refines towards analytic price");
    auto result = qbp::pde(c, fine);
    check(result.values.gamma.has_value() && result.values.theta.has_value(), "Bumped PDE Greeks are available");
    check(result.heat_min >= -1e-10L && result.heat_max <= c.cap + 1e-10L, "Heat maximum principle");
    check(result.base_domain_bound < 1e-15L, "Wide-domain tail bound is small");
    qbp::HeatGrid small{-0.6L, 0.7L, 512, 1};
    c.volatility = 0.7L;
    c.maturity = 2;
    const auto heat = qbp::heat_backward_euler(c, small);
    check(heat.lambda > 10000, "Large-lambda case exceeds explicit Euler stability range");
    for (qbp::Real value : heat.values) check(value >= -1e-10L && value <= c.cap + 1e-10L, "Large-lambda bounded heat solution");
    c = qbp::Contract{};
    qbp::HeatGrid trace_grid{-1, 1, 16, 3};
    const auto trace = qbp::heat_backward_euler(c, trace_grid, true);
    check(trace.trace.size() == 4 && trace.trace.front().size() == 17, "Trace retains initial plus every time step");
    for (const auto& row : trace.trace) check(row.front() == 0 && row.back() == c.cap, "Trace has exact supplied boundary values");
    rejects([&] { qbp::heat_backward_euler(c, qbp::HeatGrid{-1, 1, 999999, 1}, true); },
            "Trace allocation limit includes its initial row");
    rejects([&] { qbp::pde(c, qbp::HeatGrid{-0.1L, 0.1L, 32, 32}); }, "PDE rejects a domain omitting the upper kink");
    c.spot = 1e10L;
    rejects([&] { qbp::pde(c, fine); }, "PDE does not silently clip a query outside its domain");
}

} // namespace

int main() {
    try {
        analytic_checks();
        degenerate_checks();
        pde_checks();
        qbp::Contract invalid;
        invalid.spot = -1;
        rejects([&] { qbp::analytic(invalid); }, "Invalid spot rejected");
        std::cout << "{\"native_assertions\":" << assertions << ",\"status\":\"passed\"}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Native check failed after " << assertions << " assertions: " << error.what() << '\n';
        return 1;
    }
}
