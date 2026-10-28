"""Outward price images for a capped Black--Scholes claim.

Every input denotes an exact rational number: a float denotes its binary64
value, whereas a decimal string denotes that decimal rational. Arb supplies
outward transcendental evaluation. The price-image primitive certifies a
supplied volatility interval. The optional Gaussian bridge also encloses
chi-squared quantiles and rounding-cell scatter, under an explicit latent
Gaussian sampling premise; neither procedure validates the market model.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
import math
from typing import Any

import flint
from flint import arb, ctx

F = Fraction


def rational(value: Any, name: str = "value") -> Fraction:
    """Parse finite scalar data without a decimal-to-binary64 round trip."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite rational scalar.")
    try:
        result = F(value)
    except (ValueError, TypeError, OverflowError, ZeroDivisionError) as error:
        raise ValueError(f"{name} must be a finite rational scalar.") from error
    return result


@dataclass(frozen=True)
class ExactCappedGBM:
    """Known positive spot/strike/cap/maturity and real continuous rates."""
    spot: Fraction = F(100)
    strike: Fraction = F(100)
    cap: Fraction = F(20)
    maturity: Fraction = F(1)
    rate: Fraction = F(3, 100)
    dividend: Fraction = F(1, 100)

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, rational(getattr(self, name), name))
        for name in ("spot", "strike", "cap", "maturity"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be strictly positive.")

    def exact_inputs(self) -> dict[str, str]:
        return {name: str(getattr(self, name)) for name in self.__dataclass_fields__}


def _arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def _ball(lower: Fraction, upper: Fraction) -> arb:
    if lower > upper:
        raise ValueError("Reversed rational interval.")
    return _arb((lower + upper) / 2) + arb(0, _arb((upper - lower) / 2).upper())


def _endpoint(value: arb, upwards: bool, bits: int) -> Fraction:
    """Coarsen an exact dyadic Arb endpoint outwards to denominator 2**bits.

    Avoid materialising enormous denominators for negligible Gaussian tails.
    The floor/ceiling operation uses only exact integers, including negatives.
    """
    if not value.is_finite():
        raise ArithmeticError("Arb did not return a finite price enclosure.")
    endpoint = value.upper() if upwards else value.lower()
    mantissa, exponent = map(int, endpoint.man_exp())
    shift = exponent + bits
    if shift >= 0:
        if shift > 100000:
            raise ArithmeticError("Endpoint exceeds the declared serialisation range.")
        integer = mantissa << shift
    else:
        width = -shift
        if width > abs(mantissa).bit_length() + 1:
            floor = -int(mantissa < 0)
            exact = mantissa == 0
        else:
            denominator = 1 << width
            floor, remainder = divmod(mantissa, denominator)
            exact = remainder == 0
        integer = floor + int(upwards and not exact)
    return F(integer, 1 << bits)


def _bounds(value: arb, bits: int) -> tuple[Fraction, Fraction]:
    return _endpoint(value, False, bits), _endpoint(value, True, bits)


def _serial(interval: tuple[Fraction, Fraction]) -> dict[str, str]:
    lower, upper = interval
    return {"lower_exact_fraction": str(lower), "upper_exact_fraction": str(upper),
            "width_exact_fraction": str(upper - lower)}


def _normal_cdf(value: arb) -> arb:
    return (-value / arb(2).sqrt()).erfc() / 2


def _price(contract: ExactCappedGBM, volatility: arb) -> arb:
    """Capped call via a tail-aware normal CDF and an exact strike spread.

    Precision refinement resolves cancellation; no floating clipping or
    heuristic error estimate is substituted for an enclosure.
    """
    if not volatility.is_finite() or not volatility > 0:
        raise ArithmeticError("A strictly positive finite volatility ball is required.")
    c = contract
    s, k, h, t, r, q = map(_arb, (c.spot, c.strike, c.cap, c.maturity,
                                c.rate, c.dividend))
    root_time = t.sqrt()
    v = volatility * root_time
    drift = (r - q + volatility * volatility / 2) * t
    d1k = ((s / k).log() + drift) / v
    d1j = d1k - (h / k).log1p() / v
    d2k, d2j = d1k - v, d1j - v
    # An erfc difference avoids cancellation of two CDFs near one. Each
    # branch computes the identical real expression on the whole input ball.
    def difference(first: arb, second: arb) -> arb:
        if second > 0:
            return (_normal_cdf(-second) - _normal_cdf(-first))
        return _normal_cdf(first) - _normal_cdf(second)
    return (s * (-q * t).exp() * difference(d1k, d1j)
            - k * (-r * t).exp() * difference(d2k, d2j)
            + h * (-r * t).exp() * _normal_cdf(d2j))


def _critical_square(contract: ExactCappedGBM) -> arb:
    c = contract
    return ((2 * (_arb(c.strike) / _arb(c.spot)).log()
             + (_arb(c.cap) / _arb(c.strike)).log1p()) / _arb(c.maturity)
            - 2 * _arb(c.rate - c.dividend))


def _options(absolute_tolerance: Any, start_bits: int,
             max_bits: int) -> tuple[Fraction, int]:
    tolerance = rational(absolute_tolerance, "absolute_tolerance")
    if tolerance <= 0:
        raise ValueError("absolute_tolerance must be strictly positive.")
    if (isinstance(start_bits, bool) or isinstance(max_bits, bool)
            or not isinstance(start_bits, int) or not isinstance(max_bits, int)
            or start_bits < 64 or max_bits < start_bits or max_bits > 16384):
        raise ValueError("Precision limits must be integers with 64 <= start <= max <= 16384.")
    # Endpoint coarsening contributes at most 2**(-bits) per boundary.
    endpoint_bits = max(128, tolerance.denominator.bit_length()
                        - tolerance.numerator.bit_length() + 16)
    return tolerance, endpoint_bits


def certified_price(contract: ExactCappedGBM, volatility: Any, *,
                    absolute_tolerance: Any = "1e-24", start_bits: int = 128,
                    max_bits: int = 4096) -> dict[str, Any]:
    """Return rational endpoints enclosing the exact model price.

    Failure to meet the requested width raises, rather than returning an
    uncertified answer. Width is an achieved enclosure width, not an estimate.
    """
    sigma = rational(volatility, "volatility")
    if sigma <= 0:
        raise ValueError("volatility must be strictly positive.")
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    bits = start_bits
    while True:
        with ctx.workprec(bits):
            enclosure = _bounds(_price(contract, _arb(sigma)), endpoint_bits)
        if enclosure[1] - enclosure[0] <= tolerance:
            return {"schema": "qbp-validated-gbm-price-v1",
                    "contract_exact": contract.exact_inputs(),
                    "volatility_exact": str(sigma), "price": _serial(enclosure),
                    "precision_bits": bits, "endpoint_dyadic_bits": endpoint_bits,
                    "requested_width_exact": str(tolerance),
                    "python_flint_version": flint.__version__,
                    "target": "Exact capped GBM price at declared rational inputs."}
        if bits == max_bits:
            raise ArithmeticError("Requested price enclosure width was not attained.")
        bits = min(2 * bits, max_bits)


def certified_price_image(contract: ExactCappedGBM, volatility_lower: Any,
                          volatility_upper: Any, *,
                          absolute_tolerance: Any = "1e-24", start_bits: int = 128,
                          max_bits: int = 4096) -> dict[str, Any]:
    """Enclose both exact global extrema on a supplied positive interval.

    No sampled optimisation is used. The exact vega sign is sign(c*-sigma²),
    so the minimum is at an endpoint and at most one interior maximum exists.
    The maximiser is sqrt(clamp(c*, lower², upper²)); boundary and nearly
    stationary cases require no unresolved branch decision.
    """
    lower, upper = (rational(volatility_lower, "volatility_lower"),
                    rational(volatility_upper, "volatility_upper"))
    if not 0 < lower <= upper:
        raise ValueError("Volatility endpoints must be positive and ordered.")
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    bits = start_bits
    while True:
        with ctx.workprec(bits):
            first = _bounds(_price(contract, _arb(lower)), endpoint_bits)
            second = _bounds(_price(contract, _arb(upper)), endpoint_bits)
            critical_bits = max(endpoint_bits, bits + 16)
            critical = _bounds(_critical_square(contract), critical_bits)
            minimum = min(first[0], second[0]), min(first[1], second[1])
            # The exact maximiser is sqrt(clamp(c*, lower², upper²)).
            # Clamping is monotone, including equal/boundary cases; it needs
            # no numerical root-classification decision.
            clamp = lambda value: min(max(value, lower * lower), upper * upper)
            peak_square = clamp(critical[0]), clamp(critical[1])
            peak_ball = _ball(*peak_square).sqrt()
            if not peak_ball.is_finite() or not peak_ball > 0:
                raise ArithmeticError("The clamped maximiser has no positive finite Arb enclosure.")
            peak_location = _bounds(peak_ball, endpoint_bits)
            peak_value = _bounds(_price(contract, peak_ball), endpoint_bits)
            maximum = peak_value
            if lower == upper:
                classification = "singleton"
            elif critical[1] <= lower * lower:
                classification = "decreasing_on_interval"
            elif critical[0] >= upper * upper:
                classification = "increasing_on_interval"
            elif lower * lower < critical[0] and critical[1] < upper * upper:
                classification = "strict_interior_maximum"
            else:
                classification = "boundary_stationary_point_enclosed"
        if max(minimum[1] - minimum[0], maximum[1] - maximum[0]) <= tolerance:
            return {"schema": "qbp-validated-gbm-image-v1",
                    "contract_exact": contract.exact_inputs(),
                    "volatility_interval_exact": [str(lower), str(upper)],
                    "minimum_price": _serial(minimum), "maximum_price": _serial(maximum),
                    "image_outer": _serial((minimum[0], maximum[1])),
                    "endpoint_price_enclosures": [_serial(first), _serial(second)],
                    "critical_square": _serial(critical),
                    "internal_clamp_dyadic_bits": critical_bits,
                    "classification": classification,
                    "possible_peak_location": None if peak_location is None else _serial(peak_location),
                    "possible_peak_price": None if peak_value is None else _serial(peak_value),
                    "precision_bits": bits, "endpoint_dyadic_bits": endpoint_bits,
                    "requested_extremum_width_exact": str(tolerance),
                    "python_flint_version": flint.__version__,
                    "statistical_coverage_claim": False,
                    "target": "Exact model price image of the supplied rational volatility interval.",
                    "scope": "Known contract/rates, stipulated constant-volatility GBM; no empirical or quantum-advantage claim."}
        if bits == max_bits:
            raise ArithmeticError("Requested extremum enclosure widths were not attained.")
        bits = min(2 * bits, max_bits)


def _integer(value: Any, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}.")
    return value


def _chi_square_cdf(value: Fraction, degrees_freedom: int) -> arb:
    if value == 0:
        return arb(0)
    return (_arb(value) / 2).gamma_lower(_arb(F(degrees_freedom, 2)), regularized=1)


@lru_cache(maxsize=128)
def _quantile_core(degrees_freedom: int, probability: Fraction,
                   tolerance: Fraction, start_bits: int, max_bits: int,
                   max_iterations: int, endpoint_bits: int) -> tuple:
    """Exact bisection with every ordering decision proved by an Arb ball."""
    precision_used = start_bits
    calls = 0

    def sign(value: Fraction) -> int:
        nonlocal precision_used, calls
        bits = start_bits
        while True:
            if calls >= max_iterations:
                raise ArithmeticError("Chi-squared CDF evaluation budget exhausted.")
            calls += 1
            with ctx.workprec(bits):
                cdf = _chi_square_cdf(value, degrees_freedom)
                if not cdf.is_finite():
                    raise ArithmeticError("Nonfinite chi-squared CDF enclosure.")
                p = _arb(probability)
                decision = -1 if cdf < p else 1 if cdf > p else 0
            precision_used = max(precision_used, bits)
            if decision:
                return decision
            if bits == max_bits:
                raise ArithmeticError("Chi-squared CDF ordering remained unresolved.")
            bits = min(2 * bits, max_bits)

    lower, upper = F(0), F(max(1, degrees_freedom))
    while sign(upper) < 0:
        lower, upper = upper, 2 * upper
    while upper - lower > tolerance:
        midpoint = (lower + upper) / 2
        if sign(midpoint) < 0:
            lower = midpoint
        else:
            upper = midpoint
    # The returned rational witnesses must themselves prove the bracket.
    # A valid internal comparison is insufficient if endpoint coarsening
    # subsequently erases its sign. Grow witness resolution with precision,
    # check the exact rational inequalities, and refine or abstain explicitly.
    witness_bits = max(precision_used, start_bits)
    while True:
        if calls + 2 > max_iterations:
            raise ArithmeticError("Chi-squared CDF witness budget exhausted.")
        calls += 2
        with ctx.workprec(witness_bits):
            cdf_bits = max(endpoint_bits, witness_bits + 16)
            lo_cdf = _bounds(_chi_square_cdf(lower, degrees_freedom), cdf_bits)
            hi_cdf = _bounds(_chi_square_cdf(upper, degrees_freedom), cdf_bits)
        precision_used = max(precision_used, witness_bits)
        if lo_cdf[1] < probability < hi_cdf[0]:
            break
        if witness_bits == max_bits:
            raise ArithmeticError("Retained chi-squared CDF witnesses did not prove the strict bracket.")
        witness_bits = min(2 * witness_bits, max_bits)
    return lower, upper, lo_cdf, hi_cdf, precision_used, calls


def certified_chi_square_quantile(degrees_freedom: int, probability: Any, *,
                                  absolute_tolerance: Any = "1e-24",
                                  start_bits: int = 128, max_bits: int = 4096,
                                  max_iterations: int = 4096) -> dict[str, Any]:
    """Certify the unique positive chi-squared quantile by CDF bracketing.

    Bisection invariants are F(lower) < p < F(upper); the lower boundary may
    be exactly zero. No SciPy or ordinary floating-point quantile is used.
    The evaluation budget includes precision refinements, upper expansion and
    the final rational CDF witnesses, whose strict inequalities are checked.
    """
    nu = _integer(degrees_freedom, "degrees_freedom", 1)
    budget = _integer(max_iterations, "max_iterations", 1)
    p = rational(probability, "probability")
    if not 0 < p < 1:
        raise ValueError("probability must lie strictly between zero and one.")
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    lower, upper, lo_cdf, hi_cdf, bits, calls = _quantile_core(
        nu, p, tolerance, start_bits, max_bits, budget, endpoint_bits)
    return {"schema": "qbp-validated-chi-square-quantile-v1",
            "degrees_freedom": nu, "probability_exact": str(p),
            "quantile": _serial((lower, upper)),
            "cdf_at_lower": _serial(lo_cdf), "cdf_at_upper": _serial(hi_cdf),
            "strict_bracket_invariant": "F(lower) < probability < F(upper)",
            "precision_bits": bits, "cdf_evaluations": calls,
            "requested_width_exact": str(tolerance),
            "python_flint_version": flint.__version__}


def _rational_sqrt_bounds(value: Fraction, bits: int) -> tuple[Fraction, Fraction]:
    """Outward dyadic square-root bounds using only integer arithmetic."""
    if value < 0:
        raise ValueError("Square-root input must be nonnegative.")
    scale = 1 << bits
    floor = math.isqrt((value.numerator << (2 * bits)) // value.denominator)
    lower = F(floor, scale)
    upper = lower if lower * lower == value else F(floor + 1, scale)
    return lower, upper


def binary64_rounding_radius(value: float) -> Fraction:
    """Conservative nearest-even binary64 rounding-cell radius, exactly.

    For a finite binary64 result z, the radius is the larger adjacent half
    spacing. At the largest finite magnitude the missing finite neighbour
    is continued by one spacing, covering the finite overflow threshold.
    This concerns one final rounding operation; it does not certify random
    generation, preceding arithmetic, data-provider rounding or model fit.
    """
    if isinstance(value, bool) or not isinstance(value, float) or not math.isfinite(value):
        raise ValueError("A finite Python float is required.")
    centre = F.from_float(value)
    previous, following = math.nextafter(value, -math.inf), math.nextafter(value, math.inf)
    left = F.from_float(previous) if math.isfinite(previous) else 2 * centre - F.from_float(following)
    right = F.from_float(following) if math.isfinite(following) else 2 * centre - F.from_float(previous)
    return max(centre - left, right - centre) / 2


def exact_scatter_enclosure(observations: Any, rounding_radii: Any = None, *,
                            sqrt_bits: int = 192) -> dict[str, Any]:
    """Enclose latent centred scatter when each supplied |R_i-z_i| <= eps_i.

    The Euclidean orthogonal projection P=I-11'/n is a contraction, so
    |sqrt(C_R)-sqrt(C_z)| <= ||R-z||_2 <= sqrt(sum eps_i²).
    All observation arithmetic and reported bounds use exact rationals.
    """
    _integer(sqrt_bits, "sqrt_bits", 64)
    try:
        values = [rational(x, "observation") for x in observations]
    except TypeError as error:
        raise ValueError("observations must be an iterable of finite rational scalars.") from error
    if len(values) < 2:
        raise ValueError("At least two observations are required.")
    if rounding_radii is None:
        radii = [F(0)] * len(values)
    else:
        try:
            radii = [rational(x, "rounding_radius") for x in rounding_radii]
        except TypeError as error:
            raise ValueError("rounding_radii must match observations.") from error
        if len(radii) != len(values) or any(x < 0 for x in radii):
            raise ValueError("Nonnegative rounding_radii must match observations.")
    count = len(values)
    total = sum(values, F(0))
    scatter = sum((x * x for x in values), F(0)) - total * total / count
    radius_square = sum((x * x for x in radii), F(0))
    if radius_square == 0:
        lower = upper = scatter
    else:
        scatter_root = _rational_sqrt_bounds(scatter, sqrt_bits)
        error_root = _rational_sqrt_bounds(radius_square, sqrt_bits)
        lower = max(F(0), scatter_root[0] - error_root[1]) ** 2
        upper = (scatter_root[1] + error_root[1]) ** 2
    return {"observations": count, "observed_scatter_exact": str(scatter),
            "rounding_l2_radius_squared_exact": str(radius_square),
            "latent_scatter": _serial((lower, upper)), "sqrt_bits": sqrt_bits,
            "input_relation": "Each latent observation R_i satisfies |R_i-z_i| <= supplied epsilon_i.",
            "zero_rounding_radii": radius_square == 0}


def certified_gaussian_price_interval(contract: ExactCappedGBM, observations: Any,
                                      observation_interval: Any,
                                      alpha: Any = "1/20", rounding_radii: Any = None, *,
                                      absolute_tolerance: Any = "1e-24",
                                      quantile_tolerance: Any = "1e-24",
                                      start_bits: int = 128, max_bits: int = 4096,
                                      max_iterations: int = 4096) -> dict[str, Any]:
    """Outward scalar Gaussian-volatility interval and exact capped price image.

    The statistical premise concerns latent iid N(m*Delta,sigma²*Delta)
    observations with unknown mean, not their rounded discrete outputs.
    quantile_tolerance has dimensionless chi-squared units; absolute_tolerance
    is the monetary width requested for each price-extremum enclosure.
    Each latent observation must belong to its supplied deterministic cell.
    On a successful numerical return, containment follows whenever the exact
    two-sided chi-squared pivot interval covers sigma. Thus
    P(success and returned interval excludes the true model price) <= alpha.
    This is not P(exclusion | success) <= alpha; abstentions remain outcomes.
    Malformed inputs raise ValueError. Exhausted numerical budgets and a
    nonpositive scatter lower bound return an explicit abstention record.
    """
    delta = rational(observation_interval, "observation_interval")
    a = rational(alpha, "alpha")
    if delta <= 0 or not 0 < a < 1:
        raise ValueError("Observation interval must be positive and alpha in (0,1).")
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    _integer(max_iterations, "max_iterations", 1)
    q_tolerance, _ = _options(quantile_tolerance, start_bits, max_bits)
    scatter = exact_scatter_enclosure(observations, rounding_radii,
                                      sqrt_bits=max(192, endpoint_bits))
    result = {"schema": "qbp-validated-gaussian-price-interval-v1",
              "status": "abstained", "alpha_exact": str(a),
              "observation_interval_exact": str(delta),
              "requested_price_extremum_width_exact": str(tolerance),
              "requested_quantile_width_exact": str(q_tolerance),
              "contract_exact": contract.exact_inputs(), "scatter": scatter,
              "coverage_statement": "P(success and returned interval excludes true model price) <= alpha, under stated premises.",
              "premises": "Latent independent exactly Gaussian log returns at equal exact spacing, constant positive diffusion shared under P and Q, fixed known contract/rates, and all latent observations in the supplied deterministic rounding cells.",
              "scope": "Scalar model uncertainty only; no market-validation, Greek, basket or quantum-advantage claim."}
    lo_scatter = F(scatter["latent_scatter"]["lower_exact_fraction"])
    hi_scatter = F(scatter["latent_scatter"]["upper_exact_fraction"])
    if lo_scatter <= 0:
        result["abstention_reason"] = "Latent scatter lower bound is not strictly positive."
        return result
    try:
        nu = scatter["observations"] - 1
        lower_q = certified_chi_square_quantile(nu, a / 2, absolute_tolerance=q_tolerance,
                                               start_bits=start_bits, max_bits=max_bits,
                                               max_iterations=max_iterations)
        upper_q = certified_chi_square_quantile(nu, 1 - a / 2, absolute_tolerance=q_tolerance,
                                               start_bits=start_bits, max_bits=max_bits,
                                               max_iterations=max_iterations)
        result["lower_tail_quantile"] = lower_q
        result["upper_tail_quantile"] = upper_q
        q_low = F(lower_q["quantile"]["lower_exact_fraction"])
        q_high = F(upper_q["quantile"]["upper_exact_fraction"])
        if q_low <= 0:
            raise ArithmeticError("Lower-tail quantile lower bound is not strictly positive.")
        variance = lo_scatter / (delta * q_high), hi_scatter / (delta * q_low)
        volatility = (_rational_sqrt_bounds(variance[0], max(192, endpoint_bits))[0],
                      _rational_sqrt_bounds(variance[1], max(192, endpoint_bits))[1])
        if volatility[0] <= 0:
            raise ArithmeticError("Outward volatility lower bound is not strictly positive.")
        result["variance_interval"] = _serial(variance)
        result["volatility_interval_exact"] = [str(x) for x in volatility]
        result["price_image"] = certified_price_image(contract, *volatility,
                                                       absolute_tolerance=tolerance,
                                                       start_bits=start_bits,
                                                       max_bits=max_bits)
        result["status"] = "certified"
        result["numerical_failure_probability"] = 0
    except (ArithmeticError, OverflowError) as error:
        result["abstention_reason"] = str(error)
    return result
