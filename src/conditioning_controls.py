"""Exact GBM conditioning and independently fitted stock control variates.

The financial model and conditional expectations are classical. Numerical
moment checks use floating-point Gaussian quadrature, not validated bounds.
Only NumPy and the standard library are required, alongside capped_gbm.py.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from typing import Iterable

import numpy as np

from capped_gbm import Contract, analytic, payoff, normal_cdf


METHODS = (
    "raw_payoff", "exact_conditional", "raw_pilot_stock_control",
    "conditional_raw_pilot_coefficient", "conditional_own_pilot_coefficient",
    "conditional_unconditioned_stock_control",
)


def _split_check(u: float) -> None:
    if not math.isfinite(u) or not 0.0 <= u <= 1.0:
        raise ValueError("The conditioning fraction must be in [0, 1].")


def _cdf(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.fromiter((0.5 * math.erfc(-v / math.sqrt(2.0)) for v in x.flat),
                       dtype=float, count=x.size).reshape(x.shape)


def _cdf_interval(lower: np.ndarray, upper: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Phi(upper)−Phi(lower) and Phi(lower), sharing two CDF evaluations.

    Subtract survival functions in the upper tail to avoid subtracting two
    values rounded to one. Both arrays must have matching/broadcastable shapes.
    """
    upper_tail = lower > 0.0
    first = _cdf(np.where(upper_tail, -lower, lower))
    second = _cdf(np.where(upper_tail, -upper, upper))
    return (np.where(upper_tail, first - second, second - first),
            np.where(upper_tail, 1.0 - first, first))


def conditional_value(c: Contract, u: float, state: np.ndarray,
                      diagnostics: dict | None = None) -> np.ndarray:
    """Discount to time zero the analytic residual-life call spread.

    u=0 returns the known analytic price. u=1 evaluates the terminal payoff.
    The interior calculation requires four normal-CDF evaluations per state.
    """
    _split_check(u)
    state = np.asarray(state, dtype=float)
    if np.any(~np.isfinite(state)) or np.any(state <= 0.0):
        raise ValueError("Stock states must be finite and positive.")
    if u == 0.0:
        return np.full_like(state, analytic(c)["price"])
    if u == 1.0:
        return math.exp(-c.rate * c.maturity) * payoff(state, c)
    tau = (1.0 - u) * c.maturity
    sd = c.volatility * math.sqrt(tau)
    d1_low = (np.log(state / c.strike)
              + (c.rate - c.dividend + 0.5 * c.volatility ** 2) * tau) / sd
    d1_high = d1_low - math.log1p(c.cap / c.strike) / sd
    difference_d1, _ = _cdf_interval(d1_high, d1_low)
    difference_d2, cdf_d2_high = _cdf_interval(d1_high - sd, d1_low - sd)
    remaining_price = (state * math.exp(-c.dividend * tau) * difference_d1
                       - c.strike * math.exp(-c.rate * tau) * difference_d2
                       + c.cap * math.exp(-c.rate * tau) * cdf_d2_high)
    upper_bound = c.cap * math.exp(-c.rate * tau)
    bounded = np.clip(remaining_price, 0.0, upper_bound)
    time_zero_discount = math.exp(-c.rate * u * c.maturity)
    # This numerical fallback is recorded, not presented as exact unbiasedness.
    if diagnostics is not None:
        adjustment = np.abs(bounded - remaining_price) * time_zero_discount
        diagnostics["conditional_state_evaluations"] = diagnostics.get("conditional_state_evaluations", 0) + state.size
        diagnostics["lower_range_violations"] = diagnostics.get("lower_range_violations", 0) + int(np.count_nonzero(remaining_price < 0.0))
        diagnostics["upper_range_violations"] = diagnostics.get("upper_range_violations", 0) + int(np.count_nonzero(remaining_price > upper_bound))
        diagnostics["maximum_time_zero_clipping_adjustment"] = max(
            diagnostics.get("maximum_time_zero_clipping_adjustment", 0.0), float(np.max(adjustment, initial=0.0)))
        diagnostics["sum_absolute_time_zero_clipping_adjustments"] = diagnostics.get(
            "sum_absolute_time_zero_clipping_adjustments", 0.0) + float(np.sum(adjustment))
    return time_zero_discount * bounded


def split_observations(c: Contract, u: float, z1: np.ndarray,
                       z2: np.ndarray, diagnostics: dict | None = None) -> dict[str, np.ndarray]:
    """Couple a terminal payoff and its conditional expectation on the same path."""
    _split_check(u)
    z1, z2 = np.broadcast_arrays(np.asarray(z1, dtype=float), np.asarray(z2, dtype=float))
    if np.any(~np.isfinite(z1)) or np.any(~np.isfinite(z2)):
        raise ValueError("Gaussian inputs must be finite.")
    t = u * c.maturity
    state = c.spot * np.exp(c.log_drift * t + c.volatility * math.sqrt(t) * z1)
    terminal = c.spot * np.exp(c.log_drift * c.maturity
                              + c.volatility * math.sqrt(c.maturity)
                              * (math.sqrt(u) * z1 + math.sqrt(1.0 - u) * z2))
    discount = math.exp(-c.rate * c.maturity)
    mean_stock = c.spot * math.exp(-c.dividend * c.maturity)
    return {
        "raw": discount * payoff(terminal, c),
        "conditional": conditional_value(c, u, state, diagnostics),
        "stock_control": discount * terminal - mean_stock,
        "conditional_stock_control": (math.exp(-c.rate * t - c.dividend
                                                   * (c.maturity - t)) * state - mean_stock),
    }


def _unclipped_fit(values: np.ndarray, controls: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    controls = np.asarray(controls, dtype=float)
    if controls.ndim == 1:
        controls = controls[:, None]
    if (values.ndim != 1 or controls.ndim != 2 or len(values) != len(controls)
            or len(values) < 2
            or not np.all(np.isfinite(values)) or not np.all(np.isfinite(controls))):
        raise ValueError("A finite pilot of at least two paired observations is required.")
    centred = controls - controls.mean(axis=0)
    coefficient = np.linalg.lstsq(centred, values - values.mean(), rcond=None)[0]
    if not np.all(np.isfinite(coefficient)):
        raise ArithmeticError("The numerical pilot fit returned a non-finite coefficient.")
    return coefficient


def fit_controls(values: np.ndarray, controls: np.ndarray,
                 bound: float = 1e6) -> np.ndarray:
    """Centred Moore-Penrose least squares followed by Euclidean-norm clipping.

    Fit on independent pilot observations. Duplicate/zero controls are allowed.
    Clipping keeps the fitted coefficient norm bounded; it is not a claim that
    the resulting random coefficients are population-optimal. The scalar case
    is clipping to [−bound,bound].
    """
    if not math.isfinite(bound) or bound <= 0.0:
        raise ValueError("The coefficient bound must be finite and positive.")
    coefficient = _unclipped_fit(values, controls)
    norm = float(np.linalg.norm(coefficient))
    return coefficient if norm <= bound else coefficient * (bound / norm)


@dataclass(frozen=True)
class PilotCoefficients:
    raw: float
    conditional: float
    size: int
    bound: float = 1e6
    raw_unclipped: float | None = None
    conditional_unclipped: float | None = None
    raw_clipped: bool = False
    conditional_clipped: bool = False


def fit_pilot(c: Contract, u: float, size: int, seed: int,
              bound: float = 1e6, diagnostics: dict | None = None) -> PilotCoefficients:
    if not isinstance(size, int) or size < 2:
        raise ValueError("The pilot requires at least two observations.")
    if not math.isfinite(bound) or bound <= 0.0:
        raise ValueError("The coefficient bound must be finite and positive.")
    rng = np.random.default_rng(seed)
    obs = split_observations(c, u, rng.standard_normal(size), rng.standard_normal(size), diagnostics)
    raw = float(_unclipped_fit(obs["raw"], obs["stock_control"])[0])
    conditional = float(_unclipped_fit(obs["conditional"], obs["conditional_stock_control"])[0])
    return PilotCoefficients(float(np.clip(raw, -bound, bound)),
                             float(np.clip(conditional, -bound, bound)),
                             size, bound, raw, conditional,
                             abs(raw) > bound, abs(conditional) > bound)


def estimator_values(observations: dict[str, np.ndarray],
                     coefficients: PilotCoefficients) -> dict[str, np.ndarray]:
    y, m = observations["raw"], observations["conditional"]
    control, projected = observations["stock_control"], observations["conditional_stock_control"]
    return {
        "raw_payoff": y,
        "exact_conditional": m,
        "raw_pilot_stock_control": y - coefficients.raw * control,
        "conditional_raw_pilot_coefficient": m - coefficients.raw * projected,
        "conditional_own_pilot_coefficient": m - coefficients.conditional * projected,
        "conditional_unconditioned_stock_control": m - coefficients.raw * control,
    }


def truncated_stock_moment(c: Contract, power: int, lower: float,
                           upper: float = math.inf) -> float:
    """E[S_T**power 1(lower<S_T<upper)] under the specified risk-neutral GBM."""
    if power not in (0, 1, 2) or lower < 0.0 or upper <= lower:
        raise ValueError("Use powers 0, 1, 2 and a non-empty non-negative stock interval.")
    sd = c.volatility * math.sqrt(c.maturity)
    mu = math.log(c.spot) + c.log_drift * c.maturity
    shifted = mu + power * sd * sd
    left = -math.inf if lower == 0.0 else (math.log(lower) - shifted) / sd
    right = math.inf if math.isinf(upper) else (math.log(upper) - shifted) / sd
    # Survival-function differences avoid cancellation in the upper tail.
    probability = (normal_cdf(-left) - normal_cdf(-right) if left > 0.0
                   else normal_cdf(right) - normal_cdf(left))
    return math.exp(power * mu + 0.5 * power * power * sd * sd) * probability


def raw_population_moments(c: Contract) -> dict[str, float]:
    low, high = c.strike, c.strike + c.cap
    p0 = truncated_stock_moment(c, 0, low, high)
    p1 = truncated_stock_moment(c, 1, low, high)
    p2 = truncated_stock_moment(c, 2, low, high)
    tail = truncated_stock_moment(c, 0, high)
    tail_stock = truncated_stock_moment(c, 1, high)
    discount = math.exp(-c.rate * c.maturity)
    mean = analytic(c)["price"]
    mean_stock = c.spot * math.exp(-c.dividend * c.maturity)
    second = discount ** 2 * (p2 - 2.0 * low * p1 + low ** 2 * p0 + c.cap ** 2 * tail)
    variance = max(0.0, second - mean ** 2)
    covariance = discount ** 2 * (p2 - low * p1 + c.cap * tail_stock) - mean * mean_stock
    control_variance = mean_stock ** 2 * math.expm1(c.volatility ** 2 * c.maturity)
    beta = covariance / control_variance
    return {"mean": mean, "variance": variance,
            "stock_control_variance": control_variance,
            "stock_covariance": covariance, "oracle_beta": beta,
            "oracle_control_variance": max(0.0, variance - covariance ** 2 / control_variance)}


def gaussian_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Hermite nodes/weights for N(0,1), suitable for moderate order."""
    if not isinstance(order, int) or not 2 <= order <= 320:
        raise ValueError("Use a Gauss-Hermite order between 2 and 320.")
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    return math.sqrt(2.0) * nodes, weights / math.sqrt(math.pi)


def conditional_population_moments(c: Contract, u: float,
                                   order: int = 320, diagnostics: dict | None = None) -> dict[str, float]:
    """Smooth 1D quadrature diagnostics, not an outward-rounded certificate."""
    _split_check(u)
    raw = raw_population_moments(c)
    mean = raw["mean"]
    if u == 0.0:
        variance = covariance = control_variance = beta = optimal = 0.0
        quadrature_mean = mean
        quadrature_control_mean = 0.0
    elif u == 1.0:
        variance, covariance = raw["variance"], raw["stock_covariance"]
        control_variance, beta = raw["stock_control_variance"], raw["oracle_beta"]
        optimal, quadrature_mean = raw["oracle_control_variance"], mean
        quadrature_control_mean = 0.0
    else:
        nodes, weights = gaussian_rule(order)
        t = u * c.maturity
        state = c.spot * np.exp(c.log_drift * t + c.volatility * math.sqrt(t) * nodes)
        m = conditional_value(c, u, state, diagnostics)
        control = (math.exp(-c.rate * t - c.dividend * (c.maturity - t)) * state
                   - c.spot * math.exp(-c.dividend * c.maturity))
        # Centre using the exact analytic price/control means, not the rule mean.
        variance = float(weights @ ((m - mean) ** 2))
        covariance = float(weights @ ((m - mean) * control))
        control_variance = (c.spot * math.exp(-c.dividend * c.maturity)) ** 2 * math.expm1(
            c.volatility ** 2 * t)
        beta = covariance / control_variance
        optimal = max(0.0, variance - covariance ** 2 / control_variance)
        quadrature_mean, quadrature_control_mean = float(weights @ m), float(weights @ control)
    same_beta_variance = (variance - 2.0 * raw["oracle_beta"] * covariance
                          + raw["oracle_beta"] ** 2 * control_variance)
    naive_added = raw["oracle_beta"] ** 2 * (raw["stock_control_variance"] - control_variance)
    return {"mean": mean, "quadrature_mean": quadrature_mean,
            "quadrature_control_mean": quadrature_control_mean,
            "variance": variance, "stock_control_variance": control_variance,
            "stock_covariance": covariance, "oracle_beta": beta,
            "oracle_control_variance": optimal,
            "raw_oracle_coefficient_conditional_variance": same_beta_variance,
            "unconditioned_control_variance": same_beta_variance + naive_added,
            "unconditioned_control_added_variance": naive_added,
            "projection_removed_variance": raw["variance"] - variance,
            "same_beta_projection_removed_variance": raw["oracle_control_variance"] - same_beta_variance,
            "coefficient_reoptimisation_reduction": (raw["oracle_beta"] - beta) ** 2 * control_variance}


def residual_population_second_moment(c: Contract, u: float, beta: float,
                                      order: int = 320) -> float:
    """Integrate exact residual-life lognormal conditional variances in 1D.

    This is a separate numerical check of the projection variance identity.
    Inner expectations are analytic truncated moments; outer quadrature is not
    a rigorous certificate and cancellation remains floating-point arithmetic.
    """
    _split_check(u)
    if not math.isfinite(beta):
        raise ValueError("The control coefficient must be finite.")
    if u == 1.0:
        return 0.0
    if u == 0.0:
        raw = raw_population_moments(c)
        return (raw["variance"] - 2.0 * beta * raw["stock_covariance"]
                + beta ** 2 * raw["stock_control_variance"])
    nodes, weights = gaussian_rule(order)
    t = u * c.maturity
    states = c.spot * np.exp(c.log_drift * t + c.volatility * math.sqrt(t) * nodes)
    conditional_variances = []
    for state in states:
        local = raw_population_moments(replace(c, spot=float(state), maturity=c.maturity - t))
        conditional_variances.append(local["variance"] - 2.0 * beta * local["stock_covariance"]
                                     + beta ** 2 * local["stock_control_variance"])
    return math.exp(-2.0 * c.rate * t) * float(weights @ np.asarray(conditional_variances))


def cost_counts(method: str, production_size: int, pilot_size: int) -> dict[str, int]:
    """Interior-split stand-alone kernel counts; shared study costs separate.

    Raw controls need one normal/terminal stock per draw. Conditional controls
    need one normal/intermediate stock, but four CDF calls per draw. Naively
    unconditioned controls additionally simulate a second normal/terminal stock.
    Independent conditional and raw pilots are charged to their own methods.
    """
    if method not in METHODS or min(production_size, pilot_size) < 0:
        raise ValueError("Unknown estimator or negative count.")
    conditional = method.startswith("conditional") or method == "exact_conditional"
    pilot = pilot_size if "pilot" in method or method == "conditional_unconditioned_stock_control" else 0
    raw_pilot = pilot if "raw_pilot" in method or method == "conditional_unconditioned_stock_control" else 0
    conditional_pilot = pilot - raw_pilot
    naive = method == "conditional_unconditioned_stock_control"
    return {
        "production_draws": production_size, "independent_pilot_draws": pilot,
        "pilot_least_squares_fits": int(pilot > 0),
        "normal_variates": production_size * (2 if naive else 1) + pilot,
        "terminal_payoff_evaluations": (0 if conditional else production_size) + raw_pilot,
        "terminal_stock_exponentials": (production_size if not conditional or naive else 0) + raw_pilot,
        "intermediate_stock_exponentials": (production_size if conditional else 0) + conditional_pilot,
        "conditional_call_spread_evaluations": (production_size if conditional else 0) + conditional_pilot,
        "normal_cdf_evaluations": 4 * ((production_size if conditional else 0) + conditional_pilot),
    }

