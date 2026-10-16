"""Synthetic capped-GBM pricing companion, using only NumPy and the standard library.

The finite-difference solution is backward Euler on a fixed uniform heat grid.
Its constant tridiagonal matrix is diagonalised using a DST-I; the FFT is only
an implementation of the tridiagonal recurrence, not a spectral space scheme.
Exact-arithmetic stability does not certify floating-point arithmetic errors.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Iterable

import numpy as np


GREEKS = ("price", "delta", "gamma", "vega", "rho", "theta")


@dataclass(frozen=True)
class Contract:
    spot: float = 100.0
    strike: float = 100.0
    cap: float = 20.0
    maturity: float = 1.0
    rate: float = 0.03
    dividend: float = 0.01
    volatility: float = 0.2

    def __post_init__(self) -> None:
        vals = asdict(self)
        if not all(math.isfinite(x) for x in vals.values()):
            raise ValueError("All contract parameters must be finite.")
        if min(self.spot, self.strike, self.cap, self.maturity,
               self.volatility) <= 0.0:
            raise ValueError("spot, strike, cap, maturity and volatility must be positive.")

    @property
    def log_drift(self) -> float:
        return self.rate - self.dividend - 0.5 * self.volatility ** 2


def normal_cdf(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def normal_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def call_analytic(c: Contract, strike: float) -> dict[str, float]:
    """Black-Scholes call price and Greeks; theta is calendar-time theta."""
    s, t, r, q, sig = c.spot, c.maturity, c.rate, c.dividend, c.volatility
    rt = math.sqrt(t)
    d1 = (math.log(s / strike) + (r - q + 0.5 * sig * sig) * t) / (sig * rt)
    d2 = d1 - sig * rt
    dq, dr = math.exp(-q * t), math.exp(-r * t)
    n1, n2, density = normal_cdf(d1), normal_cdf(d2), normal_pdf(d1)
    return {
        "price": s * dq * n1 - strike * dr * n2,
        "delta": dq * n1,
        "gamma": dq * density / (s * sig * rt),
        "vega": s * dq * density * rt,
        "rho": strike * t * dr * n2,
        "theta": (-s * dq * density * sig / (2.0 * rt)
                  + q * s * dq * n1 - r * strike * dr * n2),
    }


def analytic(c: Contract) -> dict[str, float]:
    """Capped call = call(K) - call(K+H); rates/volatility use unit bumps."""
    low, high = call_analytic(c, c.strike), call_analytic(c, c.strike + c.cap)
    return {key: low[key] - high[key] for key in GREEKS}


def payoff(terminal: np.ndarray | float, c: Contract) -> np.ndarray:
    return np.minimum(np.maximum(np.asarray(terminal) - c.strike, 0.0), c.cap)


@dataclass(frozen=True)
class HeatGrid:
    lower: float
    upper: float
    intervals: int
    time_steps: int

    def __post_init__(self) -> None:
        if not math.isfinite(self.lower) or not math.isfinite(self.upper):
            raise ValueError("Heat boundaries must be finite.")
        if self.upper <= self.lower:
            raise ValueError("Heat upper bound must exceed its lower bound.")
        if not isinstance(self.intervals, int) or self.intervals < 2:
            raise ValueError("At least two integer space intervals are required.")
        if not isinstance(self.time_steps, int) or self.time_steps < 1:
            raise ValueError("At least one integer time step is required.")

    @property
    def spacing(self) -> float:
        return (self.upper - self.lower) / self.intervals

    def validate_for(self, c: Contract) -> None:
        if not self.lower < 0.0 < math.log1p(c.cap / c.strike) < self.upper:
            raise ValueError("The domain must contain both payoff kinks strictly.")
        y = math.log(c.spot / c.strike) + c.log_drift * c.maturity
        if not self.lower <= y <= self.upper:
            raise ValueError("The evaluation coordinate lies outside the heat domain.")


def _dst_one(values: np.ndarray) -> np.ndarray:
    """Unnormalised DST-I; applying it twice gives 2*(n+1)*values."""
    n = len(values)
    extended = np.zeros(2 * (n + 1), dtype=float)
    extended[1:n + 1] = values
    extended[n + 2:] = -values[::-1]
    return -np.fft.fft(extended).imag[1:n + 1]


def heat_initial(c: Contract, grid: HeatGrid) -> tuple[np.ndarray, np.ndarray]:
    grid.validate_for(c)
    nodes = np.linspace(grid.lower, grid.upper, grid.intervals + 1)
    # The capped expression avoids evaluating exp(y) for very large upper bounds.
    f = np.zeros_like(nodes)
    capped = nodes >= math.log1p(c.cap / c.strike)
    uncapped = (nodes > 0.0) & ~capped
    f[capped] = c.cap
    f[uncapped] = c.strike * np.expm1(nodes[uncapped])
    f[0], f[-1] = 0.0, c.cap
    return nodes, f


def heat_backward_euler(c: Contract, grid: HeatGrid) -> tuple[np.ndarray, np.ndarray]:
    """Apply exactly N backward-Euler steps algebraically, up to float round-off."""
    nodes, initial = heat_initial(c, grid)
    lift = c.cap * (nodes - grid.lower) / (grid.upper - grid.lower)
    z = initial[1:-1] - lift[1:-1]
    lam = 0.5 * c.volatility ** 2 * (c.maturity / grid.time_steps) / grid.spacing ** 2
    modes = np.arange(1, grid.intervals, dtype=float)
    eig = 4.0 * np.sin(math.pi * modes / (2.0 * grid.intervals)) ** 2
    amplification = np.exp(-grid.time_steps * np.log1p(lam * eig))
    coefficients = _dst_one(z)
    evolved = _dst_one(coefficients * amplification) / (2.0 * grid.intervals)
    values = lift.copy()
    values[1:-1] += evolved
    values[0], values[-1] = 0.0, c.cap
    return nodes, values


def heat_backward_euler_thomas(c: Contract, grid: HeatGrid) -> tuple[np.ndarray, np.ndarray]:
    """Direct tridiagonal recurrence used as an independent implementation check."""
    nodes, values = heat_initial(c, grid)
    n = grid.intervals - 1
    lam = 0.5 * c.volatility ** 2 * (c.maturity / grid.time_steps) / grid.spacing ** 2
    cp = np.zeros(n)
    inv = np.zeros(n)
    inv[0] = 1.0 / (1.0 + 2.0 * lam)
    if n > 1:
        cp[0] = -lam * inv[0]
    for j in range(1, n):
        inv[j] = 1.0 / (1.0 + 2.0 * lam + lam * cp[j - 1])
        if j < n - 1:
            cp[j] = -lam * inv[j]
    for _ in range(grid.time_steps):
        rhs = values[1:-1].copy()
        rhs[-1] += lam * c.cap
        rhs[0] *= inv[0]
        for j in range(1, n):
            rhs[j] = (rhs[j] + lam * rhs[j - 1]) * inv[j]
        for j in range(n - 2, -1, -1):
            rhs[j] -= cp[j] * rhs[j + 1]
        values[1:-1] = rhs
    return nodes, values


def pde_price(c: Contract, grid: HeatGrid) -> float:
    """Monotone linear interpolation in heat coordinates followed by discounting."""
    nodes, values = heat_backward_euler(c, grid)
    y = math.log(c.spot / c.strike) + c.log_drift * c.maturity
    return math.exp(-c.rate * c.maturity) * float(np.interp(y, nodes, values))


def domain_price_bound(c: Contract, grid: HeatGrid) -> float:
    """Bound at both interpolation nodes, replacing whole-line heat by 0/H.

    This bounds the domain component of the interpolated numerical price.
    Interpolation error against the continuous whole-line solution is a
    separate discretisation component. Constants are not outward rounded.
    """
    grid.validate_for(c)
    y = math.log(c.spot / c.strike) + c.log_drift * c.maturity
    sd = c.volatility * math.sqrt(c.maturity)
    left = min(grid.intervals - 1, max(0, math.floor((y - grid.lower) / grid.spacing)))
    bracket_nodes = (grid.lower + left * grid.spacing,
                     grid.lower + (left + 1) * grid.spacing)
    exit_bound = max(min(1.0, 2.0 * normal_cdf(-(node - grid.lower) / sd)
                         + 2.0 * normal_cdf(-(grid.upper - node) / sd))
                     for node in bracket_nodes)
    return math.exp(-c.rate * c.maturity) * c.cap * exit_bound


def default_bumps(c: Contract, grid: HeatGrid) -> dict[str, float]:
    """Shrinking bumps compatible with a conservative O(h**(1/2)) price bound.

    Spot bump d~h**(1/8) makes E/d**2 vanish; the other first-derivative bumps
    d~h**(1/4) make E/d vanish. The added spot h term resolves the interpolation
    cells on coarse grids. h is dimensionless. Explicit numerical constants
    are design choices; parameter/domain/arithmetic errors remain separate.
    """
    factor = grid.spacing ** 0.25
    return {"spot": min(0.10 * c.spot, c.spot * (0.006 * grid.spacing ** 0.125 + grid.spacing)),
            "volatility": min(0.25 * c.volatility, 0.10 * c.volatility * factor),
            "rate": 0.01 * factor,
            "maturity": min(0.25 * c.maturity, 0.05 * c.maturity * factor)}


def bumped_greeks(c: Contract, grid: HeatGrid,
                  bumps: dict[str, float] | None = None,
                  price_function=pde_price) -> dict[str, float]:
    """Central bumped prices, not derivatives of a piecewise-linear interpolant."""
    bumps = default_bumps(c, grid) if bumps is None else bumps
    for key in ("spot", "volatility", "rate", "maturity"):
        if key not in bumps or not math.isfinite(bumps[key]) or bumps[key] <= 0.0:
            raise ValueError("All four parameter bumps must be positive and finite.")
    base = price_function(c, grid)
    derivatives = {"price": base}
    for parameter, greek in (("spot", "delta"), ("volatility", "vega"),
                             ("rate", "rho"), ("maturity", "theta")):
        step = bumps[parameter]
        low = replace(c, **{parameter: getattr(c, parameter) - step})
        high = replace(c, **{parameter: getattr(c, parameter) + step})
        lp, hp = price_function(low, grid), price_function(high, grid)
        derivatives[greek] = (hp - lp) / (2.0 * step)
        if parameter == "maturity":
            derivatives[greek] *= -1.0
        if parameter == "spot":
            derivatives["gamma"] = (hp - 2.0 * base + lp) / step ** 2
    return {key: derivatives[key] for key in GREEKS}


def analytic_bumped_greeks(c: Contract, grid: HeatGrid,
                           bumps: dict[str, float] | None = None) -> dict[str, float]:
    return bumped_greeks(c, grid, bumps,
                         price_function=lambda contract, unused: analytic(contract)["price"])


def pde_greek_domain_bounds(c: Contract, grid: HeatGrid,
                            bumps: dict[str, float]) -> dict[str, float]:
    """Only domain-error amplification; no grid or round-off certificate is implied."""
    base = domain_price_bound(c, grid)
    out = {"price": base}
    for parameter, greek in (("spot", "delta"), ("volatility", "vega"),
                             ("rate", "rho"), ("maturity", "theta")):
        step = bumps[parameter]
        lo = domain_price_bound(replace(c, **{parameter: getattr(c, parameter) - step}), grid)
        hi = domain_price_bound(replace(c, **{parameter: getattr(c, parameter) + step}), grid)
        out[greek] = (lo + hi) / (2.0 * step)
        if parameter == "spot":
            out["gamma"] = (lo + 2.0 * base + hi) / step ** 2
    return out


def monte_carlo_observations(c: Contract, normals: np.ndarray) -> np.ndarray:
    """Unbiased six-metric observations plus a bounded mixed Gamma rival.

    Columns follow GREEKS, with mixed pathwise/likelihood-ratio Gamma as the
    seventh column. All outputs reuse the same exact terminal evaluation.
    """
    z = np.asarray(normals, dtype=float)
    t, sig, s = c.maturity, c.volatility, c.spot
    root_t = math.sqrt(t)
    terminal = s * np.exp(c.log_drift * t + sig * root_t * z)
    pay = payoff(terminal, c)
    active = (terminal > c.strike) & (terminal < c.strike + c.cap)
    d = math.exp(-c.rate * t)
    path = np.where(active, terminal, 0.0)
    v = sig * root_t
    gamma_weight = ((z * z - 1.0) / (s * s * v * v) - z / (s * s * v))
    mixed_gamma = d * path / (s * s) * (z / v - 1.0)
    return np.column_stack((d * pay,
                            d * path / s,
                            d * pay * gamma_weight,
                            d * path * (root_t * z - sig * t),
                            d * t * (path - pay),
                            d * (c.rate * pay - path * (c.log_drift + sig * z / (2.0 * root_t))),
                            mixed_gamma))


def mixed_gamma_absolute_bound(c: Contract) -> float:
    """Bound D*J/S**2*(max active |Z|/v+1) for the mixed Gamma observation."""
    v = c.volatility * math.sqrt(c.maturity)
    za = (math.log(c.strike / c.spot) - c.log_drift * c.maturity) / v
    zb = (math.log((c.strike + c.cap) / c.spot) - c.log_drift * c.maturity) / v
    return (math.exp(-c.rate * c.maturity) * (c.strike + c.cap) / c.spot ** 2
            * (max(abs(za), abs(zb)) / v + 1.0))


def greek_second_moment_bounds(c: Contract) -> dict[str, float]:
    """Conservative bounds on raw estimator second moments, not variances."""
    d2 = math.exp(-2.0 * c.rate * c.maturity)
    j, s, t, sig, h = c.strike + c.cap, c.spot, c.maturity, c.volatility, c.cap
    v = sig * math.sqrt(t)
    return {"price": d2 * h * h,
            "delta": d2 * j * j / (s * s),
            "gamma": d2 * h * h / s ** 4 * (2.0 / v ** 4 + 1.0 / v ** 2),
            "vega": 2.0 * d2 * j * j * (t + sig * sig * t * t),
            "rho": d2 * t * t * (j + h) ** 2,
            "theta": 3.0 * d2 * ((c.rate * h) ** 2 + (j * c.log_drift) ** 2
                                    + j * j * sig * sig / (4.0 * t))}


def raw_observation_ranges(c: Contract) -> dict[str, tuple[float, float]]:
    """Mathematical bounds; floating-point evaluation is not outward rounded."""
    s, j, h, t, sig = c.spot, c.strike + c.cap, c.cap, c.maturity, c.volatility
    d, rt = math.exp(-c.rate * t), math.sqrt(t)
    v = sig * rt
    za = (math.log(c.strike / s) - c.log_drift * t) / v
    zb = (math.log(j / s) - c.log_drift * t) / v
    maxz = max(abs(za), abs(zb))
    vmax = d * j * (rt * maxz + sig * t)
    thetamax = d * (abs(c.rate) * h + j * (abs(c.log_drift) + sig * maxz / (2.0 * rt)))
    return {"price": (0.0, d * h), "delta": (0.0, d * j / s),
            "vega": (-vmax, vmax), "rho": (-t * d * h, t * d * c.strike),
            "theta": (-thetamax, thetamax)}


def _mc_summary(observations: np.ndarray, c: Contract, alpha: float,
                method: str, counts: dict, coefficients=None,
                mixed_gamma: bool = False) -> list[dict]:
    n = len(observations)
    if n < 2:
        raise ValueError("At least two independent observations are required.")
    means = observations.mean(axis=0)
    ses = observations.std(axis=0, ddof=1) / math.sqrt(n)
    ref = analytic(c)
    moment = greek_second_moment_bounds(c)
    ranges = raw_observation_ranges(c)
    if not 0.0 < alpha < 1.0:
        raise ValueError("The finite-sample error probability must be in (0,1).")
    rows = []
    metrics = [(6, "gamma")] if mixed_gamma else list(enumerate(GREEKS))
    for i, key in metrics:
        row = {"method": method, "quantity": key, "estimate": float(means[i]),
               "gamma_estimator": ("mixed_pathwise_likelihood_ratio" if mixed_gamma else "payoff_likelihood_ratio") if key == "gamma" else "not_applicable",
               "analytic": ref[key], "error": float(means[i] - ref[key]),
               "sample_standard_error": float(ses[i]),
               "asymptotic_normal_low": float(means[i] - 1.959963984540054 * ses[i]),
               "asymptotic_normal_high": float(means[i] + 1.959963984540054 * ses[i]),
               "interval_status": "pointwise 95 per cent asymptotic normal diagnostic, not certified finite-sample coverage",
               "independent_observations": n, **counts}
        if method.startswith("raw_exact_terminal"):
            gamma_bound = mixed_gamma_absolute_bound(c) if mixed_gamma else None
            second_moment_bound = gamma_bound ** 2 if mixed_gamma else moment[key]
            row["second_moment_upper_bound"] = second_moment_bound
            row["chebyshev_half_width"] = math.sqrt(second_moment_bound / (n * alpha))
            row["finite_bound_error_probability"] = alpha
            row["chebyshev_status"] = "finite-sample bound in ideal exact sampling/arithmetic; constants not outward rounded"
            if key in ranges or mixed_gamma:
                low, high = (-gamma_bound, gamma_bound) if mixed_gamma else ranges[key]
                row["observation_lower_bound"] = low
                row["observation_upper_bound"] = high
                row["hoeffding_half_width"] = ((high - low) * math.sqrt(math.log(2.0 / alpha) / (2.0 * n)))
                row["finite_bound_half_width"] = row["hoeffding_half_width"]
                row["finite_bound_type"] = "Hoeffding"
            else:
                row["finite_bound_half_width"] = row["chebyshev_half_width"]
                row["finite_bound_type"] = "Chebyshev"
        if coefficients is not None:
            row["independent_pilot_control_coefficient"] = float(coefficients[i])
        rows.append(row)
    return rows


def monte_carlo(c: Contract, observations: int = 131072, pilot_pairs: int = 16384,
                seed: int = 732451, alpha: float = 0.05) -> list[dict]:
    """Raw MC and independent-pilot antithetic terminal-stock controls.

    Raw uses 2*N independent normal draws. Antithetic uses N independent
    pairs (Z,-Z). The controlled method allocates N-pilot_pairs production
    pairs plus a separately generated pilot, matching the same total 2*N
    terminal evaluation budget. Price itself is analytically known and is the strongest
    baseline. The centred stock control has known mean S*exp(-q*T).
    """
    if observations < 2 or pilot_pairs < 2 or observations - pilot_pairs < 2:
        raise ValueError("At least two production observations and pilot pairs are required.")
    children = np.random.SeedSequence(seed).spawn(3)
    raw_rng, pilot_rng, pair_rng = [np.random.default_rng(child) for child in children]
    raw = monte_carlo_observations(c, raw_rng.standard_normal(2 * observations))
    raw_counts = {"production_normals": 2 * observations,
                  "production_terminal_evaluations": 2 * observations,
                  "pilot_normals": 0, "pilot_terminal_evaluations": 0,
                  "total_terminal_evaluations": 2 * observations}
    raw_rows = _mc_summary(raw, c, alpha, "raw_exact_terminal", raw_counts)
    raw_rows += _mc_summary(raw, c, alpha, "raw_exact_terminal_mixed_gamma", raw_counts,
                            mixed_gamma=True)

    def paired(z):
        vals = 0.5 * (monte_carlo_observations(c, z) + monte_carlo_observations(c, -z))
        dt = math.exp(-c.rate * c.maturity)
        st = c.spot * math.exp(c.log_drift * c.maturity) * np.cosh(c.volatility * math.sqrt(c.maturity) * z)
        stock = dt * st - c.spot * math.exp(-c.dividend * c.maturity)
        return vals, stock

    pilot, pilot_stock = paired(pilot_rng.standard_normal(pilot_pairs))
    centred = pilot_stock - pilot_stock.mean()
    denom = float(centred @ centred)
    if denom <= np.finfo(float).tiny:
        coefficients = np.zeros(len(GREEKS) + 1)
    else:
        coefficients = centred @ (pilot - pilot.mean(axis=0)) / denom
    # A declared deterministic cap makes coefficients bounded even on an
    # unusually ill-conditioned pilot. Independence still gives unbiasedness.
    coefficients = np.clip(coefficients, -1e6, 1e6)
    z = pair_rng.standard_normal(observations)
    pairs, pair_stock = paired(z)
    controlled_count = observations - pilot_pairs
    # Re-evaluate only the declared controlled production budget. The shared
    # normal prefix deliberately couples the two comparisons but does not
    # remove terminal-evaluation charges for either estimator.
    controlled_pairs, controlled_stock = paired(z[:controlled_count])
    control = controlled_pairs - controlled_stock[:, None] * coefficients[None, :]
    common_counts = {"production_normals": controlled_count,
                     "production_terminal_evaluations": 2 * controlled_count,
                     "pilot_normals": pilot_pairs,
                     "pilot_terminal_evaluations": 2 * pilot_pairs,
                     "total_terminal_evaluations": 2 * observations}
    pair_counts = {**common_counts, "production_normals": observations,
                   "production_terminal_evaluations": 2 * observations,
                   "pilot_normals": 0, "pilot_terminal_evaluations": 0}
    pair_rows = _mc_summary(pairs, c, alpha, "antithetic_exact_terminal", pair_counts)
    pair_rows += _mc_summary(pairs, c, alpha, "antithetic_exact_terminal_mixed_gamma",
                             pair_counts, mixed_gamma=True)
    control_rows = _mc_summary(control, c, alpha, "antithetic_stock_control",
                               common_counts, coefficients)
    control_rows += _mc_summary(control, c, alpha, "antithetic_stock_control_mixed_gamma",
                                common_counts, coefficients, mixed_gamma=True)
    raw_ses = {(row["quantity"], row["gamma_estimator"]): row["sample_standard_error"] for row in raw_rows}
    raw_gamma_se = raw_ses[("gamma", "payoff_likelihood_ratio")]
    for row in pair_rows + control_rows:
        se = row["sample_standard_error"]
        row["estimated_variance_reduction_equal_total_terminal_evaluations"] = (raw_ses[(row["quantity"], row["gamma_estimator"])] / se) ** 2 if se else None
    for row in raw_rows + pair_rows + control_rows:
        if row["quantity"] == "gamma":
            se = row["sample_standard_error"]
            row["estimated_variance_reduction_against_raw_payoff_lr_gamma"] = (raw_gamma_se / se) ** 2 if se else None
    return raw_rows + pair_rows + control_rows


def synthetic_contracts() -> Iterable[tuple[str, Contract]]:
    """Disclosure-safe designed test cases, not estimated market parameters."""
    yield "lower_kink", Contract()
    yield "upper_kink", Contract(spot=120.0)
    yield "inside_spread", Contract(spot=110.0)
    yield "deep_cap", Contract(spot=160.0)
    yield "low_vol_short_dated", Contract(spot=100.0, maturity=0.05, volatility=0.06)
    yield "negative_rate_dividend", Contract(spot=105.0, rate=-0.015, dividend=0.045, volatility=0.25)
    yield "high_volatility", Contract(spot=95.0, volatility=0.65, maturity=2.0, cap=30.0)


def study_domain(c: Contract, sigmas: float = 9.0) -> tuple[float, float]:
    """Fixed per-contract domain, reused unchanged for all Greek bumps."""
    y = math.log(c.spot / c.strike) + c.log_drift * c.maturity
    kink = math.log1p(c.cap / c.strike)
    margin = max(0.30, sigmas * c.volatility * math.sqrt(c.maturity))
    return min(0.0, y) - margin, max(kink, y) + margin
