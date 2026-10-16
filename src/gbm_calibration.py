"""Synthetic Gaussian-return inference for a classical GBM pricing model.

The probability statements use independent exactly Gaussian returns and exact
arithmetic. SciPy quantiles and floating-point evaluation are numerical
implementations, not outward-rounded certificates. No market data are used.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.stats import beta, chi2, t

from capped_gbm import Contract, analytic


@dataclass(frozen=True)
class GaussianEstimate:
    observations: int
    dimensions: int
    interval: float
    log_drift: np.ndarray
    scatter: np.ndarray
    covariance_ml_candidate: np.ndarray
    covariance_unbiased: np.ndarray
    scatter_rank: int
    covariance_ml_status: str


@dataclass(frozen=True)
class ParameterIntervals:
    log_drift_lower: np.ndarray
    log_drift_upper: np.ndarray
    variance_lower: np.ndarray
    variance_upper: np.ndarray
    volatility_lower: np.ndarray
    volatility_upper: np.ndarray
    marginal_alpha: float
    familywise: bool


def _returns_array(returns: np.ndarray, interval: float) -> np.ndarray:
    data = np.asarray(returns, dtype=float)
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] < 1:
        raise ValueError("At least two return observations and one dimension are required.")
    if not math.isfinite(interval) or interval <= 0.0 or not np.isfinite(data).all():
        raise ValueError("The observation interval must be positive and all returns finite.")
    return data


def estimate_gaussian_returns(returns: np.ndarray, interval: float) -> GaussianEstimate:
    """Estimate annual log drift and covariance, flagging singular ML candidates.

    A positive-definite unrestricted covariance MLE exists only for full-rank
    scatter, which requires n>d. Reported rank uses NumPy's floating-point
    tolerance and is a numerical diagnostic, not an exact rank certificate.
    Marginal inference remains available when sample covariance is
    singular but each marginal scatter is positive; that is not a full-rank ML
    estimate of the joint law.
    """
    data = _returns_array(returns, interval)
    n, d = data.shape
    centre = data.mean(axis=0)
    residual = data - centre
    scatter = residual.T @ residual
    rank = int(np.linalg.matrix_rank(scatter))
    status = ("full_rank_covariance_mle" if n > d and rank == d
              else "singular_scatter_no_positive_definite_covariance_mle")
    return GaussianEstimate(n, d, interval, centre / interval, scatter,
                            scatter / (n * interval), scatter / ((n - 1) * interval),
                            rank, status)


def negative_log_likelihood(returns: np.ndarray, interval: float,
                            log_drift: np.ndarray, covariance: np.ndarray) -> float:
    """Negative Gaussian log likelihood, including normalising constants."""
    data = _returns_array(returns, interval)
    n, d = data.shape
    mean = np.asarray(log_drift, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    if mean.shape != (d,) or cov.shape != (d, d) or not np.isfinite(mean).all():
        raise ValueError("Parameter dimensions must match returns.")
    if not np.isfinite(cov).all() or not np.allclose(cov, cov.T, rtol=1e-12, atol=1e-15):
        raise ValueError("The covariance must be finite and symmetric.")
    try:
        chol = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError as error:
        raise ValueError("The likelihood requires a positive-definite covariance.") from error
    residual = data - interval * mean
    whitened = np.linalg.solve(chol, residual.T)
    log_det = 2.0 * float(np.log(np.diag(chol)).sum())
    return 0.5 * (n * (d * math.log(2.0 * math.pi * interval) + log_det)
                  + float(np.sum(whitened ** 2)) / interval)


def trace_target_shrinkage(estimate: GaussianEstimate, strength: float = 0.1) -> np.ndarray:
    """Return fixed-strength sample-trace regularisation, never labelled ML.

    Positive trace and strength in (0,1] make this matrix positive definite even
    if scatter is singular. The target depends on observations; only strength
    is fixed. Zero trace cannot be remedied by this target. Exact intervals
    cannot be transferred to this shrunken covariance.
    """
    if not math.isfinite(strength) or not 0.0 < strength <= 1.0:
        raise ValueError("Shrinkage strength must lie in (0,1].")
    cov = estimate.covariance_ml_candidate
    target = float(np.trace(cov)) / estimate.dimensions
    if target <= 0.0 or not math.isfinite(target):
        raise ValueError("Trace-target shrinkage requires strictly positive scatter trace.")
    return (1.0 - strength) * cov + strength * target * np.eye(estimate.dimensions)


def gaussian_parameter_intervals(estimate: GaussianEstimate, alpha: float = 0.05,
                                 familywise: bool = False) -> ParameterIntervals:
    """Marginal exact t/chi-square intervals, or Bonferroni for 2d quantities.

    Intervals always use the original unbiased covariance diagonal, never the
    shrunken matrix. The marginal observation law must be non-degenerate normal.
    """
    if not math.isfinite(alpha) or not 0.0 < alpha < 1.0:
        raise ValueError("Alpha must lie in (0,1).")
    variance = np.diag(estimate.covariance_unbiased)
    if np.any(variance <= 0.0):
        raise ValueError("Normal marginal confidence intervals require positive marginal scatter.")
    nu = estimate.observations - 1
    marginal_alpha = alpha / (2 * estimate.dimensions) if familywise else alpha
    radius = float(t.ppf(1.0 - marginal_alpha / 2.0, nu)) * np.sqrt(
        variance / (estimate.observations * estimate.interval))
    lo = nu * variance / float(chi2.ppf(1.0 - marginal_alpha / 2.0, nu))
    hi = nu * variance / float(chi2.ppf(marginal_alpha / 2.0, nu))
    return ParameterIntervals(estimate.log_drift - radius, estimate.log_drift + radius,
                              lo, hi, np.sqrt(lo), np.sqrt(hi), marginal_alpha, familywise)


def capped_volatility_peak(c: Contract) -> float | None:
    """Unique positive interior stationary volatility, when it exists.

    The log expression is evaluated as a sum to avoid forming large products.
    At this volatility the capped call has a strict global volatility maximum.
    """
    square = ((math.log(c.strike / c.spot) + math.log((c.strike + c.cap) / c.spot))
              / c.maturity - 2.0 * (c.rate - c.dividend))
    return math.sqrt(square) if square > 0.0 else None


def capped_price_interval(c: Contract, volatility_lower: float,
                          volatility_upper: float) -> dict[str, float | bool | None]:
    """Exact mathematical image of a positive volatility interval, numerically evaluated.

    The lower price is the smaller endpoint value; the upper also evaluates the
    interior volatility maximum if present. Under fixed S,K,H,T,r,q, coverage
    of the true volatility implies coverage of its true risk-neutral price.
    """
    if (not math.isfinite(volatility_lower) or not math.isfinite(volatility_upper)
            or not 0.0 < volatility_lower <= volatility_upper):
        raise ValueError("Volatility endpoints must be finite, positive and ordered.")
    first = analytic(replace(c, volatility=volatility_lower))["price"]
    second = analytic(replace(c, volatility=volatility_upper))["price"]
    peak = capped_volatility_peak(c)
    interior = peak is not None and volatility_lower < peak < volatility_upper
    upper = max(first, second)
    if interior:
        upper = max(upper, analytic(replace(c, volatility=peak))["price"])
    return {"lower": min(first, second), "upper": upper,
            "endpoint_only_upper": max(first, second),
            "interior_peak_included": interior, "volatility_peak": peak}


def binomial_coverage_interval(successes: int, repetitions: int,
                               alpha: float = 0.05) -> tuple[float, float]:
    """Clopper-Pearson coverage diagnostic for repeated independent datasets."""
    if (not isinstance(successes, int) or not isinstance(repetitions, int)
            or repetitions < 1 or not 0 <= successes <= repetitions
            or not 0.0 < alpha < 1.0):
        raise ValueError("Valid integer binomial counts and alpha are required.")
    lo = 0.0 if successes == 0 else float(beta.ppf(alpha / 2.0, successes,
                                                 repetitions - successes + 1))
    hi = 1.0 if successes == repetitions else float(beta.ppf(1.0 - alpha / 2.0,
                                                            successes + 1,
                                                            repetitions - successes))
    return lo, hi
