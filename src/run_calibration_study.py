"""Reproduce synthetic Gaussian calibration and capped-price uncertainty checks."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy

import capped_gbm
from capped_gbm import Contract, analytic, synthetic_contracts
from gbm_calibration import (binomial_coverage_interval, capped_price_interval,
                             capped_volatility_peak, estimate_gaussian_returns,
                             gaussian_parameter_intervals, trace_target_shrinkage)


CORRELATION = np.array([[1.0, 0.4, -0.2], [0.4, 1.0, 0.25], [-0.2, 0.25, 1.0]])
PHYSICAL_TOTAL_DRIFT = np.array([0.06, 0.04, 0.08])
FAMILY_ALPHA = 0.05
STUDY_RESULT_FILENAMES = ("parameter_estimation.csv", "coverage.csv", "price_uncertainty.csv",
                          "frequency_horizon.csv", "frequency_horizon_coverage.csv",
                          "counterexamples.json")


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def case_parameters(c: Contract) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vols = c.volatility * np.array([1.0, 0.8, 1.2])
    covariance = np.outer(vols, vols) * CORRELATION
    dividend = np.array([c.dividend, 0.015, 0.02])
    log_drift = PHYSICAL_TOTAL_DRIFT - dividend - 0.5 * np.diag(covariance)
    return vols, covariance, dividend, log_drift


def simulate_case(c: Contract, observations: int, interval: float,
                  repetitions: int, seed: int) -> tuple[list[dict], list[dict], dict]:
    """Aggregate independent datasets; no return arrays are retained in output."""
    if not isinstance(observations, int) or observations < 4:
        raise ValueError("At least four observations are required for finite Student-pivot variance diagnostics.")
    vols, covariance, dividend, log_drift = case_parameters(c)
    rng = np.random.default_rng(seed)
    chol = np.linalg.cholesky(covariance)
    n, d, nu = observations, 3, observations - 1
    mean_estimates, variance_estimates, mu_estimates = [], [], []
    student_pivots, variance_pivots = [], []
    mean_cover = np.zeros(d, dtype=int)
    variance_cover = np.zeros(d, dtype=int)
    family_mean_cover = np.zeros(d, dtype=int)
    family_variance_cover = np.zeros(d, dtype=int)
    family_mu_cover = np.zeros(d, dtype=int)
    family_joint = 0
    proper_price_cover = naive_price_cover = volatility_cover = peak_count = 0
    price_lo, price_hi, endpoint_hi, sigma_lo, sigma_hi = [], [], [], [], []
    plugin_prices, ml_ranks = [], []
    frobenius_errors, correlations = [], []
    truth_price = analytic(c)["price"]
    for _ in range(repetitions):
        # These are independent physical-measure returns. Their mean is not
        # inserted into risk-neutral pricing; only the first marginal sigma is.
        returns = interval * log_drift + math.sqrt(interval) * rng.standard_normal((n, d)) @ chol.T
        estimate = estimate_gaussian_returns(returns, interval)
        ordinary = gaussian_parameter_intervals(estimate, FAMILY_ALPHA)
        simultaneous = gaussian_parameter_intervals(estimate, FAMILY_ALPHA, familywise=True)
        vhat = np.diag(estimate.covariance_unbiased)
        mean_estimates.append(estimate.log_drift)
        variance_estimates.append(vhat)
        mu_hat = estimate.log_drift + dividend + 0.5 * vhat
        mu_estimates.append(mu_hat)
        student_pivots.append((estimate.log_drift - log_drift) / np.sqrt(vhat / (n * interval)))
        variance_pivots.append(nu * vhat / np.diag(covariance))
        mean_cover += ((ordinary.log_drift_lower <= log_drift)
                       & (log_drift <= ordinary.log_drift_upper))
        variance_cover += ((ordinary.variance_lower <= np.diag(covariance))
                           & (np.diag(covariance) <= ordinary.variance_upper))
        cm = ((simultaneous.log_drift_lower <= log_drift)
              & (log_drift <= simultaneous.log_drift_upper))
        cv = ((simultaneous.variance_lower <= np.diag(covariance))
              & (np.diag(covariance) <= simultaneous.variance_upper))
        family_mean_cover += cm
        family_variance_cover += cv
        family_joint += int(np.all(cm & cv))
        mu_lo = simultaneous.log_drift_lower + dividend + 0.5 * simultaneous.variance_lower
        mu_hi = simultaneous.log_drift_upper + dividend + 0.5 * simultaneous.variance_upper
        family_mu_cover += (mu_lo <= PHYSICAL_TOTAL_DRIFT) & (PHYSICAL_TOTAL_DRIFT <= mu_hi)
        low, high = float(ordinary.volatility_lower[0]), float(ordinary.volatility_upper[0])
        image = capped_price_interval(c, low, high)
        # Tiny tolerance guards floating comparison only; it is disclosed and
        # is not a mathematical enlargement by an arithmetic error bound.
        price_tolerance = 1e-12 * max(1.0, abs(truth_price))
        proper_price_cover += int(image["lower"] - price_tolerance <= truth_price
                                  <= image["upper"] + price_tolerance)
        naive_price_cover += int(image["lower"] - price_tolerance <= truth_price
                                 <= image["endpoint_only_upper"] + price_tolerance)
        volatility_cover += int(low <= c.volatility <= high)
        peak_count += int(image["interior_peak_included"])
        price_lo.append(image["lower"])
        price_hi.append(image["upper"])
        endpoint_hi.append(image["endpoint_only_upper"])
        sigma_lo.append(low)
        sigma_hi.append(high)
        plugin_prices.append(analytic(replace(c, volatility=math.sqrt(vhat[0])))["price"])
        ml_ranks.append(estimate.scatter_rank)
        frobenius_errors.append(float(np.linalg.norm(estimate.covariance_ml_candidate - covariance)))
        est_vol = np.sqrt(vhat)
        correlations.append(estimate.covariance_unbiased / np.outer(est_vol, est_vol))
    means, variances, mus = np.asarray(mean_estimates), np.asarray(variance_estimates), np.asarray(mu_estimates)
    student, variance_pivot = np.asarray(student_pivots), np.asarray(variance_pivots)
    parameter_rows = []
    for j in range(d):
        parameter_rows.append({
            "dimension": j + 1, "observations": n, "interval": interval,
            "observation_horizon": n * interval, "repetitions": repetitions, "seed": seed,
            "true_log_drift": float(log_drift[j]), "true_physical_total_drift": float(PHYSICAL_TOTAL_DRIFT[j]),
            "true_variance": float(covariance[j, j]), "true_volatility": float(vols[j]),
            "log_drift_estimate_mean": float(means[:, j].mean()),
            "log_drift_bias": float(means[:, j].mean() - log_drift[j]),
            "log_drift_empirical_standard_deviation": float(means[:, j].std(ddof=1)),
            "log_drift_theoretical_standard_deviation": float(vols[j] / math.sqrt(n * interval)),
            "variance_unbiased_estimate_mean": float(variances[:, j].mean()),
            "variance_bias": float(variances[:, j].mean() - covariance[j, j]),
            "variance_empirical_standard_deviation": float(variances[:, j].std(ddof=1)),
            "variance_theoretical_standard_deviation": float(covariance[j, j] * math.sqrt(2.0 / nu)),
            "physical_total_drift_estimate_mean": float(mus[:, j].mean()),
            "physical_total_drift_estimate_standard_deviation": float(mus[:, j].std(ddof=1)),
            "student_pivot_mean": float(student[:, j].mean()),
            "student_pivot_variance": float(student[:, j].var(ddof=1)),
            "student_pivot_theoretical_variance": nu / (nu - 2),
            "chi_square_pivot_mean": float(variance_pivot[:, j].mean()),
            "chi_square_pivot_variance": float(variance_pivot[:, j].var(ddof=1)),
            "chi_square_pivot_theoretical_mean": nu,
            "chi_square_pivot_theoretical_variance": 2 * nu,
            "volatility_interval_mean_width": (float(np.mean(np.asarray(sigma_hi) - np.asarray(sigma_lo)))
                                                if j == 0 else None),
            "covariance_ml_candidate_frobenius_error_mean": float(np.mean(frobenius_errors)),
            "correlation_estimate_12_mean": float(np.mean(np.asarray(correlations)[:, 0, 1])),
            "correlation_estimate_13_mean": float(np.mean(np.asarray(correlations)[:, 0, 2])),
            "correlation_estimate_23_mean": float(np.mean(np.asarray(correlations)[:, 1, 2])),
            "full_rank_ml_count": sum(rank == d for rank in ml_ranks),
        })
    coverage_rows = []
    events = []
    for j in range(d):
        events.extend([(f"log_drift_{j+1}_marginal", int(mean_cover[j]), FAMILY_ALPHA),
                       (f"variance_{j+1}_marginal", int(variance_cover[j]), FAMILY_ALPHA),
                       (f"log_drift_{j+1}_bonferroni", int(family_mean_cover[j]), FAMILY_ALPHA / (2*d)),
                       (f"variance_{j+1}_bonferroni", int(family_variance_cover[j]), FAMILY_ALPHA / (2*d)),
                       (f"physical_total_drift_{j+1}_bonferroni_image", int(family_mu_cover[j]), FAMILY_ALPHA)])
    events.extend([("joint_2d_parameters_bonferroni", family_joint, FAMILY_ALPHA),
                   ("first_volatility_marginal", volatility_cover, FAMILY_ALPHA),
                   ("proper_capped_price_image", proper_price_cover, FAMILY_ALPHA),
                   ("endpoint_only_price_image_no_coverage_guarantee", naive_price_cover, None)])
    for event, count, alpha in events:
        lo, hi = binomial_coverage_interval(count, repetitions)
        coverage_rows.append({"event": event, "observations": n, "interval": interval,
                              "observation_horizon": n*interval, "repetitions": repetitions,
                              "seed": seed, "covered": count, "coverage": count/repetitions,
                              "coverage_binomial_lower_95": lo, "coverage_binomial_upper_95": hi,
                              "construction_alpha": alpha,
                              "status": "synthetic independent-dataset coverage diagnostic; quantiles not certified"})
    price_summary = {"observations": n, "interval": interval, "repetitions": repetitions, "seed": seed,
                     "true_price": truth_price, "true_volatility": c.volatility,
                     "volatility_peak": capped_volatility_peak(c),
                     "interior_peak_included_count": peak_count,
                     "plugin_price_mean": float(np.mean(plugin_prices)),
                     "plugin_price_bias": float(np.mean(plugin_prices) - truth_price),
                     "plugin_price_rmse": float(np.sqrt(np.mean((np.asarray(plugin_prices)-truth_price)**2))),
                     "mean_volatility_interval_lower": float(np.mean(sigma_lo)),
                     "mean_volatility_interval_upper": float(np.mean(sigma_hi)),
                     "mean_proper_price_interval_lower": float(np.mean(price_lo)),
                     "mean_proper_price_interval_upper": float(np.mean(price_hi)),
                     "mean_proper_price_interval_width": float(np.mean(np.asarray(price_hi)-np.asarray(price_lo))),
                     "mean_endpoint_only_price_upper": float(np.mean(endpoint_hi)),
                     "proper_price_covered": proper_price_cover, "endpoint_only_price_covered": naive_price_cover,
                     "arithmetic_comparison_tolerance": 1e-12*max(1.0, abs(truth_price))}
    return parameter_rows, coverage_rows, price_summary


def run(output: Path, repetitions: int = 1000, seed: int = 926341) -> dict:
    if not isinstance(repetitions, int) or repetitions < 2:
        raise ValueError("At least two repeated datasets are required.")
    output.mkdir(parents=True, exist_ok=True)
    cases = list(synthetic_contracts())
    base = Contract()
    peak = capped_volatility_peak(base)
    assert peak is not None
    cases.append(("volatility_peak", replace(base, volatility=peak)))
    parameter_rows, coverage_rows, price_rows, case_meta = [], [], [], []
    for index, (name, c) in enumerate(cases):
        parameters, coverage, prices = simulate_case(c, 252, 1.0/252.0, repetitions, seed+1009*index)
        parameter_rows.extend({"case": name, **row} for row in parameters)
        coverage_rows.extend({"case": name, **row} for row in coverage)
        price_rows.append({"case": name, **prices})
        vols, covariance, dividend, drift = case_parameters(c)
        case_meta.append({"case": name, "contract": asdict(c), "volatilities": vols.tolist(),
                          "covariance": covariance.tolist(), "dividends": dividend.tolist(),
                          "physical_total_drift": PHYSICAL_TOTAL_DRIFT.tolist(),
                          "physical_log_drift": drift.tolist()})
    frequency_rows, frequency_coverage = [], []
    regimes = [("fixed_one_year_horizon", n, 1.0/n) for n in (63, 252, 1008)]
    regimes += [("fixed_daily_frequency_longer_horizon", n, 1.0/252.0) for n in (252, 1008)]
    for index, (regime, n, interval) in enumerate(regimes):
        parameters, coverage, _ = simulate_case(base, n, interval, repetitions, seed+100000+1009*index)
        frequency_rows.extend({"regime": regime, **row} for row in parameters)
        frequency_coverage.extend({"regime": regime, **row} for row in coverage)
    deterministic_image = capped_price_interval(replace(base, volatility=peak), 0.9*peak, 1.1*peak)
    counterexample = {"contract": asdict(replace(base, volatility=peak)),
                      "volatility_lower": 0.9*peak, "volatility_upper": 1.1*peak,
                      "true_price_at_peak": analytic(replace(base, volatility=peak))["price"],
                      **deterministic_image}
    singular_returns = np.array([[0.01, -0.01, 0.005], [0.02, 0.01, -0.005],
                                 [-0.01, 0.0, 0.02]])
    deficient = estimate_gaussian_returns(singular_returns, 1.0/252.0)
    shrunk = trace_target_shrinkage(deficient, 0.1)
    pathology = {"observations": 3, "dimensions": 3, "scatter_rank": deficient.scatter_rank,
                 "covariance_ml_status": deficient.covariance_ml_status,
                 "covariance_candidate_eigenvalues": np.linalg.eigvalsh(deficient.covariance_ml_candidate).tolist(),
                 "shrinkage_strength": 0.1, "shrunken_eigenvalues": np.linalg.eigvalsh(shrunk).tolist(),
                 "remedy_status": "fixed-strength sample-trace target regularisation; not maximum likelihood; no exact CI transplanted"}
    for filename, rows in (("parameter_estimation.csv", parameter_rows), ("coverage.csv", coverage_rows),
                           ("price_uncertainty.csv", price_rows), ("frequency_horizon.csv", frequency_rows),
                           ("frequency_horizon_coverage.csv", frequency_coverage)):
        write_csv(output / filename, rows)
    (output / "counterexamples.json").write_text(json.dumps({"nonmonotone_price": counterexample,
                                                             "singular_covariance": pathology}, indent=2)+"\n")
    manifest = {
        "data_status": "entirely synthetic Gaussian returns; no licensed, downloaded or empirical observations",
        "raw_return_storage": "not retained; only code, seeds and aggregate synthetic summaries are saved",
        "study_design": "predeclared eight cases; 1000 independent datasets per case by default; no post-hoc exclusions",
        "model": "independent equally spaced multivariate physical Gaussian log returns of constant-coefficient GBM",
        "physical_pricing_distinction": "physical total drift estimates are not inserted into risk-neutral pricing; r-q remains pricing drift",
        "pricing_scope": "one-asset capped call using calibrated first marginal volatility; other marginals demonstrate joint inference only",
        "correlation_matrix": CORRELATION.tolist(),
        "correlation_eigenvalues": np.linalg.eigvalsh(CORRELATION).tolist(),
        "observations_per_case": 252, "interval_per_case": 1.0/252.0,
        "repetitions": repetitions, "seed": seed,
        "case_seed_rule": "seed+1009*case_index", "regime_seed_rule": "seed+100000+1009*regime_index",
        "family_alpha": FAMILY_ALPHA,
        "marginal_intervals": "exact Student-t log-drift and chi-square variance under independent Gaussian returns",
        "family_intervals": "Bonferroni split alpha/(2d) across d log drifts and d variances; union-bound coverage at least 1-alpha",
        "covariance_estimation": "C/(n*delta) ML only for n>d and full-rank scatter; C/((n-1)*delta) unbiased distinct",
        "rank_status": "NumPy floating-point tolerance diagnostic; not an exact rank certificate; n<=d always flagged as no positive-definite covariance MLE",
        "regularisation": "fixed lambda=.1 with data-dependent sample-trace target in singular demonstration only; not ML or exact inferred coverage",
        "price_interval": "image of marginal volatility CI under fixed contract/rate/dividend; endpoints for minimum and unique positive interior volatility maximum for maximum",
        "price_uncertainty_exclusions": "no uncertainty in S,K,H,T,r,q, model form, market microstructure or numerical discretisation is covered",
        "drift_information": "at fixed observation horizon, log-drift standard deviation sigma/sqrt(n*delta) does not vanish with sampling frequency",
        "coverage_diagnostic": "Clopper-Pearson 95 per cent pointwise intervals for counts over independent synthetic datasets; no joint coverage assertion across diagnostic rows",
        "arithmetic_status": "float64 normal generation, Student/chi-square/beta quantiles and prices; not outward rounded or formally certified",
        "price_coverage_comparison_tolerance": "1e-12*max(1,absolute true price), not a certified arithmetic error bound",
        "originality_status": "established Gaussian inference plus explicitly derived project-specific capped-price uncertainty propagation; no literature priority claimed",
        "cases": case_meta, "frequency_horizon_regimes": [{"regime": label, "observations": n, "interval": delta} for label,n,delta in regimes],
        "counterexamples": "counterexamples.json",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "code_sha256": {"run_calibration_study.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                          "gbm_calibration.py": hashlib.sha256(Path(__file__).with_name("gbm_calibration.py").read_bytes()).hexdigest(),
                          "capped_gbm.py": hashlib.sha256(Path(capped_gbm.__file__).read_bytes()).hexdigest()},
        "results": {name: {"sha256": hashlib.sha256((output/name).read_bytes()).hexdigest(),
                            "bytes": (output/name).stat().st_size} for name in STUDY_RESULT_FILENAMES},
    }
    (output / "study_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/capped_gbm/calibration"))
    parser.add_argument("--repetitions", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=926341)
    args = parser.parse_args()
    manifest = run(args.output, args.repetitions, args.seed)
    print(json.dumps({"output": str(args.output.resolve()), "cases": len(manifest["cases"]),
                      "repetitions": manifest["repetitions"], "data_status": manifest["data_status"]}, indent=2))
