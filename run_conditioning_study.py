"""Reproduce the synthetic projection/control study without acquiring data."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import capped_gbm

from capped_gbm import analytic, synthetic_contracts
from conditioning_controls import (METHODS, conditional_population_moments,
                                    cost_counts, estimator_values, fit_pilot,
                                    raw_population_moments, residual_population_second_moment,
                                    split_observations)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def controlled_variance(variance: float, covariance: float,
                        control_variance: float, beta: float) -> float:
    return variance - 2.0 * beta * covariance + beta ** 2 * control_variance


def run(output: Path, pilot_size: int = 8192, production_size: int = 65536,
        seed: int = 840731, repeated_pilots: int = 32) -> dict:
    if min(pilot_size, production_size) < 2 or repeated_pilots < 1:
        raise ValueError("Use at least two pilot/production draws and one repeated pilot.")
    output.mkdir(parents=True, exist_ok=True)
    fractions = (0.25, 0.5, 0.9)
    cases = list(synthetic_contracts())
    population_rows, estimates, pilot_rows, quadrature_rows, identities = [], [], [], [], []
    meta, main_pilots = [], []
    range_diagnostics = {phase: {} for phase in (
        "production", "main_pilots", "repeated_pilots", "gh160_moments", "gh320_moments", "kinked_joint_diagnostic")}
    for case_index, (name, c) in enumerate(cases):
        raw = raw_population_moments(c)
        meta.append({"case": name, "parameters": asdict(c), "analytic_price": raw["mean"]})
        for split_index, u in enumerate(fractions):
            stream_base = seed + case_index * 10007 + split_index * 1009
            pilot = fit_pilot(c, u, pilot_size, stream_base, diagnostics=range_diagnostics["main_pilots"])
            main_pilots.append({"case": name, "conditioning_fraction": u, **asdict(pilot)})
            rng = np.random.default_rng(stream_base + 4000001)
            obs = split_observations(c, u, rng.standard_normal(production_size),
                                     rng.standard_normal(production_size), range_diagnostics["production"])
            coarse = conditional_population_moments(c, u, 160, range_diagnostics["gh160_moments"])
            cond = conditional_population_moments(c, u, 320, range_diagnostics["gh320_moments"])
            variance_by_method = {
                "raw_payoff": raw["variance"], "exact_conditional": cond["variance"],
                "raw_pilot_stock_control": controlled_variance(raw["variance"], raw["stock_covariance"],
                                                                raw["stock_control_variance"], pilot.raw),
                "conditional_raw_pilot_coefficient": controlled_variance(cond["variance"], cond["stock_covariance"],
                                                                         cond["stock_control_variance"], pilot.raw),
                "conditional_own_pilot_coefficient": controlled_variance(cond["variance"], cond["stock_covariance"],
                                                                         cond["stock_control_variance"], pilot.conditional),
            }
            variance_by_method["conditional_unconditioned_stock_control"] = (
                variance_by_method["conditional_raw_pilot_coefficient"]
                + pilot.raw ** 2 * (raw["stock_control_variance"] - cond["stock_control_variance"]))
            # Shared study streams couple rows. They are not independent replications.
            for method, values in estimator_values(obs, pilot).items():
                mean = float(values.mean())
                variance = float(values.var(ddof=1))
                beta = (pilot.conditional if method == "conditional_own_pilot_coefficient"
                        else pilot.raw if "control" in method or "coefficient" in method else 0.0)
                estimates.append({
                    "case": name, "conditioning_fraction": u, "method": method,
                    "estimate": mean, "analytic_price": raw["mean"], "error": mean - raw["mean"],
                    "sample_variance": variance, "sample_standard_error": math.sqrt(variance / production_size),
                    "population_variance_diagnostic_given_pilot": variance_by_method[method],
                    "fitted_coefficient": beta,
                    "population_variance_diagnostic_to_raw": variance_by_method[method] / raw["variance"],
                    "raw_variance_to_population_diagnostic": raw["variance"] / variance_by_method[method],
                    **cost_counts(method, production_size, pilot_size),
                    "source_arithmetic": "float64; moment quadrature is not a validated certificate",
                    "confidence_status": "reported standard error is descriptive; no Gaussian or finite-sample coverage claim",
                })
            exact_counts = {key: 0 for key in cost_counts("raw_payoff", 0, 0)}
            estimates.append({"case": name, "conditioning_fraction": u,
                              "method": "analytic_call_spread", "estimate": analytic(c)["price"],
                              "analytic_price": raw["mean"], "error": 0.0,
                              "sample_variance": 0.0, "sample_standard_error": 0.0,
                              "population_variance_diagnostic_given_pilot": 0.0,
                              "fitted_coefficient": 0.0, **exact_counts,
                              "normal_cdf_evaluations": 4,
                              "source_arithmetic": "float64 analytic evaluation; zero sampling error, round-off not certified",
                              "confidence_status": "no sampling"})
            population_rows.append({
                "case": name, "conditioning_fraction": u,
                "raw_variance": raw["variance"], "raw_oracle_beta": raw["oracle_beta"],
                "raw_oracle_control_variance": raw["oracle_control_variance"],
                "conditional_variance": cond["variance"], "conditional_oracle_beta": cond["oracle_beta"],
                "conditional_oracle_control_variance": cond["oracle_control_variance"],
                "conditional_with_raw_oracle_beta_variance": cond["raw_oracle_coefficient_conditional_variance"],
                "conditional_unconditioned_control_variance": cond["unconditioned_control_variance"],
                "unconditioned_control_added_variance": cond["unconditioned_control_added_variance"],
                "raw_oracle_to_conditional_oracle_variance_ratio": raw["oracle_control_variance"] / cond["oracle_control_variance"],
                "production_kernel_cost_break_even_ratio": raw["oracle_control_variance"] / cond["oracle_control_variance"],
                "break_even_scope": "oracle asymptotic production variance only; independent pilot and kernel costs excluded",
                "raw_control_variance": raw["stock_control_variance"],
                "conditional_control_variance": cond["stock_control_variance"],
                "raw_covariance": raw["stock_covariance"], "conditional_covariance": cond["stock_covariance"],
                "projection_removed_variance": cond["projection_removed_variance"],
                "same_beta_projection_removed_variance": cond["same_beta_projection_removed_variance"],
                "coefficient_reoptimisation_reduction": cond["coefficient_reoptimisation_reduction"],
                "arithmetic_status": "raw truncated moments analytic in float64; conditional moments GH320, unvalidated",
            })
            for field in ("quadrature_mean", "quadrature_control_mean", "variance", "stock_covariance",
                          "oracle_beta", "oracle_control_variance"):
                quadrature_rows.append({"case": name, "conditioning_fraction": u, "quantity": field,
                                        "gh160": coarse[field], "gh320": cond[field],
                                        "refinement_difference": cond[field] - coarse[field],
                                        "status": "quadrature refinement diagnostic; not a rigorous remainder bound"})
            # Numerical Pythagoras check using an independent conditional inner rule.
            from conditioning_controls import conditional_value, gaussian_rule
            nodes, weights = gaussian_rule(80)
            outer = nodes[:, None]
            inner = nodes[None, :]
            grid_obs = split_observations(c, u, outer, inner, range_diagnostics["kinked_joint_diagnostic"])
            m, cg = grid_obs["conditional"], grid_obs["conditional_stock_control"]
            residual = ((grid_obs["raw"] - m) - raw["oracle_beta"]
                        * (grid_obs["stock_control"] - cg))
            residual_second = float(weights @ (residual ** 2 @ weights))
            difference = raw["oracle_control_variance"] - cond["raw_oracle_coefficient_conditional_variance"]
            residual_smooth_160 = residual_population_second_moment(c, u, raw["oracle_beta"], 160)
            residual_smooth_320 = residual_population_second_moment(c, u, raw["oracle_beta"], 320)
            identities.append({"case": name, "conditioning_fraction": u,
                               "projection_removed_variance_by_moments": difference,
                               "residual_second_moment_analytic_inner_gh160": residual_smooth_160,
                               "residual_second_moment_analytic_inner_gh320": residual_smooth_320,
                               "smooth_identity_discrepancy": residual_smooth_320 - difference,
                               "smooth_residual_quadrature_difference": residual_smooth_320 - residual_smooth_160,
                               "residual_second_moment_2d_gh80": residual_second,
                               "identity_discrepancy": residual_second - difference,
                               "orthogonality_raw_residual_times_conditional": float(
                                   weights @ ((grid_obs["raw"] - m) * (m - raw["mean"]) @ weights)),
                               "status": "smooth analytic-inner check and deliberately weaker kinked 2D quadrature; neither is a certificate"})
            for replicate in range(repeated_pilots):
                fitted = fit_pilot(c, u, pilot_size, stream_base + 8000003 + replicate * 1543,
                                   diagnostics=range_diagnostics["repeated_pilots"])
                raw_fit_cond_penalty = (fitted.raw - cond["oracle_beta"]) ** 2 * cond["stock_control_variance"]
                own_fit_penalty = (fitted.conditional - cond["oracle_beta"]) ** 2 * cond["stock_control_variance"]
                pilot_rows.append({"case": name, "conditioning_fraction": u, "replicate": replicate,
                                   "pilot_size": pilot_size, "raw_fitted_beta": fitted.raw,
                                   "conditional_fitted_beta": fitted.conditional,
                                   "raw_unclipped_beta": fitted.raw_unclipped,
                                   "conditional_unclipped_beta": fitted.conditional_unclipped,
                                   "raw_clipped": fitted.raw_clipped,
                                   "conditional_clipped": fitted.conditional_clipped,
                                   "raw_fit_conditional_excess_variance": raw_fit_cond_penalty,
                                   "own_fit_conditional_excess_variance": own_fit_penalty,
                                   "own_fit_worse_than_raw_fit_for_this_pilot": own_fit_penalty > raw_fit_cond_penalty,
                                   "status": "32 independent pilot diagnostics; not a universal ranking or expected-cost theorem"})
    write_csv(output / "conditioning_price_comparison.csv", estimates)
    write_csv(output / "population_variance_comparison.csv", population_rows)
    write_csv(output / "quadrature_refinement.csv", quadrature_rows)
    write_csv(output / "projection_identity_diagnostics.csv", identities)
    write_csv(output / "independent_pilot_risk.csv", pilot_rows)
    write_csv(output / "main_pilot_coefficients.csv", main_pilots)
    positive_gains = [row["raw_oracle_to_conditional_oracle_variance_ratio"] for row in population_rows]
    adverse = [row for row in pilot_rows if row["own_fit_worse_than_raw_fit_for_this_pilot"]]
    manifest = {
        "data_status": "entirely synthetic designed contracts; no market-data acquisition or redistribution",
        "purpose": "test L2 conditioning/control geometry and estimator cost limitations within classical GBM",
        "underlying_law": "classical risk-neutral GBM, continuous dividend yield, exact terminal simulation",
        "price_only_scope": "no conditioned Greek formula is claimed",
        "cases": meta, "conditioning_fractions": list(fractions),
        "main_independent_pilot_size": pilot_size, "production_draws_per_case_split": production_size,
        "coefficient_absolute_bound": 1e6, "repeated_independent_pilots": repeated_pilots,
        "clipping_rule": "Euclidean coefficient norm at most 1e6; scalar [−1e6,1e6]",
        "main_pilot_clipping_events": sum(r["raw_clipped"] + r["conditional_clipped"] for r in main_pilots),
        "repeated_pilot_clipping_events": sum(r["raw_clipped"] + r["conditional_clipped"] for r in pilot_rows),
        "seed": seed,
        "analytic_baseline": "call spread; strongest exact-price baseline here, no Monte Carlo or quantum pricing advantage",
        "stream_design": "independent pilot and production; methods share production draws; each split has its own independent stream",
        "time_coupling_limitation": "different fractions are marginal split simulations, not a jointly sampled Brownian path across times",
        "coefficient_sign": "estimate Y−beta C or m−beta c, with controls centred at their known zero means",
        "pilot_scope": "bounded pilot coefficients preserve ideal-arithmetic conditional unbiasedness; no universal finite-pilot superiority",
        "cost_scope": "per-estimator stand-alone normal/exponential/payoff/CDF counts; no equal-wall-clock claim",
        "actual_shared_study_generation_counts": {
            "main_normal_variates": 2 * len(population_rows) * (pilot_size + production_size),
            "repeated_pilot_normal_variates": 2 * len(population_rows) * repeated_pilots * pilot_size,
            "pilot_least_squares_fits": 2 * len(population_rows) * (1 + repeated_pilots),
            "note": "shared joint observations compute both raw and conditional targets; stand-alone kernel counts differ"},
        "quadrature": "smooth conditional moments and analytic-inner residual checks GH160/GH320; kinked joint identity GH80 by GH80 diagnostic only",
        "variance_ratio_raw_oracle_to_conditional_oracle_min": min(positive_gains),
        "variance_ratio_raw_oracle_to_conditional_oracle_max": max(positive_gains),
        "pilot_adverse_rankings_retained": len(adverse), "repeated_pilot_rows": len(pilot_rows),
        "no_practical_quantum_advantage": True,
        "floating_point_status": "NumPy float64, pseudo-random normal generation, stdlib erfc; no outward-rounded certificate",
        "conditional_price_formula": "stable CDF-interval call-spread expression, four CDF evaluations/state; upper-tail intervals use survival differences",
        "conditional_range_fallback_diagnostics": range_diagnostics,
        "numerical_unbiasedness_limitation": "ideal exact-arithmetic projections are unbiased; floating-point CDF evaluation and range clipping are not a certified unbiased implementation",
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
    }
    hashes = {}
    for path in sorted(output.glob("*.csv")):
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name in ("conditioning_controls.py", "run_conditioning_study.py"):
        path = Path(__file__).with_name(name)
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    hashes["capped_gbm.py"] = hashlib.sha256(Path(capped_gbm.__file__).read_bytes()).hexdigest()
    test_name = "test_conditioning_controls.py"
    test_paths = (Path(__file__).with_name(test_name), Path(__file__).parent / "checks" / test_name)
    test_path = next((path for path in test_paths if path.is_file()), None)
    if test_path is not None:
        hashes[test_name] = hashlib.sha256(test_path.read_bytes()).hexdigest()
    manifest["sha256"] = hashes
    (output / "study_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/capped_gbm/conditioning"))
    parser.add_argument("--pilot-size", type=int, default=8192)
    parser.add_argument("--production-size", type=int, default=65536)
    parser.add_argument("--repeated-pilots", type=int, default=32)
    parser.add_argument("--seed", type=int, default=840731)
    args = parser.parse_args()
    print(json.dumps(run(args.output, args.pilot_size, args.production_size,
                         args.seed, args.repeated_pilots), indent=2))
