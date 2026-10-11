"""Reproduce the wholly synthetic capped-GBM price/Greek comparison study."""

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

from capped_gbm import (Contract, GREEKS, HeatGrid, analytic,
                        analytic_bumped_greeks, bumped_greeks, default_bumps,
                        domain_price_bound, heat_backward_euler,
                        heat_backward_euler_thomas, monte_carlo,
                        pde_greek_domain_bounds, pde_price, study_domain,
                        synthetic_contracts)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path, mc_pairs: int = 131072, seed: int = 732451) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    cases = list(synthetic_contracts())
    meshes = (128, 256, 512, 1024, 2048, 4096, 8192, 16384)
    family_alpha = 0.05
    raw_alpha = family_alpha / (len(cases) * (len(GREEKS) + 1))
    pde_rows, mc_rows, domain_rows, solver_rows, refinement_rows = [], [], [], [], []
    case_meta = []
    for case_index, (name, c) in enumerate(cases):
        lower, upper = study_domain(c)
        reference = analytic(c)
        case_meta.append({"case": name, "parameters": asdict(c),
                          "heat_lower": lower, "heat_upper": upper,
                          "analytic": reference})
        for n in meshes:
            grid = HeatGrid(lower, upper, n, math.ceil(n * n / 8))
            bumps = default_bumps(c, grid)
            estimate = bumped_greeks(c, grid, bumps)
            exact_bump = analytic_bumped_greeks(c, grid, bumps)
            domain = pde_greek_domain_bounds(c, grid, bumps)
            for quantity in GREEKS:
                pde_rows.append({"case": name, "quantity": quantity,
                                 "space_intervals": n, "time_steps": grid.time_steps,
                                 "heat_lower": lower, "heat_upper": upper,
                                 "space_spacing": grid.spacing,
                                 "time_spacing": c.maturity / grid.time_steps,
                                 "spot_bump": bumps["spot"],
                                 "volatility_bump": bumps["volatility"],
                                 "rate_bump": bumps["rate"],
                                 "maturity_bump": bumps["maturity"],
                                 "estimate": estimate[quantity], "analytic": reference[quantity],
                                 "total_error": estimate[quantity] - reference[quantity],
                                 "analytic_bumped_quantity": exact_bump[quantity],
                                 "analytic_bump_bias": exact_bump[quantity] - reference[quantity],
                                 "pde_bumped_minus_analytic_bumped": estimate[quantity] - exact_bump[quantity],
                                 "domain_error_upper_bound": domain[quantity],
                                 "arithmetic_status": "float64; no outward-rounded arithmetic certificate"})
        fine_grid = HeatGrid(lower, upper, meshes[-1], math.ceil(meshes[-1] ** 2 / 8))
        for nt in (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
            grid = HeatGrid(lower, upper, meshes[-1], nt)
            price = pde_price(c, grid)
            refinement_rows.append({"case": name, "refinement": "time_at_fixed_space",
                                    "space_intervals": grid.intervals, "time_steps": nt,
                                    "space_spacing": grid.spacing, "time_spacing": c.maturity / nt,
                                    "estimate": price, "analytic": reference["price"],
                                    "error": price - reference["price"]})
        for n in meshes:
            grid = HeatGrid(lower, upper, n, fine_grid.time_steps)
            price = pde_price(c, grid)
            refinement_rows.append({"case": name, "refinement": "space_at_fixed_time",
                                    "space_intervals": n, "time_steps": grid.time_steps,
                                    "space_spacing": grid.spacing, "time_spacing": c.maturity / grid.time_steps,
                                    "estimate": price, "analytic": reference["price"],
                                    "error": price - reference["price"]})
        # Match h and k exactly when enlarging the domain, rather than changing
        # the space discretisation at the same time as the boundary truncation.
        base = fine_grid
        extra = math.ceil(0.5 * (upper - lower) / base.spacing)
        wide = HeatGrid(lower - extra * base.spacing, upper + extra * base.spacing,
                        base.intervals + 2 * extra, base.time_steps)
        p0, p1 = pde_price(c, base), pde_price(c, wide)
        domain_rows.append({"case": name, "base_lower": lower, "base_upper": upper,
                            "enlarged_lower": wide.lower, "enlarged_upper": wide.upper,
                            "base_intervals": base.intervals, "enlarged_intervals": wide.intervals,
                            "base_spacing": base.spacing, "enlarged_spacing": wide.spacing,
                            "time_steps": base.time_steps, "base_price": p0, "enlarged_price": p1,
                            "difference": p1 - p0, "base_domain_bound": domain_price_bound(c, base),
                            "enlarged_domain_bound": domain_price_bound(c, wide)})
        for nt in (1, 3, 64):
            grid = HeatGrid(lower, upper, 64, nt)
            _, spectral = heat_backward_euler(c, grid)
            _, direct = heat_backward_euler_thomas(c, grid)
            solver_rows.append({"case": name, "space_intervals": 64, "time_steps": nt,
                                "lambda": 0.5 * c.volatility ** 2 * (c.maturity / nt) / grid.spacing ** 2,
                                "maximum_absolute_dst_thomas_difference": float(np.max(np.abs(spectral - direct))),
                                "minimum_heat_value": float(spectral.min()),
                                "maximum_heat_value": float(spectral.max()), "cap": c.cap,
                                "status": "numerical cross-check; exact-arithmetic stability proved separately"})
        for row in monte_carlo(c, mc_pairs, max(2048, mc_pairs // 8), seed + case_index * 1009,
                              alpha=raw_alpha):
            mc_rows.append({"case": name, **row})
    write_csv(output / "pde_price_greek_comparison.csv", pde_rows)
    write_csv(output / "independent_mesh_refinement.csv", refinement_rows)
    write_csv(output / "monte_carlo_comparison.csv", mc_rows)
    write_csv(output / "domain_enlargement_check.csv", domain_rows)
    write_csv(output / "backward_euler_solver_check.csv", solver_rows)
    final_rows = [row for row in pde_rows if row["space_intervals"] == meshes[-1]]
    final_errors = {q: max(abs(row["total_error"]) for row in final_rows if row["quantity"] == q) for q in GREEKS}
    failed_raw = [f"{row['case']}:{row['quantity']}" for row in mc_rows
                  if row["method"].startswith("raw_exact_terminal")
                  and abs(row["error"]) > row["finite_bound_half_width"]]
    manifest = {"data_status": "entirely synthetic designed contracts; no LSEG or other proprietary data",
                "model": "risk-neutral GBM with continuous dividend yield and capped call payoff",
                "greek_units": {"delta": "per unit spot", "gamma": "per unit spot squared",
                                "vega": "per unit volatility", "rho": "per unit continuously compounded rate",
                                "theta": "per year of calendar time; minus derivative in time to maturity"},
                "analytic_baseline": "closed-form call spread; dominant baseline for this model",
                "pde_method": "uniform central heat space difference/backward Euler; DST-I implements same tridiagonal recurrence",
                "time_space_design": "coupled N=ceil(M squared/8), plus independent refinements",
                "domain_status": "fixed per-case nine-sigma boxes; finite domain bounds separate from grid/bump error",
                "convergence_scope": "exact arithmetic theorem with growing domains; fixed-domain numerical study with explicit tail bounds",
                "bump_design": "spot S*(0.006*h**(1/8)+h) capped at0.1*S; first parameter bumps proportional h**(1/4); explicit per-row steps",
                "initial_payoff": "actual kinked payoff; no numerical smoothing or Rannacher start-up",
                "raw_mc_familywise_error_probability": family_alpha,
                "raw_mc_number_of_outputs": len(cases) * (len(GREEKS) + 1),
                "raw_mc_individual_error_probability": raw_alpha,
                "finite_confidence_scope": "Hoeffding for bounded raw price/pathwise Greeks/mixed PW-LR Gamma, Chebyshev payoff-LR Gamma; ideal independent exact normals/arithmetic",
                "gamma_estimators": ["payoff_likelihood_ratio", "mixed_pathwise_likelihood_ratio"],
                "floating_point_confidence_limitation": "bound constants and normal generation are not formally certified or outward rounded",
                "normal_interval_scope": "pointwise asymptotic diagnostics only; no familywise finite coverage claim for antithetic/control output",
                "controlled_mc_scope": "independent pilot fitted coefficients clipped to plus/minus10^6; production observations independent conditionally on pilot",
                "cross_method_coupling": "antithetic and controlled production use a shared normal prefix; terminal evaluations charged separately; no confidence claim for paired method differences",
                "mc_pairs": mc_pairs, "raw_production_normals_per_case": 2 * mc_pairs,
                "antithetic_production_pairs_per_case": mc_pairs,
                "controlled_production_pairs_per_case": mc_pairs - max(2048, mc_pairs // 8),
                "control_pilot_pairs_per_case": max(2048, mc_pairs // 8), "seed": seed,
                "space_intervals": list(meshes), "cases": case_meta,
                "finest_maximum_absolute_errors": final_errors,
                "raw_outputs_outside_conservative_bounds": failed_raw,
                "solver_maximum_absolute_crosscheck_difference": max(row["maximum_absolute_dst_thomas_difference"] for row in solver_rows),
                "domain_maximum_absolute_matched_grid_difference": max(abs(row["difference"]) for row in domain_rows),
                "environment": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()},
                "code_sha256": {Path(__file__).name: hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                 "capped_gbm.py": hashlib.sha256(Path(__file__).with_name("capped_gbm.py").read_bytes()).hexdigest()},
                "results": {file.name: {"sha256": hashlib.sha256(file.read_bytes()).hexdigest(), "bytes": file.stat().st_size}
                            for file in sorted(output.glob("*.csv"))},
                "originality_status": "established classical methods and project-specific transparent verification; no novel quantum or pricing theory claimed"}
    (output / "study_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/capped_gbm"))
    parser.add_argument("--mc-pairs", type=int, default=131072)
    parser.add_argument("--seed", type=int, default=732451)
    args = parser.parse_args()
    result = run(args.output, args.mc_pairs, args.seed)
    print(json.dumps({"output": str(args.output.resolve()),
                      "finest_maximum_absolute_errors": result["finest_maximum_absolute_errors"],
                      "raw_outputs_outside_conservative_bounds": result["raw_outputs_outside_conservative_bounds"],
                      "solver_crosscheck": result["solver_maximum_absolute_crosscheck_difference"]}, indent=2))
