"""Preregistered synthetic budget-constrained static GBM hedging comparison."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import warnings

import numpy as np
import scipy

import capped_gbm
from capped_gbm import call_analytic, synthetic_contracts
from gbm_hedging import (hedge_financing, hedge_moments, quadratic_loss,
                         quadrature_moments, residual_square_integral, simulated_gains,
                         slsqp_hedge, solve_constrained_hedge)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def risk_fields(c, mu, theta, moments):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        positive, error = residual_square_integral(c, mu, theta)
    raw = quadratic_loss(theta, moments.matrix, moments.cross, moments.constant)
    term_scale = (abs(moments.constant)+2*float(np.dot(np.abs(moments.cross), np.abs(theta)))
                  +float(np.abs(theta)@np.abs(moments.matrix)@np.abs(theta)))
    arithmetic_floor = 128*np.finfo(float).eps*max(term_scale, np.finfo(float).tiny)
    return {"population_loss_nonnegative_integral": positive, "quadrature_error_estimate": error,
            "population_loss_moment_quadratic": raw, "moment_minus_nonnegative_integral": raw-positive,
            "quadratic_cancellation_scale": term_scale, "heuristic_roundoff_floor": arithmetic_floor,
            "risk_below_heuristic_roundoff_floor": positive <= arithmetic_floor,
            "negative_moment_quadratic": raw < 0,
            "quadrature_warnings": " | ".join(str(item.message) for item in caught),
            "expected_discounted_residual": float(moments.residual_mean-moments.gain_mean@theta),
            "status": "non-certified floating-point residual-square quadrature; raw quadratic retained"}


def solution_fields(solution):
    return {"stock_holding": float(solution.theta[0]), "call_holding": float(solution.theta[1]),
            "active_constraints": ";".join(solution.active_labels),
            "candidates_considered": solution.candidates_considered,
            "kkt_primal_violation": solution.primal_violation, "kkt_dual_violation": solution.dual_violation,
            "kkt_stationarity_residual": solution.stationarity_residual,
            "kkt_complementarity_residual": solution.complementarity_residual,
            "kkt_multipliers": json.dumps(solution.multipliers.tolist()),
            "matrix_minimum_eigenvalue": float(solution.matrix_eigenvalues[0]),
            "matrix_condition_number": solution.matrix_condition_number}


def run(output: Path, protocol_path: Path) -> dict:
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use a fresh output directory; prior outcomes must not be overwritten.")
    output.mkdir(parents=True, exist_ok=True)
    protocol_bytes = protocol_path.read_bytes()
    protocol = json.loads(protocol_bytes)
    source_paths = {"run_hedging_study.py": Path(__file__),
                    "gbm_hedging.py": Path(__file__).with_name("gbm_hedging.py"),
                    "capped_gbm.py": Path(capped_gbm.__file__)}
    registration = {"protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
                    "protocol": protocol, "code_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                                            for name,path in source_paths.items()},
                    "registration_status": "written before any study outcome is computed"}
    (output/"registration.json").write_text(json.dumps(registration, indent=2)+"\n")
    cases = list(synthetic_contracts())
    if [name for name,_ in cases] != protocol["cases"]:
        raise ValueError("The reused contract menu differs from the registered protocol.")
    constrained_rows, centred_rows, projection_rows = [], [], []
    unconstrained_rows, nohedge_rows, quadrature_rows, mc_rows, external_rows = [], [], [], [], []
    failure_rows = []
    seed, paths = protocol["simulation_seed"], protocol["simulation_paths_per_regime"]
    for case_index, (name,c) in enumerate(cases):
        for mu_index, mu in enumerate(protocol["physical_total_returns"]):
            base = {"case": name, "physical_total_return": mu, "pricing_rate": c.rate,
                    "physical_equals_risk_neutral_total_drift": mu == c.rate,
                    "data_status": "synthetic stipulated GBM law"}
            moments = hedge_moments(c, mu)
            eig, cov_eig = np.linalg.eigvalsh(moments.matrix), np.linalg.eigvalsh(moments.gain_covariance)
            if eig[0] <= 0 or cov_eig[0] <= 0:
                failure_rows.append({**base, "stage": "moment_positive_definiteness",
                                     "matrix_eigenvalues": eig.tolist(), "covariance_eigenvalues": cov_eig.tolist(),
                                     "resolution": "not repaired; regime cannot be solved"})
                continue
            theta_mse = np.linalg.solve(moments.matrix, moments.cross)
            theta_variance = np.linalg.solve(moments.gain_covariance, moments.centred_cross)
            inverse_mean = np.linalg.solve(moments.gain_covariance, moments.gain_mean)
            zeta = float(moments.gain_mean@inverse_mean)
            eta = float(moments.residual_mean-moments.gain_mean@theta_variance)
            audit_theta = theta_variance+inverse_mean*eta/(1+zeta)
            gap_formula = eta*eta*zeta/(1+zeta)
            zero = np.zeros(2)
            no_fields = risk_fields(c, mu, zero, moments)
            un_fields = risk_fields(c, mu, theta_mse, moments)
            var_fields = risk_fields(c, mu, theta_variance, moments)
            nohedge_rows.append({**base, **hedge_financing(c,moments,zero), **no_fields})
            unconstrained_rows.append({**base, **hedge_financing(c,moments,theta_mse), **un_fields,
                                       "variance_stock_holding": float(theta_variance[0]),
                                       "variance_call_holding": float(theta_variance[1]),
                                       "variance_hedge_loss_integral": var_fields["population_loss_nonnegative_integral"],
                                       "variance_hedge_loss_quadratic": var_fields["population_loss_moment_quadratic"],
                                       "variance_hedge_expected_residual": var_fields["expected_discounted_residual"],
                                       "variance_mse_gap_formula": gap_formula,
                                       "variance_mse_gap_integral": var_fields["population_loss_nonnegative_integral"]-un_fields["population_loss_nonnegative_integral"],
                                       "sherman_morrison_position_error": float(np.max(np.abs(audit_theta-theta_mse))),
                                       "matrix_condition_number": float(eig[-1]/eig[0]),
                                       "centred_matrix_condition_number": float(cov_eig[-1]/cov_eig[0]),
                                       "scope": "formal unconstrained L2 relaxation; not necessarily an admissible trading strategy"})
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                quad_values, quad_error = quadrature_moments(c, mu)
            exact_values = np.array([*moments.gain_mean, moments.matrix[0,0], moments.matrix[0,1],
                                     moments.matrix[1,1], *moments.cross, moments.constant])
            for quantity, closed, numerical in zip(("gain_mean_stock","gain_mean_call","M11","M12","M22","h1","h2","c"),exact_values,quad_values):
                quadrature_rows.append({**base,"quantity":quantity,"closed_moment":float(closed),
                                        "independent_quadrature":float(numerical),"difference":float(closed-numerical),
                                        "vector_quadrature_error_estimate":quad_error,
                                        "warnings":" | ".join(str(item.message) for item in caught),
                                        "arithmetic_status":"float64 closed moments and independent adaptive quadrature; not certified"})
            gains, residual = simulated_gains(c,mu,paths,seed+1009*(3*case_index+mu_index))
            choices = [("no_hedge", None, zero, no_fields), ("unconstrained_mse",None,theta_mse,un_fields)]
            for kappa in protocol["budget_multipliers"]:
                budget = kappa*moments.premium
                row_base = {**base,"budget_multiplier":kappa,"budget":budget,"premium":moments.premium}
                try:
                    solution = solve_constrained_hedge(moments.matrix,moments.cross,moments.constant,moments.costs,budget)
                    variance_solution = solve_constrained_hedge(moments.gain_covariance,moments.centred_cross,
                                                                 moments.residual_variance,moments.costs,budget)
                    projection = solve_constrained_hedge(np.eye(2),theta_variance,float(theta_variance@theta_variance),moments.costs,budget)
                    generic = slsqp_hedge(moments,budget)
                    generic_risk = risk_fields(c,mu,generic["theta"],moments)
                    risk = risk_fields(c,mu,solution.theta,moments)
                    centred_risk = risk_fields(c,mu,variance_solution.theta,moments)
                    projected_risk = risk_fields(c,mu,projection.theta,moments)
                except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
                    failure_rows.append({**row_base,"stage":"optimisation","error":str(error),"resolution":"retained; no fallback or clipping"})
                    continue
                constrained_rows.append({**row_base,**solution_fields(solution),**hedge_financing(c,moments,solution.theta),**risk,
                                         "population_loss_reduction_from_nohedge":no_fields["population_loss_nonnegative_integral"]-risk["population_loss_nonnegative_integral"],
                                         "constraint_penalty_above_unconstrained_integral":risk["population_loss_nonnegative_integral"]-un_fields["population_loss_nonnegative_integral"],
                                         "constraint_penalty_quadratic_formula":float((solution.theta-theta_mse)@moments.matrix@(solution.theta-theta_mse)),
                                         "slsqp_stock_holding":float(generic["theta"][0]),"slsqp_call_holding":float(generic["theta"][1]),
                                         "slsqp_success":generic["success"],"slsqp_status":generic["status"],"slsqp_message":generic["message"],
                                         "slsqp_iterations":generic["iterations"],"slsqp_primal_violation":generic["primal_violation"],
                                         "slsqp_loss_integral":generic_risk["population_loss_nonnegative_integral"],
                                         "slsqp_minus_active_set_loss_integral":generic_risk["population_loss_nonnegative_integral"]-risk["population_loss_nonnegative_integral"],
                                         "slsqp_position_maximum_difference":float(np.max(np.abs(generic["theta"]-solution.theta))),
                                         "constrained_variance_mse_gap_integral":centred_risk["population_loss_nonnegative_integral"]-risk["population_loss_nonnegative_integral"],
                                         "projected_variance_mse_gap_integral":projected_risk["population_loss_nonnegative_integral"]-risk["population_loss_nonnegative_integral"]})
                centred_rows.append({**row_base,**solution_fields(variance_solution),**hedge_financing(c,moments,variance_solution.theta),**centred_risk,
                                     "minimised_objective":"centred physical residual variance; not mean-square risk"})
                projection_rows.append({**row_base,**solution_fields(projection),**hedge_financing(c,moments,projection.theta),**projected_risk,
                                        "minimised_objective":"Euclidean distance from unconstrained centred-variance holdings; MSE evaluated separately"})
                choices.extend([("constrained_mse",kappa,solution.theta,risk),
                                ("constrained_variance",kappa,variance_solution.theta,centred_risk),
                                ("projected_unconstrained_variance",kappa,projection.theta,projected_risk)])
            for method,kappa,theta,population in choices:
                realised = residual-gains@theta
                squared = realised**2
                mean = float(squared.mean())
                standard_error = float(squared.std(ddof=1)/math.sqrt(paths))
                mc_rows.append({**base,"method":method,"budget_multiplier":kappa,"paths":paths,
                                "seed":seed+1009*(3*case_index+mu_index),"sample_mse":mean,"sample_mse_standard_error":standard_error,
                                "population_loss_nonnegative_integral":population["population_loss_nonnegative_integral"],
                                "difference":mean-population["population_loss_nonnegative_integral"],
                                "standardised_difference":((mean-population["population_loss_nonnegative_integral"])/standard_error if standard_error>0 else None),
                                "sample_expected_residual":float(realised.mean()),"sample_residual_standard_deviation":float(realised.std(ddof=1)),
                                "sample_maximum_absolute_residual":float(np.max(np.abs(realised))),
                                "diagnostic_status":"matched independent validation paths; sample standard errors are asymptotic diagnostics, not finite coverage certificates"})
            upper_cost = call_analytic(c,c.strike+c.cap)["price"]
            external_rows.append({**base,"lower_call_holding":1,"upper_call_holding":-1,
                                  "lower_call_cost":float(moments.costs[1]),"upper_call_cost":upper_cost,
                                  "net_cost":float(moments.costs[1])-upper_cost,"premium":moments.premium,
                                  "exact_arithmetic_replication_loss":0,
                                  "scope":"external instrument-set identity; upper call absent from constrained universe; not a same-space baseline"})
    named_rows = {"constrained_mse.csv":constrained_rows,"constrained_variance.csv":centred_rows,
                  "projected_variance.csv":projection_rows,"unconstrained.csv":unconstrained_rows,
                  "nohedge.csv":nohedge_rows,"moment_quadrature.csv":quadrature_rows,
                  "monte_carlo.csv":mc_rows,"external_spread_reference.csv":external_rows}
    for filename,rows in named_rows.items():
        write_csv(output/filename,rows)
    (output/"failures.json").write_text(json.dumps(failure_rows,indent=2)+"\n")
    manifest = {**registration,"data_status":"entirely synthetic; no raw paths retained",
                "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"platform":platform.platform()},
                "contracts":[{"case":name,"parameters":asdict(c)} for name,c in cases],
                "row_counts":{name:len(rows) for name,rows in named_rows.items()},"failure_count":len(failure_rows),
                "risk_reporting":"non-negative residual-square Gaussian quadrature is primary numerical risk; cancellation-prone moment quadratic and all discrepancies retained",
                "roundoff_floor":"128*machine_epsilon*(sum absolute quadratic terms), heuristic diagnostic only; not an error certificate",
                "constraint_scope":"long-only holdings <=1, initial instrument-cost budget; finite initial cash/borrowing retained; unconstrained relaxation may be unbounded below",
                "objective_scope":"uncentred physical-law MSE includes bias; centred-variance and projected-variance comparators are separately labelled",
                "mc_scope":"same independent normal paths per case/physical law across holdings; no training data or path archive; no finite interval guarantee",
                "arithmetic_status":"all moment formulas, KKT diagnostics, SciPy quadrature/SLSQP and simulated normals are floating-point implementations, not certified",
                "originality_status":"established quadratic optimisation and Gaussian moment methods with project-specific financial comparison; no theoretical priority claim",
                "results":{name:{"sha256":hashlib.sha256((output/name).read_bytes()).hexdigest(),"bytes":(output/name).stat().st_size}
                           for name in (*named_rows,"failures.json","registration.json")}}
    (output/"study_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path("results/capped_gbm/hedging"))
    parser.add_argument("--protocol",type=Path,default=Path(__file__).with_name("hedging_protocol.json"))
    args=parser.parse_args()
    manifest=run(args.output,args.protocol)
    print(json.dumps({"output":str(args.output.resolve()),"row_counts":manifest["row_counts"],
                      "failure_count":manifest["failure_count"]},indent=2))
