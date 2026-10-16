"""Public feasibility screen for conditional interval encodings of basket spreads.

Requires NumPy. This is a classical numerical screen and a small state-vector
correctness check, not a compiled reversible arithmetic circuit or speedup test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
SQRT2 = math.sqrt(2.0)
NORMAL = NormalDist()


def cdf(x):
    x = np.asarray(x, dtype=float)
    out = np.fromiter((0.5 * math.erfc(-v / SQRT2) for v in x.ravel()),
                      dtype=float, count=x.size)
    return out.reshape(x.shape)


def interval_probability(lo, hi):
    """Stable normal interval probability, avoiding subtraction near one."""
    lo, hi = np.broadcast_arrays(lo, hi)
    result = np.empty_like(lo, dtype=float)
    positive = lo >= 0
    result[positive] = cdf(-lo[positive]) - cdf(-hi[positive])
    result[~positive] = cdf(hi[~positive]) - cdf(lo[~positive])
    return np.clip(result, 0.0, 1.0)


def payoff(a, strike, cap):
    return np.clip(a - strike, 0.0, cap)


def call(log_b, loading, strike):
    z = (math.log(strike) - log_b) / loading
    return np.exp(log_b + loading**2 / 2) * cdf(loading - z) - strike * cdf(-z)


def common_conditional(log_b, loading, strike, cap, order=48):
    """Positive truncated-payoff integral; a stable established comparator.

    Quadrature is only a numerical reference. Its error is checked by order
    refinement, rather than described as a rigorous quadrature certificate.
    """
    log_b = np.asarray(log_b)
    z0 = (math.log(strike) - log_b) / loading
    z1 = (math.log(strike + cap) - log_b) / loading
    nodes, weights = np.polynomial.legendre.leggauss(order)
    z = (z0[:, None] + z1[:, None]) / 2 + (z1 - z0)[:, None] * nodes / 2
    integrand = np.clip(np.exp(log_b[:, None] + loading * z) - strike, 0, cap)
    integrand *= np.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    return (z1 - z0) * (integrand @ weights) / 2 + cap * cdf(-z1)


def model(d, correlation, volatility=0.25, maturity=1.0):
    # Stipulated ex-dividend GBM inputs. No observed market data are used.
    sigmas = volatility * np.linspace(0.8, 1.2, d)
    corr = (1 - correlation) * np.eye(d) + correlation * np.ones((d, d))
    cov = np.outer(sigmas, sigmas) * corr * maturity
    chol = np.linalg.cholesky(cov)
    direction = np.linalg.solve(chol, np.ones(d))
    direction /= np.linalg.norm(direction)
    loading = float((chol @ direction)[0])
    _, _, vh = np.linalg.svd(direction[None, :], full_matrices=True)
    residual = chol @ vh[1:].T
    assert np.max(np.abs(residual @ residual.T + loading**2 - cov)) < 1e-12
    r, dividend = 0.03, 0.02
    mu = math.log(100) + (r - dividend - sigmas**2 / 2) * maturity
    weights = np.full(d, 1 / d)
    return mu, cov, residual, loading, weights, math.exp(-r * maturity)


def summarise_batches(values, batches):
    a = np.asarray(values).reshape(batches, -1)
    means = a.mean(axis=1)
    return {"mean": float(means.mean()),
            "sd": float(np.std(a.ravel(), ddof=1)),
            "batch_mean_se": float(np.std(means, ddof=1) / math.sqrt(batches))}


def scenario(d, correlation, strike, cap, paths, batches, seed):
    mu, cov, residual, loading, w, discount = model(d, correlation)
    rng = np.random.default_rng(seed)
    n = paths * batches
    y = rng.standard_normal((n, d - 1))
    z = rng.standard_normal(n)
    u = rng.random(n)
    logs = mu + y @ residual.T
    log_a = np.log(np.exp(logs) @ w)
    log_g = logs @ w
    assert np.min(log_a - log_g) >= -1e-12
    a, g = np.exp(log_a + loading * z), np.exp(log_g + loading * z)
    direct_a, direct_g = payoff(a, strike, cap), payoff(g, strike, cap)
    full_a = common_conditional(log_a, loading, strike, cap)
    full_g = common_conditional(log_g, loading, strike, cap)
    k = strike + cap * u
    root_a, root_g = (np.log(k) - log_a) / loading, (np.log(k) - log_g) / loading
    random_a = cap * cdf(-root_a)
    random_residual = cap * interval_probability(root_a, root_g)
    geo_mu = float(w @ mu)
    geo_sd = math.sqrt(float(w @ cov @ w))
    geo_mean = float(call(np.array(geo_mu), geo_sd, strike)
                     - call(np.array(geo_mu), geo_sd, strike + cap))
    outputs = {
        "direct": direct_a,
        "direct_geometric_residual": direct_a - direct_g,
        "full_conditional": full_a,
        "full_conditional_geometric_residual": full_a - full_g,
        "random_strike_conditional": random_a,
        "random_strike_conditional_geometric_residual": random_residual,
    }
    statistics = {name: summarise_batches(value, batches) for name, value in outputs.items()}
    for name, stat in statistics.items():
        correction = geo_mean if "residual" in name else 0.0
        stat["price_estimate"] = discount * (stat["mean"] + correction)
        stat["price_batch_se"] = discount * stat["batch_mean_se"]
    ratios = {}
    target_sd = statistics["random_strike_conditional_geometric_residual"]["sd"]
    for baseline in ["direct_geometric_residual", "full_conditional_geometric_residual"]:
        sd = statistics[baseline]["sd"]
        ratios[baseline] = {
            "candidate_sd_ratio": target_sd / sd,
            "maximum_oracle_overhead_for_leading_variance_sensitive_break_even": sd / target_sd,
        }
    # An analytic conservative variance bound, not the observed sample variance.
    mean_a = float(w @ np.exp(mu + np.diag(cov) / 2))
    mean_g = math.exp(geo_mu + geo_sd**2 / 2)
    residual_sd_bound = math.sqrt(min(cap**2 / 4, cap * max(mean_a - mean_g, 0)))
    # Independent within-input checks of full conditioning and strike integration.
    check_idx = np.linspace(0, n - 1, 16, dtype=int)
    refined = common_conditional(log_a[check_idx], loading, strike, cap, order=96)
    refinement_error = float(np.max(np.abs(refined - full_a[check_idx])))
    grid_n = 4096
    grid = strike + cap * (np.arange(grid_n) + 0.5) / grid_n
    integrated = cap * cdf((log_a[check_idx, None] - np.log(grid)) / loading).mean(axis=1)
    grid_error = float(np.max(np.abs(integrated - refined)))
    grid_bound = cap / (2 * grid_n)
    assert grid_error <= grid_bound + refinement_error + 1e-10
    return {"dimension": d, "correlation": correlation, "strike": strike, "cap": cap,
            "paths_per_batch": paths, "independent_batches": batches, "seed": seed,
            "common_loading": loading, "statistics": statistics, "break_even_screen": ratios,
            "analytic_residual_sd_upper_bound": residual_sd_bound,
            "checks": {"quadrature_order_48_96_max_difference": refinement_error,
                       "midpoint_strike_grid_max_difference": grid_error,
                       "midpoint_strike_grid_bound": grid_bound}}


def scalar_stress():
    results = []
    nodes, weights = np.polynomial.legendre.leggauss(128)
    for loading in [0.0001, 0.01, 0.1, 0.8]:
        # Split around the narrow normal transition so the stress case is
        # resolved even for loading=1e-4; a global polynomial grid misses it.
        cuts = sorted(set([100.0, 110.0] +
                          [float(np.clip(105 * math.exp(sign * 10 * loading), 100, 110))
                           for sign in [-1, 1]]))
        mean, second = 0.0, 0.0
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            strikes = (lo + hi) / 2 + (hi - lo) * nodes / 2
            h = 10 * cdf((math.log(105) - np.log(strikes)) / loading)
            mean += float(weights @ h) * (hi - lo) / 20
            second += float(weights @ (h * h)) * (hi - lo) / 20
        var = max(0.0, second - mean * mean)
        zcuts = sorted(set([-12.0, 12.0] +
                          [float(np.clip(math.log(k / 105) / loading, -12, 12))
                           for k in [100, 110]]))
        direct_mean, direct_second = 0.0, 0.0
        for lo, hi in zip(zcuts[:-1], zcuts[1:]):
            z = (lo + hi) / 2 + (hi - lo) * nodes / 2
            f = payoff(105 * np.exp(loading * z), 100, 10)
            density = np.exp(-z * z / 2) / math.sqrt(2 * math.pi)
            direct_mean += float(weights @ (f * density)) * (hi - lo) / 2
            direct_second += float(weights @ (f * f * density)) * (hi - lo) / 2
        direct_var = max(0.0, direct_second - direct_mean * direct_mean)
        results.append({"conditional_basket_base": 105, "strike": 100, "cap": 10,
                        "loading": loading, "random_strike_variance": var,
                        "direct_payoff_variance": direct_var,
                        "mean_identity_difference": abs(mean - direct_mean)})
    return results


def state_vector_check():
    """Exact small comparison-register oracle, using host-evaluated probabilities.

    The truth table is deliberately labelled. It is NOT a scalable arithmetic
    circuit. It verifies Born probabilities and Grover interference only.
    """
    y = np.array([NORMAL.inv_cdf((j + 0.5) / 4) for j in range(4)])
    u = (np.arange(4) + 0.5) / 4
    mu, _, residual, loading, w, _ = model(2, 0.6)
    logs = mu + y[:, None] * residual[:, 0]
    log_a, log_g = np.log(np.exp(logs) @ w), logs @ w
    strikes = 100 + 10 * u
    za = (np.log(strikes)[None, :] - log_a[:, None]) / loading
    zg = (np.log(strikes)[None, :] - log_g[:, None]) / loading
    probabilities = interval_probability(za, zg)
    comparison_levels = 256
    thresholds = np.floor(comparison_levels * probabilities).astype(int)
    inputs = 4 * 4 * comparison_levels
    psi = np.zeros(inputs * 2, dtype=complex)
    for iy in range(4):
        for iu in range(4):
            for ir in range(comparison_levels):
                inp = (iy * 4 + iu) * comparison_levels + ir
                flag = int(ir < thresholds[iy, iu])
                psi[2 * inp + flag] = 1 / math.sqrt(inputs)
    finite_a = float(thresholds.mean() / comparison_levels)
    probability = float(np.sum(np.abs(psi[1::2])**2))
    assert abs(probability - finite_a) < 1e-14
    assert abs(np.vdot(psi, psi).real - 1) < 1e-14
    theta = math.asin(math.sqrt(finite_a))
    state = psi.copy()
    errors = []
    for iterations in range(8):
        actual = float(np.sum(np.abs(state[1::2])**2))
        predicted = math.sin((2 * iterations + 1) * theta)**2
        errors.append(abs(actual - predicted))
        phase = state.copy()
        phase[1::2] *= -1
        state = 2 * psi * np.vdot(psi, phase) - phase
    assert max(errors) < 1e-12
    assert 0 < probability < 1
    return {"state_vector_qubits": int(math.log2(inputs * 2)), "host_evaluated_probability_truth_table": True,
            "finite_output_success_probability": probability,
            "pre_quantisation_success_probability": float(probabilities.mean()),
            "comparison_quantisation_bound": 1 / comparison_levels,
            "maximum_grover_probability_error": max(errors),
            "interpretation": "Quantum-mechanical correctness only; no gate cost or speedup claim."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paths", type=int, default=4096)
    parser.add_argument("--batches", type=int, default=8)
    args = parser.parse_args()
    start = time.perf_counter()
    # Fixed before seeing results; twelve cases include narrow and wider caps.
    cases = [(4, rho, strike, cap) for rho in [0.1, 0.6, 0.9]
             for strike, cap in [(100, 1), (100, 10), (100, 25), (115, 10)]]
    scenarios = [scenario(*case, args.paths, args.batches, 20261005 + i)
                 for i, case in enumerate(cases)]
    result = {"status": "classical feasibility screen; not compiled quantum resources",
              "data": "fully stipulated synthetic GBM; no LSEG observations",
              "confidence": "batch standard errors are diagnostics, not certified pricing intervals or variance bounds",
              "numpy_version": np.__version__,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scenarios": scenarios, "scalar_stress": scalar_stress(),
              "quantum_correctness": state_vector_check(),
              "elapsed_cpu_workflow_seconds": time.perf_counter() - start}
    output = ROOT / "results" / "pilot.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(output), "scenarios": len(scenarios),
                      "quantum_correctness": result["quantum_correctness"],
                      "screen": [{"rho": s["correlation"], "strike": s["strike"], "cap": s["cap"],
                                  "maximum_candidate_oracle_overhead_vs_direct_control":
                                  s["break_even_screen"]["direct_geometric_residual"]
                                  ["maximum_oracle_overhead_for_leading_variance_sensitive_break_even"]}
                                 for s in scenarios]}, indent=2))


if __name__ == "__main__":
    main()
