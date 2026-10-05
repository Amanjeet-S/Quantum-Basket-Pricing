"""Independent numerical diagnostics of the retained monetary certificates."""
from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np

from pilot import cdf

ROOT = Path(__file__).resolve().parent


def root(coefficients, loadings, strike):
    lo, hi = -40.0, 40.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if np.sum(coefficients * np.exp(loadings * mid)) < strike:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def positive_reference(coefficients, loadings, strike, cap):
    lo = root(coefficients, loadings, strike)
    hi = root(coefficients, loadings, strike + cap)
    nodes, weights = np.polynomial.legendre.leggauss(96)
    z = (lo + hi) / 2 + (hi - lo) * nodes / 2
    f = np.clip(np.exp(z[:, None] * loadings) @ coefficients - strike, 0, cap)
    density = np.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    return float((hi - lo) * (weights @ (f * density)) / 2 + cap * cdf(-hi))


def main():
    root_records = []
    for coefficients, loadings in [([40, 40], [.1, .3]),
                                   ([50, 50], [.2, .4]),
                                   ([60, 60], [.15, .25])]:
        c, b = np.array(coefficients), np.array(loadings)
        reference = positive_reference(c, b, 100, 10)
        for n in [8, 32, 128]:
            strikes = 100 + 10 * (np.arange(n) + .5) / n
            for rho in [0, .001, .01, .1, 1]:
                targets = strikes + rho * np.where(np.arange(n) % 2, 1, -1)
                approximate_roots = np.array([root(c, b, k) for k in targets])
                actual_residual = float(np.max(np.abs(
                    np.exp(approximate_roots[:, None] * b) @ c - strikes)))
                mean = float(10 * cdf(-approximate_roots).mean())
                error, bound = abs(mean - reference), rho + 10 / (2 * n)
                assert actual_residual <= rho + 1e-10
                assert error <= bound + 1e-10
                root_records.append({"coefficients": coefficients, "loadings": loadings,
                                     "midpoint_strikes": n, "residual_allowance": rho,
                                     "observed_residual": actual_residual,
                                     "observed_mean_error": error, "bound": bound})
    normal = NormalDist()
    nodes, weights = np.polynomial.legendre.leggauss(96)
    quantile_records = []
    for n in [1, 2, 4, 16, 256]:
        mids = np.array([normal.inv_cdf((j + .5) / n) for j in range(n)])
        inner = np.array([normal.inv_cdf(j / n) for j in range(1, n)])
        phi = lambda x: np.exp(-x * x / 2) / math.sqrt(2 * math.pi)
        formula = float(2 * (phi(mids).sum() - phi(inner).sum()))
        edges = np.concatenate(([-12.0], inner, [12.0]))
        integral = 0.0
        for j, mid in enumerate(mids):
            for lo, hi in [(edges[j], mid), (mid, edges[j + 1])]:
                z = (lo + hi) / 2 + (hi - lo) * nodes / 2
                integral += float(weights @ (np.abs(z - mid) * phi(z))) * (hi - lo) / 2
        assert abs(formula - integral) < 1e-11
        quantile_records.append({"midpoint_quantiles": n, "formula": formula,
                                 "independent_integral": integral,
                                 "difference": abs(formula - integral)})
    output = {"interpretation": "Numerical diagnostics; proofs are in mathematical_note.md. No rounding enclosure or gate advantage is certified.",
              "root_grid_checks": root_records, "quantile_coupling_checks": quantile_records}
    (ROOT / "results" / "certificate_checks.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"root_grid_checks": len(root_records), "quantile_checks": len(quantile_records),
                      "maximum_quantile_integral_difference": max(x["difference"] for x in quantile_records)}))


if __name__ == "__main__":
    main()
