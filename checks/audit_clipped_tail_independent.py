#!/usr/bin/env python3
"""Independent exact-rational clipped-tail audit; standard library only.

This retains the exact parameter settings and checks executed in the original
481,241-configuration audit. It does not import clipping_certificate.py or any
project certificate implementation. Square-root regimes are explicitly skipped
by the rational enumeration, rather than silently checked with floating point.
The analytic constant-law proof is retained in clip_review.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path


F = Fraction
PARAMETER_SETTINGS = {
    "rho": {"numerators": [1, 16], "denominator": 4},
    "r": {"numerators": [0, 16], "denominator": 8},
    "u": {"numerators": [0, 32], "denominator": 8},
    "v": {"numerators": [0, 64], "denominator": 8},
    "majorant_grid": {"numerators": [0, 16], "denominator": 16,
                       "definition": "x=rho*k/16"},
}
EXPECTED_COUNTS = {
    "zero": 157456,
    "r0-mean": 22138,
    "second-only": 92128,
    "both-active": 94831,
    "endpoint": 114688,
}


def audit() -> dict:
    counts: dict[str, int] = {}
    skipped: dict[str, int] = {}
    checked = 0
    total = 0
    comparisons = 0
    for rho, r, u, v in product(
        [F(k, 4) for k in range(1, 17)],
        [F(k, 8) for k in range(0, 17)],
        [F(k, 8) for k in range(0, 33)],
        [F(k, 8) for k in range(0, 65)],
    ):
        total += 1
        up = min(u, rho)
        vp = min(v, rho * up, rho * rho)
        if rho <= r or up == 0 or vp == 0:
            val = F(0)
            case = "zero"
            pts = [(F(0), F(1))]
            coef = None
        elif r == 0:
            if up * up <= vp:
                val = up
                pts = [(up, F(1))]
                coef = (F(0), F(1), F(0))
                case = "r0-mean"
            else:
                skipped["r0-sqrt"] = skipped.get("r0-sqrt", 0) + 1
                continue
        elif rho < 2 * r:
            val = (rho - r) * vp / (rho * rho)
            mass = vp / (rho * rho)
            pts = [(F(0), 1 - mass), (rho, mass)]
            coef = (F(0), F(0), (rho - r) / (rho * rho))
            case = "endpoint"
        elif vp <= 4 * r * r and up >= vp / (2 * r):
            val = vp / (4 * r)
            mass = vp / (4 * r * r)
            pts = [(F(0), 1 - mass), (2 * r, mass)]
            coef = (F(0), F(0), 1 / (4 * r))
            case = "second-only"
        elif vp > 4 * r * r and up * up >= vp:
            skipped["constant-sqrt"] = skipped.get("constant-sqrt", 0) + 1
            continue
        else:
            val = up - r * up * up / vp
            mass = up * up / vp
            point = vp / up
            pts = [(F(0), 1 - mass), (point, mass)]
            coef = (F(0), 1 - 2 * r * up / vp,
                    r * up * up / (vp * vp))
            case = "both-active"

        assert sum(p for x, p in pts) == 1
        assert all(F(0) <= x <= rho and p >= 0 for x, p in pts)
        assert sum(x * p for x, p in pts) <= up
        assert sum(x * x * p for x, p in pts) <= vp
        assert sum(max(x - r, F(0)) * p for x, p in pts) == val
        if coef:
            c0, c1, c2 = coef
            assert c1 >= 0 and c2 >= 0
            for x in [rho * F(k, 16) for k in range(17)]:
                assert c0 + c1 * x + c2 * x * x >= max(x - r, F(0))
                comparisons += 1
            assert c0 + c1 * up + c2 * vp == val
        checked += 1
        counts[case] = counts.get(case, 0) + 1

    assert counts == EXPECTED_COUNTS
    assert checked == 481241
    assert total == 583440
    assert checked + sum(skipped.values()) == total
    return {
        "schema": "qbp-independent-clipped-tail-rational-audit-v1",
        "status": "passed",
        "arithmetic": "exact fractions.Fraction, no binary floating point",
        "imports_project_certificate": False,
        "parameter_settings": PARAMETER_SETTINGS,
        "parameter_configurations_total": total,
        "exact_rational_regimes_checked": checked,
        "counts": counts,
        "skipped_counts": skipped,
        "majorant_point_comparisons": comparisons,
        "checks": [
            "attaining masses sum to one",
            "attaining masses nonnegative and support lies in [0,rho]",
            "attaining mean is within effective upper moment bound",
            "attaining second moment is within effective upper moment bound",
            "attaining clipped-tail expectation equals proposed bound exactly",
            "coefficients for upper moment bounds are nonnegative",
            "quadratic majorant dominates payoff at 17 rational grid points",
            "dual moment objective equals proposed bound exactly",
            "configuration and branch counts match the original executed audit",
        ],
        "limitations": [
            "Constant square-root branches are proved analytically in the memo "
            "and are explicitly excluded from this rational enumeration.",
            "A finite set of majorant evaluations is diagnostic; the analytic "
            "pointwise proof supplies the universal majorant statement.",
            "Attainment is for arbitrary probability laws on the continuous "
            "support interval, not necessarily a prescribed finite grid.",
            "This mathematical implementation audit is not an empirical "
            "market-data validation or a theoretical novelty claim.",
        ],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default_output = (Path(__file__).resolve().parents[1] / "results" /
                      "clipped_tail_independent_audit.json")
    parser.add_argument("--output", type=Path, default=default_output)
    args = parser.parse_args()
    result = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"],
                      "checked": result["exact_rational_regimes_checked"],
                      "counts": result["counts"],
                      "skipped_counts": result["skipped_counts"],
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
