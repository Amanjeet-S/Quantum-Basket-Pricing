"""Narrow exact-rational check of monetary admission after serialization.

This deliberately does not replay either full mathematical audit.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
DEFAULT_STAGE = PROJECT / "cost_boundary" if (PROJECT / "cost_boundary").is_dir() else PROJECT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path, default=DEFAULT_STAGE)
    ap.add_argument("--out", type=Path, default=ROOT / "saved_price_endpoint_audit.json")
    args = ap.parse_args()
    records = []
    for name in ("residual_same_information_taylor.json", "caldana_conditioning_benchmark.json"):
        path = args.stage / "results" / name
        source = json.loads(path.read_text())
        rows = []
        distinct_saved = set()
        for record in source["records"]:
            for order in record.get("orders", [record]):
                saved = order["continuous_price_interval"]
                discounted = order["discounted_price_interval"]
                radius = (F(saved["upper_exact_fraction"]) - F(saved["lower_exact_fraction"])) / 2
                discounted_radius = (F(discounted["upper_exact_fraction"]) - F(discounted["lower_exact_fraction"])) / 2
                admitted = radius <= F(1, 100)
                discounted_admitted = discounted_radius <= F(1, 100)
                assert admitted == order["certifies_undiscounted_cent"]
                assert admitted == discounted_admitted
                if admitted:
                    distinct_saved.add((record["model_id"], record["cap"]))
                rows.append({"model_id": record["model_id"], "cap": record["cap"],
                             "order": order.get("order"),
                             "saved_undiscounted_halfwidth": str(radius),
                             "saved_discounted_halfwidth": str(discounted_radius),
                             "saved_undiscounted_admitted": admitted,
                             "saved_discounted_admitted": discounted_admitted,
                             "recorded_admitted": order["certifies_undiscounted_cent"]})
        records.append({"file": name, "result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "order_rows": len(rows), "passing_order_rows": sum(r["recorded_admitted"] for r in rows),
                        "passing_distinct_contracts": len(distinct_saved), "rows": rows})
    assert records[0]["passing_distinct_contracts"] == 33
    assert records[1]["passing_distinct_contracts"] == 0
    out = {"schema": "qbp-saved-price-endpoint-admission-v1", "status": "passed",
           "method": "Exact Fraction half-widths of saved interval endpoints, at threshold 1/100; no floating-point threshold or full producer replay",
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "records": records,
           "scope": "No classification mismatch in this cohort. This does not equate the pre-serialization radius field with a saved Arb interval half-width for arbitrary inputs."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "results": [
        {k: r[k] for k in ("file", "order_rows", "passing_order_rows", "passing_distinct_contracts")}
        for r in records]}, indent=2))


if __name__ == "__main__":
    main()
