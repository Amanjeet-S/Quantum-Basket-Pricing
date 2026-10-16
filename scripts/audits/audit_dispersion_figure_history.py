"""Independent exact data/coordinate and preserved-history audit.

Only the coordinate check uses transcendental arithmetic. Model/payoff
selection, allowances, flags, bars, historical hashes and admission counts
are reconstructed from exact fractions and the full retained results.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json
from flint import arb, ctx
import independent_dispersion_math as im

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path, default=ROOT.parent)
    ap.add_argument("--out", type=Path, default=ROOT / "dispersion_figure_history_audit.json")
    args = ap.parse_args()
    stage = args.stage
    ctx.prec = 1024
    fullpath = stage / "results/dispersion_separation_results.json"
    full = json.loads(fullpath.read_text())
    summarypath = stage / "results/dispersion_separation_summary.json"
    summary = json.loads(summarypath.read_text())
    protocolpath = stage / "results/dispersion_separation_protocol.json"
    selectionpath = stage / "figure1_protocol.json"
    selection = json.loads(selectionpath.read_text())
    figdir = stage / "figures"
    datapath = figdir / "dispersion_separation_plot_data.json"
    data = json.loads(datapath.read_text())
    coords = json.loads((figdir / "dispersion_separation_coordinates.json").read_text())
    assert summary["full_result_sha256"] == digest(fullpath)
    assert data["source_full_result_sha256"] == digest(fullpath)
    assert data["source_summary_sha256"] == digest(summarypath)
    assert data["experiment_protocol_sha256"] == digest(protocolpath)
    assert data["selection_sha256"] == digest(selectionpath)
    assert selection["model_family"] == data["representative_model"] == "n8_lambda0.2"
    assert selection["payoff"] == "vanilla" and data["representative_payoff"] == "ordinary call"
    assert data["scale_values"] == selection["scales"] == ["1", "3/2", "2", "3"]
    assert data["threshold_exact_fraction"] == "1"
    assert F(data["transfer_allocation_exact_fraction"]) == F(1, 200)
    assert F(data["price_allocation_exact_fraction"]) == F(1, 100)
    counters = {"historical_file_hashes": 0, "historical_input_result_links": 0,
                "historical_immutable_degree3_snapshots": 0, "unchanged_historical_raw_radii": 0,
                "representative_exact_curve_points": 0, "full_panel_exact_bar_counts": 0,
                "native_coordinate_checks": 0, "native_mark_checks": 0}
    histories = []
    finals = {(r["model_id"], r["cap"]): r for r in full["records"]}
    for name in ("registered_run1_before_intersection", "registered_run2_before_exact_interval_storage"):
        directory = stage / "history" / name
        manifest = json.loads((directory / "snapshot_manifest.json").read_text())
        for file, expected in manifest.items():
            assert digest(directory / file) == expected
            counters["historical_file_hashes"] += 1
        prior = json.loads((directory / "dispersion_separation_results.json").read_text())
        assert prior["source_sha256"]["dispersion_separation_study.py"] == digest(directory / "dispersion_separation_study.py")
        assert prior["protocol_sha256"] == digest(directory / "dispersion_separation_protocol.json")
        counters["historical_input_result_links"] += 2
        for row in prior["records"]:
            assert [o["order"] for o in row["degree3_continuous_price"]["orders"]] == [1, 2, 3]
            assert len(row["degree3_continuous_price"]["phase_work"]) == 1
            counters["historical_immutable_degree3_snapshots"] += 1
            final = finals[(row["model_id"], row["cap"])]
            assert F(row["aggregate_R2_transport"]["financial_transfer_upper"]) == F(final["aggregate_R2_transport"]["financial_transfer_upper"])
            for old, new in zip(row["degree3_continuous_price"]["orders"], final["degree3_continuous_price"]["orders"]):
                assert F(old["undiscounted_price_radius_upper"]) == F(new["undiscounted_price_radius_upper"])
                counters["unchanged_historical_raw_radii"] += 1
        histories.append({"history": name, "result_sha256": digest(directory / "dispersion_separation_results.json"),
                          "protocol_sha256": prior["protocol_sha256"],
                          "source_sha256": prior["source_sha256"]["dispersion_separation_study.py"]})
    # The initial figure family is retained in the original registered
    # protocol. Later comparator strengthening is explicitly post-outcome.
    firstprotocol = json.loads((stage / "history/registered_run1_before_intersection/dispersion_separation_protocol.json").read_text())
    registered = firstprotocol["figure_fixed_selection"]
    assert registered["model_family"] == "n8"
    assert F(registered["common_lambda"]) == F(1, 5)
    assert registered["payoff"] == "vanilla"
    assert registered["scales"] == selection["scales"]
    assert selection["model_selection_fixed_before_outcomes"] is True
    assert selection["fixed_before_outcomes"] is False
    assert selection["no_selection_from_observed_separators"] is True

    pointlookup = {}
    for curve in data["curves"]:
        source = next(r for r in full["records"] if r["original_model_id"] == "n8_lambda0.2"
                      and r["cap"] is None and r["dispersion_scale"] == curve["dispersion_scale"])
        ratio_t = F(source["aggregate_R2_transport"]["financial_transfer_upper"]) / F(1, 200)
        ratio_p = F(source["degree3_continuous_price"]["intersection_radius_upper"]) / F(1, 100)
        assert ratio_t == F(curve["transfer_ratio_exact_fraction"])
        assert ratio_p == F(curve["classical_ratio_exact_fraction"])
        assert curve["transfer_pass"] == (ratio_t < 1)
        assert curve["classical_price_pass"] == (ratio_p <= 1)
        assert curve["model_id"] == source["model_id"]
        assert F(str(curve["x_float"])) == F(curve["dispersion_scale"])
        pointlookup[(curve["model_id"], "transfer_ratio_exact_fraction")] = (F(curve["dispersion_scale"]), ratio_t)
        pointlookup[(curve["model_id"], "classical_ratio_exact_fraction")] = (F(curve["dispersion_scale"]), ratio_p)
        counters["representative_exact_curve_points"] += 2
    assert len(pointlookup) == 8
    countlookup = {}
    for row in data["counts"]:
        cohort = [r for r in full["records"] if r["dispersion_scale"] == row["dispersion_scale"]]
        assert len(cohort) == row["contracts"] == 36
        provisional = sum(F(r["aggregate_R2_transport"]["financial_transfer_upper"]) < F(1, 200)
                          and not r["degree3_continuous_price"]["intersection_certifies_cent"] for r in cohort)
        survivors = sum(r["separates_after_degree5"] for r in cohort)
        assert row["through_degree3"] == provisional
        assert row["after_targeted_degree5"] == survivors <= provisional
        countlookup[(row["dispersion_scale"], "through_degree3")] = provisional
        countlookup[(row["dispersion_scale"], "after_targeted_degree5")] = survivors
        counters["full_panel_exact_bar_counts"] += 2
    assert len(countlookup) == 8

    mapping = {k: F(str(v)) for k, v in coords["mapping"].items()}
    native = (figdir / "dispersion_separation_figure.tex").read_text()
    assert "not observed errors" in native and "complete quantum prices" in native
    assert len(coords["points"]) == 16
    for point in coords["points"]:
        if point["panel"] == "a":
            c, ratio = pointlookup[(point["model_id"], point["series"])]
            assert F(point["exact_fraction"]) == ratio
            x = mapping["left"] + mapping["width"] * (c - mapping["xmin"]) / (mapping["xmax"] - mapping["xmin"])
            logratio = im.ball(ratio).log() / arb(10).log()
            y = im.ball(mapping["bottom"]) + im.ball(mapping["height"] / (mapping["log10_ymax"] - mapping["log10_ymin"])) * (logratio - im.ball(mapping["log10_ymin"]))
            assert abs(F(str(point["x_mm"])) - x) <= F(1, 10 ** 9)
            assert abs(im.ball(F(str(point["y_mm"]))) - y).upper() < im.ball(F(1, 10 ** 9))
            if point["series"] == "transfer_ratio_exact_fraction":
                mark = rf'\put({point["x_mm"]:.4f},{point["y_mm"]:.4f}){{\circle*{{1.2}}}}'
            else:
                mark = rf'\put({point["x_mm"]-.55:.4f},{point["y_mm"]-.55:.4f}){{\rule{{1.1mm}}{{1.1mm}}}}'
            assert mark in native
            counters["native_coordinate_checks"] += 2
            counters["native_mark_checks"] += 1
        else:
            count = countlookup[(point["dispersion_scale"], point["series"])]
            assert point["count"] == count
            height = mapping["height"] * count / mapping["count_ymax"]
            assert abs(F(str(point["height_mm"])) - height) <= F(1, 10 ** 9)
            assert F(str(point["y_mm"])) == mapping["bottom"]
            i = ["1", "3/2", "2", "3"].index(point["dispersion_scale"])
            offset = F(-19, 5) if point["series"] == "through_degree3" else F(3, 10)
            x = mapping["right"] + mapping["bar_width"] * F(2 * i + 1, 8) + offset
            assert abs(F(str(point["x_mm"])) - x) <= F(1, 10 ** 9)
            if count:
                mark = rf'\put({point["x_mm"]:.4f},{point["y_mm"]:.4f}){{\rule{{3.5mm}}{{{point["height_mm"]:.4f}mm}}}}'
                assert mark in native
                counters["native_mark_checks"] += 1
            counters["native_coordinate_checks"] += 3
    out = {"schema": "qbp-dispersion-figure-history-independent-audit-v1", "status": "passed",
           "audit_counts": counters, "history_records": histories,
           "representative_point_count": 8, "full_panel_bar_count": 8,
           "source_sha256": digest(__file__), "result_sha256": digest(fullpath),
           "plot_data_sha256": digest(datapath), "selection_sha256": digest(selectionpath),
           "native_figure_sha256": digest(figdir / "dispersion_separation_figure.tex"),
           "coordinate_sha256": digest(figdir / "dispersion_separation_coordinates.json"),
           "scope": "Exact allowance/count selection and native mark geometry. Neither certificate is an observed price error. Original family is pre-outcome; strengthened opponent and aggregate presentation are transparently subsequent."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "audit_counts": counters}, indent=2))


if __name__ == "__main__":
    main()
