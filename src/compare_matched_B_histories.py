"""Retain the executed paired B/clean diagnostic comparison and scope checks."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
sys.dont_write_bytecode = True


def read(path):
    path = Path(path)
    return json.load(gzip.open(path, "rt")) if path.suffix == ".gz" else json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def without_resource_counts(value):
    if isinstance(value, dict):
        return {k: without_resource_counts(v) for k, v in value.items()
                if k not in ("complete_batch_T", "cumulative_pilot_T", "pilot_T_count")}
    if isinstance(value, list):
        return [without_resource_counts(v) for v in value]
    return value


def main(args):
    repo, out = Path(args.repo), Path(args.output)
    data = read(out / "matched_B/matched_chebyshev_resources.json.gz")
    protocol = read(out / "matched_B/matched_chebyshev_protocol.json")
    old = read(repo / "results/matched_chebyshev/amended/matched_chebyshev_resources.json.gz")
    inputs = read(out / "matched_B_inputs.json")
    deferred = read(out / "deferred_cleanup_resources.json.gz")
    assert data["protocol_sha256"] == sha(out / "matched_B/matched_chebyshev_protocol.json")
    index = {(r["record_index"], r["architecture"]): r for r in deferred["records"]}
    cases, histories, looks, profiles = [], 0, 0, 0
    same_degrees = same_stops = same_numeric_controller = 0
    for ci, (fresh, previous) in enumerate(zip(data["cases"], old["cases"])):
        assert fresh["model_id"] == previous["model_id"] == inputs["cases"][ci]["model_id"]
        assert fresh["common_signal_interval"] == previous["common_signal_interval"]
        methods = {}
        for method in protocol["methods"]:
            new, earlier = fresh["methods"][method], previous["methods"][method]
            info = inputs["cases"][ci]["methods"][method]
            record = index[info["matched_variants_proof"]["deferred_record_index"], info["architecture"]]
            half = record["deferred_B_half_charge"]
            fallback = new["half_alpha_QPE_fallback"]
            for field in ("metrics", "per_repetition_metrics", "peak_logical_qubits", "complete_oracle_T", "controlled_Q_T", "controlled_reflection", "repetitions", "controlled_Q_calls", "phase_register_measurements", "complete_reinitialisations"):
                assert fallback["complete_charge"][field] == half[field]
            assert new["full_alpha_canonical_T"] == record["deferred_B"]["T_full_alpha"]
            assert new["whole_history_T_cap"] == record["deferred_B"]["T_worst_history_upper"]
            for sample in new["sampler_profiles"].values():
                assert sample["preparation_symbol"] == "B"
                assert sample["preparation_conjugacy_proof_sha256"] == deferred["checked_proof_sha256"]
                assert "dirty" in sample["financial_cleanup_scope"]
                assert "already includes freshly emitted Boolean-body cleanup" not in sample["financial_cleanup_scope"]
                profiles += 1
            for h, h_old in zip(new["diagnostics"], earlier["diagnostics"]):
                for key in ("mesh_index", "mesh_fraction", "truth_probability", "seed_index", "seed"):
                    assert h[key] == h_old[key]
                same_stops += h["controller_result"]["status"] == h_old["controller_result"]["status"]
                same_numeric_controller += without_resource_counts(h["controller_result"]) == without_resource_counts(h_old["controller_result"])
                degrees = [b["degree"] for b in h["raw_batches"]]
                old_degrees = [b["degree"] for b in h_old["raw_batches"]]
                same_degrees += degrees == old_degrees
                if degrees == old_degrees:
                    for batch, b_old in zip(h["raw_batches"], h_old["raw_batches"]):
                        for key in ("look", "degree", "shots", "successes", "paired_seeded_words_sha256", "signal"):
                            assert batch[key] == b_old[key]
                assert h["actual_history_T_count"] <= new["whole_history_T_cap"]
                histories += 1
                looks += len(h["raw_batches"])
            summary = new["diagnostic_summary"]
            old_summary = earlier["diagnostic_summary"]
            methods[method] = {"full_alpha_canonical_T": new["full_alpha_canonical_T"],
                "half_alpha_fallback_T": fallback["T_count"], "whole_history_T_cap": new["whole_history_T_cap"],
                "pilot_certified": summary["pilot_certified"], "previous_clean_pilot_certified": old_summary["pilot_certified"],
                "common_mesh_mean_complete_T_exact": summary["observed_history_T_mean"],
                "common_mesh_median_complete_T_exact": summary["observed_history_T_median"],
                "previous_clean_common_mesh_mean_T_exact": old_summary["observed_history_T_mean"],
                "previous_clean_common_mesh_median_T_exact": old_summary["observed_history_T_median"],
                "coverage_diagnostic_failures": summary["coverage_diagnostic_failures"]}
        cases.append({"model_id": fresh["model_id"], "methods": methods,
            "common_mesh_mean_winner": min(methods, key=lambda m: (F(methods[m]["common_mesh_mean_complete_T_exact"]), m)),
            "common_mesh_median_winner": min(methods, key=lambda m: (F(methods[m]["common_mesh_median_complete_T_exact"]), m)),
            "complete_worst_history_cap_winner": min(methods, key=lambda m: (methods[m]["whole_history_T_cap"], m))})
    assert histories == 1200
    result = {"schema": "qbp-matched-deferred-B-history-comparison-v1", "status": "passed trace/resource binding comparison",
              "source_sha256": sha(__file__), "matched_B_resources_sha256": sha(out / "matched_B/matched_chebyshev_resources.json.gz"),
              "matched_B_protocol_sha256": sha(out / "matched_B/matched_chebyshev_protocol.json"),
              "matched_B_inputs_sha256": sha(out / "matched_B_inputs.json"),
              "previous_clean_amended_resources_sha256": sha(repo / "results/matched_chebyshev/amended/matched_chebyshev_resources.json.gz"),
              "deferred_resources_sha256": sha(out / "deferred_cleanup_resources.json.gz"),
              "histories": histories, "paired_histories": 400, "fresh_looks": looks,
              "histories_with_unchanged_degree_sequences": same_degrees,
              "histories_with_unchanged_stopping_status": same_stops,
              "histories_with_unchanged_nonresource_controller_records": same_numeric_controller,
              "sampler_profiles_with_corrected_B_metadata": profiles, "cases": cases,
              "interpretation": ("The clean/B paired histories were executed and compared rather than inferred to be unchanged. "
                  "All financial information and deterministic word batches agree. The numerical conditional rival retains "
                  "the lowest descriptive common-mesh mean in all four cases; shared-threshold retains the lowest median "
                  "and guaranteed complete worst-history cap. Mesh means and medians are not expected resources at an "
                  "actual financial flag, and no adaptive estimator dominance or new query law is established."),
              "scope": "Trace, metadata and registered cost binding comparison. Independent confidence/circuit formula verification is a separate audit."}
    target = out / "matched_B_history_comparison.json"
    if target.exists():
        raise FileExistsError("Preserve executed comparison")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "histories", "fresh_looks", "histories_with_unchanged_degree_sequences", "histories_with_unchanged_stopping_status", "histories_with_unchanged_nonresource_controller_records", "sampler_profiles_with_corrected_B_metadata")}))
    for case in cases:
        print(case["model_id"], case["common_mesh_mean_winner"], case["common_mesh_median_winner"],
              {m: r["common_mesh_mean_complete_T_exact"] for m, r in case["methods"].items()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--output", required=True)
    main(parser.parse_args())
