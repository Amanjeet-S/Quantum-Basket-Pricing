"""Extract the registered B canonical winners for a matched diagnostic panel."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil
import sys
sys.dont_write_bytecode = True


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def read(path):
    path = Path(path)
    return json.load(gzip.open(path, "rt")) if path.suffix == ".gz" else json.loads(path.read_text())


def main(args):
    repo, out, proof = Path(args.repo), Path(args.output), Path(args.proof)
    resources = read(out / "deferred_cleanup_resources.json.gz")
    protocol = read(out / "deferred_cleanup_protocol.json")
    assert resources["protocol_sha256"] == sha(out / "deferred_cleanup_protocol.json")
    assert sha(proof) == protocol["checked_proof_sha256"]
    rows = {(r["record_index"], r["architecture"]): r for r in resources["records"]}
    methods = (("shared_threshold", "shared_threshold"), ("raw_payoff", "raw_payoff"),
               ("numerical_conditional", "numerical_conditional_price"))
    cases = []
    disclosure = ("All previous clean canonical and matched Chebyshev diagnostics, the cleanup-deferral identities, "
                  "and every B/optimised-clean canonical and half-alpha cost are already known. No fresh B diagnostic "
                  "history has yet run. This is a disclosed sequential matched resource challenge using the same "
                  "prior mesh, seeds, controller, confidence allocation, limits and pilot-cap ratio.")
    for case in resources["selected_comparisons"]:
        case_methods = {}
        for output_method, source_method in methods:
            choice = case["variants"]["deferred_B"]["methods"][source_method]
            row = rows[choice["canonical_record_index"], choice["canonical_architecture"]]
            assert digest(row) == choice["canonical_record_sha256"]
            assert row["method"] == source_method and row["model_id"] == case["model_id"]
            assert row["architecture"] == "unitary8t"
            result = row["deferred_B_full_result"]
            flag, prior = result["flags"][0], result["correction_prior"]
            profile = resources["profiles"][row["deferred_profile_sha256"]]
            assert digest(profile) == row["deferred_profile_sha256"] == flag["financial_profile_sha256"]
            assert profile["work_dirty_after_preparation"]
            case_methods[output_method] = {
                "architecture": row["architecture"],
                "financial_preparation_profile": profile,
                "financial_preparation_profile_sha256": digest(profile),
                "probability_lower": prior["probability_lower"], "probability_upper": prior["probability_upper"],
                "finite_ideal_mean_lower": prior["finite_ideal_mean_lower"], "finite_ideal_mean_upper": prior["finite_ideal_mean_upper"],
                "probability_tolerance": flag["probability_tolerance"],
                "monetary_signal_factor_upper": result["monetary_signal_factor_upper"],
                "deterministic_price_error_upper": result["deterministic_error_upper"],
                "full_alpha_canonical_T": result["T_count"],
                "selected_strategy_sha256": digest({"architecture": row["architecture"], "result": result}),
                "target_monetary_radius": "1/100", "alpha": "1/3600",
                "matched_variants_proof": {
                    "financial_preparation_variant": "deferred_B=P U H",
                    "proof_filename": proof.name, "proof_sha256": sha(proof),
                    "deferred_resources_filename": "deferred_cleanup_resources.json.gz",
                    "deferred_resources_sha256": sha(out / "deferred_cleanup_resources.json.gz"),
                    "deferred_record_index": row["record_index"], "deferred_record_sha256": digest(row),
                    "source_clean_profile_sha256": profile["source_clean_profile_sha256"],
                    "source_clean_profile_location": {"results_relative_path": row["source_path"],
                        "source_row_index": row["source_row_index"], "source_row_sha256": row["source_row_sha256"],
                        "architecture": row["architecture"], "flag": "centred"},
                    "probability_certificate_sha256": row["probability_certificate_sha256"],
                    "complete_deterministic_budget_sha256": row["complete_deterministic_budget_sha256"],
                    "coefficient_informed_prior_sha256": row["coefficient_informed_prior_sha256"],
                    "financial_information_unchanged": True,
                    "legacy_adapter_symbol_A_means_supplied_B": True,
                    "sampler_metadata_annotated_by_deferred_wrapper": True,
                    "physical_cleanup_scope": ("The supplied B retains correlated arithmetic work. Every B-adjoint "
                        "clears the work before the same actual zero-reflection tape. Odd-degree terminal work is "
                        "dirty and all allocated wires are measured/reset; even-degree core-zero observation is "
                        "made in the restored clean-work invariant. The legacy numerical calculator's A "
                        "symbol denotes the supplied B. The staged wrapper explicitly corrects sampler "
                        "cleanup metadata without changing counts, numerical logic or controller decisions."),
                    "fallback_scope": ("The legacy complete-charge arithmetic is reused on the exact B primitive "
                        "Counter. Physical validity follows from the actual controlled-tape conjugacy proof, "
                        "not the legacy charge adapter's clean-work docstring. Complete charge equality is "
                        "checked against the registered deferred half-alpha result."),
                    "known_prior_outcomes": disclosure}}
        cases.append({"model_id": case["model_id"], "methods": case_methods,
                      "matched_variants_proof": {"proof_filename": proof.name, "proof_sha256": sha(proof),
                          "preparation_variant": "deferred_B", "prior_outcome_disclosure": disclosure}})
    inputs = {"schema": "qbp-matched-chebyshev-inputs-v1", "status": "checked registered canonical B winner extraction",
              "builder_sha256": sha(__file__), "parent_evidence_sha256": {
                  "deferred_cleanup_resources.json.gz": sha(out / "deferred_cleanup_resources.json.gz"),
                  "deferred_cleanup_protocol.json": sha(out / "deferred_cleanup_protocol.json"),
                  "results/common_shift/common_shift_resources.json.gz": sha(repo / "results/common_shift/common_shift_resources.json.gz"),
                  "results/matched_chebyshev/amended/matched_chebyshev_resources.json.gz": sha(repo / "results/matched_chebyshev/amended/matched_chebyshev_resources.json.gz"),
                  proof.name: sha(proof)},
              "prior_outcome_disclosure": disclosure, "cases": cases,
              "scope": "Same financial priors and exact B winner profiles. Stipulated common-signal diagnostic mesh only, not actual financial probabilities, expected hardware costs or new amplitude-estimation theory."}
    target = out / "matched_B_inputs.json"
    if target.exists():
        # No matched protocol or history has been registered/executed. Retain
        # the first extraction before this solely semantic clarification.
        previous = out / "before_matched_metadata_clarification"
        previous.mkdir(exist_ok=True)
        old_target = previous / target.name
        if old_target.exists():
            raise FileExistsError("Preserve prior input clarification")
        shutil.copyfile(target, old_target)
    target.write_text(json.dumps(inputs, indent=2) + "\n")
    adapter = out / "matched_chebyshev.py"
    if not adapter.exists():
        shutil.copyfile(repo / "matched_chebyshev.py", adapter)
    assert sha(adapter) == sha(repo / "matched_chebyshev.py")
    wrapper = out / "matched_deferred_chebyshev.py"
    if wrapper.exists():
        raise FileExistsError("Preserve executed metadata wrapper")
    text = adapter.read_text()
    old_key = "result = {'matched_chebyshev.py': sha(__file__)}"
    assert text.count(old_key) == 1
    text = text.replace(old_key, "result = {Path(__file__).name: sha(__file__)}")
    marker = "def register(repo, output, inputs):"
    helper = '''def annotate_deferred_sampler_metadata(resources):
    """Correct historical clean-work wording; preserve all numerical logic."""
    original = resources.build_profile
    def annotated(preparation, reflection, degree):
        assert preparation['work_dirty_after_preparation']
        profile = original(preparation, reflection, degree)
        profile['financial_cleanup_scope'] = (
            'The supplied preparation is B=P U H, with correlated arithmetic work. '
            'Its exact adjoint restores zero work before the same actual reflection tape. '
            'Odd terminal work remains dirty until full-allocation observation/reset; '
            'even core-zero observation occurs in the restored clean-work invariant. '
            'The full original allocated-wire measurement/reset charge is retained.')
        profile['preparation_symbol'] = 'B'
        profile['preparation_conjugacy_proof_sha256'] = preparation['proof_sha256']
        profile['financial_preparation_profiles_reused'] = True
        return profile
    resources.build_profile = annotated


'''
    assert text.count(marker) == 1
    text = text.replace(marker, helper + marker)
    call_marker = "    protocol = read(output/'matched_chebyshev_protocol.json')\n    verify(repo, protocol)"
    assert text.count(call_marker) == 1
    text = text.replace(call_marker, "    protocol = read(output/'matched_chebyshev_protocol.json')\n    annotate_deferred_sampler_metadata(resources)\n    verify(repo, protocol)")
    wrapper.write_text(text)
    print(json.dumps({"status": "prepared", "matched_B_inputs_sha256": sha(target),
                      "untouched_legacy_adapter_sha256": sha(adapter),
                      "metadata_wrapper_sha256": sha(wrapper), "cases": len(cases), "method_profiles": 12}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True)
    parser.add_argument("--proof", required=True)
    parser.add_argument("--output", required=True)
    main(parser.parse_args())
