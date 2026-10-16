"""Registered fair cleanup-deferral and cancellation-optimised clean comparison.

The published circuit identity removes redundant compute/uncompute pairs.
This producer reuses retained exact body/predicate profiles. It does not
claim fresh emission, a new quantum primitive or an estimator speed-up.
"""
from collections import Counter
from copy import deepcopy
from dataclasses import asdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import sys
sys.dont_write_bytecode = True

ALPHA = F(1, 3600)
HALF_ALPHA = F(1, 7200)
EPS = F(1, 100)
ARCHITECTURES = ("exact7t", "unitary8t")
METHODS = ("shared_threshold", "raw_payoff", "numerical_conditional_price")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def read(path):
    path = Path(path)
    return json.load(gzip.open(path, "rt")) if path.suffix == ".gz" else json.loads(path.read_text())


def write(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError("Preserve previous registered comparison: " + str(path))
    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
                stream.write(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode())
    else:
        path.write_text(json.dumps(value, indent=2, default=str) + "\n")


def source_files(repo):
    old = read(repo / "results/common_shift/common_shift_protocol.json")
    files = {"results/common_shift/common_shift_protocol.json": sha(repo / "results/common_shift/common_shift_protocol.json"),
             "results/common_shift/common_shift_resources.json.gz": sha(repo / "results/common_shift/common_shift_resources.json.gz"),
             "results/common_shift/common_shift_poolwise_fallback.json.gz": sha(repo / "results/common_shift/common_shift_poolwise_fallback.json.gz")}
    for path, h in old["resources_sha256"].items():
        assert sha(repo / "results" / path) == h
        files["results/" + path] = h
    for path in ("integral_pricing_budget.py", "qpe_estimator.py", "qft_synthesis.py", "correction_information_priors.py",
                 "boundary_precision_selection.py", "check_expanded_precision_resources.py", "quantum_gates.py",
                 "results/qft_rotation_tapes.json", "results/qft_rotation_certificates.json"):
        files[path] = sha(repo / path)
    return files


def bindings(row, index):
    return {"record_index": index, "record_sha256": digest(row), "model_id": row["model_id"], "method": row["method"],
            "source_path": row["source_path"], "source_row_index": row["source_row_index"],
            "source_row_sha256": row["source_row_sha256"], "profiles_sha256": row["profiles_sha256"],
            "probability_certificate_sha256": digest(row["probability_certificate"]),
            "complete_deterministic_budget_sha256": digest(row["complete_deterministic_budget"]),
            "coefficient_informed_prior_sha256": digest(row["coefficient_informed_prior"]),
            "configuration": {k: row[k] for k in ("normal_bits", "PG_exp_bits", "PA_exp_bits", "PG_bits", "PA_bits") if k in row}}


def inverse_counts(counts):
    inverse = {"T": "Tdg", "Tdg": "T", "S": "Sdg", "Sdg": "S"}
    result = Counter()
    for op, value in counts.items():
        result[inverse.get(op, op)] += value
    return result


def deferred_profile(clean, proof_sha):
    """Exact retained count transformation, with dirty work disclosed."""
    forward = Counter(clean["primitive_counts_forward_body"])
    field = "predicate_primitive_counts" if "predicate_primitive_counts" in clean else "comparator_primitive_counts"
    predicate = Counter(clean[field])
    uniform = Counter({"H": clean["uniform_H"]})
    assert Counter(clean["primitive_counts"]) == forward + forward + predicate + uniform
    counts = forward + predicate + uniform
    terminal = inverse_counts(forward)
    assert sum(counts.values()) + sum(forward.values()) == sum(clean["primitive_counts"].values())
    # These are retained serial block bounds for one forward body. Each
    # block is already a complete coherent primitive with clean private work.
    depth = sum(b["explicit_T_depth_serial_upper"] for b in clean["blocks"])
    five = sum(b["conservative_5layer_T_depth_serial_upper"] for b in clean["blocks"])
    assert depth <= clean["explicit_T_depth_serial_upper"] and five <= clean["conservative_5layer_T_depth_serial_upper"]
    return {"schema": "qbp-retained-deferred-financial-profile-v1", "architecture": clean["architecture"],
            "primitive_counts": dict(counts), "primitive_counts_forward_body": dict(forward),
            "predicate_primitive_counts": dict(predicate), "terminal_inverse_body_counts": dict(terminal),
            "logical_T_count": counts["T"] + counts["Tdg"],
            "terminal_inverse_body_T": terminal["T"] + terminal["Tdg"],
            "reflection_core_qubits": clean["reflection_core_qubits"], "uniform_H": clean["uniform_H"],
            "logical_qubits_upper": clean["logical_qubits_upper"],
            "explicit_T_depth_serial_upper": clean["explicit_T_depth_serial_upper"],
            "conservative_5layer_T_depth_serial_upper": clean["conservative_5layer_T_depth_serial_upper"],
            "terminal_inverse_body_explicit_T_depth_serial_upper": depth,
            "terminal_inverse_body_conservative_5layer_T_depth_serial_upper": five,
            "source_clean_profile_sha256": digest(clean), "source_blocks_sha256": digest(clean["blocks"]),
            "source_emitter_sha256": clean["source_sha256"],
            "proof_sha256": proof_sha, "frozen_body_profiles_reused": True,
            "fresh_whole_preparation_emission_claim": False, "work_dirty_after_preparation": True,
            "private_work_in_zero_reflection_core": False,
            "actual_reflection_tape_borrows_clean_work": True,
            "depth_scope": "Old clean serial upper bounds are conservatively retained for B; no new depth improvement or scheduling claim.",
            "scope": "B=P U H with retained coherent body U and predicate P. The stated conjugacy proof validates the same reflection core and good signal. Private work remains correlated and dirty until all-wire terminal measurement/reset, or one final U-adjoint for the optimised clean alternative."}


def charge_deferred_plan(setup, profile, qft):
    """Justified dirty-work extension of the retained complete-charge formula."""
    import integral_pricing_budget as budget
    import qpe_estimator as established
    if setup.zero_shots:
        return {"repetitions": 0, "T": 0, "scope": "probability-prior zero-shot decision"}
    budget.validate_oracle_profile(profile)
    assert profile["work_dirty_after_preparation"] and not profile["private_work_in_zero_reflection_core"]
    bp = Counter(profile["primitive_counts"])
    reflection = budget.reflection_metrics(profile["reflection_core_qubits"])
    bt = bp["T"] + bp["Tdg"]
    bm = {"T": bt, "measurements": 0, "feedback_CZ": 0, "reset_X": 0,
          "expanded_primitive_operations": sum(bp.values()), "uniform_H": profile["uniform_H"],
          "explicit_T_depth_serial_upper": profile["explicit_T_depth_serial_upper"],
          "conservative_5layer_T_depth_serial_upper": profile["conservative_5layer_T_depth_serial_upper"]}
    bm.update({"gate_" + op: value for op, value in bp.items()})
    qp = bp + inverse_counts(bp) + Counter(reflection["gate_counts"])
    qm = {"T": 2 * bt + reflection["T"], "measurements": 0, "feedback_CZ": 0, "reset_X": 0,
          "expanded_primitive_operations": sum(qp.values()), "uniform_H": 2 * profile["uniform_H"],
          "explicit_T_depth_serial_upper": 2 * bm["explicit_T_depth_serial_upper"] + reflection["explicit_T_depth_serial_upper"],
          "conservative_5layer_T_depth_serial_upper": 2 * bm["conservative_5layer_T_depth_serial_upper"] + reflection["conservative_5layer_T_depth_serial_upper"]}
    qm.update({"gate_" + op: value for op, value in qp.items()})
    peak = max(profile["logical_qubits_upper"], reflection["logical_qubits"])
    charge = established.charge(setup, {"metrics": bm, "logical_qubits": peak}, {"metrics": qm, "logical_qubits": peak}, qft)
    wires = charge["peak_logical_qubits"]
    for key, per in (("measurements", wires - setup.phase_register_bits), ("reset_X", wires), ("expanded_primitive_operations", 2 * wires)):
        charge["per_repetition_metrics"][key] += per
        charge["metrics"][key] += setup.repetitions * per
    charge.update(complete_oracle_T=bt, controlled_Q_T=qm["T"], controlled_reflection=reflection,
                  work_dirty_before_terminal_reset=True,
                  source="one retained deferred B, M-1 controlled Q_B and one emitted inverse QFT per fresh repetition",
                  reset_scope="All allocated wires, including correlated private work, are measured and conditionally X reset. Phase measurements are counted once. Terminal dirty work does not persist between repetitions.")
    assert all(charge["metrics"][k] == setup.repetitions * v for k, v in charge["per_repetition_metrics"].items())
    return charge


def charge_optimised_clean(deferred_charge, profile):
    """Add exactly the terminal U-adjoint omitted by deferred B."""
    charge = deepcopy(deferred_charge)
    repetitions = charge["repetitions"]
    per = {"T": profile["terminal_inverse_body_T"],
           "expanded_primitive_operations": sum(profile["terminal_inverse_body_counts"].values()),
           "explicit_T_depth_serial_upper": profile["terminal_inverse_body_explicit_T_depth_serial_upper"],
           "conservative_5layer_T_depth_serial_upper": profile["terminal_inverse_body_conservative_5layer_T_depth_serial_upper"]}
    per.update({"gate_" + op: count for op, count in profile["terminal_inverse_body_counts"].items()})
    for key, value in per.items():
        charge["per_repetition_metrics"][key] = charge["per_repetition_metrics"].get(key, 0) + value
        charge["metrics"][key] = charge["metrics"].get(key, 0) + repetitions * value
    charge.update(terminal_cleanup_T_per_repetition=profile["terminal_inverse_body_T"],
                  terminal_cleanup_inverse_body_counts=profile["terminal_inverse_body_counts"],
                  work_dirty_before_terminal_reset=False,
                  source="Cancellation-optimised clean QPE: B-QPE followed by exactly one terminal U-adjoint per fresh repetition. Published conjugacy cancels all other redundant U/U-adjoint pairs.")
    assert all(charge["metrics"][k] == repetitions * v for k, v in charge["per_repetition_metrics"].items())
    return charge


def fixed_inputs(args):
    repo = Path(args.repo)
    return {"producer_source_sha256": sha(__file__), "checked_proof_sha256": sha(args.proof),
            "checked_proof_filename": Path(args.proof).name, "inputs_sha256": source_files(repo)}


def register(args):
    repo, out = Path(args.repo), Path(args.output)
    data = read(repo / "results/common_shift/common_shift_resources.json.gz")
    assert len(data["records"]) == 416 and data["uncosted_admissions"] == 0
    protocol = {"schema": "qbp-deferred-cleanup-protocol-v1", "registration_utc": datetime.now(timezone.utc).isoformat(),
                **fixed_inputs(args), "registered_before_deferred_or_optimised_clean_cost_outcomes": True,
                "known_outcomes": "All retained clean common-shift financial certificates and complete costs are known. The published compute/uncompute cancellation identity is known. This is a disclosed sequential cost challenge, not an initially blinded experiment.",
                "records": [bindings(r, i) for i, r in enumerate(data["records"])],
                "configurations": data["configurations"], "architectures": list(ARCHITECTURES),
                "variants": ["deferred_B", "cancellation_optimised_clean"],
                "profile_rule": "Retained clean primitive counts must equal2F+P+H. Set B counts=F+P+H, preserve old core, allocated-wire and conservative depth bounds, and remove the active-after-clean-preparation claim. Add exactly inverse F once per repetition for the cancellation-optimised clean comparator.",
                "mathematical_rule": "B=U A and [U,Z_good]=0 give Q_B=U Q_A U-adjoint. Controlled-query/QFT signal laws and all financial certificates, probability priors, monetary budgets and estimation tolerances are unchanged. Cancellation-optimised clean QPE is B-QPE followed by one terminal U-adjoint. No private work wire is added to or removed from the retained S0 reflection core.",
                "retained_financial_information": "Exact common_shift eta, finite ideal mean interval, implemented prior, coefficient/mass/transfer interval, monetary radius1/100 and full alpha1/3600. No price input. Every rejection is retained, including12 numerical rows with no frozen profiles.",
                "fallback_rule": "Rebuild every admitted architecture with half alpha1/7200 and identical target, upper prior, QFT tape and B profile. For each variant T_worst=T_fallback+floor(T_fallback/10). Optimise this guaranteed upper cap over the full retained pool independently for each method/variant/case.",
                "selection": "Minimise each complete cap over all admitted records/architectures separately by method/case/variant. Tie-break by normal bits, PG bits, PA bits, architecture, record index. Retain naive clean, deferred and optimised-clean comparisons.",
                "scope": "Published cleanup deferral/circuit cancellation applied fairly to all implemented financial encodings. Exact retained primitive-resource transformations and complete confidence/reset charge caps; not a fresh emitted whole circuit, new primitive, globally optimal estimator, actual finite-price observation, hardware runtime or practical quantum advantage."}
    write(out / "deferred_cleanup_protocol.json", protocol)
    print("Registered 416 financial settings/two architectures and both cleanup variants before fresh costs", flush=True)


def tie(row):
    c = row["configuration"]
    return (c.get("normal_bits", 0), c.get("PG_exp_bits", c.get("PG_bits", 0)),
            c.get("PA_exp_bits", c.get("PA_bits", 0)), row["architecture"], row["record_index"])


def select(rows, old_data, old_fallback):
    old_min = {r["model_id"]: r for r in old_data["selected_comparisons"]}
    old_half = {r["model_id"]: r for r in old_fallback["selected_comparisons"]}
    comparisons = []
    for mid in dict.fromkeys(r["model_id"] for r in rows):
        variants = {}
        for variant in ("deferred_B", "cancellation_optimised_clean"):
            methods = {}
            for method in METHODS:
                group = [r for r in rows if r["model_id"] == mid and r["method"] == method]
                canonical = min(group, key=lambda r: (r[variant]["T_full_alpha"], tie(r)))
                worst = min(group, key=lambda r: (r[variant]["T_worst_history_upper"], tie(r)))
                old_T = old_min[mid]["methods"][method]["T_count"]
                old_worst = old_half[mid]["methods"][method]["T_worst_history_upper"]
                methods[method] = {"canonical_record_index": canonical["record_index"], "canonical_architecture": canonical["architecture"],
                                   "canonical_configuration": canonical["configuration"], "T_full_alpha": canonical[variant]["T_full_alpha"],
                                   "canonical_record_sha256": digest(canonical),
                                   "worst_history_record_index": worst["record_index"], "worst_history_architecture": worst["architecture"],
                                   "worst_history_configuration": worst["configuration"], "T_half_alpha_fallback": worst[variant]["T_half_alpha_fallback"],
                                   "T_worst_history_upper": worst[variant]["T_worst_history_upper"],
                                   "previous_naive_clean_T_full_alpha": old_T, "previous_naive_clean_T_worst_history_upper": old_worst,
                                   "saving_vs_naive_clean_canonical_exact": str(F(old_T - canonical[variant]["T_full_alpha"], old_T)),
                                   "saving_vs_naive_clean_worst_upper_exact": str(F(old_worst - worst[variant]["T_worst_history_upper"], old_worst))}
            rivals = {m: methods[m]["T_full_alpha"] for m in ("raw_payoff", "numerical_conditional_price")}
            strong = min(rivals, key=lambda m: (rivals[m], m))
            T, E = rivals[strong], methods["shared_threshold"]["T_full_alpha"]
            modern_rivals = {m: methods[m]["T_worst_history_upper"] for m in ("raw_payoff", "numerical_conditional_price")}
            modern = min(modern_rivals, key=lambda m: (modern_rivals[m], m))
            TW, EW = modern_rivals[modern], methods["shared_threshold"]["T_worst_history_upper"]
            variants[variant] = {"methods": methods, "strongest_rival": strong,
                                 "saving_against_strongest_canonical_exact": str(F(T - E, T)),
                                 "saving_against_strongest_canonical_per_cent_diagnostic": float(100 * F(T - E, T)),
                                 "strongest_worst_upper_rival": modern,
                                 "saving_against_strongest_worst_upper_exact": str(F(TW - EW, TW)),
                                 "saving_against_strongest_worst_upper_per_cent_diagnostic": float(100 * F(TW - EW, TW))}
        variants["terminal_cleanup_difference"] = {
            m: {"canonical_deferred_saving_vs_optimised_clean_exact": str(F(
                variants["cancellation_optimised_clean"]["methods"][m]["T_full_alpha"] - variants["deferred_B"]["methods"][m]["T_full_alpha"],
                variants["cancellation_optimised_clean"]["methods"][m]["T_full_alpha"]))} for m in METHODS}
        comparisons.append({"model_id": mid, "variants": variants})
    return comparisons


def run(args):
    repo, out = Path(args.repo), Path(args.output)
    protocol = read(out / "deferred_cleanup_protocol.json")
    assert all(protocol[k] == v for k, v in fixed_inputs(args).items())
    old_data = read(repo / "results/common_shift/common_shift_resources.json.gz")
    old_half = read(repo / "results/common_shift/common_shift_poolwise_fallback.json.gz")
    half_index = {(r["record_index"], r["architecture"]): r for r in old_half["records"]}
    sys.path.insert(0, str(repo))
    import integral_pricing_budget as budget
    import correction_information_priors as frontend
    original_charge = budget.charge_full_oracle
    budget.charge_full_oracle = charge_deferred_plan
    sources, profiles, rows, rejections = {}, {}, [], []
    try:
        for i, (row, binding) in enumerate(zip(old_data["records"], protocol["records"])):
            assert bindings(row, i) == binding
            path = row["source_path"]
            if path not in sources:
                sources[path] = read(repo / "results" / path)
            base = sources[path]["rows" if row["method"] == "shared_threshold" else "records"][row["source_row_index"]]
            assert digest(base) == row["source_row_sha256"] and digest(base["profiles"]) == row["profiles_sha256"]
            if not row["strategies"]:
                assert row.get("rejection") and not base["profiles"]
                rejections.append({**binding, "retained_rejection": row["rejection"]})
                continue
            for architecture in ARCHITECTURES:
                clean = base["profiles"][architecture]["centred"]
                B = deferred_profile(clean, protocol["checked_proof_sha256"])
                h = digest(B)
                assert h not in profiles
                profiles[h] = B
                result = frontend.informed_complete_plan("positive", row["complete_deterministic_budget"],
                            row["probability_certificate"], row["coefficient_informed_prior"], {"centred": B})
                assert result["kind"] == "complete_quantum_resource_plan" and len(result["flags"]) == 1
                previous = next(s["result"] for s in row["strategies"] if s["architecture"] == architecture)
                flag, old_flag = result["flags"][0], previous["flags"][0]
                for field in ("plan", "fourier_profile", "probability_tolerance", "prior_upper"):
                    assert digest(flag[field]) == digest(old_flag[field])
                for field in ("deterministic_error_upper", "statistical_price_error_upper", "total_price_error_upper", "failure_probability_upper", "coefficient_informed_probability_interval"):
                    assert digest(result[field]) == digest(previous[field])
                bcharge = flag["complete_charge"]
                ccharge = charge_optimised_clean(bcharge, B)
                half_plan, half_qft = budget.compiled_extended_plan(F(flag["probability_tolerance"]), HALF_ALPHA, F(flag["prior_upper"]))
                half_b = charge_deferred_plan(half_plan, B, half_qft)
                half_c = charge_optimised_clean(half_b, B)
                hp = asdict(half_plan)
                assert digest(half_qft) == digest(flag["fourier_profile"])
                assert hp["M"] == flag["plan"]["M"] and (flag["plan"]["repetitions"], hp["repetitions"]) == (25, 27)
                previous_T = previous["T_count"]
                repetitions, M = flag["plan"]["repetitions"], flag["plan"]["M"]
                assert previous_T - result["T_count"] == repetitions * (2 * M - 1) * B["terminal_inverse_body_T"]
                assert ccharge["metrics"]["T"] - bcharge["metrics"]["T"] == repetitions * B["terminal_inverse_body_T"]
                for full, half in ((bcharge, half_b), (ccharge, half_c)):
                    assert F(half["metrics"]["T"], full["metrics"]["T"]) == F(27, 25)
                rec = {**binding, "architecture": architecture, "deferred_profile_sha256": h,
                       "previous_current_clean_result_sha256": digest(previous),
                       "previous_naive_clean_T_full_alpha": previous_T,
                       "previous_naive_clean_T_half_alpha_fallback": half_index[i, architecture]["T_half_alpha_fallback"],
                       "deferred_B_full_result": result, "half_plan": hp,
                       "half_fourier_profile_sha256": digest(half_qft),
                       "half_probability_tolerance": flag["probability_tolerance"], "half_prior_upper": flag["prior_upper"],
                       "deferred_B_half_charge": half_b,
                       "cancellation_optimised_clean_full_charge": ccharge,
                       "cancellation_optimised_clean_half_charge": half_c,
                       "financial_target_prior_budget_phaseplan_unchanged": True}
                for variant, full, half in (("deferred_B", bcharge, half_b), ("cancellation_optimised_clean", ccharge, half_c)):
                    fullT, halfT = full["metrics"]["T"], half["metrics"]["T"]
                    rec[variant] = {"T_full_alpha": fullT, "T_half_alpha_fallback": halfT,
                                    "pilot_T_budget": halfT // 10, "T_worst_history_upper": halfT + halfT // 10}
                rows.append(rec)
            if (i + 1) % 52 == 0:
                print("Rebuilt", i + 1, "financial settings", len(rows), "paired architecture plans", flush=True)
    finally:
        budget.charge_full_oracle = original_charge
    assert len(rows) == 808 and len(rejections) == 12
    assert all(protocol[k] == v for k, v in fixed_inputs(args).items())
    comparisons = select(rows, old_data, old_half)
    resources = {"schema": "qbp-deferred-cleanup-resources-v1", "status": "completed retaining every registered outcome",
                 "protocol_sha256": sha(out / "deferred_cleanup_protocol.json"), **fixed_inputs(args),
                 "profiles": profiles, "records": rows, "rejections": rejections, "selected_comparisons": comparisons,
                 "registered_configurations": 416, "admitted_architecture_plans": 808,
                 "variant_canonical_plans": 1616, "variant_half_alpha_fallback_plans": 1616,
                 "rejected_configurations": 12, "frozen_mathematical_profiles_reused": True,
                 "data_origin": "Synthetic rational models and retained certified coherent primitives; no licensed or private observations.",
                 "scope": protocol["scope"]}
    write(out / "deferred_cleanup_resources.json.gz", resources)
    write(out / "deferred_cleanup_summary.json", {k: v for k, v in resources.items() if k not in ("records", "profiles", "rejections")})
    print("Completed 808 paired canonical and half-alpha plans; 12 retained rejections", flush=True)
    for case in comparisons:
        print(case["model_id"], {v: {m: x["T_full_alpha"] for m, x in case["variants"][v]["methods"].items()}
                               for v in ("deferred_B", "cancellation_optimised_clean")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("register", "run"))
    parser.add_argument("--repo", required=True)
    parser.add_argument("--proof", required=True)
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args()
    register(arguments) if arguments.mode == "register" else run(arguments)
