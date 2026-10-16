"""Independent audit of stress transfer and continuous Taylor certificates.

The audit uses 1,024-bit factorial transport series and direct rational U
polynomials. It checks order-intersection arithmetic and immutable degree-three
work snapshots, as well as the retained scale-one mathematical reproduction.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
import argparse, hashlib, json, time
from flint import arb, ctx
from audit_dispersion_models import parameters
import independent_dispersion_math as im

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_case(case, model, geometry, norm, maximum, counters):
    assert [o["order"] for o in case["orders"]] == list(range(1, maximum + 1))
    assert len(case["phase_work"]) == (1 if maximum == 3 else 2)
    assert case["phase_work"][0]["maximum_coefficient_degree"] == 3
    if maximum == 5:
        assert case["phase_work"][1]["starting_coefficient_degree"] == 3
        assert case["phase_work"][1]["maximum_coefficient_degree"] == 5
    w, R, mu, lam, rg, mg, d, s, Q, V = parameters(model)
    cap = None if case["cap"] is None else F(case["cap"])
    coeff = im.coefficients(geometry, mg, lam, cap, maximum)
    assert im.contains_record(case["c0_c1_certificate"]["c0"], coeff[0])
    assert im.contains_record(case["c0_c1_certificate"]["c1"], coeff[1])
    counters["coefficient_enclosures"] += 2
    partial = coeff[0]
    first = None
    for order in case["orders"]:
        k = order["order"]
        assert im.contains_record(order["coefficient_interval"], coeff[k]), (case["model_id"], cap, k, "coefficient")
        partial += coeff[k]
        assert im.contains_record(order["partial_sum_interval"], partial), (case["model_id"], cap, k, "partial")
        counters["coefficient_enclosures"] += 1
        M = F(norm["C_r_upper"][str(k + 1)])
        moment = order["remainder_U_moment_interval"]
        assert im.contains_record(moment, geometry.full_moments[k + 1])
        assert F(order["global_C_r_norm_upper"]) == M
        radius = M * F(moment["upper_exact_fraction"]) / factorial(k + 1)
        assert F(order["remainder_upper"]) == radius
        savedpartial = order["partial_sum_interval"]
        rawlo = max(F(0), F(savedpartial["lower_exact_fraction"]) -
                    (F(0) if cap is None and k == 1 else radius))
        rawhi = F(savedpartial["upper_exact_fraction"]) + radius
        if cap is not None:
            rawhi = min(cap, rawhi)
        assert F(order["undiscounted_price_radius_upper"]) == (rawhi - rawlo) / 2
        price = order["continuous_price_interval"]
        assert F(price["lower_exact_fraction"]) <= rawlo <= rawhi <= F(price["upper_exact_fraction"])
        assert ((rawhi - rawlo) / 2 <= F(1, 100)) == order["certifies_undiscounted_cent"]
        if order["certifies_undiscounted_cent"] and first is None:
            first = k
        counters["Taylor_remainders_and_raw_halfwidths"] += 1
    assert first == case["first_successful_order"]
    assert F(case["best_price_radius_upper"]) == min(F(o["undiscounted_price_radius_upper"]) for o in case["orders"])

    # Strengthened classical screen: exact rational intersections of saved
    # valid intervals and the correctly directed intrinsic forward cap.
    EA = 100 * im.ball(F(model["risk_free_rate"]) - F(model["dividend_yield"])).exp()
    assert im.contains_record(case["EA_interval"], EA)
    intrinsic = case["intrinsic_price_interval"]
    assert F(intrinsic["lower_exact_fraction"]) == 0
    hi = F(intrinsic["upper_exact_fraction"])
    assert hi == cap if cap is not None else hi >= im.endpoint(EA, True)
    lo = F(0)
    intersection_first = None
    for order, progress in zip(case["orders"], case["intersection_progress"]):
        saved = order["continuous_price_interval"]
        lo = max(lo, F(saved["lower_exact_fraction"]))
        hi = min(hi, F(saved["upper_exact_fraction"]))
        assert lo <= hi
        interval = progress["price_interval"]
        assert F(interval["lower_exact_fraction"]) == lo
        assert F(interval["upper_exact_fraction"]) == hi
        radius = (hi - lo) / 2
        assert F(progress["price_radius_upper"]) == radius
        assert progress["through_degree"] == order["order"]
        assert progress["certifies_undiscounted_cent"] == (radius <= F(1, 100))
        if radius <= F(1, 100) and intersection_first is None:
            intersection_first = order["order"]
        counters["exact_interval_intersections"] += 1
    assert case["intersection_price_interval"] == case["intersection_progress"][-1]["price_interval"]
    assert F(case["intersection_radius_upper"]) == (hi - lo) / 2
    assert case["intersection_certifies_cent"] == ((hi - lo) / 2 <= F(1, 100))
    assert case["intersection_first_successful_order"] == intersection_first

    # Count actual raw-T visits and distinct mean-shift points separately;
    # the two are not interchangeable work quantities.
    n = len(w)
    phases = [(2, 3)] if maximum == 3 else [(2, 3), (4, 5)]
    for work, (firstk, lastk) in zip(case["phase_work"], phases):
        count = sum(sum(comb(n + j - 1, j) for j in range(k + 1)) for k in range(firstk, lastk + 1))
        assert work["weighted_kernel_term_visits"] == count
        counters["phase_visit_counts"] += 1
    shifts = {geometry.basis_value(a)[3] for poly in geometry.powers[:maximum + 1] for a in poly}
    bc = case["boundary_operation_counts"]
    assert bc["distinct_tilted_boundary_points"] == len(shifts)
    assert bc["Gaussian_density_evaluations"] == len(shifts) * (1 if cap is None else 2)
    assert bc["falling_boundary_derivative_evaluations"] == len(shifts) * (maximum - 1)
    expected_updates = sum(2 + 3 * (r - 2) for r in range(2, maximum))
    assert bc["polynomial_coefficient_recurrence_updates"] == expected_updates
    counters["snapshot_operation_counts"] += 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path, default=ROOT.parent)
    ap.add_argument("--out", type=Path, default=ROOT / "dispersion_certificate_audit.json")
    args = ap.parse_args()
    stage = args.stage
    ctx.prec = im.PRECISION
    started = time.perf_counter()
    resultpath = stage / "results/dispersion_separation_results.json"
    result = json.loads(resultpath.read_text())
    protocol = stage / "results/dispersion_separation_protocol.json"
    models = {m["model_id"]: m for m in json.loads(protocol.read_text())["models"]}
    assert result["protocol_sha256"] == digest(protocol)
    for name, value in result["source_sha256"].items():
        assert value == digest(stage / name)
    for name, value in result["norm_input_sha256"].items():
        assert value == digest(stage / "results" / name)
    assert len(models) == 36 and len(result["records"]) == 144
    normrows = json.loads((stage / "results/boundary_derivative_norms.json").read_text())["records"]
    counters = {"transport_input_enclosures": 0, "hybrid_moment_enclosures": 0,
                "continuous_moment_enclosures": 0, "outside_moment_enclosures": 0,
                "derivative_envelope_components": 0, "complete_transfer_upper_checks": 0,
                "coefficient_enclosures": 0, "Taylor_remainders_and_raw_halfwidths": 0,
                "exact_interval_intersections": 0, "phase_visit_counts": 0,
                "snapshot_operation_counts": 0, "snapshot_immutability_checks": 0,
                "scale_one_bound_and_order_matches": 0}
    transport = im.IndependentTransport(result["transport_actual_inputs"])
    counters["transport_input_enclosures"] = transport.input_checks
    assert im.contains_record(result["transport_profile"]["mass_interval"], transport.mass)
    groups = {}
    moments = {}
    tailvalues = {}
    for saved in result["continuous_geometry_records"]:
        key = (saved["assets"], saved["dispersion_scale"])
        model = next(m for m in models.values() if (m["assets"], m["dispersion_scale"]) == key)
        w, R, mu, lam, rg, mg, d, s, Q, V = parameters(model)
        degree = max(map(int, saved["continuous_U_moment_intervals"]))
        geometry = im.IndependentGeometry(w, d, s, rg, degree)
        groups[key] = geometry
        for k, moment in enumerate(geometry.full_moments):
            assert im.contains_record(saved["continuous_U_moment_intervals"][str(k)], moment), (key, k)
            counters["continuous_moment_enclosures"] += 1
        assert saved["monomial_counts_by_degree"] == [comb(len(w) + k - 1, k) for k in range(degree + 1)]
        assert saved["phase_work"][0]["maximum_moment_degree"] == 4
        assert len(saved["phase_work"]) == (2 if degree == 6 else 1)
        for phase in saved["phase_work"]:
            first, last = phase["first_new_moment_degree"], phase["maximum_moment_degree"]
            assert phase["new_Gaussian_geometry_exponential_evaluations"] == sum(comb(len(w) + k - 1, k) for k in range(first, last + 1))
    for saved in result["tail_records"]:
        key = (saved["assets"], saved["dispersion_scale"])
        tailvalues[key] = groups[key].outside(groups[key].powers[2])
        assert im.contains_record(saved["U2_outside_interval"], tailvalues[key]), key
        counters["outside_moment_enclosures"] += 1
    for saved in result["hybrid_geometry_records"]:
        key = (saved["assets"], saved["dispersion_scale"])
        geometry = groups[key]
        rows = {}
        for row in saved["coordinates"]:
            j = row["coordinate"]
            assert F(row["geometric_loading"]) == geometry.rg[j]
            polys = im.base_polynomials(geometry.w, geometry.s, j)
            rows[j] = {}
            for kind in (1, 2):
                rows[j][kind] = {}
                for name, poly in polys.items():
                    value = geometry.hybrid(poly, j, kind, transport)
                    assert im.contains_record(row["moments"][str(kind)][name], value), (key, j, kind, name)
                    rows[j][kind][name] = value
                    counters["hybrid_moment_enclosures"] += 1
        moments[key] = rows
    oldtransfer = json.loads((stage / "results/tensor_transfer.json").read_text())["records"]
    oldprice = json.loads((stage / "results/original_degree3_price_intervals.json").read_text())["records"]
    outcomes = []
    for row in result["records"]:
        model = models[row["model_id"]]
        key = (row["assets"], row["dispersion_scale"])
        g = groups[key]
        cap = None if row["cap"] is None else F(row["cap"])
        norm = next(n for n in normrows if F(n["lambda"]) == F(model["lambda"])
                    and n["payoff"] == row["payoff"]
                    and (None if n["cap"] is None else F(n["cap"])) == cap)
        transfer = row["aggregate_R2_transport"]
        sum_upper = F(0)
        for coordinate in transfer["coordinates"]:
            j = coordinate["coordinate"]
            for kind, label in [(1, "first_transport"), (2, "second_transport")]:
                components, constants = im.envelope_components(norm, g.rg[j], moments[key][j][kind], kind)
                saved = coordinate[label]
                for name, value in components.items():
                    assert im.endpoint(value, True) <= F(saved["components_upper"][name]), (row["model_id"], cap, j, kind, name)
                    counters["derivative_envelope_components"] += 1
                assert F(saved["error_upper"]) == sum(map(F, saved["components_upper"].values()))
                sum_upper += F(saved["error_upper"])
        assert {k: F(v) for k, v in transfer["norms"].items()} == constants
        tail = F(norm["M2_upper"]) * im.endpoint(tailvalues[key], True) / 2
        assert tail <= F(transfer["tail_upper"])
        assert F(transfer["financial_transfer_upper"]) == sum_upper + F(transfer["tail_upper"])
        passed = F(transfer["financial_transfer_upper"]) < F(1, 200)
        assert passed == (transfer["status"] == "passed")
        counters["complete_transfer_upper_checks"] += 1
        check_case(row["degree3_continuous_price"], model, g, norm, 3, counters)
        expected_selection = passed and not row["degree3_continuous_price"]["intersection_certifies_cent"]
        assert row["selected_for_degree5"] == row["separates_after_degree3"] == expected_selection
        later = row["degree5_continuous_price"]
        if expected_selection:
            assert later is not None
            check_case(later, model, g, norm, 5, counters)
            assert later["orders"][:3] == row["degree3_continuous_price"]["orders"]
            assert later["phase_work"][:1] == row["degree3_continuous_price"]["phase_work"]
            assert later["boundary_operation_counts"] != row["degree3_continuous_price"]["boundary_operation_counts"]
            counters["snapshot_immutability_checks"] += 3
        else:
            assert later is None
        assert row["separates_after_degree5"] == (expected_selection and not later["intersection_certifies_cent"])
        if row["dispersion_scale"] == "1":
            old = next(r for r in oldtransfer if r["model_id"] == row["original_model_id"] and r["cap"] == row["cap"])
            assert abs(F(old["financial_transfer_upper"]) - F(transfer["financial_transfer_upper"])) < F(1, 10 ** 60)
            prior = next(r for r in oldprice if r["model_id"] == row["original_model_id"] and r["cap"] == row["cap"])
            assert len(row["degree3_continuous_price"]["orders"]) == len(prior["orders"]) == 3
            for neworder, oldorder in zip(row["degree3_continuous_price"]["orders"], prior["orders"]):
                assert abs(F(neworder["undiscounted_price_radius_upper"]) - F(oldorder["undiscounted_radius_upper"])) < F(1, 10 ** 60)
                counters["scale_one_bound_and_order_matches"] += 1
            counters["scale_one_bound_and_order_matches"] += 1
        outcomes.append({"scale": row["dispersion_scale"], "transfer_pass": passed,
                         "degree3_pass": row["degree3_continuous_price"]["intersection_certifies_cent"],
                         "degree3_individual_pass": row["degree3_continuous_price"]["first_successful_order"] is not None,
                         "selected": expected_selection,
                         "degree5_rescued": expected_selection and later["intersection_certifies_cent"],
                         "remaining": row["separates_after_degree5"]})
    summaries = []
    for scale in ("1", "3/2", "2", "3"):
        selected = [o for o in outcomes if o["scale"] == scale]
        summaries.append({"dispersion_scale": scale, "cases": len(selected),
                          "transport_passes": sum(o["transfer_pass"] for o in selected),
                          "degree3_price_passes": sum(o["degree3_pass"] for o in selected),
                          "degree3_individual_interval_price_passes": sum(o["degree3_individual_pass"] for o in selected),
                          "degree3_separators": sum(o["selected"] for o in selected),
                          "degree5_evaluated": sum(o["selected"] for o in selected),
                          "degree5_rescued": sum(o["degree5_rescued"] for o in selected),
                          "remaining_degree5_separators": sum(o["remaining"] for o in selected)})
    assert summaries == result["summary_by_scale"]
    out = {"schema": "qbp-dispersion-separation-independent-certificate-audit-v1", "status": "passed",
           "precision_bits": im.PRECISION, "transport_series_order": im.SERIES,
           "audit_counts": counters, "summary_by_scale": summaries,
           "independent_scalar_transform_queries": len(transport.cache),
           "elapsed_seconds_diagnostic": time.perf_counter() - started,
           "method": "Direct exact rational U polynomials, cancellation-safe factorial Peano transforms, independent falling-operator/Hermite kernels, all full tails and rational interval intersections; no production helper imported",
           "source_sha256": digest(__file__), "independent_math_source_sha256": digest(ROOT / "independent_dispersion_math.py"),
           "model_audit_source_sha256": digest(ROOT / "audit_dispersion_models.py"),
           "result_sha256": digest(resultpath), "protocol_sha256": digest(protocol),
           "scope": "Financial mathematical certificates and declared comparator admission, not a complete quantum pricing algorithm or classical impossibility result"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "audit_counts": counters,
                      "summary_by_scale": summaries}, indent=2))


if __name__ == "__main__":
    main()
