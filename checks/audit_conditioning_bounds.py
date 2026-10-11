"""Independent audit of geometric exercise-family and AG bounds.

The audit reconstructs all Gaussian covariances, events and geometric calls,
checks every root-history sign and endpoint comparison, and independently
propagates all 36 call-spread intervals. No production helper is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
from flint import arb, ctx

PRECISION = 1024


def ball(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def endpoint(x, top):
    y = x.upper() if top else x.lower()
    m, e = map(int, y.man_exp())
    return F(m * (1 << e)) if e >= 0 else F(m, 1 << -e)


def interval(lo, hi):
    lo, hi = F(lo), F(hi)
    return ball((lo + hi) / 2) + arb(0, ball((hi - lo) / 2).upper())


def contained(rec, value):
    return (F(rec["lower_exact_fraction"]) <= endpoint(value, False)
            <= endpoint(value, True) <= F(rec["upper_exact_fraction"]))


def sign(x):
    return 1 if endpoint(x, False) > 0 else -1 if endpoint(x, True) < 0 else 0


def phi(x):
    return (-x * x / 2).exp() / (2 * arb.pi()).sqrt()


def Phi(x):
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def local_range(function, derivative, lo, hi):
    # Independent mean-value enclosure is tighter than simply rebuilding
    # an Arb natural interval. The latter need not be contained in another
    # valid natural interval because radius magnitudes have finite precision.
    lo, hi = F(lo), F(hi)
    midpoint = (lo + hi) / 2
    value = function(ball(midpoint))
    error = ball((hi - lo) / 2) * abs(derivative(interval(lo, hi)))
    radius = endpoint(error, True)
    return interval(endpoint(value, False) - radius, endpoint(value, True) + radius)


class IndependentFamily:
    def __init__(self, model):
        self.model = model
        self.w = tuple(map(F, model["weights"]))
        self.mu = tuple(map(F, model["log_drifts"]))
        self.r = tuple(tuple(map(F, r)) for r in model["loadings"])
        self.lam = F(model["lambda"])
        assert self.lam > 0 and min(self.w) > 0 and sum(self.w) == 1
        self.rg = tuple(sum(self.w[i] * self.r[i][j] for i in range(len(self.w)))
                        for j in range(model["outer_dimension"]))
        self.mg = dot(self.w, self.mu)
        self.var = self.lam ** 2 + dot(self.rg, self.rg)
        self.sigma = ball(self.var).sqrt()
        self.cov = tuple(self.lam ** 2 + dot(row, self.rg) for row in self.r)
        assert dot(self.w, self.cov) == self.var > 0
        self.a = tuple(ball(c) / self.sigma for c in self.cov)
        # The model's martingale condition gives a common expectation; I use
        # that equality instead of the producer's separate exponent rows.
        drift = F(model["risk_free_rate"]) - F(model["dividend_yield"])
        assert all(self.mu[i] + (dot(self.r[i], self.r[i]) + self.lam ** 2) / 2 == drift
                   for i in range(len(self.w)))
        self.forward = tuple(100 * ball(w) * ball(drift).exp() for w in self.w)
        self.EA = 100 * ball(drift).exp()
        self.EG = 100 * ball(self.mg + self.var / 2).exp()
        self.classification = "mixed" if min(self.cov) < 0 < max(self.cov) else "nonnegative"

    def S(self, z, order=0):
        return sum((f * a ** order * (a * z - a * a / 2).exp()
                    for f, a in zip(self.forward, self.a)), arb(0))

    def L(self, z, strike):
        return sum((f * Phi(a - z) for f, a in zip(self.forward, self.a)), arb(0)) - ball(strike) * Phi(-z)

    def geometric(self, strike):
        x = (ball(self.mg) - (ball(strike) / 100).log()) / self.sigma
        return self.EG * Phi(x + self.sigma) - ball(strike) * Phi(x)


def history_check(history, function, counters):
    for item in history:
        value = function(ball(F(item["mid"])))
        assert contained(item["score_interval"], value)
        recorded = item["certified_sign"]
        if recorded:
            assert recorded == sign(value)
        else:
            assert F(item["score_interval"]["lower_exact_fraction"]) <= 0 <= F(item["score_interval"]["upper_exact_fraction"])
        for key in ("quarter_left", "quarter_right"):
            if key in item:
                r = item[key]
                assert contained(r["score"], function(ball(F(r["x"]))))
        counters["bisection_history_entries"] += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path,
                    default=Path("/private/tmp/qbp-residual-priority-step/cost_boundary"))
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "conditioning_bounds_independent_audit.json")
    args = ap.parse_args()
    ctx.prec = PRECISION
    stage = args.stage
    protocol = stage / "results/generalisation_protocol.json"
    results = stage / "results/caldana_conditioning_benchmark.json"
    producer = json.loads(results.read_text())
    models = {m["model_id"]: m for m in json.loads(protocol.read_text())["models"]}
    assert producer["source_sha256"] == hashlib.sha256(
        (stage / "caldana_conditioning_benchmark.py").read_bytes()).hexdigest()
    assert producer["protocol_sha256"] == hashlib.sha256(protocol.read_bytes()).hexdigest()
    for name, digest in producer["dependencies"].items():
        assert digest == hashlib.sha256((stage / name).read_bytes()).hexdigest()
    counters = {"family_parameter_intervals": 0, "minimum_brackets": 0,
                "critical_root_brackets": 0, "bisection_history_entries": 0,
                "Lprime_identity_checks": 0, "candidate_enclosures": 0,
                "geometric_AG_intervals": 0, "call_intervals": 0,
                "contract_propagations": 0, "discounted_intervals": 0}
    families = {}
    family_records = {}
    for rec in producer["families"]:
        f = IndependentFamily(models[rec["model_id"]])
        families[rec["model_id"]] = f
        family_records[rec["model_id"]] = rec
        assert f.classification == rec["classification"]
        assert tuple(map(F, rec["geometric_covariances_exact"])) == f.cov
        for r, value in [(rec["sigma_G_interval"], f.sigma),
                         (rec["EA_interval"], f.EA), (rec["EG_interval"], f.EG)]:
            assert contained(r, value)
            counters["family_parameter_intervals"] += 1
        for r, value in zip(rec["standardised_covariance_intervals"], f.a):
            assert contained(r, value)
            counters["family_parameter_intervals"] += 1
        for r, value in zip(rec["weighted_forward_intervals"], f.forward):
            assert contained(r, value)
            counters["family_parameter_intervals"] += 1
        assert endpoint(f.EA - f.EG, False) >= 0
        minimum = rec["minimum_certificate"]
        if minimum is not None:
            assert f.classification == "mixed"
            lo, hi = F(minimum["lo"]), F(minimum["hi"])
            assert sign(f.S(ball(lo), 1)) < 0 < sign(f.S(ball(hi), 1))
            assert sign(f.S(interval(lo, hi), 2)) > 0
            assert contained(minimum["S_minimum_interval"], local_range(f.S, lambda z: f.S(z, 1), lo, hi))
            assert contained(minimum["left_derivative_interval"], f.S(ball(lo), 1))
            assert contained(minimum["right_derivative_interval"], f.S(ball(hi), 1))
            history_check(minimum["bisection"], lambda z: f.S(z, 1), counters)
            counters["minimum_brackets"] += 1

    calls = {}
    for rec in producer["call_bounds"]:
        f = families[rec["model_id"]]
        strike = F(rec["strike"])
        assert rec["optimisation_complete"] is True
        root = rec["critical_root"]
        assert root is not None
        lo, hi = F(root["lo"]), F(root["hi"])
        score = lambda z: f.S(z) - ball(strike)
        assert sign(score(ball(lo))) < 0 < sign(score(ball(hi)))
        assert sign(f.S(interval(lo, hi), 1)) > 0
        assert contained(root["left_score"], score(ball(lo)))
        assert contained(root["right_score"], score(ball(hi)))
        if f.classification == "mixed":
            minimum = family_records[rec["model_id"]]["minimum_certificate"]
            assert lo > F(minimum["hi"])
            assert F(minimum["S_minimum_interval"]["upper_exact_fraction"]) < strike
        history_check(root["bisection"], score, counters)
        counters["critical_root_brackets"] += 1
        for z in map(ball, (F(-3), F(-1, 4), F(0), F(1, 2), F(2))):
            direct = ball(strike) * phi(z) - sum(
                (w * phi(a - z) for w, a in zip(f.forward, f.a)), arb(0))
            factored = phi(z) * (ball(strike) - f.S(z))
            assert (direct - factored).contains(0)
            counters["Lprime_identity_checks"] += 1
        for candidate in rec["candidate_values"]:
            kind = candidate["kind"]
            if kind == "positive_part_of_forward_endpoint":
                value = f.EA - ball(strike)
            elif kind == "zero_endpoint":
                value = arb(0)
            elif kind == "fixed_z_zero":
                value = f.L(arb(0), strike)
            elif kind == "upper_critical_root":
                value = local_range(lambda z: f.L(z, strike),
                                    lambda z: phi(z) * (ball(strike) - f.S(z)), lo, hi)
            elif kind == "root_midpoint_fixed_threshold":
                value = f.L(ball(F(candidate["threshold"])), strike)
            else:
                raise ValueError(kind)
            assert contained(candidate["value"], value)
            counters["candidate_enclosures"] += 1
        family_lo = max(F(0), *(F(c["value"]["lower_exact_fraction"])
                               for c in rec["candidate_values"]))
        family_hi = max(F(0), *(F(c["value"]["upper_exact_fraction"])
                               for c in rec["candidate_values"]))
        family_interval = rec["exercise_family_supremum_interval"]
        assert F(family_interval["lower_exact_fraction"]) <= family_lo <= family_hi <= F(family_interval["upper_exact_fraction"])
        geometric = f.geometric(strike)
        AG = geometric + f.EA - f.EG
        assert contained(rec["geometric_call_interval"], geometric)
        assert contained(rec["AG_gap_interval"], f.EA - f.EG)
        assert contained(rec["AG_upper_interval"], AG)
        counters["geometric_AG_intervals"] += 3
        call = rec["call_price_interval"]
        assert F(call["lower_exact_fraction"]) <= family_lo
        assert endpoint(AG, True) <= F(call["upper_exact_fraction"])
        counters["call_intervals"] += 1
        calls[(rec["model_id"], strike)] = call

    for rec in producer["records"]:
        base = calls[(rec["model_id"], F(100))]
        blo, bhi = F(base["lower_exact_fraction"]), F(base["upper_exact_fraction"])
        cap = None if rec["cap"] is None else F(rec["cap"])
        if cap is None:
            lo, hi = blo, bhi
        else:
            highcall = calls[(rec["model_id"], F(100) + cap)]
            alo, ahi = F(highcall["lower_exact_fraction"]), F(highcall["upper_exact_fraction"])
            lo, hi = max(F(0), blo - ahi), min(cap, bhi - alo)
        radius = (hi - lo) / 2
        assert F(rec["undiscounted_radius_upper"]) == radius
        assert (radius <= F(1, 100)) == rec["certifies_undiscounted_cent"]
        target = rec["continuous_price_interval"]
        assert F(target["lower_exact_fraction"]) <= lo <= hi <= F(target["upper_exact_fraction"])
        counters["contract_propagations"] += 1
        model = models[rec["model_id"]]
        discount = ball(-F(model["risk_free_rate"]) * F(model["maturity"])).exp()
        dlo = endpoint(ball(lo) * discount, False)
        dhi = endpoint(ball(hi) * discount, True)
        target = rec["discounted_price_interval"]
        assert F(target["lower_exact_fraction"]) <= dlo <= dhi <= F(target["upper_exact_fraction"])
        counters["discounted_intervals"] += 1
    assert len(producer["records"]) == len(producer["call_bounds"]) == 36
    assert sum(r["certifies_undiscounted_cent"] for r in producer["records"]) == 0
    out = {"schema": "qbp-geometric-conditioning-independent-audit-v1",
           "status": "passed", "precision_bits": PRECISION,
           "family_count": len(families), "strike_count": 36,
           "contract_count": 36, "certified_cent_count": 0,
           "audit_counts": counters,
           "method": "Independent exact covariances and martingale forwards, Gaussian tilting, all root-history signs, endpoint/critical candidates, AG upper and exact rational spread propagation; no production helper",
           "proof_scope": "Mixed-sign S is strictly convex and coercive; its upper crossing is the only finite maximum of L. Endpoint values complete the supremum comparison. Cap is an independently bounded call spread, not a native optimal cap bound.",
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "result_sha256": hashlib.sha256(results.read_bytes()).hexdigest(),
           "producer_sha256": producer["source_sha256"]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "audit_counts": counters}, indent=2))


if __name__ == "__main__":
    main()
