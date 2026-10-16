"""Independent polynomial reconstruction of the continuous degree-three test.

This audit directly expands (T-1)**k using exact rational polynomial products,
instead of taking the producer's binomial differences of raw T**j moments.
It independently integrates the signed polynomials using Gaussian tilting and
uses no production pricing, parameter, geometry or interval helpers.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
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


def serial(x):
    return {"lower_exact_fraction": str(endpoint(x, False)),
            "upper_exact_fraction": str(endpoint(x, True))}


def contains_record(record, value):
    return (F(record["lower_exact_fraction"]) <= endpoint(value, False)
            <= endpoint(value, True) <= F(record["upper_exact_fraction"]))


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def phi(x):
    return (-x * x / 2).exp() / (2 * arb.pi()).sqrt()


def Phi(x):
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            alpha = tuple(x + y for x, y in zip(a, b))
            out[alpha] = out.get(alpha, F(0)) + ca * cb
    return {a: c for a, c in out.items() if c}


def geometry(model):
    w = tuple(map(F, model["weights"]))
    mu = tuple(map(F, model["log_drifts"]))
    r = tuple(tuple(map(F, row)) for row in model["loadings"])
    lam = F(model["lambda"])
    assert sum(w) == 1 and min(w) > 0 and lam > 0
    assert F(model["maturity"]) == 1 and model["spot"] == 100
    drift = F(model["risk_free_rate"]) - F(model["dividend_yield"])
    assert all(mu[i] + (dot(r[i], r[i]) + lam ** 2) / 2 == drift
               for i in range(len(w)))
    mg = dot(w, mu)
    rg = tuple(sum(w[i] * r[i][j] for i in range(len(w)))
               for j in range(model["outer_dimension"]))
    d = tuple(x - mg for x in mu)
    s = tuple(tuple(x - y for x, y in zip(row, rg)) for row in r)
    assert dot(w, d) == 0
    assert all(sum(w[i] * s[i][j] for i in range(len(w))) == 0
               for j in range(len(rg)))
    zero = (0,) * len(w)
    u = {zero: F(-1)}
    for i, wi in enumerate(w):
        alpha = tuple(int(k == i) for k in range(len(w)))
        u[alpha] = wi
    powers = [{zero: F(1)}]
    for _ in range(4):
        powers.append(mul(powers[-1], u))
    basis = {}
    for poly in powers:
        for alpha in poly:
            if alpha not in basis:
                v = tuple(sum(alpha[i] * s[i][j] for i in range(len(w)))
                          for j in range(len(rg)))
                basis[alpha] = (ball(dot(alpha, d) + dot(v, v) / 2).exp(),
                                dot(rg, v))
    moments = [sum((ball(c) * basis[a][0] for a, c in p.items()), arb(0))
               for p in powers]
    return mg, rg, lam, powers, basis, moments


def call_kernel(m, sigma, strike, order):
    strike = ball(strike)
    ell = (strike / 100).log()
    z = (m - ell) / sigma
    growth = 100 * (m + sigma * sigma / 2).exp()
    if order == 0:
        return growth * Phi(z + sigma) - strike * Phi(z)
    if order == 1:
        return growth * Phi(z + sigma)
    boundary = strike * phi(z) / sigma
    if order == 2:
        return boundary
    if order == 3:
        return boundary * (-z / sigma - 2)
    raise ValueError("Unsupported order")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path,
                    default=Path("/private/tmp/qbp-residual-priority-step/cost_boundary"))
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "continuous_degree3_independent_audit.json")
    args = ap.parse_args()
    ctx.prec = PRECISION
    stage = args.stage
    protocol = stage / "results/generalisation_protocol.json"
    result_file = stage / "results/residual_same_information_taylor.json"
    producer = json.loads(result_file.read_text())
    models = {m["model_id"]: m for m in json.loads(protocol.read_text())["models"]}
    norms = json.loads((stage / "results/boundary_derivative_norms.json").read_text())["records"]
    geometry_records = {r["assets"]: r for r in producer["geometry_evidence"]}
    assert producer["source_sha256"] == hashlib.sha256(
        (stage / "residual_same_information_taylor.py").read_bytes()).hexdigest()
    assert producer["protocol_sha256"] == hashlib.sha256(protocol.read_bytes()).hexdigest()
    assert producer["norm_input_sha256"] == hashlib.sha256(
        (stage / "results/boundary_derivative_norms.json").read_bytes()).hexdigest()
    for name, digest in producer["dependencies"].items():
        assert digest == hashlib.sha256((stage / name).read_bytes()).hexdigest()
    counters = {"model_premises": 0, "continuous_moment_enclosures": 0,
                "coefficient_enclosures": 0, "exact_remainder_arithmetic": 0,
                "price_interval_containment": 0, "discounted_interval_containment": 0,
                "exact_constructed_halfwidth": 0,
                "cent_status_and_first_success": 0}
    cache = {}
    successes = []
    for record in producer["records"]:
        model = models[record["model_id"]]
        if model["model_id"] not in cache:
            cache[model["model_id"]] = geometry(model)
            counters["model_premises"] += 1
        mg, rg, lam, powers, basis, moments = cache[model["model_id"]]
        for k in range(5):
            assert contains_record(geometry_records[model["assets"]]["continuous_U_moments"][str(k)], moments[k])
            counters["continuous_moment_enclosures"] += 1
        cap = None if record["cap"] is None else F(record["cap"])
        norm = next(n for n in norms if n["payoff"] == record["payoff"]
                    and F(n["lambda"]) == lam
                    and (None if n["cap"] is None else F(n["cap"])) == cap)
        sigma = ball(lam ** 2 + dot(rg, rg)).sqrt()
        def payoff(m, order):
            value = call_kernel(m, sigma, 100, order)
            if cap is not None:
                value -= call_kernel(m, sigma, 100 + cap, order)
            return value
        coefficients = [payoff(ball(mg), 0)]
        for k in (1, 2, 3):
            coefficients.append(sum((ball(c) * basis[a][0] *
                payoff(ball(mg + basis[a][1]), k)
                for a, c in powers[k].items()), arb(0)) / factorial(k))
        assert contains_record(record["c0_c1_input_certificate"]["c0"], coefficients[0])
        counters["coefficient_enclosures"] += 1
        assert contains_record(record["c0_c1_input_certificate"]["c1"], coefficients[1])
        counters["coefficient_enclosures"] += 1
        partial = coefficients[0]
        first = None
        for order in record["orders"]:
            k = order["order"]
            assert contains_record(order["coefficient_interval"], coefficients[k])
            counters["coefficient_enclosures"] += 1
            partial += coefficients[k]
            assert contains_record(order["partial_sum_interval"], partial)
            M = F(norm["C_r_upper"][str(k + 1)])
            mu_upper = F(order["continuous_ratio_moment_interval"]["upper_exact_fraction"])
            radius = M * mu_upper / factorial(k + 1)
            assert F(order["remainder_M_upper"]) == M
            assert F(order["remainder_upper"]) == radius
            counters["exact_remainder_arithmetic"] += 1
            lo = max(F(0), endpoint(partial, False) -
                     (F(0) if cap is None and k == 1 else radius))
            hi = endpoint(partial, True) + radius
            if cap is not None:
                hi = min(cap, hi)
            target = order["continuous_price_interval"]
            assert F(target["lower_exact_fraction"]) <= lo <= hi <= F(target["upper_exact_fraction"])
            counters["price_interval_containment"] += 1
            halfwidth = F(order["undiscounted_radius_upper"])
            producer_low = max(F(0), F(order["partial_sum_interval"]["lower_exact_fraction"]) -
                               (F(0) if cap is None and k == 1 else radius))
            producer_high = F(order["partial_sum_interval"]["upper_exact_fraction"]) + radius
            if cap is not None:
                producer_high = min(cap, producer_high)
            assert halfwidth == (producer_high - producer_low) / 2
            counters["exact_constructed_halfwidth"] += 1
            certified = halfwidth <= F(1, 100)
            assert certified == order["certifies_undiscounted_cent"]
            if certified and first is None:
                first = k
            discount = ball(-F(model["risk_free_rate"]) * F(model["maturity"])).exp()
            # The producer's price interval is itself saved with an extra Arb
            # serialization radius. Checking its exact raw construction would
            # need those original endpoints; use the independent narrower range.
            disc_lo = endpoint(ball(lo) * discount, False)
            disc_hi = endpoint(ball(hi) * discount, True)
            disc_target = order["discounted_price_interval"]
            assert F(disc_target["lower_exact_fraction"]) <= disc_lo <= disc_hi <= F(disc_target["upper_exact_fraction"])
            counters["discounted_interval_containment"] += 1
        assert first == record["first_successful_order"]
        counters["cent_status_and_first_success"] += 1
        if first is not None:
            successes.append([record["model_id"], record["cap"], first])
    assert len(producer["records"]) == 36 and len(successes) == 33
    out = {"schema": "qbp-degree3-independent-polynomial-audit-v1",
           "status": "passed", "precision_bits": PRECISION,
           "record_count": len(producer["records"]), "successes": len(successes),
           "audit_counts": counters,
           "method": "Direct exact rational polynomial multiplication of (T-1)^k, full Gaussian transforms and independently coded convolved call boundaries; no production helpers",
           "mathematical_scope": "The finite Taylor remainder requires U>=0, positive conditional smoothing and the recorded certified C_(k+1) bound. The audit verifies coefficients and propagated intervals, not literature priority.",
           "result_sha256": hashlib.sha256(result_file.read_bytes()).hexdigest(),
           "producer_sha256": producer["source_sha256"],
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "successes": len(successes),
                      "audit_counts": counters}, indent=2))


if __name__ == "__main__":
    main()
