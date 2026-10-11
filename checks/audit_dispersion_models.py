"""Exact independent model audit for the controlled dispersion stress.

The audit checks every scaled primitive and derived field against the
unchanged original protocol. Covariance positivity follows structurally
from the retained factorisation, rather than a floating eigenvalue test.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

ROOT = Path(__file__).resolve().parent


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def exact_equal(a, b):
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(exact_equal(x, y) for x, y in zip(a, b))
    try:
        return F(a) == F(b)
    except (ValueError, TypeError):
        return a == b


def parameters(model):
    w = tuple(map(F, model["weights"]))
    R = tuple(tuple(map(F, row)) for row in model["loadings"])
    mu = tuple(map(F, model["log_drifts"]))
    lam = F(model["lambda"])
    RG = tuple(sum(w[i] * R[i][j] for i in range(len(w))) for j in range(len(R[0])))
    mg = dot(w, mu)
    d = tuple(x - mg for x in mu)
    s = tuple(tuple(x - y for x, y in zip(row, RG)) for row in R)
    Q = sum(w[i] * dot(R[i], R[i]) for i in range(len(w)))
    V = Q - dot(RG, RG)
    return w, R, mu, lam, RG, mg, d, s, Q, V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=Path, default=ROOT.parent)
    ap.add_argument("--out", type=Path, default=ROOT / "dispersion_model_audit.json")
    args = ap.parse_args()
    original_path = args.stage / "results/generalisation_protocol.json"
    protocol_path = args.stage / "results/dispersion_separation_protocol.json"
    original = json.loads(original_path.read_text())
    protocol = json.loads(protocol_path.read_text())
    base = {m["model_id"]: m for m in original["models"]}
    assert len(protocol["models"]) == 36 and len(base) == 9
    counters = {"models": 0, "drift_entries": 0, "covariance_entries": 0,
                "scaled_loading_entries": 0, "centred_loading_entries": 0,
                "derived_dispersion_fields": 0, "weighted_centring_checks": 0,
                "scale_one_original_fields": 0}
    model_checks = []
    for m in protocol["models"]:
        c = F(m["dispersion_scale"])
        assert c in (F(1), F(3, 2), F(2), F(3))
        b = base[m["original_model_id"]]
        bw, bR, bmu, blam, bRG, bmg, bd, bs, bQ, bV = parameters(b)
        w, R, mu, lam, RG, mg, d, s, Q, V = parameters(m)
        assert w == bw and lam == blam > 0 and sum(w) == 1 and min(w) > 0
        assert m["assets"] == b["assets"] and m["outer_dimension"] == b["outer_dimension"]
        carry = F(m["risk_free_rate"]) - F(m["dividend_yield"])
        assert carry == F(b["risk_free_rate"]) - F(b["dividend_yield"])
        assert F(m["maturity"]) == 1 and m["spot"] == m["strike"] == 100
        for i in range(len(w)):
            assert R[i] == tuple(c * v for v in bR[i])
            counters["scaled_loading_entries"] += len(R[i])
            assert mu[i] == carry - (dot(R[i], R[i]) + lam * lam) / 2
            assert d[i] == c * c * bd[i]
            counters["drift_entries"] += 1
            assert s[i] == tuple(c * v for v in bs[i])
            counters["centred_loading_entries"] += len(s[i])
            for j in range(len(w)):
                assert F(m["covariance"][i][j]) == dot(R[i], R[j]) + lam * lam
                counters["covariance_entries"] += 1
        assert RG == tuple(c * v for v in bRG)
        assert mg == bmg - (c * c - 1) * bQ / 2
        assert Q == c * c * bQ and V == c * c * bV > 0
        assert V == sum(wi * dot(si, si) for wi, si in zip(w, s))
        assert F(m["weighted_residual_dispersion"]) == V
        assert F(m["geometric_outer_variance"]) == dot(RG, RG)
        counters["derived_dispersion_fields"] += 2
        assert dot(w, d) == 0
        assert all(sum(w[i] * s[i][j] for i in range(len(w))) == 0 for j in range(len(RG)))
        counters["weighted_centring_checks"] += 1 + len(RG)
        if c == 1:
            for key, value in b.items():
                if key == "model_id":
                    assert m["original_model_id"] == value
                else:
                    assert exact_equal(m[key], value), (m["model_id"], key)
                counters["scale_one_original_fields"] += 1
        counters["models"] += 1
        model_checks.append({"model_id": m["model_id"], "scale": str(c),
                             "mean_asset_forward": "100*exp(r_f-q)",
                             "weighted_residual_dispersion_exact": str(V),
                             "geometric_outer_variance_exact": str(dot(RG, RG)),
                             "effective_geometric_variance_exact": str(dot(RG, RG) + lam * lam),
                             "covariance_PSD_reason": "c^2 RR^T + lambda^2 11^T; exact factorisation",
                             "Jensen_reason": "Positive weights of sum one and exact weighted-centred d,s"})
    out = {"schema": "qbp-dispersion-model-independent-audit-v1", "status": "passed",
           "audit_counts": counters, "models": model_checks,
           "original_protocol_sha256": hashlib.sha256(original_path.read_bytes()).hexdigest(),
           "stress_protocol_sha256": hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "scope": "Exact rational model invariants and structural covariance validity; standard pricing-measure assumption remains stipulated. No empirical sample or pricing outcome used."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "audit_counts": counters}, indent=2))


if __name__ == "__main__":
    main()
