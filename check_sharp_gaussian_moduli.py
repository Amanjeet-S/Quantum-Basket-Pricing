"""Independent interval checks of the sharp universal Gaussian moduli.

The proof fragment supplies the global optimisation and attainment proof.
This script checks scalar roots, directed formula evaluation and strict
comparisons without importing pricing, graph or confidence code.
"""
from fractions import Fraction
from pathlib import Path
import json
import argparse
import hashlib
import sys
sys.dont_write_bytecode=True
from flint import arb, fmpq, ctx

ctx.prec = 320


def aa(value):
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def interval(lo, hi):
    return arb(aa((lo + hi) / 2), aa((hi - lo) / 2))


def phi(value):
    return (-(value * value) / 2).exp() / (2 * arb.pi()).sqrt()


def Phi(value):
    return (1 + (value / arb(2).sqrt()).erf()) / 2


def bisect(fun, lo, hi, bits=120):
    assert fun(aa(lo)) < 0
    assert fun(aa(hi)) > 0
    while hi - lo > Fraction(1, 2**bits):
        mid = (lo + hi) / 2
        result = fun(aa(mid))
        if result < 0:
            lo = mid
        elif result > 0:
            hi = mid
        else:
            raise ArithmeticError("Increase interval precision for an undecided sign")
    assert fun(aa(lo)) < 0
    assert fun(aa(hi)) > 0
    return lo, hi


def k(t):
    t = Fraction(t)
    if not t:
        return arb(0), (Fraction(0), Fraction(0))
    half = t / 2
    lo, hi = bisect(lambda u: u * (aa(half) * u).tanh() - aa(half),
                    max(Fraction(1), half), 1 + half)
    u = interval(lo, hi)
    return phi(u - aa(half)) - phi(u + aa(half)), (lo, hi)


def omega(t, e):
    t, e = Fraction(t), Fraction(e)
    if not t or not e:
        return arb(0), (Fraction(0), Fraction(0))
    b, c = (t + e) / 2, abs(t - e) / 2
    fun = lambda u: (aa(b) * u).cosh().log() - (aa(c) * u).cosh().log() - aa(t * e / 2)
    lo, hi = bisect(fun, max(t, e) / 2, max(t, e) / 2 + 1 / min(t, e))
    u = interval(lo, hi)
    r = lambda x: Phi(x + aa(t / 2)) - Phi(x - aa(t / 2))
    return r(u - aa(e / 2)) - r(u + aa(e / 2)), (lo, hi)


def A(t):
    return 2 * Phi(aa(Fraction(t) / 2)) - 1


def root_record(root):
    return {"lower": str(root[0]), "upper": str(root[1])}


def main(output):
    t_values = [Fraction(1, 100), Fraction(1, 10), Fraction(1, 2),
                Fraction(1), Fraction(3, 2), Fraction(2), Fraction(4), Fraction(10)]
    e_values = [Fraction(1, 1024), Fraction(1, 10), Fraction(1), Fraction(4)]
    rows = []
    for t in t_values:
        kval, kroot = k(t)
        old = min(phi(arb(0)), aa(t) * phi(arb(1)))
        assert kval > 0 and kval < old
        examples = []
        for e in e_values:
            oval, oroot = omega(t, e)
            reverse, _ = omega(e, t)
            ke, _ = k(e)
            assert oval.overlaps(reverse)
            assert oval > 0 and oval < kval * aa(e)
            assert oval < ke * aa(t)
            assert oval < A(t) and oval < A(e)
            assert A(e) < phi(arb(0)) * aa(e)
            examples.append({"e": str(e), "Omega": str(oval),
                             "root": root_record(oroot),
                             "linear_sharp_bound": str(kval * aa(e)),
                             "symmetric_root_check": True})
        rows.append({"t": str(t), "k": str(kval), "root": root_record(kroot),
                     "previous_min_coefficient": str(old),
                     "relative_coefficient_reduction": str(1 - kval / old),
                     "finite_shift_examples": examples})
    example_k, _ = k(Fraction(1))
    assert example_k / aa(Fraction(1, 1000)) > 101
    result = {"status": "passed", "precision_bits": ctx.prec,
              "derivative_coefficients_checked": len(rows),
              "finite_moduli_checked": len(rows) * len(e_values),
              "proof_scope": "Exact scalar roots and strict inequalities; global optimality is proved in the adjacent TeX fragment",
              "financial_price_input": False,
              "generic_novelty_claim": False,
              "cap_example": {"K": 100, "H": 1, "lambda": "1/1000", "d": "1/1000",
                              "universal_kappa": str(example_k / aa(Fraction(1, 1000))),
                              "cap_kappa_certified_upper": 101,
                              "relative_reduction_at_least": str(1 - aa(101) / (example_k / aa(Fraction(1, 1000))))},
              "rows": rows}
    destination=Path(output)
    if destination.exists():raise FileExistsError("Preserve earlier scalar checks: "+str(destination))
    result["checker_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["proof_scope"]="Global optimisation and attainment are proved in proof_sources/sharp_common_information_fragment.tex; finite scalar checks are complementary."
    destination.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "derivative_coefficients": len(rows),
                      "finite_moduli": result["finite_moduli_checked"],
                      "t_one_reduction": rows[3]["relative_coefficient_reduction"],
                      "t_three_halves_reduction": rows[4]["relative_coefficient_reduction"],
                      "cap_example_reduction": result["cap_example"]["relative_reduction_at_least"]}))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    main(parser.parse_args().output)
