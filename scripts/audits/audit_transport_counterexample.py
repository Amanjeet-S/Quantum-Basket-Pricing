"""Independent, outward audit of a non-uniform transport-order example.

The functions here are generic C2 test functions, not financial residuals.
The certificate uses exact primitive calculations and independent midpoint
quadrature with analytic Lipschitz remainders. No production helper is used.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json
from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
PRECISION = 512
MIDPOINTS = 4096


def ball(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def endpoint(x, top):
    y = x.upper() if top else x.lower()
    mantissa, exponent = map(int, y.man_exp())
    return (F(mantissa * (1 << exponent)) if exponent >= 0
            else F(mantissa, 1 << -exponent))


def serial(x):
    if not x.is_finite():
        raise ArithmeticError("Non-finite enclosure")
    return {"lower_exact_fraction": str(endpoint(x, False)),
            "upper_exact_fraction": str(endpoint(x, True))}


def interval(lo, hi):
    lo, hi = F(lo), F(hi)
    return ball((lo + hi) / 2) + arb(0, ball((hi - lo) / 2).upper())


def padded(x, error):
    top = endpoint(error, True)
    return interval(endpoint(x, False) - top, endpoint(x, True) + top)


def contains(x, y):
    return endpoint(x, False) <= endpoint(y, False) <= endpoint(y, True) <= endpoint(x, True)


def phi(x):
    return (-x * x / 2).exp() / (2 * arb.pi()).sqrt()


def Phi(x):
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def midpoint(func, lipschitz):
    total = sum((func(ball(F(2 * i + 1, 2 * MIDPOINTS)))
                 for i in range(MIDPOINTS)), arb(0)) / MIDPOINTS
    # On each interval: integral |x-mid| = width**2/4.
    return padded(total, lipschitz / (4 * MIDPOINTS))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "transport_counterexample_audit.json")
    args = ap.parse_args()
    ctx.prec = PRECISION
    pi = arb.pi()
    k = 2 * pi
    # Integral primitive for (1-z) cos(kz), evaluated at exact quarter
    # points using exact sin/cos values at 0, pi/2, 3pi/2 and 2pi.
    g0 = -1 / (k * k)
    gquarter = arb(3) / (4 * k)
    gthreequarter = -1 / (4 * k)
    gone = -1 / (k * k)
    ftc_primitive = k * (gquarter - g0 - gthreequarter + gquarter
                         + gone - gthreequarter)
    # Integral primitive for (1-z)**2 sin(kz), at 0, 1/2 and 1.
    h0 = -1 / k + 2 / (k ** 3)
    hhalf = 1 / (4 * k) - 2 / (k ** 3)
    hone = 2 / (k ** 3)
    peano_primitive = k * k / 2 * (2 * hhalf - h0 - hone)
    ftc_formula = arb(2)
    peano_formula = 3 * pi / 2 - 2 / pi
    if not (ftc_primitive - ftc_formula).contains(0):
        raise ArithmeticError("FTC primitive discrepancy")
    if not (peano_primitive - peano_formula).contains(0):
        raise ArithmeticError("Peano primitive discrepancy")

    ftc_quad = k * midpoint(lambda x: (1 - x) * abs((k * x).cos()), 1 + k)
    peano_quad = k * k / 2 * midpoint(
        lambda x: (1 - x) ** 2 * abs((k * x).sin()), 2 + k)
    if not contains(ftc_quad, ftc_formula) or not contains(peano_quad, peano_formula):
        raise ArithmeticError("Uniform independent quadrature discrepancy")

    h = ball(F(1, 10))
    ph0, phh, cdfh = phi(arb(0)), phi(h), Phi(h)
    gauss_ftc = 2 * k * midpoint(
        lambda x: abs((k * x).cos()) * (cdfh - Phi(h * x)),
        h * ph0 * (1 + k))
    gauss_peano = 2 * k * k / h * midpoint(
        lambda x: abs((k * x).sin()) *
                  (phi(h * x) - phh + h * x * (Phi(h * x) - cdfh)),
        h * h * ph0 * (1 + k / 2))
    ftc_upper = 4 * h * ph0
    peano_lower = 2 * h * phh * peano_formula
    analytic_ratio_lower = (-h * h / 2).exp() * (3 * pi / 4 - 1 / pi)
    if endpoint(analytic_ratio_lower, False) <= 2:
        raise ArithmeticError("Analytic ratio failed to exceed two")
    if endpoint(gauss_peano, False) <= endpoint(gauss_ftc, True):
        raise ArithmeticError("Independent Gaussian quadrature failed strict ordering")
    # The primitive-quadrature intervals should be compatible with the
    # coarser analytic density bounds. They need not be contained in them.
    if endpoint(gauss_ftc, False) > endpoint(ftc_upper, True):
        raise ArithmeticError("FTC density bound inconsistency")
    if endpoint(gauss_peano, True) < endpoint(peano_lower, False):
        raise ArithmeticError("Peano density bound inconsistency")

    out = {"schema": "qbp-generic-transport-order-counterexample-audit-v1",
           "status": "passed", "precision_bits": PRECISION,
           "midpoints_per_integral": MIDPOINTS,
           "uniform_FTC_exact": "2",
           "uniform_Peano_exact": "3*pi/2-2/pi",
           "uniform_FTC_primitive_interval": serial(ftc_primitive),
           "uniform_Peano_primitive_interval": serial(peano_primitive),
           "uniform_FTC_independent_quadrature": serial(ftc_quad),
           "uniform_Peano_independent_quadrature": serial(peano_quad),
           "Gaussian_h": "1/10", "Gaussian_f": "sin(2*pi*y/h)",
           "Gaussian_node": "0", "Gaussian_centred_first_moment": "0 by symmetry",
           "Gaussian_FTC_independent_quadrature": serial(gauss_ftc),
           "Gaussian_Peano_independent_quadrature": serial(gauss_peano),
           "Gaussian_FTC_analytic_upper": serial(ftc_upper),
           "Gaussian_Peano_analytic_lower": serial(peano_lower),
           "Gaussian_ratio_analytic_lower": serial(analytic_ratio_lower),
           "quadrature_remainder": "Lipschitz_constant/(4*midpoints) on [0,1]; positive prefactor applied afterwards",
           "scope": "Generic smooth functions. This rejects a uniform transport-order theorem; it does not exhibit a frozen financial-residual failure.",
           "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"],
                      "uniform_FTC": float(ftc_formula),
                      "uniform_Peano": float(peano_formula),
                      "Gaussian_FTC_quadrature": str(gauss_ftc),
                      "Gaussian_Peano_quadrature": str(gauss_peano),
                      "Gaussian_ratio_analytic_lower": str(analytic_ratio_lower)}, indent=2))


if __name__ == "__main__":
    main()
