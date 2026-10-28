"""Independent small-matrix checks of the exact trace-bound implementation."""

from fractions import Fraction
import importlib.util
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cpp_residual_certificate", ROOT / "src/cpp_residual_certificate.py")
certificate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(certificate)


def exact_step(previous, lam, height):
    n = len(previous)
    a = [[Fraction(1 + 2 * lam) if i == j else (-lam if abs(i-j) == 1 else Fraction(0))
          for j in range(n)] for i in range(n)]
    b = list(previous)
    b[-1] += lam * height
    for i in range(n):
        for j in range(i+1, n):
            factor = a[j][i] / a[i][i]
            for k in range(i, n):
                a[j][k] -= factor * a[i][k]
            b[j] -= factor * b[i]
    answer = [Fraction(0)] * n
    for i in reversed(range(n)):
        answer[i] = (b[i] - sum(a[i][j]*answer[j] for j in range(i+1,n))) / a[i][i]
    return answer


def fraction(item):
    return Fraction(int(item['numerator']), int(item['denominator']))


class ResidualCertificateTests(unittest.TestCase):
    def test_exact_one_node_step(self):
        result = certificate.certify_trace([[0,0,2], [0,.5,2]], .5, 2)
        self.assertEqual(fraction(result['final_grid_error_bound']), 0)

    def test_one_node_perturbation_is_bounded(self):
        result = certificate.certify_trace([[0,0,2], [0,.75,2]], .5, 2)
        self.assertGreaterEqual(fraction(result['final_grid_error_bound']), Fraction(1,4))

    def test_accumulated_bound_matches_independent_rational_solves(self):
        for n in (1,2,5):
            for lam in (Fraction(0),Fraction(1,32),Fraction(1),Fraction(10000)):
                with self.subTest(n=n,lam=lam):
                    height = Fraction(2)
                    target = [Fraction.from_float(float(Fraction(j, n+1))) for j in range(1,n+1)]
                    trace = [[0,*map(float,target),2]]
                    for step in range(3):
                        target = exact_step(target,lam,height)
                        trace.append([0,*[float(v)+(-1)**(j+step)*2**(-10-step) for j,v in enumerate(target)],2])
                    result = certificate.certify_trace(trace,float(lam),2)
                    actual = max(abs(Fraction.from_float(trace[-1][j+1])-v) for j,v in enumerate(target))
                    self.assertLessEqual(actual,fraction(result['final_grid_error_bound']))

    def test_scalar_bound_covers_interpolation_and_evaluation(self):
        trace = [[0,0,2],[0,.75,2]]
        result = certificate.certify_trace(trace,.5,2,price=.4,left_index=0,right_weight=.5,discount=.9)
        target = Fraction.from_float(.9)*Fraction(1,2)*Fraction(1,2)
        self.assertLessEqual(abs(Fraction.from_float(.4)-target),fraction(result['scalar_price_error_bound']))

    def test_overshoots_are_retained(self):
        result = certificate.certify_trace([[0,0,1],[0,2,1]],.25,1)
        self.assertTrue(result['observed_range_violation'])
        self.assertEqual(fraction(result['max_observed_range_excess']),1)
        self.assertTrue(result['range_excess_within_residual_bound'])
        self.assertGreater(fraction(result['final_grid_error_bound']),0)

    def test_outside_initial_range_does_not_imply_a_payoff_range_certificate(self):
        result = certificate.certify_trace([[0,2,1],[0,2,1]],0,1)
        self.assertFalse(result['initial_range_valid'])
        self.assertIsNone(result['range_excess_within_residual_bound'])
        self.assertEqual(fraction(result['final_grid_error_bound']),0)

    def test_upward_conversion_encloses_exact_bound(self):
        for value in (Fraction(0),Fraction(1,3),Fraction(1,10),Fraction(1,2**1075)):
            upper = certificate.upper_binary64(value)
            self.assertGreaterEqual(Fraction.from_float(upper),value)
        self.assertIsNone(certificate.upper_binary64(Fraction(2**2048)))

    def test_invalid_trace_and_coefficients_are_rejected(self):
        cases = [([[0,0,1]],.2,1), ([[0,1],[0,1]],.2,1),
                 ([[0,0,1],[0,0,0]],.2,1), ([[0,0,1],[0,1,1,1]],.2,1),
                 ([[0,0,1],[0,math.nan,1]],.2,1), ([[0,0,1],[0,0,1]],-.2,1),
                 ([[0,0,1],[0,0,1]],math.inf,1)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                certificate.certify_trace(*args)

    def test_partial_or_invalid_scalar_configuration_is_rejected(self):
        for kwargs in ({'price':.1}, {'price':.1,'left_index':True,'right_weight':.5,'discount':1},
                       {'price':.1,'left_index':0,'right_weight':1.1,'discount':1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                certificate.certify_trace([[0,0,1],[0,.1,1]],.1,1,**kwargs)


if __name__ == '__main__':
    unittest.main()
