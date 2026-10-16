"""Independent polynomial and root certificates for the stronger remainder."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import json
import os
import sys
import unittest

HERE=Path(__file__).resolve().parent
AUDITOR=HERE/'audit_gaussian_derivative_norm.py'
if not AUDITOR.exists():AUDITOR=HERE.parent/'audit_gaussian_derivative_norm.py'
spec=importlib.util.spec_from_file_location('independent_norm_audit',AUDITOR)
independent=importlib.util.module_from_spec(spec);spec.loader.exec_module(independent)
PRODUCER=Path(os.environ.get('QBP_NORM_PRODUCER',str(HERE.parent/'gaussian_derivative_norm.py')))
if not PRODUCER.exists():PRODUCER=Path('/private/tmp/qbp-gaussian-norm-step/gaussian_derivative_norm.py')
pspec=importlib.util.spec_from_file_location('norm_producer_under_test',PRODUCER)
producer=importlib.util.module_from_spec(pspec);pspec.loader.exec_module(producer)


class GaussianDerivativeNormTests(unittest.TestCase):
    def test_closed_coefficients_and_gaussian_first_derivative(self):
        self.assertEqual(independent.coefficients(independent.independent_norm_polynomial(1,F(1))), (F(1),))
        self.assertEqual(independent.coefficients(independent.independent_norm_polynomial(2,F(1))), (F(-1),F(-1)))
        self.assertEqual(independent.coefficients(independent.independent_norm_polynomial(3,F(1))), (F(1),F(3),F(1)))
        self.assertEqual(independent.coefficients(independent.independent_norm_polynomial(4,F(1))), (F(0),F(-8),F(-6),F(-1)))

    def test_distinct_exact_polynomial_constructions(self):
        for lam in (F(1,20),F(1,5),F(2,5)):
            for r in range(1,13):
                self.assertEqual(independent.coefficients(independent.independent_norm_polynomial(r,lam)),producer.norm_polynomial(r,lam))

    def test_rational_and_multiple_stationary_roots(self):
        cases=independent.custom_root_cases(PRODUCER)
        self.assertEqual(sum(x['distinct_real_roots'] for x in cases),4)
        self.assertEqual(cases[1]['squarefree_multiplicities'],[3])

    def test_retained_complete_norm_profiles(self):
        data=json.loads((PRODUCER.parent/'results/gaussian_derivative_norms.json').read_text())
        checked=[]
        for row in data['records']:
            p=independent.independent_norm_polynomial(row['derivative_order'],F(row['lambda']))
            checked.append(independent.audit_profile(row,p))
        for row in data['unit_norms']:
            j=row['unit_derivative_order']
            checked.append(independent.audit_profile(row,independent.independent_hermite(j)*((-1)**j),'unit_Hermite'))
        self.assertEqual(len(checked),36)
        self.assertEqual(sum(row['real_stationary_points'] for row in checked),180)
        self.assertTrue(all(row['all_Arb_root_boxes_contain_exact_rational_endpoints'] for row in checked))


if __name__=='__main__':unittest.main()
