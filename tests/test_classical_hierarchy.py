"""Tests of financially relevant early-screen and numerical-fallback behaviour."""
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch
import json
import os
import sys
import unittest
from flint import arb,ctx

ROOT=Path(os.environ.get('QBP_AUDIT_ROOT',str(Path(__file__).resolve().parents[1])))
sys.path.insert(0,str(ROOT))
import taylor_order_benchmark as tb
from taylor_pricing import model_parameters


class HierarchyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.models=json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']

    def setUp(self):
        self.saved_precision=ctx.prec;ctx.prec=224

    def tearDown(self):
        ctx.prec=self.saved_precision

    def model(self,n,lam):
        return next(m for m in self.models if m['assets']==n and F(m['lambda'])==F(lam))

    def bank(self,model):
        w,d,s,mg,rg,lam=model_parameters(model)
        return tb.Geometry((w,d,s,rg))

    def test_accepted_control_never_evaluates_taylor_coefficients(self):
        for n in (2,4,8):
            model=self.model(n,F(2,5));bank=self.bank(model)
            with patch.object(tb,'coefficients',side_effect=AssertionError('Taylor coefficients must not be paid for after an accepted control certificate')):
                result=tb.evaluate(model,1,'monotone_then_screen',bank)
            self.assertTrue(result['selected']['one_sided_monotonicity'])
            self.assertTrue(result['selected']['certified'])
            self.assertEqual(result['kernel_calls'],1)
            self.assertEqual(result['new_geometry_exponentials'],n+1)
            self.assertEqual(result['screens'],[])

    def test_failed_width_screen_skips_the_scalar_control_call(self):
        model=self.model(2,F(1,20));bank=self.bank(model)
        original=tb.kernel_derivatives;orders=[]
        def observed(*args):
            orders.append(args[-1]);return original(*args)
        with patch.object(tb,'kernel_derivatives',side_effect=observed):
            result=tb.evaluate(model,1,'monotone_then_screen',bank)
        self.assertIsNone(result['one_sided_preliminary_interval'])
        self.assertGreater(F(result['one_sided_screen_upper_exact_fraction']),F(1,50))
        self.assertNotIn(0,orders)
        self.assertEqual(len(orders),result['kernel_calls'])
        self.assertTrue(result['selected']['certified'])

    def test_valid_wide_control_enclosure_falls_back_and_is_charged(self):
        model=self.model(2,F(2,5));bank=self.bank(model)
        original=tb.kernel_derivatives;orders=[]
        def widened(*args):
            orders.append(args[-1]);values=original(*args)
            if args[-1]==0:values=[values[0]+arb(0,'0.1')]
            return values
        with patch.object(tb,'kernel_derivatives',side_effect=widened):
            result=tb.evaluate(model,1,'monotone_then_screen',bank)
        preliminary=result['one_sided_preliminary_interval']
        self.assertIsNotNone(preliminary)
        self.assertFalse(preliminary['certified'])
        self.assertLessEqual(F(result['one_sided_screen_upper_exact_fraction']),F(1,50))
        self.assertTrue(result['selected']['certified'])
        self.assertFalse(result['selected'].get('one_sided_monotonicity',False))
        self.assertEqual(orders.count(0),1)
        self.assertEqual(result['kernel_calls'],len(orders))
        self.assertEqual(result['kernel_calls'],4)
        self.assertEqual(result['new_geometry_exponentials'],6)

    def test_order_limit_retains_an_unresolved_price_enclosure(self):
        model=self.model(2,F(1,20));bank=self.bank(model)
        result=tb.evaluate(model,10,'monotone_then_screen',bank)
        self.assertFalse(result['selected']['certified'])
        self.assertEqual(result['selected']['order'],8)
        self.assertEqual([o['order']for o in result['evaluated_orders']],[8])
        self.assertEqual(result['kernel_calls'],45)
        self.assertEqual(result['new_geometry_exponentials'],55)
        self.assertIsNone(result['one_sided_preliminary_interval'])


if __name__=='__main__':unittest.main()
