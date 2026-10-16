"""Boundary regressions for the exact finite classical comparator."""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import classical_engine as ce
from encoding_resources import Graph


class ExactClassicalTests(unittest.TestCase):
    def test_product_above_uint64(self):
        q,w = 30,37
        a = np.asarray([0,(1 << w)-1,1 << (w-1)],np.uint64)
        b = np.asarray([(1 << w)-1,(1 << w)-1,(1 << w)-1],np.uint64)
        expected = [((int(x)*int(y)) >> q)&((1 << w)-1) for x,y in zip(a,b)]
        self.assertEqual(ce.fixed_multiply_u64(a,b,q,w).tolist(),expected)

    def test_full_unsigned_product_windows(self):
        a = np.asarray([(1 << 64)-1,1 << 63,(1 << 32)+1],np.uint64)
        b = np.asarray([(1 << 64)-1,(1 << 64)-1,(1 << 63)+1],np.uint64)
        for drop in (0,30,63,64,96,127):
            actual = ce.multiply_shift_u64(a,b,drop,37)
            expected = [((int(x)*int(y)) >> drop)&((1 << 37)-1) for x,y in zip(a,b)]
            self.assertEqual(actual.tolist(),expected)

    def test_signed_modular_clamp(self):
        q,w = 18,25
        g = Graph(q);g.width = w
        x = g.inp('x',w);g.output = g.clamp(g.sub(x,g.const(1)),1)
        inputs = np.asarray([0,(1 << q)-1,1 << q,(1 << q)+1,2 << q,(1 << (w-1))+1],np.uint64)
        engine = ce.FiniteEvaluator(g,{})
        actual = engine.evaluate({'x':inputs})
        expected = [g.evaluate({'x':int(v)},{}) for v in inputs]
        self.assertEqual(actual.tolist(),expected)
        self.assertEqual(actual[:3].tolist(),[0,0,0])

    def test_interpolation_fallback_above_uint64(self):
        q,w,bits = 30,37,4
        g = Graph(q);g.width = w
        x = g.inp('x',q+3);g.output = g.function('general_exp',x,bits)
        mask = (1 << w)-1
        table = [(j&mask)|(mask << w) for j in range(1 << bits)]
        engine = ce.FiniteEvaluator(g,{('general_exp',bits,q):table},require_normalised_output=False)
        inputs = np.asarray([0,(1 << 29)-1,1 << 29,(8 << q)-1],np.uint64)
        self.assertFalse(engine.plan[0]['native_product_safe'])
        actual = engine.evaluate({'x':inputs})
        expected = [g.evaluate({'x':int(v)},{('general_exp',bits,q):table}) for v in inputs]
        self.assertEqual(actual.tolist(),expected)

    def test_squared_sum_exceeds_uint64(self):
        a = np.full(20,1 << 30,np.uint64)
        moments = ce.ExactMoments(30).add(a)
        self.assertGreater(moments.total_squares,1 << 64)
        self.assertEqual(moments.total_squares,20*(1 << 60))
        self.assertEqual(moments.as_dict()['normalised_sample_variance_exact_fraction'],'0')


if __name__ == '__main__':unittest.main()
