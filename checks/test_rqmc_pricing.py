"""Exact digital-shift and integer-reduction contract checks."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import rqmc_pricing as engine


class RandomisedDigitalShiftTests(unittest.TestCase):
    def test_author_direction_fixture(self):
        expected=np.array([[0,0,0],[.5,.5,.5],[.75,.25,.25],[.25,.75,.75],
            [.375,.375,.625],[.875,.875,.125],[.625,.125,.875],
            [.125,.625,.375],[.1875,.3125,.9375],[.6875,.8125,.4375]])
        actual=engine.sobol_words(4,3)[:10]
        self.assertTrue(np.array_equal(actual,expected*(1<<32)))

    def test_every_net_point_uniform_under_all_finite_shifts(self):
        points=engine.digital_points({'x':2,'y':3},net_bits=2)
        shifts={'x':np.repeat(np.arange(4,dtype=np.uint64),8),
                'y':np.tile(np.arange(8,dtype=np.uint64),4)}
        shifted=engine.shifted_inputs(points,shifts)
        for point in range(4):
            outcomes={(int(shifted['x'][j*4+point]),int(shifted['y'][j*4+point])) for j in range(32)}
            self.assertEqual(outcomes,{(x,y) for x in range(4) for y in range(8)})

    def test_dependent_net_points_do_not_become_independent_samples(self):
        points=engine.digital_points({'x':2},net_bits=2)['x']
        self.assertEqual(set(map(int,points)),set(range(4)))
        sums=[]
        for shift in range(4):
            inputs=points^np.uint64(shift)
            sums.append(sum(map(int,inputs)))
        self.assertEqual(sums,[6]*4)

    def test_uint64_avoids_all_batch_sum_and_square_overflows(self):
        maximum=1<<36
        values=np.array([maximum,maximum-1,(1<<35)+17,65535,65536,0],dtype=np.uint64)
        self.assertEqual(engine.exact_u64_square_sum(values),sum(int(x)**2 for x in values))
        many=np.full(1<<16,maximum,dtype=np.uint64)
        self.assertEqual(engine.exact_u64_square_sum(many),(1<<16)*maximum**2)

    def test_exact_statistics_include_cross_terms_and_disjoint_pairs(self):
        values=[0,1,7,2,11,13]
        d=64*(1<<20)
        actual=engine.batch_statistics(values,20)
        self.assertEqual(int(actual['integer_batch_sum']),sum(values))
        self.assertEqual(int(actual['integer_batch_square_sum']),sum(x*x for x in values))
        self.assertEqual(int(actual['integer_paired_difference_square_sum']),1+25+4)
        self.assertEqual(Fraction(actual['paired_variance_estimate_exact_fraction']),Fraction(30,6*d*d))

    def test_prospective_variance_bound_all_two_point_laws(self):
        rho=Fraction(5,8)
        for numerator in range(0,65):
            probability=Fraction(numerator,64)
            mean=probability*rho; second=probability*rho*rho
            cert={'finite_lattice_range_fraction':str(rho),'fair_mean_upper_fraction':str(mean),
                  'fair_second_moment_upper_fraction':str(second)}
            _,_,_,upper=engine.pure_variance_bound(cert)
            self.assertGreaterEqual(upper,second-mean*mean)

    def test_prespecified_confidence_allocation(self):
        self.assertEqual(len(engine.CHECKPOINTS),11)
        self.assertEqual(2*len(engine.CHECKPOINTS)*engine.ALPHA_EVENT,Fraction(1,100))
        self.assertEqual(engine.CHECKPOINTS[0],64)
        self.assertEqual(engine.CHECKPOINTS[-1],65536)


if __name__=='__main__':unittest.main()
