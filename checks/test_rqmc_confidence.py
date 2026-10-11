"""Distribution-sensitive contract checks for the shift-batch confidence rule."""
from decimal import Decimal,localcontext
from fractions import Fraction as F
import importlib.util
import math
import os
from pathlib import Path
import unittest

MODULE=Path(os.environ.get('RQMC_CONFIDENCE_SOURCE',str(Path(__file__).resolve().parents[1]/'rqmc_confidence.py')))
spec=importlib.util.spec_from_file_location('confidence_under_review',MODULE)
confidence=importlib.util.module_from_spec(spec);spec.loader.exec_module(confidence)


def statistics(values,denominator=1):
    R=len(values);return {'replicates':R,'batch_denominator':str(denominator),
        'integer_batch_sum':str(sum(values)),'integer_batch_square_sum':str(sum(x*x for x in values)),
        'independent_pairs':R//2,'integer_paired_difference_square_sum':str(sum((values[2*i]-values[2*i+1])**2 for i in range(R//2)))}


def result(values,denominator=1,**kwargs):
    return confidence.confidence_from_batch_statistics(statistics(values,denominator),rho=kwargs.pop('rho',F(1)),
        alpha_variance=F(1,2200),alpha_mean=F(1,2200),**kwargs)


class ShiftBatchConfidenceTests(unittest.TestCase):
    def test_exact_sqrt_encloses_tiny_near_square_and_large_inputs(self):
        step=F(1,1<<confidence.ROOT_BITS)
        values=[F(0),F(1,1<<400),F(1,3),F(1),F(1)+F(1,1<<400),F(1<<512),F(1<<512)+F(1,3)]
        for a in values:
            upper=confidence.sqrt_upper(a)
            self.assertGreaterEqual(upper*upper,a)
            if upper:
                self.assertLess((upper-step)*(upper-step),a)
        with self.assertRaises(ValueError):confidence.sqrt_upper(F(-1,3))

    def test_logarithm_enclosure_has_checked_absolute_margin(self):
        for a in [F(1)+F(1,1<<500),F(2200),F(4400),F(10**40),F(10**40)-F(1,7)]:
            with localcontext() as ctx:
                ctx.prec=300;true=(Decimal(a.numerator)/Decimal(a.denominator)).ln()
                reported=F(confidence.log_upper(a));approx=Decimal(reported.numerator)/Decimal(reported.denominator)
                self.assertGreater(approx,true)
                self.assertLess(approx-true,Decimal('1.00000000000000000001e-120'))
        for a in (F(1),F(0),F(10**40)+1):
            with self.assertRaises(ValueError):confidence.log_upper(a)

    def test_scale_equivariance_for_exact_binary_ranges(self):
        values=[0,0,1,1,0,1,1,0]*8
        first=result(values)
        second=result(values,8,rho=F(1,8))
        for field in ('sample_mean_fraction','mean_lower_fraction','mean_upper_fraction','halfwidth_upper_fraction'):
            # Rational roots have an absolute 2^-160 increment, so they need
            # not be exactly homogeneous. The discrepancy is enclosed.
            self.assertLessEqual(abs(F(first[field])/8-F(second[field])),F(4,1<<confidence.ROOT_BITS))
        self.assertEqual(F(second['paired_variance_observation_fraction']),F(first['paired_variance_observation_fraction'])/64)

    def test_known_deterministic_variance_prior_is_used_without_dividing_by_net_size(self):
        values=[0,1]*32
        first=result(values,variance_prior=F(1,16))
        self.assertEqual(F(first['variance_upper_fraction']),F(1,16))
        self.assertEqual(F(first['variance_prior_upper_fraction']),F(1,16))
        self.assertEqual(first['replicates'],64)
        self.assertEqual(first['pairs'],32)

    def test_rho_zero_gives_exact_zero_interval(self):
        record=result([0]*64,rho=0,variance_prior=0,mean_prior=0)
        self.assertEqual(record['status'],'nonempty')
        self.assertEqual(F(record['mean_lower_fraction']),0)
        self.assertEqual(F(record['mean_upper_fraction']),0)
        self.assertEqual(F(record['halfwidth_upper_fraction']),0)

    def test_empty_mean_intersection_is_explicit_failure(self):
        record=result([1]*1024,mean_prior=F(1,100),variance_prior=F(99,10000))
        self.assertEqual(record['status'],'empty_intersection')
        self.assertIsNone(record['mean_lower_fraction'])
        self.assertIsNone(record['mean_upper_fraction'])
        self.assertIsNone(record['halfwidth_upper_fraction'])
        self.assertNotIn('estimate_fraction',record)

    def test_odd_number_of_batches_uses_only_disjoint_complete_pairs(self):
        values=[0,1]*32+[1]
        record=result(values)
        self.assertEqual(record['replicates'],65)
        self.assertEqual(record['pairs'],32)
        self.assertEqual(F(record['sample_mean_fraction']),F(33,65))
        self.assertEqual(F(record['paired_variance_observation_fraction']),F(1,2))

    def test_exact_large_integer_statistics_do_not_overflow(self):
        D=1<<200;values=[D,0]*32
        record=result(values,D)
        self.assertEqual(F(record['sample_mean_fraction']),F(1,2))
        self.assertEqual(F(record['paired_variance_observation_fraction']),F(1,2))
        self.assertEqual(F(record['variance_upper_fraction']),F(1,4))

    def test_invalid_range_statistics_and_priors_are_rejected(self):
        good=statistics([0,1]*32)
        mutations=[('replicates',1),('batch_denominator',0),('independent_pairs',33),
            ('integer_batch_sum',-1),('integer_batch_square_sum',0),
            ('integer_paired_difference_square_sum',65),('integer_batch_sum',65)]
        for key,value in mutations:
            changed={**good,key:value}
            with self.assertRaises(ValueError,msg=(key,value)):
                confidence.confidence_from_batch_statistics(changed,rho=1,alpha_variance=F(1,2200),alpha_mean=F(1,2200))
        for kwargs in ({'rho':F(5,4)},{'rho':F(-1)},{'variance_prior':F(1,3)},
                       {'variance_prior':-1},{'mean_prior':F(5,4)},{'mean_prior':-1}):
            with self.assertRaises(ValueError):result([0,1]*32,**kwargs)

    def test_non_integral_and_float_counters_are_rejected_without_truncation(self):
        good=statistics([0,1]*32)
        for key,value in good.items():
            for malformed in (F(int(value)*2+1,2),float(value),str(value)+'.0'):
                changed={**good,key:malformed}
                with self.assertRaises(ValueError,msg=(key,malformed)):
                    confidence.confidence_from_batch_statistics(changed,rho=1,alpha_variance=F(1,2200),alpha_mean=F(1,2200))

    def test_boolean_counters_are_rejected_in_every_field(self):
        good=statistics([0,1]*32)
        for key in good:
            for boolean in (False,True):
                with self.assertRaises(ValueError,msg=(key,boolean)):
                    confidence.confidence_from_batch_statistics({**good,key:boolean},rho=1,alpha_variance=F(1,2200),alpha_mean=F(1,2200))

    def test_impossible_aggregate_moments_are_rejected(self):
        # Each tuple passes the older mean/range/Cauchy checks, but cannot
        # arise from nonnegative observations in its declared range.
        too_large_second={'replicates':64,'batch_denominator':10,'integer_batch_sum':64,
            'integer_batch_square_sum':400,'independent_pairs':32,
            'integer_paired_difference_square_sum':0}
        too_large_pairs={**too_large_second,'integer_batch_square_sum':128,
            'integer_paired_difference_square_sum':129}
        for statistics_record,rho in ((too_large_second,F(1,2)),(too_large_pairs,F(1))):
            R=statistics_record['replicates'];D=statistics_record['batch_denominator']
            S=statistics_record['integer_batch_sum'];S2=statistics_record['integer_batch_square_sum']
            Q=statistics_record['integer_paired_difference_square_sum'];m=R//2
            self.assertLessEqual(F(S,R*D),rho)
            self.assertLessEqual(S*S,R*S2)
            self.assertLessEqual(F(S2,R*D*D),rho*rho)
            self.assertLessEqual(F(Q,m*D*D),rho*rho)
            with self.assertRaises(ValueError):
                confidence.confidence_from_batch_statistics(statistics_record,rho=rho,alpha_variance=F(1,2200),alpha_mean=F(1,2200))

    def test_digit_strings_preserve_exact_valid_statistics(self):
        values=[0,1]*32
        good=statistics(values)
        integer_result=confidence.confidence_from_batch_statistics(good,rho=1,alpha_variance=F(1,2200),alpha_mean=F(1,2200))
        string_result=confidence.confidence_from_batch_statistics({k:str(v) for k,v in good.items()},rho=1,alpha_variance=F(1,2200),alpha_mean=F(1,2200))
        self.assertEqual(integer_result,string_result)

    def test_prescribed_event_allocation_sums_exactly(self):
        checkpoints=[1<<k for k in range(6,17)]
        self.assertEqual(len(checkpoints),11)
        self.assertEqual(2*len(checkpoints)*F(1,2200),F(1,100))
        record=result([0,1]*32)
        self.assertEqual(F(record['failure_budget_fraction']),F(1,1100))


if __name__=='__main__':unittest.main()
