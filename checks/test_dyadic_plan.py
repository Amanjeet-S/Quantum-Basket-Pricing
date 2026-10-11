"""Exact dyadic identities and independent scalar-resource replay checks."""
from fractions import Fraction as F
import math
import unittest
from dyadic_plan import (partition,exact_bucket_means,build_plan,capped_allocation,
                         allocation_weights,execute,choose_plan,MIN_EPSILON)
from dyadic_plan import sharp_grid_mean_upper
import full_batch_miqae
import iterative_estimator


def resource_profiles(q):
    return {b:{'A':{'logical_T_count':31+b,'measurement_count':1},
               'Q':{'logical_T_count':65+2*b,'measurement_count':2},
               'qubits':8+b} for b in range(q+1)}


class DyadicTests(unittest.TestCase):
    def test_certified_support_maximum_removes_unavailable_buckets(self):
        q=12;maximum=31
        b,tail=partition(q,0,F(1,20),F(1,200),max_integer=maximum)
        self.assertEqual(tail,0)
        self.assertTrue(all(bucket.prior==0 for bucket in b if bucket.index>=5))
        self.assertEqual(b[4].prior,F(31,32))
        profiles={i:p for i,p in resource_profiles(q).items() if i<=4}
        p=build_plan(q,F(1,20),F(1,200),F(1,1000),F(1,100),profiles,
                     'three_five','cost_moment',max_integer=maximum,sharp_priors=True)
        self.assertEqual(p['max_integer'],maximum)
        self.assertTrue(all(bucket['maximum_shots']==0 for bucket in p['buckets'] if bucket['index']>=5))
        self.assertLessEqual(p['estimate_clipping_interval'][1],F(maximum,1 << q))

    def test_sharp_finite_moment_bound_and_attainment(self):
        for q in range(1,7):
            delta=F(1,1 << q)
            for lower in range(1,1 << q):
                upper=min(1 << q,lower+3)
                lo,hi=lower*delta,upper*delta
                for v in (F(0),lo*lo/2,lo*lo,(lo*lo+hi*hi)/2,hi*hi,F(1)):
                    m=sharp_grid_mean_upper(q,v,lower,upper)
                    self.assertLessEqual(m,hi)
                    if m<lo:
                        # Exact zero/l primal witness and second moment.
                        probability=m/lo
                        self.assertTrue(0<=probability<=1)
                        self.assertLessEqual(probability*lo*lo,v)
                    elif m==hi:self.assertLessEqual(hi*hi,v)
                    else:
                        j=m//delta;x=j*delta;y=x+delta
                        p=(m-x)/(y-x)
                        self.assertTrue(0<=p<=1)
                        self.assertEqual((1-p)*x+p*y,m)
                        self.assertLessEqual((1-p)*x*x+p*y*y,v)
                    # The convex supporting line bounds every allowed atom.
                    for a in [0]+list(range(lower,upper+1)):
                        z=a*delta
                        if m<lo:self.assertLessEqual(lo*z,z*z)
                        elif m<hi:
                            j=m//delta;x=j*delta;y=x+delta
                            self.assertLessEqual((x+y)*z-x*y,z*z)

    def test_exact_identity_all_small_atoms(self):
        for q in range(1,10):
            buckets,_=partition(q,0,1,1)
            for a in range((1 << q)+1):
                ps=exact_bucket_means([a],q)
                self.assertEqual(sum((b.weight*ps[b.index] for b in buckets),F(0)),F(a,1 << q))

    def test_moment_priors_and_omitted_prefix_all_atoms(self):
        for q in range(1,8):
            for a in range((1 << q)+1):
                x=F(a,1 << q);moments=exact_bucket_means([a],q)
                for s in range(q+1):
                    buckets,tail=partition(q,s,x,x*x)
                    omitted=x if a<1 << s else F(0)
                    self.assertLessEqual(omitted,tail)
                    for b in buckets:self.assertLessEqual(moments[b.index],b.prior)
                    reconstructed=omitted+sum((b.weight*moments[b.index] for b in buckets),F(0))
                    self.assertEqual(reconstructed,x)

    def test_aggregate_second_moment_relation(self):
        for q in range(2,9):
            values=list(range((1 << q)+1));p=exact_bucket_means(values,q)
            v=sum((F(a*a,1 << (2*q)) for a in values),F(0))/len(values)
            b,_=partition(q,0,1,1)
            self.assertLessEqual(sum((a.weight*a.weight*p[a.index] for a in b),F(0)),2*v)

    def test_exact_capped_waterfill_and_declared_range(self):
        b,_=partition(12,5,F(1,20),F(1,500))
        profiles=resource_profiles(12)
        for style in ('equal','dyadic','cost_moment'):
            weights=allocation_weights(b,profiles,style)
            errors=capped_allocation(b,F(1,500),weights)
            self.assertLessEqual(sum(errors,F(0)),F(1,500))
            for bucket,error in zip(b,errors):
                self.assertLessEqual(error,bucket.weight*bucket.prior/2)
                if error<bucket.weight*bucket.prior/2:self.assertGreaterEqual(error/bucket.weight,MIN_EPSILON)
        with self.assertRaises(ValueError):capped_allocation(b,F(1,10**12),tuple(F(1) for _ in b))

    def test_whole_prior_zero_shot_and_top_atom(self):
        profiles=resource_profiles(12)
        p=build_plan(12,F(1,100),F(1,10**6),F(1,1000),F(1,100),profiles)
        self.assertEqual(p['status'],'whole_prior_zero_shot')
        self.assertLessEqual(p['total_error_upper'],F(1,1000))
        b,_=partition(12,0,1,F(1,10))
        self.assertEqual(b[-1].index,12);self.assertTrue(b[-1].top_atom)
        self.assertEqual(b[-1].prior,F(1,10));self.assertEqual(b[-1].weight,F(1))

    def test_metrics_against_independent_scalar_plan(self):
        q=9;profiles=resource_profiles(q)
        for controller in ('miqae','three_five'):
            p=build_plan(q,F(1,20),F(1,200),F(1,1000),F(1,100),profiles,controller)
            expected={name:0 for name in ('logical_T_count','measurement_count')}
            for b in p['buckets']:
                if b['status']=='prior_zero_shot':continue
                scalar=(full_batch_miqae if controller=='miqae' else iterative_estimator).plan(
                    b['probability_tolerance'],b['alpha'],b['prior_upper'])
                for metric in expected:
                    ca,cq=profiles[b['index']]['A'][metric],profiles[b['index']]['Q'][metric]
                    expected[metric]+=(scalar.cost_upper(ca,cq)['cost_upper'] if controller=='miqae'
                                     else scalar.total_shots_upper*ca+scalar.total_grover_calls_upper*cq)
            self.assertEqual(p['aggregate_metrics_upper'],expected)
            self.assertLessEqual(sum((b['alpha'] for b in p['buckets']),F(0)),p['alpha'])

    def test_expected_born_count_full_bucket_execution(self):
        q=8;values=[0,1,7,15,37,127,256];profiles=resource_profiles(q)
        mean=sum((F(a,1 << q) for a in values),F(0))/len(values)
        second=sum((F(a*a,1 << (2*q)) for a in values),F(0))/len(values)
        probabilities=exact_bucket_means(values,q)
        for controller in ('miqae','three_five'):
            p=build_plan(q,mean,second,F(1,100),F(1,100),profiles,controller)
            def backend(index,k,n):
                theta=math.asin(math.sqrt(float(probabilities[index])))
                return int(round(n*math.sin((2*k+1)*theta)**2))
            result=execute(p,backend,profiles)
            self.assertEqual(result['status'],'certified')
            lo,hi=result['mean_interval']
            self.assertLessEqual(float(lo),float(mean)+1e-13)
            self.assertGreaterEqual(float(hi),float(mean)-1e-13)
            self.assertLessEqual(abs(result['estimate']-mean),p['tolerance'])

    def test_finite_candidate_selection_retains_adverse_peers(self):
        p=choose_plan(8,F(1,20),F(1,200),F(1,100),F(1,100),resource_profiles(8),
                      controller='three_five',prefixes=(0,1),allocations=('equal','cost_moment'))
        self.assertEqual(len(p['candidates']),4)
        winner=p['winner']['aggregate_metrics_upper']['logical_T_count']
        self.assertTrue(all(winner<=r['aggregate_metrics_upper']['logical_T_count'] for r in p['candidates']))


if __name__=='__main__':
    unittest.main(verbosity=2)
