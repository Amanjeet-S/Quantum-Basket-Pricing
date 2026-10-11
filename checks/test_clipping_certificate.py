"""Independent feasible-law and dual-majorant clipping checks."""
from fractions import Fraction as F
import unittest
from clipping_certificate import clipped_tail_bound,basic_bounds,polynomial,sqrt_interval


class ClippingTests(unittest.TestCase):
    def test_all_branch_witnesses_and_dual_majorants(self):
        cases=set()
        for rho in (F(1,4),F(1,2),F(1),F(2)):
            for r in (F(0),rho/8,rho/4,rho/2,3*rho/4,rho,2*rho):
                for u in (F(0),rho/100,rho/8,rho/2,rho,2*rho):
                    for v in (F(0),rho*rho/1000,rho*rho/32,rho*rho/4,rho*rho,2*rho*rho):
                        record=clipped_tail_bound(rho,u,v,r);cases.add(record['case'])
                        self.assertLessEqual(record['lower'],record['upper'])
                        self.assertGreaterEqual(record['upper'],0)
                        self.assertLessEqual(record['upper']-record['lower'],F(1,1 << 100))
                        if record['dual_quadratic'] is not None:
                            for index in range(101):
                                x=rho*F(index,100)
                                self.assertGreaterEqual(polynomial(record,x),max(F(0),x-r))
                            a,b,c=record['dual_quadratic']
                            self.assertTrue(b>=0 and c>=0)
                            dual_cost=a+b*record['effective_mean_upper']+c*record['effective_second_moment_upper']
                            self.assertLessEqual(dual_cost,record['upper'])
                        witness=record['attaining_law']
                        if isinstance(witness,tuple):
                            self.assertEqual(sum((p for x,p in witness),F(0)),1)
                            self.assertTrue(all(0<=x<=rho and p>=0 for x,p in witness))
                            mean=sum((x*p for x,p in witness),F(0))
                            second=sum((x*x*p for x,p in witness),F(0))
                            value=sum((max(F(0),x-r)*p for x,p in witness),F(0))
                            self.assertLessEqual(mean,u);self.assertLessEqual(second,v)
                            self.assertEqual(value,record['upper'])
        self.assertEqual(cases,{'zero','mean_only','constant_root','range_second_moment',
                               'second_moment_tangent','mean_second_moment'})

    def test_arbitrary_finite_feasible_distributions(self):
        for rho in (F(1,4),F(1),F(3)):
            values=[rho*F(i,5) for i in range(6)]
            for n in range(1,5):
                for start in range(6):
                    law=[values[(start+i)%6] for i in range(n)]
                    mean=sum(law,F(0))/n;second=sum((x*x for x in law),F(0))/n
                    for r in (F(0),rho/8,rho/4,rho/2,3*rho/4,rho):
                        bound=clipped_tail_bound(rho,mean,second,r)
                        actual=sum((max(F(0),x-r) for x in law),F(0))/n
                        self.assertLessEqual(actual,bound['upper'])

    def test_basic_bound_dominance_and_boundaries(self):
        for rho in (F(1,4),F(1),F(2)):
            for r in (rho/8,rho/4,rho/2,3*rho/4):
                for u in (rho/100,rho/8,rho/2,rho):
                    for v in (rho*rho/1000,rho*rho/32,rho*rho/4,rho*rho):
                        sharp=clipped_tail_bound(rho,u,v,r)
                        simple=basic_bounds(rho,u,v,r)
                        self.assertLessEqual(sharp['upper'],simple['minimum']+F(1,1 << 100))
        # All neighbouring formulas agree at their shared boundaries.
        rho=F(1);r=F(1,4);v=4*r*r;u=2*r
        self.assertEqual(clipped_tail_bound(rho,u,v,r)['upper'],r)
        self.assertEqual(clipped_tail_bound(2*r,u,v,r)['upper'],r)
        self.assertEqual(clipped_tail_bound(rho,F(1,3),F(1,9),r)['upper'],F(1,9))
        self.assertEqual(sqrt_interval(F(1,9)),(F(1,3),F(1,3)))


if __name__=='__main__':unittest.main(verbosity=2)
