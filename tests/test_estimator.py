"""Independent boundary, confidence-controller and complete quantum checks."""
from fractions import Fraction as F
import json
import math
import os
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).parent))
sys.path.insert(0,str(Path(os.environ.get('QBP_REPO',Path(__file__).resolve().parent))))
from reversible import Circuit
from iterative_estimator import (plan,estimate,quadrant,probability_interval,
    inverse_probability_interval,sqrt_upper_fraction,decimal_record)
from quantum_gates import (grover,replay,zero_state,flag_probability,resources,expand_toffoli)


def toy_preparation():
    c=Circuit();x=c.alloc(2);j=c.alloc(2);flag=c.alloc(1)[0]
    c.compare_into(j,x,flag)
    names={1:'X',2:'CX',3:'CCX'}
    boolean=[]
    for op,a,b,t in c.gates:
        boolean.append((names[op],t) if op==1 else (names[op],a,t) if op==2 else (names[op],a,b,t))
    return [('H',w) for w in x+j]+boolean,tuple(x+j)+(flag,),flag,c.wire_count


class EstimatorTests(unittest.TestCase):
    def test_exact_seven_T_expansion(self):
        gates=expand_toffoli([('CCX',0,1,2)])
        self.assertEqual(resources(gates,3)['T_count_exact_7T_Toffoli'],7)
        self.assertEqual(resources(gates,3)['T_depth_assigned_3layer_Toffoli_or_explicit_T_upper'],4)
        for basis in range(8):
            initial=[0j]*8;initial[basis]=1
            direct=replay([('CCX',0,1,2)],initial)
            expanded=replay(gates,initial)
            self.assertLess(max(abs(a-b) for a,b in zip(direct,expanded)),1e-14)
        A,core,flag,n=toy_preparation()
        direct=replay(A,zero_state(n));expanded=replay(expand_toffoli(A),zero_state(n))
        self.assertLess(max(abs(a-b) for a,b in zip(direct,expanded)),1e-12)

    def test_probability_enclosures(self):
        checks={F(0):0.,F(1):1.,F(1,2):.5,F(1,3):.25,F(2,3):.75}
        for u,reference in checks.items():
            a,b=probability_interval(u)
            self.assertLessEqual(a,b);self.assertLessEqual(float(a),reference+1e-15)
            self.assertGreaterEqual(float(b),reference-1e-15)
            self.assertLess(b-a,F(1,10**70))

    def test_inverse_contract(self):
        for p in (F(0),F(1),F(1,2),F(1,100000),F(1,4),F(3,4),F(99999,100000)):
            lo,hi=inverse_probability_interval(p)
            self.assertLessEqual(probability_interval(lo)[0],p)
            self.assertGreaterEqual(probability_interval(hi)[1],p)
            self.assertLessEqual(hi-lo,F(1,1<<69))

    def test_quadrant_endpoints(self):
        self.assertEqual(quadrant(3,F(0),F(1,3)),0)
        self.assertEqual(quadrant(3,F(1,3),F(2,3)),1)
        self.assertEqual(quadrant(3,F(2,3),F(1)),2)
        self.assertIsNone(quadrant(3,F(1,3)-F(1,10**12),F(1,3)+F(1,10**12)))

    def test_three_five_cover(self):
        # Interval centres include every internal 3/5 stationary location,
        # endpoints, and neighbours on both sides of those values.
        critical=[F(0),F(1),F(1,3),F(2,3)]+[F(j,5) for j in range(1,5)]
        points=sorted(set(critical+[min(F(1),max(F(0),x+s*F(1,10**8)))
                     for x in critical for s in (-1,1)]+[F(j,100) for j in range(101)]))
        for centre in points:
            p=probability_interval(centre);pl=max(F(0),p[0]-F(1,25));pu=min(F(1),p[1]+F(1,25))
            lo=inverse_probability_interval(pl)[0];hi=inverse_probability_interval(pu)[1]
            self.assertTrue(quadrant(3,lo,hi) is not None or quadrant(5,lo,hi) is not None)

    def test_plan_prior_and_confidence(self):
        full=plan(F(1,5000));small=plan(F(1,5000),prior_upper=F(1,50))
        self.assertLess(small.total_grover_calls_upper,full.total_grover_calls_upper)
        self.assertLessEqual(small.initial_multiplier*inverse_probability_interval(F(1,50))[1],1)
        self.assertTrue(small.initial_multiplier%2)
        self.assertEqual(plan(F(1,5000),prior_upper=F(0)).total_shots_upper,0)
        self.assertEqual(plan(F(1,5000),prior_upper=F(1,10000)).total_shots_upper,0)
        self.assertGreater(plan(F(1,5000),F(1,2700)).shots_per_round,full.shots_per_round)
        with self.assertRaises(ValueError):plan(F(1,5000),F(1,10**100))
        root=sqrt_upper_fraction(F(2,100));self.assertGreaterEqual(root*root,F(2,100))

    def test_controller_deterministic_expected_counts(self):
        cases=[(0.,F(1)),(1.,F(1)),(.00001,F(1,1000)),(.01,F(1,50)),(.25,F(1)),
               (.5,F(1)),(.75,F(1)),(.99999,F(1))]
        for p,prior in cases:
            theta=math.asin(math.sqrt(p))
            def backend(k,n):return int(round(n*math.sin((2*k+1)*theta)**2))
            rr=estimate(backend,F(1,5000),prior_upper=prior)
            self.assertEqual(rr['status'],'certified')
            a,b=rr['amplitude_interval'];self.assertLessEqual(b-a,F(1,2500))
            self.assertLessEqual(float(a),p+1e-12);self.assertGreaterEqual(float(b),p-1e-12)
            self.assertLessEqual(rr['total_shots'],rr['plan'].total_shots_upper)
            self.assertLessEqual(rr['total_grover_calls'],rr['plan'].total_grover_calls_upper)
        with self.assertRaises(ValueError):estimate(lambda k,n:1.5,F(1,5000))

    def test_empty_intersection_failure_is_declared(self):
        # Impossible backend data against a proved zero-like prior must not
        # be repaired into a falsely claimed interval. Initial prior stays
        # wider than 2epsilon, so at least one batch is actually queried.
        rr=estimate(lambda k,n:n,F(1,5000),prior_upper=F(1,100))
        self.assertTrue(rr['status'].startswith('confidence_failure'))

    def test_quantum_power_probabilities_and_work(self):
        A,core,flag,n=toy_preparation();work=tuple(w for w in range(n) if w not in core)
        Q=grover(A,core,flag,work);psi=replay(A,zero_state(n));p=6/16
        self.assertAlmostEqual(flag_probability(psi,flag),p,places=12)
        theta=math.asin(math.sqrt(p));state=psi
        for k in range(6):
            self.assertAlmostEqual(flag_probability(state,flag),math.sin((2*k+1)*theta)**2,places=11)
            self.assertLess(sum(abs(z)**2 for i,z in enumerate(state) if any((i>>w)&1 for w in work)),1e-24)
            state=replay(Q,state)
        self.assertEqual(resources(Q,n)['T_count_exact_7T_Toffoli'],
                         2*resources(A,n)['T_count_exact_7T_Toffoli']+14*(len(core)-2))

    def test_controlled_quantum_global_sign(self):
        A,core,flag,n=toy_preparation();control=n
        work=tuple(w for w in range(n) if w not in core)
        Q=grover(A,core,flag,work);CQ=grover(A,core,flag,work,control)
        psi=replay(A,zero_state(n));rot=[-z for z in replay(Q,psi)]
        # Control0 remains psi; control1 has the globally corrected Q psi.
        for c in (0,1):
            initial=[0j]*(1<<(n+1));offset=c*(1<<n)
            initial[offset:offset+(1<<n)]=psi
            actual=replay(CQ,initial);expected=psi if c==0 else rot
            self.assertLess(max(abs(actual[offset+i]-expected[i]) for i in range(1<<n)),1e-12)
        initial=[z/math.sqrt(2) for z in psi]+[z/math.sqrt(2) for z in psi]
        actual=replay(CQ,initial);expected=[z/math.sqrt(2) for z in psi]+[z/math.sqrt(2) for z in rot]
        self.assertLess(max(abs(a-b) for a,b in zip(actual,expected)),1e-12)
        self.assertEqual(resources(CQ,n+1)['T_count_exact_7T_Toffoli']-resources(Q,n)['T_count_exact_7T_Toffoli'],14)


def main():
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(EstimatorTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    A,core,flag,n=toy_preparation();work=tuple(w for w in range(n) if w not in core)
    out={'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
         'small_preparation_probability':6/16,'small_preparation_resources':resources(A,n),
         'exact_Toffoli_expansion_resources':resources(expand_toffoli([('CCX',0,1,2)]),3),
         'small_preparation_expanded_resources':resources(expand_toffoli(A),n),
         'small_grover_expanded_resources':resources(expand_toffoli(grover(A,core,flag,work)),n),
         'small_grover_resources':resources(grover(A,core,flag,work),n),
         'small_controlled_grover_resources':resources(grover(A,core,flag,work,n),n+1),
         'scope':'small complete logical state-vector and rational-controller checks; no quantum device or market data'}
    output=Path(__file__).resolve().parent/'results';output.mkdir(parents=True,exist_ok=True)
    (output/'estimator_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
