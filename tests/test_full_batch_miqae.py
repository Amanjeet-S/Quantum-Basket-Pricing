"""Boundary, exact chain-budget and ideal Born diagnostics."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import math
import unittest
import importlib.util
import os
from full_batch_miqae import (plan,estimate,cover_certificate,quadrant,find_next_multiplier,
    rational_record)
from published_envelopes import theorem_envelopes
from certified_scalar import probability_interval,inverse_probability_interval


class FullBatchTests(unittest.TestCase):
    def test_cover_exact_enclosure(self):
        certificate=cover_certificate()
        self.assertEqual(certificate['checked_triples'],48)
        self.assertGreater(certificate['strict_margin_lower'],F(7,10000))

    def test_dynamic_program_matches_independent_chain_recursion(self):
        for prior in (F(1),F(1,10),F(1,50)):
            setup=plan(F(1,100),prior_upper=prior)
            if not setup.maximum_rounds:continue
            table=setup.shot_table();maximum=setup.maximum_odd_multiplier
            for a_cost,q_cost in ((1,0),(0,1),(11,19)):
                @lru_cache(None)
                def tail(current):
                    return max([0]+[table[m]*(a_cost+q_cost*((m-1)//2))+tail(m)
                        for m in range(3*current,maximum+1,2)])
                m0=setup.initial_multiplier
                expected=table[m0]*(a_cost+q_cost*((m0-1)//2))+tail(m0)
                retained=setup.cost_upper(a_cost,q_cost)
                self.assertEqual(retained['cost_upper'],expected)
                chain=retained['witness_multipliers']
                self.assertEqual(chain[0],m0)
                self.assertTrue(all(b>=3*a for a,b in zip(chain,chain[1:])))
                self.assertEqual(retained['witness_shots']*a_cost+retained['witness_Q_calls']*q_cost,expected)

    def test_zero_shot_and_closed_boundaries(self):
        epsilon=F(1,5000)
        for prior in (F(0),2*epsilon):self.assertEqual(plan(epsilon,prior_upper=prior).maximum_shots,0)
        self.assertGreater(plan(epsilon,prior_upper=2*epsilon+F(1,10**8)).maximum_shots,0)
        for multiplier in (1,3,5,7):
            self.assertEqual(quadrant(multiplier,F(0),F(1,multiplier)),0)
        self.assertIsNone(quadrant(3,F(1,3)-F(1,10**12),F(1,3)+F(1,10**12)))

    def test_expected_count_paths(self):
        for cap in (1,10,25):
            for p,prior in ((0.,F(1)),(1.,F(1)),(.0001,F(1,1000)),(.01,F(1,50)),(.25,F(1)),(.5,F(1)),(.75,F(1)),(.9999,F(1))):
                theta=math.asin(math.sqrt(p))
                result=estimate(lambda k,n:int(round(n*math.sin((2*k+1)*theta)**2)),F(1,200*cap),prior_upper=prior)
                self.assertEqual(result['status'],'certified')
                lo,hi=result['amplitude_interval']
                self.assertLessEqual(hi-lo,F(1,100*cap))
                self.assertLessEqual(float(lo),p+1e-13);self.assertGreaterEqual(float(hi),p-1e-13)
                self.assertLessEqual(result['total_shots'],result['plan'].maximum_shots)
                self.assertLessEqual(result['total_Q_calls'],result['plan'].maximum_Q_calls)
                failures=sum(r['conditional_failure_upper'] for r in result['rounds'])
                self.assertLessEqual(failures,result['plan'].alpha)

    def test_failure_status_and_invalid_backend(self):
        rr=estimate(lambda k,n:n,F(1,5000),prior_upper=F(1,10))
        self.assertTrue(rr['status'].startswith('confidence_failure'))
        with self.assertRaises(ValueError):estimate(lambda k,n:1.5,F(1,5000))

    def test_complete_small_quantum_oracle_backend(self):
        source=Path(os.environ.get('QBP_REPO',Path(__file__).resolve().parents[1]))/'quantum_gates.py'
        spec=importlib.util.spec_from_file_location('qbp_small_quantum_gates',source)
        gates=importlib.util.module_from_spec(spec);spec.loader.exec_module(gates)
        # One finite fractional output bit x and an independent comparison bit.
        # The clean flag is1 iff j<x, giving the exact Born mean probability1/4.
        A=[('H',0),('H',1),('X',1),('CCX',0,1,2),('X',1)]
        Q=gates.grover(A,(0,1,2),2,(3,))
        state=gates.replay(A,gates.zero_state(4));power=0
        self.assertAlmostEqual(gates.flag_probability(state,2),.25,places=14)
        def backend(k,n):
            nonlocal state,power
            while power<k:state=gates.replay(Q,state);power+=1
            self.assertLess(sum(abs(z)**2 for i,z in enumerate(state) if i&8),1e-22)
            probability=gates.flag_probability(state,2)
            self.assertAlmostEqual(probability,math.sin((2*k+1)*math.pi/6)**2,places=10)
            return int(round(n*probability))
        result=estimate(backend,F(1,5000))
        self.assertEqual(result['status'],'certified')
        lo,hi=result['amplitude_interval'];self.assertLessEqual(float(lo),.25);self.assertGreaterEqual(float(hi),.25)

    def test_published_envelope_money_corollary(self):
        for cap in (1,10,25):
            for alpha in (F(1,100),F(1,2700)):
                full=theorem_envelopes(F(1,200*cap),alpha)
                prior=theorem_envelopes(F(1,200*cap),alpha,F(1,50))
                for key in ('iqae_published','miqae_full_batch_integer_corollary'):
                    self.assertLess(prior[key]['total_grover_calls_upper'],full[key]['total_grover_calls_upper'])
                    self.assertGreater(full[key]['total_shots_upper'],0)
                    self.assertGreater(full[key]['grover_power_upper'],0)


def ideal_simulation():
    import numpy as np
    generator=np.random.default_rng(20261006);records=[]
    pairs=[(0.,F(1)),(1.,F(1)),(.0001,F(1,1000)),(.01,F(1,50)),(.025,F(1,20)),(.25,F(1)),(.5,F(1)),(.75,F(1)),(.9999,F(1))]
    for cap in (1,10,25):
        for probability,prior in pairs:
            theta=math.asin(math.sqrt(probability));covered=flagged=0;shots=[];queries=[]
            for repeat in range(20):
                def backend(k,n):return int(generator.binomial(n,math.sin((2*k+1)*theta)**2))
                result=estimate(backend,F(1,200*cap),prior_upper=prior)
                lo,hi=result['amplitude_interval']
                covered+=float(lo)-1e-13<=probability<=float(hi)+1e-13
                flagged+=result['status'].startswith('confidence_failure')
                assert result['total_shots']<=result['plan'].maximum_shots
                assert result['total_Q_calls']<=result['plan'].maximum_Q_calls
                shots.append(result['total_shots']);queries.append(result['total_Q_calls'])
            records.append({'cap':cap,'synthetic_probability':probability,'proved_prior_upper':rational_record(prior),
                'repetitions':20,'covered':int(covered),'flagged_failures':int(flagged),
                'shots_range':[min(shots),max(shots)],'Q_calls_range':[min(queries),max(queries)]})
    return {'scope':'540 ideal-binomial software diagnostics on synthetic amplitudes; no device or market observations; coverage proof is independent of observed simulation coverage',
            'seed':20261006,'records':records,'covered':sum(r['covered'] for r in records),'flagged':sum(r['flagged_failures'] for r in records)}


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(FullBatchTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    output=Path(__file__).resolve().parents[1]/'results';output.mkdir(exist_ok=True)
    out={'test_groups':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
        'cover_certificate':rational_record(cover_certificate()),'ideal_simulation':ideal_simulation()}
    (output/'full_batch_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'covered':out['ideal_simulation']['covered'],'flagged':out['ideal_simulation']['flagged']}))
    if not result.wasSuccessful():raise SystemExit(1)
