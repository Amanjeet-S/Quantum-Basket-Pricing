"""Independent finite arithmetic and coherent-state checks for QPE plans."""
from fractions import Fraction as F
from pathlib import Path
import cmath
import json
import math
import random
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import qpe_estimator as qe
import quantum_gates as qg
import qft_synthesis as qs


def replay(gates,state):
    state=list(state)
    pending=[]
    def flush():
        nonlocal state,pending
        if pending:
            state=qg.replay(pending,state)
            pending=[]
    for g in gates:
        if g[0] not in ('CP','SWAP'):
            pending.append(g)
            continue
        flush()
        if g[0]=='CP':
            _,angle,a,b=g
            phase=cmath.exp(1j*math.pi*float(angle))
            for x in range(len(state)):
                if (x>>a)&1 and (x>>b)&1:
                    state[x]*=phase
        else:
            _,a,b=g
            for x in range(len(state)):
                if (x>>a)&1 and not (x>>b)&1:
                    y=x^(1<<a)^(1<<b)
                    state[x],state[y]=state[y],state[x]
    flush()
    return state


def truth_preparation(good):
    gates=[('H',0),('H',1)]
    for x in range(good):
        flips=[('X',w) for w in range(2) if not (x>>w)&1]
        gates+=flips+[('CCX',0,1,2)]+flips[::-1]
    return gates


def full_qpe(good,bits,include_controlled_sign=True,actual_qft=False,
             most_significant_first=False,omit_fourier_relabelling=False):
    core=tuple(range(3));phase=tuple(range(3,3+bits))
    work=(3+bits,4+bits)
    n=5+bits
    A=truth_preparation(good)
    gates=list(A)+[('H',w) for w in phase]
    for j,w in enumerate(phase):
        controlled=qg.grover(A,core,2,work,control=w)
        if not include_controlled_sign:
            assert controlled[-1]==('Z',w)
            controlled=controlled[:-1]
        exponent=bits-1-j if most_significant_first else j
        gates+=controlled*(1<<exponent)
    if actual_qft:
        info=qs.inverse_qft_tape(bits)
        gates+=[(g[0],*(phase[w if most_significant_first or omit_fourier_relabelling
                                  else bits-1-w] for w in g[1:]))
                for g in info['gates']]
    else:
        gates+=qe.ideal_inverse_qft_gates(bits,phase)
    state=replay(gates,qg.zero_state(n))
    distribution=[0.0]*(1<<bits)
    for x,z in enumerate(state):
        y=sum(((x>>w)&1)<<(bits-1-j if most_significant_first else j)
              for j,w in enumerate(phase))
        distribution[y]+=abs(z)**2
    return distribution


class QPETests(unittest.TestCase):
    def test_rational_power_ceiling(self):
        for epsilon in (F(1,200),F(1,2000),F(1,5000),F(1,10**6)):
            for prior in (F(1),F(1,2),F(1,100),F(1,10**6)):
                p=qe.plan(epsilon,prior_upper=prior)
                if p.zero_shots:
                    self.assertLessEqual(prior/2,epsilon)
                    continue
                self.assertEqual(p.M,1<<p.phase_register_bits)
                self.assertLessEqual(p.amplitude_error_upper+p.postprocess_allowance,epsilon)
                previous=2*qe.PI_UPPER*qe.sqrt_upper_fraction(p.variance_upper)/(p.M//2)+qe.PI_UPPER**2/(p.M//2)**2
                self.assertGreater(previous+qe.POSTPROCESS_ALLOWANCE,epsilon)
                self.assertEqual(p.total_controlled_grover_calls,p.repetitions*(p.M-1))

    def test_confidence_binomial_and_hoeffding(self):
        for alpha in (F(1,100),F(1,2700),F(1,10**8)):
            for confidence in ('exact_binomial','hoeffding'):
                n,s,tail,nh=qe.confidence_repetitions(alpha,confidence=confidence)
                self.assertEqual(n%2,1)
                self.assertLessEqual(tail,alpha)
                self.assertGreater(s,F(1,2))
                self.assertLessEqual(n,nh)
                if confidence=='exact_binomial' and n>1:
                    self.assertGreater(qe.binomial_failure_upper(n-2,s),alpha)
                if confidence=='hoeffding':
                    self.assertGreaterEqual(2*n*(s-F(1,2))**2,qe.log_upper(1/alpha))

    def test_inverse_fourier_matrix(self):
        comparisons=0
        for m in range(1,7):
            M=1<<m
            gates=qe.ideal_inverse_qft_gates(m)
            for x in range(M):
                state=[complex(j==x) for j in range(M)]
                result=replay(gates,state)
                target=[cmath.exp(-2j*math.pi*x*y/M)/math.sqrt(M) for y in range(M)]
                self.assertLess(max(abs(a-b) for a,b in zip(result,target)),2e-14)
                comparisons+=M
        self.assertEqual(comparisons,sum(4**m for m in range(1,7)))

    def test_controlled_grover_sign_and_distribution(self):
        for m in range(2,5):
            M=1<<m
            for good in range(5):
                actual=full_qpe(good,m)
                expected=qe.ideal_amplitude_distribution(good/4,M)
                self.assertLess(max(abs(a-b) for a,b in zip(actual,expected)),2e-13)
        wrong=full_qpe(4,3,False)
        self.assertGreater(wrong[0],1-1e-12)
        self.assertLess(wrong[4],1e-12)

    def test_literal_clifford_t_fourier_tape(self):
        for m in range(1,4):
            M=1<<m
            info=qs.inverse_qft_tape(m)
            gates=[(g[0],*(m-1-w for w in g[1:])) for g in info['gates']]
            scalar=cmath.exp(1j*math.pi*float(info['global_phase_correction_pi']))
            for x in range(M):
                actual=replay(gates,[complex(j==x) for j in range(M)])
                target=[cmath.exp(-2j*math.pi*x*y/M)/math.sqrt(M) for y in range(M)]
                vector_error=math.sqrt(sum(abs(scalar*a-b)**2 for a,b in zip(actual,target)))
                self.assertLessEqual(vector_error,float(info['norm_error_upper'])+1e-12)
        for m in range(2,4):
            for good in range(5):
                actual=full_qpe(good,m,actual_qft=True)
                target=qe.ideal_amplitude_distribution(good/4,1<<m)
                total_variation=sum(abs(a-b) for a,b in zip(actual,target))/2
                self.assertLessEqual(total_variation,float(qs.inverse_qft_profile(m)['norm_error_upper'])+1e-12)

    def test_all_retained_fourier_profiles(self):
        for m in range(1,33):
            profile=qs.inverse_qft_profile(m)
            self.assertLessEqual(profile['norm_error_upper'],F(1,64))
            gates=qs.inverse_qft_tape(m)['gates']
            counts={g:sum(t[0]==g for t in gates) for g in ('T','Tdg','H','CX','S','Sdg')}
            self.assertEqual(profile['T_count'],counts['T']+counts['Tdg'])
            self.assertEqual(profile['gate_counts'],{g:n for g,n in counts.items() if n})

    def test_nonpalindromic_phase_wire_convention(self):
        target=qe.ideal_amplitude_distribution(.25,8)
        lsb=full_qpe(1,3,actual_qft=True)
        msb=full_qpe(1,3,actual_qft=True,most_significant_first=True)
        self.assertLess(max(abs(a-b) for a,b in zip(lsb,msb)),1e-12)
        allowance=float(qs.inverse_qft_profile(3)['norm_error_upper'])
        self.assertLessEqual(sum(abs(a-b) for a,b in zip(msb,target))/2,allowance+1e-12)
        wrong=full_qpe(1,3,actual_qft=True,omit_fourier_relabelling=True)
        self.assertGreater(sum(abs(a-b) for a,b in zip(wrong,target))/2,.1)

    def test_actual_compiled_norm_confidence(self):
        for epsilon in (F(1,200),F(1,2000),F(1,5000)):
            actual,profile=qe.compiled_plan(epsilon)
            conservative=qe.plan(epsilon)
            self.assertEqual(actual.M,conservative.M)
            raw=F(profile['norm_error_upper'])
            rounded=F(qe.ceil_fraction(raw*(1<<80)),1<<80)
            self.assertEqual(actual.qft_norm_error_allowance,rounded)
            self.assertGreaterEqual(rounded,raw)
            self.assertLess(rounded-raw,F(1,1<<80))
            self.assertGreaterEqual(actual.good_probability_lower,conservative.good_probability_lower)
            self.assertLessEqual(actual.repetitions,conservative.repetitions)
            self.assertLessEqual(actual.binomial_failure_upper,actual.alpha)
        zero,profile=qe.compiled_plan(F(1,200),prior_upper=F(1,1000))
        self.assertTrue(zero.zero_shots)
        self.assertIsNone(profile)

    def test_compact_exact_confidence_record(self):
        # Dyadic confidence allocations can raise a high-precision norm
        # denominator to a large repetition power. Retain a conservative
        # binary upper enclosure without disabling Python's digit guard.
        setup,profile=qe.compiled_plan(F(1,10**6),F(1,2700))
        serialised=json.dumps(setup.record())
        self.assertLess(len(serialised),15000)
        self.assertLessEqual(setup.binomial_failure_upper,F(1,2700))
        self.assertGreaterEqual(setup.qft_norm_error_allowance,profile['norm_error_upper'])

    def test_postprocessing_and_median(self):
        for bits in range(1,10):
            M=1<<bits
            for y in range(M):
                lo,hi=qe.output_interval(y,M)
                expected=math.sin(math.pi*y/M)**2
                self.assertLess(abs(float((lo+hi)/2)-expected),2e-15)
                self.assertLessEqual(hi-lo,2*qe.POSTPROCESS_ALLOWANCE)
        setup=qe.plan(F(1,100),prior_upper=F(1,2))
        result=qe.certified_median(setup,[setup.M//4]*setup.repetitions)
        self.assertEqual(result['estimate'],F(1,2))
        with self.assertRaises(ValueError):
            qe.certified_median(setup,[0])

    def test_ideal_seeded_replay(self):
        setup=qe.plan(F(1,100))
        a=qe.execute_ideal(setup,F(1,4),random.Random(20261007))
        b=qe.execute_ideal(setup,F(1,4),random.Random(20261007))
        self.assertEqual(a,b)
        self.assertEqual(len(a['outcomes']),setup.repetitions)
        self.assertLess(abs(float(a['estimate'])-.25),.01)
        self.assertIn('no physical',a['execution_scope'])

    def test_complete_profile_charge(self):
        setup=qe.plan(F(1,100))
        a={'metrics':{'T':5,'H':2},'logical_qubits':21}
        q={'metrics':{'T':20,'H':4,'Z':1},'logical_qubits':22}
        f={'phase_register_bits':setup.phase_register_bits,'norm_error_upper':'1/64',
           'metrics':{'T':123,'H':setup.phase_register_bits},'extra_ancillas':3}
        result=qe.charge(setup,a,q,f)
        self.assertEqual(result['metrics']['T'],setup.repetitions*(5+(setup.M-1)*20+123))
        self.assertEqual(result['metrics']['H'],setup.repetitions*(2+(setup.M-1)*4+2*setup.phase_register_bits))
        self.assertEqual(result['peak_logical_qubits'],25+setup.phase_register_bits)
        self.assertEqual(result['phase_register_measurements'],setup.repetitions*setup.phase_register_bits)
        bad=dict(f,norm_error_upper='1/63')
        with self.assertRaises(ValueError):
            qe.charge(setup,a,q,bad)

    def test_good_outcome_probability(self):
        # This independent sweep is a diagnostic. The probability guarantee
        # follows from the published theorem and norm transfer, not the sweep.
        for M in (4,8,16,32,64):
            for j in range(101):
                p=j/100
                law=qe.ideal_amplitude_distribution(p,M)
                radius=2*math.pi*math.sqrt(p*(1-p))/M+math.pi**2/M**2
                good=sum(w for y,w in enumerate(law)
                         if abs(math.sin(math.pi*y/M)**2-p)<=radius+1e-14)
                self.assertGreaterEqual(good+1e-12,8/math.pi**2)


if __name__=='__main__':
    unittest.main()
