"""Meaningful finite gate, phase, convention and certificate regression checks."""
from pathlib import Path
import sys,unittest,cmath,math,json
from fractions import Fraction
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import qft_synthesis as q

def product(tape,n):
    N=1<<n;U=[[complex(int(x==y))for y in range(N)]for x in range(N)]
    for op in tape:
        if op[0]=='CX':
            cb=1<<(n-1-op[1]);tb=1<<(n-1-op[2])
            for x in range(N):
                if x&cb and not x&tb:U[x],U[x|tb]=U[x|tb],U[x]
        else:
            g,w=op;tb=1<<(n-1-w)
            if g=='H':
                for x in range(N):
                    if x&tb:continue
                    y=x|tb;a=U[x][:];b=U[y][:]
                    U[x]=[(a[j]+b[j])/math.sqrt(2)for j in range(N)]
                    U[y]=[(a[j]-b[j])/math.sqrt(2)for j in range(N)]
            else:
                phase=cmath.exp(1j*math.pi*float({'T':Fraction(1,4),'Tdg':Fraction(-1,4),'S':Fraction(1,2),'Sdg':Fraction(-1,2)}[g]))
                for x in range(N):
                    if x&tb:U[x]=[phase*a for a in U[x]]
    return U

class QFTSynthesis(unittest.TestCase):
    def test_width_bounds(self):
        for m in (0,33,1.5):
            with self.assertRaises(ValueError):q.inverse_qft_profile(m)
        with self.assertRaises(ValueError):q.inverse_qft_profile(32,Fraction(1,10**9))
    def test_exact_T_phase(self):
        r=q.rotation(2,1)
        self.assertEqual(r['chronological_gates'],('T',));self.assertEqual(r['phase_pi'],Fraction(-1,8))
    def test_negative_rotations_are_adjoints(self):
        adj={'H':'H','S':'Sdg','T':'Tdg'}
        for k in range(2,33):
            pos=q.rotation(k,1);neg=q.rotation(k,-1)
            self.assertEqual(neg['chronological_gates'],tuple(adj[g]for g in reversed(pos['chronological_gates'])))
            self.assertEqual(neg['phase_pi'],-pos['phase_pi'])
    def test_tiny_angles_have_charged_error(self):
        for k in range(19,33):
            self.assertFalse(q.rotation(k,1)['chronological_gates'])
            self.assertGreater(q.rotation(k,1)['norm_error_upper'],0)
            cp=q.controlled_phase(k-1,0,1)
            self.assertFalse(cp['gates']);self.assertGreater(cp['norm_error_upper'],0)
    def test_all_profile_counts_literal_tapes(self):
        for m in range(1,33):
            tape=q.inverse_qft_tape(m);p=q.inverse_qft_profile(m)
            c=Counter(g[0]for g in tape['gates'])
            self.assertEqual(dict(c),p['gate_counts'])
            self.assertEqual(c['T']+c['Tdg'],p['metrics']['T'])
            self.assertEqual(sum(c.values()),p['metrics']['expanded_primitive_operations'])
            self.assertLessEqual(p['norm_error_upper'],Fraction(1,64))
            self.assertEqual(p['extra_ancillas'],0)
    def test_exact_controlled_S_globalphase(self):
        cp=q.controlled_phase(1,0,1,-1)
        U=product(cp['gates'],2);phase=cmath.exp(1j*math.pi*float(cp['global_phase_correction_pi']))
        target=[1,1,1,-1j]
        self.assertLess(max(abs(phase*U[x][y]-(target[x]if x==y else 0))for x in range(4)for y in range(4)),1e-14)
    def test_coherent_inverse_dft_small_registers(self):
        for m in range(1,4):
            tape=q.inverse_qft_tape(m);U=product(tape['gates'],m);phase=cmath.exp(1j*math.pi*float(tape['global_phase_correction_pi']));N=1<<m
            max_entry=max(abs(phase*U[x][y]-cmath.exp(-2j*math.pi*x*y/N)/math.sqrt(N))for x in range(N)for y in range(N))
            self.assertLess(max_entry,float(tape['norm_error_upper'])+1e-13)
    def test_manifest_hash_and_rational_certificates(self):
        data,bounds=q._data()
        self.assertEqual(len(bounds),62)
        for value in bounds.values():self.assertLess(value,Fraction(1,2**17))
        self.assertEqual(data['generator']['pygridsynth'],'1.2.0')

if __name__=='__main__':unittest.main()
