"""Exact coherent auxiliary register for an integral-remainder price oracle.

Hadamards prepare r uniform m-bit registers. A reversible minimum produces
the midpoint discretisation of Beta(1,r). All private work is erased, and
the original uniform labels are retained for subsequent amplitude estimation.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json
from reversible import Circuit
ROOT=Path(__file__).resolve().parent

def constructor(r,m,store_gates=True):
    if isinstance(r,bool) or isinstance(m,bool) or not isinstance(r,int) or not isinstance(m,int) or r<1 or m<1:
        raise ValueError('Positive integer order and width required')
    c=Circuit(store_gates=store_gates)
    inputs=[c.alloc(m,f'uniform register {i}') for i in range(r)]
    snap=c.snapshot();current=inputs[0]
    for other in inputs[1:]:
        smaller=c.unsigned_less_than(current,other)[0]
        difference=c.alloc(m,'selection difference');out=c.alloc(m,'partial minimum')
        for a,b,d,o in zip(current,other,difference,out):
            c.cnot(a,d);c.cnot(b,d);c.cnot(b,o);c.toffoli(smaller,d,o)
        current=out
    result=c.copy_uncompute(snap,current)
    return c,tuple(inputs),result

def probabilities(r,m):
    N=1<<m
    return tuple(F((N-j)**r-(N-j-1)**r,N**r) for j in range(N))

def profile(r,m):
    c,inputs,out=constructor(r,m,False);resource=c.resources();ccx=resource['toffoli']
    return {'r':r,'midpoint_bits':m,'uniform_Hadamards':r*m,'minimum_register_bits':m,
      'raw_resources':resource,'T_exact_7T':7*ccx,
      'explicit_serial_T_depth_upper':4*ccx,'conservative_serial_T_depth_upper':5*ccx,
      'expanded_Clifford_T_operations':resource['x']+resource['cnot']+15*ccx+r*m,
      'coupling_radius':'1/'+str(1<<(m+1)),
      'scope':'Exact auxiliary preparation only; financial kernel, probability flag, estimator and final reset are additional costs.'}

def audit():
    checked=0;distribution_cases=[]
    for r,m in [(1,1),(1,3),(2,1),(2,2),(2,3),(3,2),(3,3),(4,1),(4,2),(4,3)]:
        c,inputs,out=constructor(r,m);hist=[0]*(1<<m)
        keep=set(out)|{w for reg in inputs for w in reg}
        for values in product(range(1<<m),repeat=r):
            for dirty in [0,(1<<m)-1]:
                state=c.run({**dict(zip(inputs,values)),out:dirty})
                assert c.value(state,out)==(min(values)^dirty)
                assert all(c.value(state,reg)==v for reg,v in zip(inputs,values))
                assert all(not state[w] for w in range(c.wire_count) if w not in keep)
                checked+=1
            hist[min(values)]+=1
        N=1<<m;expected=[(N-j)**r-(N-j-1)**r for j in range(N)]
        assert hist==expected and sum(probabilities(r,m),F(0))==1
        distribution_cases.append({'r':r,'m':m,'histogram':hist,'exact_minimum_distribution':True})
    profiles=[profile(r,m) for r in range(1,9) for m in [1,2,4,8,16]]
    result={'schema':'qbp-beta-minimum-circuit-v1','status':'passed','basis_checks_including_dirty_output':checked,
       'distribution_cases':distribution_cases,'profiles':profiles,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'dependency_sha256':hashlib.sha256((ROOT/'reversible.py').read_bytes()).hexdigest(),
       'proof_scope':'Basis permutation with exact phases implies the coherent action by linearity; uniform input labels remain retained. Full financial oracle is a separate obligation.'}
    (ROOT/'results/beta_minimum_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','basis_checks':checked,'profiles':len(profiles),
          'r2_m8_profile':profile(2,8)},indent=2))
if __name__=='__main__':audit()
