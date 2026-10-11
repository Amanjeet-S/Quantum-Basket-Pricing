"""Independent exact clean-predicate and Clifford+T phase audit.

The shared-threshold prototype is imported solely to execute emit_interval.
Its numerical-finance results and resource producer are not imported. Exact
phase replay takes place in Q(zeta_8), represented as four Fraction entries;
no floating tolerances are used. Exhaustive Boolean checks include unordered
endpoints and equality cases, arbitrary initial good flags and clean work.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from itertools import product
from pathlib import Path
sys.dont_write_bytecode=True
ZERO=(F(0),)*4
ONE=(F(1),F(0),F(0),F(0))
ROOT_HALF=(F(0),F(1,2),F(0),F(-1,2))  # 1/sqrt(2)


def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
    v=[F(0)]*7
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]+=x*y
    for k in range(6,3,-1):v[k-4]-=v[k]
    return tuple(v[:4])

def zeta(power):
    power%=8
    v=[F(0)]*4;v[power%4]=F(1 if power<4 else -1)
    return tuple(v)


def replay_phase(gates,state):
    state=dict(state)
    def put(out,key,value):
        new=add(out.get(key,ZERO),value)
        if new==ZERO:out.pop(key,None)
        else:out[key]=new
    for gate in gates:
        kind,*w=gate;out={}
        if kind=='H':
            bit=1<<w[0]
            for index,value in state.items():
                half=mul(ROOT_HALF,value)
                put(out,index & ~bit,half)
                put(out,index | bit,neg(half) if index & bit else half)
        elif kind in ('X','CX','CCX'):
            for index,value in state.items():
                target=index^(1<<w[-1]) if all((index>>c)&1 for c in w[:-1]) else index
                put(out,target,value)
        elif kind in ('T','Tdg','S','Sdg','Z'):
            power={'T':1,'Tdg':-1,'S':2,'Sdg':-2,'Z':4}[kind]
            for index,value in state.items():
                put(out,index,mul(value,zeta(power)) if (index>>w[0])&1 else value)
        elif kind=='CZ':
            for index,value in state.items():
                put(out,index,neg(value) if all((index>>c)&1 for c in w) else value)
        else:raise ValueError(kind)
        state=out
    return state


def replay_boolean(c,state):
    state=list(state)
    for op,a,b,t in c.gates:
        if op==1:state[t]^=1
        elif op==2:state[t]^=state[a]
        else:
            assert op==3;state[t]^=state[a]&state[b]
    return state


def physical(c,qg):
    gates=[]
    for op,a,b,t in c.gates:
        gates.append(('X',t) if op==1 else ('CX',a,t) if op==2 else ('CCX',a,b,t))
    return qg.expand_toffoli(gates)


def fill(state,register,value):
    for i,w in enumerate(register):state[w]=(value>>i)&1


def index(state):return sum(v<<w for w,v in enumerate(state))

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prototype',type=Path,required=True)
    ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--factor-resources',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    spec=importlib.util.spec_from_file_location('prototype_under_review',args.prototype)
    prototype=importlib.util.module_from_spec(spec);spec.loader.exec_module(prototype)
    prototype.configure(args.repo,args.factor_resources,args.output.parent/'configuration_only')
    base,qg=prototype.base,prototype.qg
    rows=[];phasecases=0;boolcases=0
    for width in (1,2):
        c=base.Circuit(store_gates=True)
        lower,threshold,upper=(c.alloc(width) for _ in range(3));flag=c.alloc(1)[0]
        original_active=c.active_wires
        prototype.emit_interval(c,lower,threshold,upper,flag)
        assert c.active_wires==original_active
        primitive=physical(c,qg)
        for a,k,b,g in product(range(1<<width),range(1<<width),range(1<<width),(0,1)):
            state=[0]*c.wire_count
            fill(state,lower,a);fill(state,threshold,k);fill(state,upper,b);state[flag]=g
            expected=list(state);expected[flag]=g ^ int(a<k<=b)
            observed=replay_boolean(c,state)
            assert observed==expected,(width,a,k,b,g)
            boolcases+=1
            if width==1:
                out=replay_phase(primitive,{index(state):ONE})
                assert out=={index(expected):ONE},(width,a,k,b,g,out)
                phasecases+=1
        if width==1:
            initial={};expected={}
            coefficient=(F(0),F(1,4),F(0),F(-1,4))  # 1/sqrt(8)
            for a,k,b in product(range(2),repeat=3):
                state=[0]*c.wire_count;fill(state,lower,a);fill(state,threshold,k);fill(state,upper,b)
                initial[index(state)]=coefficient
                state[flag]=int(a<k<=b);expected[index(state)]=coefficient
            assert replay_phase(primitive,initial)==expected
        rows.append({'width':width,'wire_count':c.wire_count,'logical_gates':c.gate_count,
                     'explicit_primitive_gates':len(primitive),
                     'explicit_T_count':sum(g[0] in ('T','Tdg') for g in primitive),
                     'chronological_primitive_sha256':hashlib.sha256(json.dumps(primitive,separators=(',',':')).encode()).hexdigest()})
    # An exact abstract-CDF check independent of the financial implementation.
    cdfchecks=0
    for n in (1,2,4,16,1024):
        for a in (F(1,2),F(3,4),F(999999,1000000),F(1)):
            tail=(1-a)/2;D=tail+a/(2*n)
            for j in range(n):
                p=tail+a*F(2*j+1,2*n)
                assert abs(F(j,n)-p)<=D and abs(F(j+1,n)-p)<=D
                cdfchecks+=2
            assert abs(tail+a/F(2*n))==D
    result={'schema':'qbp-independent-shared-threshold-clean-predicate-v1','status':'passed',
        'exhaustive_boolean_cases':boolcases,'exact_cyclotomic_phase_basis_cases':phasecases,
        'exact_coherent_superposition_cases':1,'abstract_left_right_CDF_checks':cdfchecks,
        'tested_profiles':rows,'prototype_source_sha256':sha(args.prototype),'audit_source_sha256':sha(__file__),
        'reversible_source_sha256':sha(args.repo/'reversible.py'),'quantum_gates_source_sha256':sha(args.repo/'quantum_gates.py'),
        'scope':'Exact predicate semantics, equality endpoints, work cleanup and literal seven-T primitive phases at small register widths. The finite-CDF algebra is checked abstractly. No financial endpoint certificate, full financial compiler, resource outcome or practical advantage is asserted.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','exhaustive_boolean_cases','exact_cyclotomic_phase_basis_cases','abstract_left_right_CDF_checks')}))

if __name__=='__main__':main()
