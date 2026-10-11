"""Independent basis permutations, exact PMFs and literal phase diagnostics."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import hashlib,json,math,sys
import numpy as np
SOURCE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE))
import beta_minimum as producer
ROOT=Path(__file__).resolve().parents[1]


def initial_index(registers,values,output,dirty):
    state=0
    for reg,value in list(zip(registers,values))+[(output,dirty)]:
        for k,wire in enumerate(reg):state|=((value>>k)&1)<<wire
    return state


def integer_replay(tape,state):
    for op,a,b,t in tape:
        if op==1 or (op==2 and state&(1<<a)) or (op==3 and state&(1<<a) and state&(1<<b)):
            state^=1<<t
    return state


def decode(state,reg):return sum(((state>>w)&1)<<k for k,w in enumerate(reg))


def raw_to_literal(tape):
    output=[]
    for op,a,b,t in tape:
        if op==1:output.append(('X',t))
        elif op==2:output.append(('CX',a,t))
        else:
            # Standalone chronological exact seven-T decomposition.
            output.extend([('H',t),('CX',b,t),('Tdg',t),('CX',a,t),('T',t),
                           ('CX',b,t),('Tdg',t),('CX',a,t),('T',b),('T',t),
                           ('H',t),('CX',a,b),('T',a),('Tdg',b),('CX',a,b)])
    return output


def literal_replay(tape,state):
    state=state.copy();indices=np.arange(len(state))
    for gate in tape:
        kind,*w=gate
        if kind=='H':
            a=indices[(indices&(1<<w[0]))==0];b=a|(1<<w[0]);x,y=state[a].copy(),state[b].copy()
            state[a],state[b]=(x+y)/math.sqrt(2),(x-y)/math.sqrt(2)
        elif kind in ('X','CX'):
            a=indices[(indices&(1<<w[-1]))==0]
            for control in w[:-1]:a=a[(a&(1<<control))!=0]
            b=a|(1<<w[-1]);state[a],state[b]=state[b].copy(),state[a].copy()
        else:
            phase=complex(math.cos(math.pi/4),math.sin(math.pi/4))
            if kind=='Tdg':phase=phase.conjugate()
            state[(indices&(1<<w[0]))!=0]*=phase
    return state


def main():
    checks=0;pmfs=[]
    for r,m in [(1,1),(1,3),(2,1),(2,2),(2,3),(3,2),(3,3),(4,1),(4,2),(4,3)]:
        c,inputs,out=producer.constructor(r,m);tape=list(c.gates);N=1<<m;hist=[0]*N
        keep_mask=sum(1<<w for reg in list(inputs)+[out] for w in reg)
        for values in product(range(N),repeat=r):
            for dirty in range(N):
                original=initial_index(inputs,values,out,dirty);actual=integer_replay(tape,original)
                assert decode(actual,out)==(dirty^min(values))
                assert all(decode(actual,reg)==v for reg,v in zip(inputs,values))
                assert actual&~keep_mask==0
                checks+=1
            hist[min(values)]+=1
        expected=tuple(F((N-j)**r-(N-j-1)**r,N**r) for j in range(N))
        assert tuple(F(n,N**r) for n in hist)==expected==producer.probabilities(r,m)
        pmfs.append({'r':r,'m':m,'exact_probabilities':list(map(str,expected)),'passed':True})
    profiles=[]
    for r in range(1,9):
        for m in [1,2,4,8,16]:
            c,inputs,out=producer.constructor(r,m);raw=Counter(g[0] for g in c.gates)
            expected=producer.profile(r,m);ccx=raw[3]
            assert expected['T_exact_7T']==7*ccx
            assert expected['expanded_Clifford_T_operations']==raw[1]+raw[2]+15*ccx+r*m
            assert expected['explicit_serial_T_depth_upper']==4*ccx
            assert expected['conservative_serial_T_depth_upper']==5*ccx
            assert expected['raw_resources']['wires']==c.wire_count
            assert expected['raw_resources']['active_wires_at_end']==(r+1)*m
            assert expected['coupling_radius']==str('1/'+str(2**(m+1)))
            profiles.append({'r':r,'m':m,'raw_gate_counts':dict(raw),'wires':c.wire_count,
                             'retained_input_output_qubits':(r+1)*m,'T_count':7*ccx,'passed':True})
    coherent=[]
    for r,m in [(2,1),(3,1),(2,2)]:
        c,inputs,out=producer.constructor(r,m);raw=list(c.gates);literal=raw_to_literal(raw)
        # An arbitrary complex logical superposition entangled with a reference.
        dimension=1<<(c.wire_count+1);state=np.zeros(dimension,complex);wanted=np.zeros(dimension,complex)
        logical=list(product(range(1<<m),repeat=r+1))
        amplitudes=np.array([complex(math.sin(k+1),math.cos(2*k+3)) for k in range(2*len(logical))])
        amplitudes/=np.linalg.norm(amplitudes)
        for k,values in enumerate(logical):
            for reference in (0,1):
                start=initial_index(inputs,values[:-1],out,values[-1])|(reference<<c.wire_count)
                end=initial_index(inputs,values[:-1],out,values[-1]^min(values[:-1]))|(reference<<c.wire_count)
                state[start]=amplitudes[2*k+reference];wanted[end]=amplitudes[2*k+reference]
        actual=literal_replay(literal,state);error=float(np.linalg.norm(actual-wanted))
        assert error<1e-12
        # Actual Hadamard preparation of original uniforms, followed by minimum.
        uniform=np.zeros(1<<c.wire_count,complex);uniform[0]=1
        prep=[('H',w) for reg in inputs for w in reg]+literal
        uniform=literal_replay(prep,uniform);pmf=[sum(abs(a)**2 for i,a in enumerate(uniform) if decode(i,out)==j) for j in range(1<<m)]
        error_p=max(abs(p-float(q)) for p,q in zip(pmf,producer.probabilities(r,m)))
        assert error_p<1e-12
        coherent.append({'r':r,'m':m,'full_complex_reference_state_error_diagnostic':error,
                         'actual_uniform_preparation_PMF_error_diagnostic':error_p,
                         'literal_primitive_gate_count':len(literal),'passed':True})
    result={'status':'passed','independent_dirty_basis_cases':checks,'exact_PMF_cases':pmfs,
            'emitted_profiles':profiles,'coherent_primitive_diagnostics':coherent,
            'source_sha256':hashlib.sha256((SOURCE/'beta_minimum.py').read_bytes()).hexdigest(),
            'dependency_sha256':hashlib.sha256((SOURCE/'reversible.py').read_bytes()).hexdigest(),
            'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact basis permutations and rational PMFs are exhaustive over the stated cases. Floating primitive-phase replay is a diagnostic; the standard exact CCX identity and linearity supply the mathematical coherent extension.'}
    (ROOT/'results/beta_minimum_independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','basis_cases':checks,'profiles':len(profiles),'coherent_cases':len(coherent)}))


if __name__=='__main__':main()
