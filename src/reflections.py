"""Explicit diagonal phase reflection with clean AND scratch.

The separate phase-aware gate list complements the Boolean arithmetic engine.
X, CCX, Z and CZ are exact logical gates. No arbitrary angle is synthesised.
"""
from collections import Counter
import json
from pathlib import Path


def zero_reflection(qubits):
    if qubits < 2:
        raise ValueError('This constructor requires at least two inputs')
    work=list(range(qubits,2*qubits-2))
    gates=[('X',i) for i in range(qubits)]
    conjunction=[]
    if qubits==2:
        gates.append(('CZ',0,1))
    else:
        conjunction.append(('CCX',0,1,work[0]))
        for j in range(2,qubits-1):
            conjunction.append(('CCX',work[j-2],j,work[j-1]))
        gates+=conjunction
        gates.append(('CZ',work[-1],qubits-1))
        gates+=list(reversed(conjunction))
    gates += [('X',i) for i in range(qubits)]
    return gates


def resources(qubits):
    gates=zero_reflection(qubits)
    counts=Counter(g[0] for g in gates)
    layers=[0]*(2*qubits-2)
    for gate in gates:
        wires=gate[1:]
        depth=max(layers[w] for w in wires)+(gate[0]=='CCX')
        for w in wires:layers[w]=depth
    return {'gate_counts':dict(counts),'logical_qubits':2*qubits-2,
            'toffoli_layer_depth':max(layers),
            'T_count_exact_7T':7*counts['CCX'],
            'T_depth_3layer_upper':3*max(layers),
            'explicit_phase_gate_list':True}


def check():
    checked=0
    for n in range(2,9):
        for value in range(1 << n):
            state=[(value>>j)&1 for j in range(n)]+[0]*(n-2)
            phase=1
            for gate in zero_reflection(n):
                if gate[0]=='X':state[gate[1]]^=1
                elif gate[0]=='CCX':
                    _,a,b,t=gate;state[t]^=state[a]&state[b]
                elif gate[0]=='CZ':
                    _,a,b=gate
                    if state[a]&state[b]:phase=-phase
                else:raise ValueError(gate)
            assert sum(state[j]<<j for j in range(n))==value
            assert all(v==0 for v in state[n:])
            assert phase==(-1 if value==0 else 1)
            checked+=1
    return {'status':'exact diagonal zero-reflection replay passed',
            'basis_checks':checked,'tested_input_qubits':[2,3,4,5,6,7,8],
            'all_inputs_preserved':True,'all_work_zero':True,
            'phase_negative_only_on_zero_input':True}


if __name__=='__main__':
    path=Path(__file__).resolve().parent/'results'/'reflection_checks.json'
    path.write_text(json.dumps(check(),indent=2)+'\n');print(path)
