"""Exact logical Clifford/Toffoli Grover and controlled-Grover templates.

State-vector replay is a small numerical diagnostic. The gates, their counts
and the control decomposition are explicit; no arbitrary rotation is used.
"""
from collections import Counter
import cmath
import math


def zero_phase(core,work,control=None):
    """I-2|0><0|, optionally controlled, with clean conjunction work."""
    core=tuple(core);controls=core+(() if control is None else (control,))
    if len(controls)<2 or len(work)<len(controls)-2:raise ValueError('Insufficient reflection wires')
    if len(set(core+tuple(work)+(() if control is None else (control,))))!=len(core)+len(work)+(control is not None):
        raise ValueError('Reflection wires overlap')
    gates=[('X',w) for w in core];body=[]
    if len(controls)==2:gates.append(('CZ',*controls))
    else:
        body.append(('CCX',controls[0],controls[1],work[0]))
        for j in range(2,len(controls)-1):body.append(('CCX',work[j-2],controls[j],work[j-1]))
        gates+=body;gates.append(('CZ',work[len(controls)-3],controls[-1]));gates+=body[::-1]
    gates += [('X',w) for w in core]
    return gates


def adjoint(gates):
    swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    return [(swap.get(g[0],g[0]),*g[1:]) for g in reversed(gates)]


def grover(A,core,flag,work,control=None):
    """Q=-A S0 A† S_good; the uncontrolled global minus is immaterial.

    A acts without control, so arithmetic Toffolis are never all upgraded to
    controlled gates. Controlled Q uses CZ, A†, controlled S0, A, Z(control).
    The last Z supplies the conditional global-minus sign.
    """
    if control is None:
        return [('Z',flag)]+adjoint(A)+zero_phase(core,work)+list(A)
    return [('CZ',control,flag)]+adjoint(A)+zero_phase(core,work,control)+list(A)+[('Z',control)]


def expand_toffoli(gates):
    """Expand each exact CCX as seven T/Tdg, six CX and two H gates.

    This explicit chronological sequence has T dependency depth four. The
    constructor's historical assigned three-layer profile is kept separate.
    """
    output=[]
    for gate in gates:
        if gate[0]!='CCX':output.append(gate);continue
        _,a,b,t=gate
        output += [('H',t),('CX',b,t),('Tdg',t),('CX',a,t),('T',t),
                   ('CX',b,t),('Tdg',t),('CX',a,t),('T',b),('T',t),
                   ('H',t),('CX',a,b),('T',a),('Tdg',b),('CX',a,b)]
    return output


def resources(gates,qubits):
    counts=Counter(g[0] for g in gates);layers=[0]*qubits;all_layers=[0]*qubits
    for gate in gates:
        wires=gate[1:]
        ld=max(layers[w] for w in wires)+(3 if gate[0]=='CCX' else 1 if gate[0] in ('T','Tdg') else 0)
        ad=max(all_layers[w] for w in wires)+1
        for w in wires:layers[w]=ld;all_layers[w]=ad
    return {'gate_counts':dict(counts),'logical_qubits':qubits,
            'T_count_exact_7T_Toffoli':7*counts['CCX']+counts['T']+counts['Tdg'],
            'T_depth_assigned_3layer_Toffoli_or_explicit_T_upper':max(layers),
            'all_gate_dependency_depth':max(all_layers),'arbitrary_rotations':0}


def replay(gates,state):
    """Numerical state-vector replay using only explicitly emitted gates."""
    state=list(state);size=len(state)
    if size==0 or size&(size-1):raise ValueError('State dimension is not a power of two')
    root2=math.sqrt(2)
    for gate in gates:
        kind,*wires=gate
        if kind=='H':
            bit=1<<wires[0]
            for i in range(size):
                if i&bit:continue
                a,b=state[i],state[i|bit];state[i]=(a+b)/root2;state[i|bit]=(a-b)/root2
        elif kind in ('X','CX','CCX'):
            target=1<<wires[-1];controls=wires[:-1]
            for i in range(size):
                if i&target or not all((i>>w)&1 for w in controls):continue
                state[i],state[i|target]=state[i|target],state[i]
        elif kind in ('Z','CZ'):
            for i in range(size):
                if all((i>>w)&1 for w in wires):state[i]=-state[i]
        elif kind in ('S','Sdg','T','Tdg'):
            phase={'S':1j,'Sdg':-1j,'T':cmath.exp(1j*math.pi/4),'Tdg':cmath.exp(-1j*math.pi/4)}[kind]
            for i in range(size):
                if (i>>wires[0])&1:state[i]*=phase
        else:raise ValueError(kind)
    return state


def zero_state(qubits):
    state=[0j]*(1<<qubits);state[0]=1+0j;return state


def flag_probability(state,flag):
    return sum(abs(z)**2 for i,z in enumerate(state) if(i>>flag)&1)
