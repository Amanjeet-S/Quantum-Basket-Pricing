"""Independent packed basis replay and literal gate accounting for audits."""
from collections import Counter
import hashlib
import math


def packed_values(values,bits):
    return [sum(((int(v)>>j)&1)<<i for i,v in enumerate(values)) for j in range(bits)]


def register_values(state,register,lanes):
    return [sum(((state[w]>>i)&1)<<j for j,w in enumerate(register)) for i in range(lanes)]


def packed_replay(circuit,assignments,*,reverse=False,initial_state=None):
    """Replay all basis lanes together using Python integer bit planes.

    This does not call Circuit.run or any graph evaluation method. A lookup
    is interpreted from its retained table and address, after checking that
    its private selector work is zero in every lane. Separate state-vector
    tests establish phase correctness of its expanded unitary operations.
    """
    lanes=len(assignments);mask=(1<<lanes)-1
    state=[0]*circuit.wire_count if initial_state is None else list(initial_state)
    if initial_state is None:
        for lane,assignment in enumerate(assignments):
            used=set()
            for register,value in assignment.items():
                if not isinstance(value,int) or value<0 or value>=1<<len(register):raise ValueError('Invalid basis assignment')
                if used.intersection(register):raise ValueError('Initial registers overlap')
                used.update(register)
                for j,w in enumerate(register):
                    state[w]|=((value>>j)&1)<<lane
    tokens=list(circuit.gates)
    if reverse:tokens=tokens[::-1]
    calls=0
    for op,c1,c2,target in tokens:
        if op==1:state[target]^=mask
        elif op==2:state[target]^=state[c1]
        elif op==3:state[target]^=state[c1]&state[c2]
        elif op==20:
            calls+=1;plan=circuit._macros[c1]
            assert all(state[w]==0 for w in plan.work),'Lookup selector work is dirty'
            # Address predicates are evaluated independently of the compiler.
            for address,word in enumerate(plan.table):
                if not word:continue
                selection=mask if plan.enable is None else state[plan.enable]
                for bit,wire in enumerate(plan.address):
                    selection&=state[wire] if address>>bit&1 else mask^state[wire]
                if not selection:continue
                for bit,wire in enumerate(plan.output):
                    if word>>bit&1:state[wire]^=selection
        else:raise ValueError(('Unknown logical tape operation',op))
    return state,{'lanes':lanes,'logical_tape_tokens':len(tokens),'lookup_invocations':calls}


def expanded_operations(circuit):
    for op,c1,c2,target in circuit.gates:
        if op==1:yield 'X',(target,)
        elif op==2:yield 'CX',(c1,target)
        elif op==3:yield 'CCX',(c1,c2,target)
        elif op==20:
            for gate,wires,condition in circuit._macros[c1].operations():
                if condition is not None:raise ValueError('Unitary audit does not invert measured operations')
                yield gate,wires
        else:raise ValueError(op)


def literal_seven_t(op,w):
    if op!='CCX':yield op,w;return
    a,b,t=w
    yield 'H',(t,);yield 'CX',(b,t);yield 'TD',(t,);yield 'CX',(a,t)
    yield 'T',(t,);yield 'CX',(b,t);yield 'TD',(t,);yield 'CX',(a,t)
    yield 'T',(b,);yield 'T',(t,);yield 'H',(t,);yield 'CX',(a,b)
    yield 'T',(a,);yield 'TD',(b,);yield 'CX',(a,b)


def independent_gate_accounting(circuit):
    primitive=Counter();literal=Counter();td=[0]*circuit.wire_count;ad=[0]*circuit.wire_count
    h=hashlib.sha256()
    for op,wires in expanded_operations(circuit):
        primitive[op]+=1;h.update((op+':'+','.join(map(str,wires))+';').encode())
        for exact,physical in literal_seven_t(op,wires):
            literal[exact]+=1
            phase=max((td[w] for w in physical),default=0)+(exact in ('T','TD'))
            layer=max((ad[w] for w in physical),default=0)+1
            for w in physical:td[w]=phase;ad[w]=layer
    assert literal['T']+literal['TD']==primitive['T']+primitive['TD']+7*primitive['CCX']
    return {'primitive_gate_counts':dict(primitive),'expanded_literal_gate_counts':dict(literal),
        'logical_T_count':literal['T']+literal['TD'],
        'actual_literal_T_dependency_depth':max(td,default=0),
        'actual_literal_all_gate_dependency_depth':max(ad,default=0),
        'expanded_unitary_tape_sha256':h.hexdigest(),'logical_wires':circuit.wire_count,
        'inverse_T_count':literal['T']+literal['TD'],
        'inverse_count_reason':'Literal unitary adjoint swaps T and TD, S and SD and reverses the emitted sequence'}
