"""Independent dense clean-work isometry replay for small primitive lookups."""
import cmath
import hashlib
import math
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))

from integral_audit_support import literal_seven_t


def apply_operations(operations,qubits,columns):
    state=np.zeros((1<<qubits,len(columns)),dtype=np.complex128)
    for col,index in enumerate(columns):state[index,col]=1
    index=np.arange(1<<qubits,dtype=np.uint64)
    for operation,wires in operations:
        if operation=='CCX':
            expanded=list(literal_seven_t(operation,wires))
        else:expanded=[(operation,wires)]
        for op,w in expanded:
            if op in ('X','CX'):
                target=1<<w[-1];selection=(index & target)==0
                for control in w[:-1]:selection&=((index>>np.uint64(control))&np.uint64(1)).astype(bool)
                left=np.flatnonzero(selection);right=left|target
                copy=state[left].copy();state[left]=state[right];state[right]=copy
            elif op=='H':
                target=1<<w[0];left=np.flatnonzero((index&np.uint64(target))==0);right=left|target
                a=state[left].copy();b=state[right].copy()
                state[left]=(a+b)/math.sqrt(2);state[right]=(a-b)/math.sqrt(2)
            elif op in ('T','TD','S','SD'):
                phase={'T':cmath.exp(1j*math.pi/4),'TD':cmath.exp(-1j*math.pi/4),'S':1j,'SD':-1j}[op]
                selection=((index>>np.uint64(w[0]))&np.uint64(1)).astype(bool)
                state[selection]*=phase
            elif op in ('Z','CZ'):
                selection=np.ones(len(index),dtype=bool)
                for wire in w:selection&=((index>>np.uint64(wire))&np.uint64(1)).astype(bool)
                state[selection]*=-1
            else:raise ValueError(('Unsupported unitary',op))
    return state


def lookup_clean_isometry(plan):
    registers=plan.address+plan.output+plan.work+(() if plan.enable is None else (plan.enable,))
    qubits=max(registers,default=-1)+1
    if qubits>13:raise ValueError('Retained dense audit is restricted to at most 13 physical wires')
    columns=[];expected=[]
    for enable in ((1,) if plan.enable is None else (0,1)):
        for address in range(1<<len(plan.address)):
            for target in range(1<<len(plan.output)):
                initial=sum((address>>i&1)<<w for i,w in enumerate(plan.address))
                initial|=sum((target>>i&1)<<w for i,w in enumerate(plan.output))
                if plan.enable is not None:initial|=enable<<plan.enable
                value=plan.table[address] if enable and address<len(plan.table) else 0
                final=initial
                for bit,wire in enumerate(plan.output):
                    if value>>bit&1:final^=1<<wire
                columns.append(initial);expected.append(final)
    operations=[(op,w) for op,w,condition in plan.operations()]
    assert all(condition is None for op,w,condition in plan.operations())
    matrix=apply_operations(operations,qubits,columns)
    reference=np.zeros_like(matrix)
    reference[expected,np.arange(len(expected))]=1
    error=float(np.max(np.abs(matrix-reference)))
    assert error<2e-12,(plan.mode,plan.table,error)
    # A general coherent superposition is an additional linear-combination
    # diagnostic, including relative phases across dirty target values.
    coefficients=np.array([cmath.exp(1j*(j+1)*.173)/(j+1) for j in range(len(columns))])
    coefficients/=np.sqrt(np.vdot(coefficients,coefficients).real)
    coherent_error=float(np.max(np.abs(matrix@coefficients-reference@coefficients)))
    assert coherent_error<2e-12
    inverse={'T':'TD','TD':'T','S':'SD','SD':'S'}
    adjoint=[(inverse.get(op,op),w) for op,w in reversed(operations)]
    reverse_matrix=apply_operations(adjoint,qubits,expected)
    reverse_reference=np.zeros_like(reverse_matrix);reverse_reference[columns,np.arange(len(columns))]=1
    reverse_error=float(np.max(np.abs(reverse_matrix-reverse_reference)))
    assert reverse_error<2e-12
    return {'mode':plan.mode,'address_bits':len(plan.address),'dirty_output_bits':len(plan.output),
        'selector_work_bits':len(plan.work),'enabled':plan.enable is not None,
        'basis_columns':len(columns),'maximum_forward_amplitude_error':error,
        'maximum_coherent_amplitude_error':coherent_error,
        'maximum_literal_adjoint_amplitude_error':reverse_error,
        'scope':'Numerical dense primitive replay; exact logical interpretation separately proved'}
