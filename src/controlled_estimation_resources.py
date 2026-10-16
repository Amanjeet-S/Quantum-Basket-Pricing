"""Emitted controlled Grover reflection and complete QPE resource charges.

The conditional global minus is a physical Z on the phase control. A/A-dagger
run unconditionally around controlled reflections; corrected instruments
remain coherent logical channels. No controlled financial measurement tape
is assumed. Reflection work starts and ends in zero.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import value_oracle as vo
import reflections
import qpe_estimator
import qft_synthesis

def controlled_zero_reflection(original_qubits):
    if isinstance(original_qubits,bool) or not isinstance(original_qubits,int) or original_qubits<2:raise ValueError('At least two original wires required')
    control=original_qubits
    # The extra phase control is tested for one, so it must not be negated.
    return tuple(g for g in reflections.zero_reflection(original_qubits+1) if not(g[0]=='X' and g[1]==control))

def expansion(counts):
    out=Counter()
    for gate,n in counts.items():
        if gate=='CCX':out.update({'H':2*n,'CX':6*n,'T':4*n,'Tdg':3*n})
        else:out[gate]+=n
    return out

@lru_cache(None)
def reflection_metrics(original_qubits):
    gates=controlled_zero_reflection(original_qubits);raw=Counter(g[0] for g in gates);expanded=expansion(raw)
    assert raw['CCX']==2*(original_qubits-1)
    return {'gate_counts':dict(raw),'expanded_gate_counts':dict(expanded),'T':7*raw['CCX'],
        'logical_qubits':2*original_qubits,'private_scratch':original_qubits-1,
        'exact_7T_CCX_count':raw['CCX'],'explicit_serial_T_depth_upper':4*raw['CCX'],
        'conservative_serial_T_depth_upper':5*raw['CCX'],'control_not_negated':True}

def emitted_profiles(flag):
    base=vo.resource_api(flag);active=flag['input_qubits']+flag['comparison_qubits']+1
    old=vo._reflection_profile(active);new=reflection_metrics(active)
    q=base['Q'].copy()
    counts=Counter({k.removeprefix('gate_'):v for k,v in q.items() if k.startswith('gate_')})
    old_counts=Counter(old['clifford_T_expanded_gate_counts']);delta=Counter(new['expanded_gate_counts'])
    for g,n in old_counts.items():delta[g]-=n
    # Good Z is replaced by controlled CZ; conditional global-minus Z added.
    delta['Z']-=1;delta['CZ']+=1;delta['Z']+=1
    counts.update(delta);assert min(counts.values())>=0
    q['T']=base['Q']['T']+14
    q['explicit_T_depth_serial_upper']=base['Q']['explicit_T_depth_serial_upper']+8
    q['conservative_5layer_T_depth_serial_upper']=base['Q']['conservative_5layer_T_depth_serial_upper']+10
    q['expanded_primitive_operations']=sum(counts.values())
    for g,n in counts.items():q['gate_'+g]=n
    assert q['T']==counts['T']+counts['Tdg']
    return {'A':{'metrics':base['A'],'logical_qubits':base['qubits']},
            'CQ':{'metrics':q,'logical_qubits':base['qubits']+2},
            'active_financial_qubits':active,
            'qubit_scope':'conservative body/reflection bound plus control and additional scratch; phase allocation adds its entire width, possibly overcounting one reusable control'}

def complete_scalar(setup,flag,qft_profile):
    if setup.zero_shots:
        return {'T':0,'runs':0,'controlled_Q_calls':0,'peak_logical_qubits':0,'metrics':{},'phase_register_bits':0}
    profiles=emitted_profiles(flag);out=qpe_estimator.charge(setup,profiles['A'],profiles['CQ'],qft_profile)
    R=setup.repetitions;m=setup.phase_register_bits;n=profiles['active_financial_qubits'];metrics=out['metrics']
    metrics['reinitialisation_additional_Z_measurements_upper']=R*n
    metrics['reinitialisation_conditional_X_upper']=R*(n+m)
    # 'measurements' already contains R*m phase observations from charge().
    metrics['all_measurements_with_reinitialisation_upper']=metrics.get('measurements',0)+R*n
    metrics['complete_primitive_operations_with_readout_and_reset_upper']=metrics.get('expanded_primitive_operations',0)+R*m+R*n+R*(n+m)
    out.update(T=metrics['T'],runs=R,phase_register_bits=m,qft_T_per_run=qft_profile['T_count'],
               active_financial_qubits=n,reset_scope='phase observations already included; financial register measured for reset, including flag; private work clean')
    return out
