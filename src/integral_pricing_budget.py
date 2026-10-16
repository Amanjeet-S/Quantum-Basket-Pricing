"""Certified QPE budget for new financial residuals, retaining complete tapes.

Existing foundational estimator and single-qubit certificates are unchanged.
This local extension permits smaller monetary-scaled probability tolerances
and wider phase registers. Unsupported distant controlled phases are replaced
by identity with a charged operator-norm bound. Financial core costs are
supplied by actual emitted oracle profiles, never a free-query assumption.
"""
from collections import Counter
from dataclasses import asdict
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json
import qft_synthesis as synthesis
import qpe_estimator as established
from certified_scalar import ceil_fraction,sqrt_upper_fraction
import quantum_gates as gates
ROOT=Path(__file__).resolve().parent


def extended_inverse_qft_tape(bits):
    if isinstance(bits,bool) or not isinstance(bits,int) or not 1<=bits<=64:
        raise ValueError('The new certified phase-width domain is 1 through 64')
    tape=[];phase=F(0);norm=F(0);omitted=[];retained=[]
    for j in range(bits//2):
        k=bits-1-j;tape.extend([('CX',j,k),('CX',k,j),('CX',j,k)])
    for j in range(bits-1,-1,-1):
        for k in range(bits-1,j,-1):
            distance=k-j
            if distance<=31:
                item=synthesis.controlled_phase(distance,k,j,-1)
                tape.extend(item['gates']);phase+=item['global_phase_correction_pi'];norm+=item['norm_error_upper']
                retained.append({'distance':distance,'control':k,'target':j,'norm_error_upper':str(item['norm_error_upper'])})
            else:
                bound=established.PI_UPPER/F(1<<distance)
                norm+=bound
                omitted.append({'distance':distance,'control':k,'target':j,'norm_error_upper':str(bound)})
        tape.append(('H',j))
    return {'gates':tuple(tape),'global_phase_correction_pi':phase,'norm_error_upper':norm,'retained':retained,'omitted':omitted}


@lru_cache(None)
def extended_inverse_qft_profile(bits):
    info=extended_inverse_qft_tape(bits);bound=info['norm_error_upper']
    if bound>established.QFT_NORM_ERROR:raise ValueError('Compiled width exceeds the certified Fourier norm budget')
    counts=Counter(g[0] for g in info['gates']);T=counts['T']+counts['Tdg']
    metrics={'T':T,'measurements':0,'feedback_CZ':0,'reset_X':0,
      'explicit_T_depth_serial_upper':T,'conservative_5layer_T_depth_serial_upper':T,
      'expanded_primitive_operations':sum(counts.values()),'uniform_H':0}
    metrics.update({f'gate_{op}':counts[op] for op in ('H','S','Sdg','T','Tdg','CX')})
    return {'phase_register_bits':bits,'norm_error_upper':bound,'global_phase_correction_pi':str(info['global_phase_correction_pi']),
      'gate_counts':dict(counts),'T_count':T,'metrics':metrics,'extra_ancillas':0,
      'controlled_phases_retained':len(info['retained']),'controlled_phases_omitted':len(info['omitted']),
      'omission_norm_upper':sum((F(x['norm_error_upper']) for x in info['omitted']),F(0)),
      'chronological_tape_sha256':hashlib.sha256(json.dumps(info['gates'],separators=(',',':')).encode()).hexdigest(),
      'proof':'Retained certified rotation sequences and charged distant-phase identity replacement; final QFT scalar phase is physically irrelevant to measurement.'}


def extended_plan(epsilon,alpha=F(1,100),prior_upper=F(1),norm_error=established.QFT_NORM_ERROR):
    epsilon,alpha,prior,norm=map(F,(epsilon,alpha,prior_upper,norm_error))
    if not established.POSTPROCESS_ALLOWANCE<epsilon<F(1,2):raise ValueError('Probability tolerance must exceed the retained postprocessing allowance and be below one half')
    if not 0<alpha<1 or not 0<=prior<=1 or not 0<=norm<=established.QFT_NORM_ERROR:raise ValueError('Invalid confidence, prior or norm budget')
    variance=prior*(1-prior) if prior<=F(1,2) else F(1,4)
    if prior<=2*epsilon:
        return established.Plan(epsilon,alpha,prior,0,1,0,variance,prior/2,F(0),F(0),F(1),F(0),0,'exact_binomial')
    root=sqrt_upper_fraction(variance);bits,M=1,2
    while 2*established.PI_UPPER*root/M+established.PI_UPPER**2/M**2+established.POSTPROCESS_ALLOWANCE>epsilon:
        bits+=1;M*=2
        if bits>64:raise ValueError('Required phase width exceeds the emitted 64-bit domain')
    error=2*established.PI_UPPER*root/M+established.PI_UPPER**2/M**2
    repeats,success,failure,hoeffding=established.confidence_repetitions(alpha,norm,'exact_binomial')
    return established.Plan(epsilon,alpha,prior,bits,M,repeats,variance,error,established.POSTPROCESS_ALLOWANCE,norm,success,failure,hoeffding,'exact_binomial')


def compiled_extended_plan(epsilon,alpha=F(1,100),prior_upper=F(1)):
    preliminary=extended_plan(epsilon,alpha,prior_upper)
    if preliminary.zero_shots:return preliminary,None
    profile=extended_inverse_qft_profile(preliminary.phase_register_bits)
    grid=1<<established.CONFIDENCE_NORM_GRID_BITS
    rounded=F(ceil_fraction(F(profile['norm_error_upper'])*grid),grid)
    final=extended_plan(epsilon,alpha,prior_upper,rounded)
    assert final.M==preliminary.M
    return final,profile


def reflection_metrics(core_qubits):
    """Actually emit controlled zero reflection and conditional global minus."""
    if isinstance(core_qubits,bool) or not isinstance(core_qubits,int) or core_qubits<2:raise ValueError('At least two core wires required')
    core=tuple(range(core_qubits));control=core_qubits
    work=tuple(range(core_qubits+1,2*core_qubits))
    ideal=[('CZ',control,core[-1])]+gates.zero_phase(core,work,control)+[('Z',control)]
    expanded=gates.expand_toffoli(ideal);counts=Counter(g[0] for g in expanded)
    T=counts['T']+counts['Tdg'];ccx=sum(g[0]=='CCX' for g in ideal)
    return {'gate_counts':dict(counts),'T':T,'explicit_T_depth_serial_upper':4*ccx,'conservative_5layer_T_depth_serial_upper':5*ccx,
      'expanded_primitive_operations':len(expanded),'logical_qubits':2*core_qubits,
      'global_minus_Z_charged':True,'controlled_good_CZ_charged':True,
      'chronological_tape_sha256':hashlib.sha256(json.dumps(expanded,separators=(',',':')).encode()).hexdigest()}


def validate_oracle_profile(profile):
    allowed={'X','H','CX','CZ','Z','S','Sdg','T','Tdg'}
    counts=profile['primitive_counts']
    if not isinstance(counts,dict) or set(counts)-allowed:
        raise ValueError('Full oracle requires expanded coherent Clifford and T atoms, without measurements or CCX placeholders')
    if any(isinstance(n,bool) or not isinstance(n,int) or n<0 for n in counts.values()):
        raise ValueError('Primitive operation counts must be nonnegative exact integers')
    fields=['reflection_core_qubits','uniform_H','logical_qubits_upper','explicit_T_depth_serial_upper','conservative_5layer_T_depth_serial_upper']
    if any(isinstance(profile[k],bool) or not isinstance(profile[k],int) or profile[k]<0 for k in fields):
        raise ValueError('Quantum widths and depth caps must be nonnegative exact integers')
    if profile['reflection_core_qubits']<2 or profile['logical_qubits_upper']<profile['reflection_core_qubits']:
        raise ValueError('Oracle allocation must contain its invariant reflection core')
    if profile['uniform_H']>counts.get('H',0):raise ValueError('Uniform Hadamards exceed emitted Hadamards')
    if profile['conservative_5layer_T_depth_serial_upper']<profile['explicit_T_depth_serial_upper']:
        raise ValueError('Conservative depth cap is smaller than its explicit serial cap')
    return True


def charge_full_oracle(setup,profile,qft):
    """Adapter contract: profile is a complete clean comparator preparation A.

    primitive_counts and both serial depth bounds include its financial body,
    flag, inverse, auxiliary law and every input/comparator Hadamard. Scratch
    starts and ends clean. Only the invariant original inputs, comparator and
    flag belong to S0; excluded private work is clean in that invariant space.
    """
    if setup.zero_shots:return {'repetitions':0,'T':0,'scope':'probability-prior zero-shot decision'}
    validate_oracle_profile(profile)
    ap=Counter(profile['primitive_counts']);core=int(profile['reflection_core_qubits']);reflection=reflection_metrics(core)
    at=ap['T']+ap['Tdg'];am={'T':at,'measurements':0,'feedback_CZ':0,'reset_X':0,
      'expanded_primitive_operations':sum(ap.values()),'uniform_H':int(profile['uniform_H']),
      'explicit_T_depth_serial_upper':int(profile['explicit_T_depth_serial_upper']),
      'conservative_5layer_T_depth_serial_upper':int(profile['conservative_5layer_T_depth_serial_upper'])}
    am.update({f'gate_{op}':ap[op] for op in ap})
    inverse_names={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    adjoint=Counter()
    for op,count in ap.items():adjoint[inverse_names.get(op,op)]+=count
    qp=ap+adjoint;qp.update(reflection['gate_counts'])
    qm={'T':2*at+reflection['T'],'measurements':0,'feedback_CZ':0,'reset_X':0,
      'expanded_primitive_operations':sum(qp.values()),'uniform_H':2*int(profile['uniform_H']),
      'explicit_T_depth_serial_upper':2*am['explicit_T_depth_serial_upper']+reflection['explicit_T_depth_serial_upper'],
      'conservative_5layer_T_depth_serial_upper':2*am['conservative_5layer_T_depth_serial_upper']+reflection['conservative_5layer_T_depth_serial_upper']}
    qm.update({f'gate_{op}':qp[op] for op in qp})
    peak=max(int(profile['logical_qubits_upper']),reflection['logical_qubits'])
    result=established.charge(setup,{'metrics':am,'logical_qubits':peak},{'metrics':qm,'logical_qubits':peak},qft)
    # Each reused ideal shot observes/reinitialises every allocated quantum wire.
    # Terminal phase measurements are present in its measurements metric,
    # but absent from its expanded-operation metric; add every observation there.
    wires=result['peak_logical_qubits']
    for key,per_shot in [('measurements',wires-setup.phase_register_bits),('reset_X',wires),('expanded_primitive_operations',2*wires)]:
        result['per_repetition_metrics'][key]+=per_shot
        result['metrics'][key]+=setup.repetitions*per_shot
    assert all(result['metrics'][key]==setup.repetitions*value for key,value in result['per_repetition_metrics'].items())
    result['reset_scope']='Each allocated wire is measured then conditionally X reset; phase observations are not counted twice in the measurement metric and are included once in expanded operations.'
    result['complete_oracle_T']=at;result['controlled_Q_T']=qm['T'];result['controlled_reflection']=reflection
    return result


def main():
    profiles=[extended_inverse_qft_profile(bits) for bits in range(1,65)]
    for bits in range(1,33):
        old=synthesis.inverse_qft_profile(bits);new=profiles[bits-1]
        assert old['T_count']==new['T_count'] and old['gate_counts']==new['gate_counts'] and old['norm_error_upper']==new['norm_error_upper']
    plans=[]
    for eps in [F(1,10**k) for k in (3,6,9,12,15,18)]:
        plan,profile=compiled_extended_plan(eps);plans.append({'epsilon':str(eps),'plan':asdict(plan),'compiled_fourier_norm':str(profile['norm_error_upper'])})
    output={'schema':'qbp-integral-extended-estimator-v1','status':'passed','profiles':profiles,'plans':plans,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'dependency_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['qpe_estimator.py','qft_synthesis.py','quantum_gates.py','certified_scalar.py']},
      'input_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['results/qft_rotation_tapes.json','results/qft_rotation_certificates.json']},
      'scope':'Actual Fourier/compiler and finite-confidence components only; financial oracle and transfer must be supplied independently.'}
    (ROOT/'results/integral_extended_estimator.json').write_text(json.dumps(output,default=str,indent=2)+'\n')
    print(json.dumps({'status':'passed','compiled_widths':64,'old_widths_preserved':32,'max_norm':float(max(F(p['norm_error_upper']) for p in profiles)),'largest_T':profiles[-1]['T_count']}))
if __name__=='__main__':main()
