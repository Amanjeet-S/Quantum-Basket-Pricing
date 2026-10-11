"""Independent extended inverse-QFT, confidence and complete-shot charge audit."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def independent_fourier(bits,rotations,bounds):
    tape=[];phase=F(0);error=F(0);omission=F(0);kept=lost=0
    for left in range(bits//2):
        right=bits-1-left;tape.extend([('CX',left,right),('CX',right,left),('CX',left,right)])
    for target in range(bits-1,-1,-1):
        for control in range(bits-1,target,-1):
            d=control-target
            if d>=32:
                error+=F(355,113*(1<<d));omission+=F(355,113*(1<<d));lost+=1;continue
            k=d+1;minus=rotations[f'-1:{k}'];plus=rotations[f'1:{k}']
            r=tuple(minus['chronological_gates']);s=tuple(plus['chronological_gates'])
            local=[(g,control) for g in r]+[(g,target) for g in r]+[('CX',control,target)]+[(g,target) for g in s]+[('CX',control,target)]
            if not r and not s:local=[]
            tape.extend(local)
            phase+=-F(1,1<<(d+2))+2*F(minus['phase_pi'])+F(plus['phase_pi'])
            error+=2*bounds[f'-1:{k}']+bounds[f'1:{k}'];kept+=1
        tape.append(('H',target))
    return tape,phase,error,omission,kept,lost


def exact_binomial(repetitions,success):
    return sum((F(math.comb(repetitions,k))*(1-success)**k*success**(repetitions-k)
                for k in range((repetitions+1)//2,repetitions+1)),F(0))


def numerical_fourier_diagnostic(extension):
    try:
        import cmath
        import numpy as np
        from integral_phase_support import apply_operations
    except ImportError:
        return {'status':'unavailable without NumPy','scope':'The exact tape/count/norm audit remains available'}
    records=[]
    for bits in range(1,6):
        info=extension.extended_inverse_qft_tape(bits)
        operations=[(('TD' if g[0]=='Tdg' else 'SD' if g[0]=='Sdg' else g[0]),g[1:]) for g in info['gates']]
        actual=apply_operations(operations,bits,list(range(1<<bits)))
        actual*=cmath.exp(1j*math.pi*float(info['global_phase_correction_pi']))
        permutation=np.array([int(f'{j:0{bits}b}'[::-1],2) for j in range(1<<bits)])
        reference=np.exp(-2j*math.pi*permutation[:,None]*permutation[None,:]/(1<<bits))/math.sqrt(1<<bits)
        error=float(np.linalg.norm(actual-reference,ord=2))
        bound=float(info['norm_error_upper'])
        assert error<=bound+2e-12
        records.append({'phase_bits':bits,'operator_norm_diagnostic':error,'certified_norm_upper':bound})
    return {'status':'all five direct Fourier matrix diagnostics passed','records':records,
        'scope':'Floating-point diagnostics of retained scalar phases and endian order, not the certified norm proof'}


def evaluate(root):
    root=Path(root);sys.path.insert(0,str(root))
    import integral_pricing_budget as extension
    import qpe_estimator as old
    tapes_path=root/'results/qft_rotation_tapes.json';certificate_path=root/'results/qft_rotation_certificates.json'
    raw=json.loads(tapes_path.read_text());certificate=json.loads(certificate_path.read_text())
    assert certificate['tapes_sha256']==digest(tapes_path)
    rotations=raw['rotations'];bounds={r['key']:F(r['operator_norm_upper_rational']) for r in certificate['certificates']}
    assert set(bounds)==set(rotations)
    assert all(r['passed'] and r['precision_bits']==384 for r in certificate['certificates'])
    records=[]
    for bits in range(1,65):
        tape,phase,norm,omission,kept,lost=independent_fourier(bits,rotations,bounds)
        actual=extension.extended_inverse_qft_tape(bits);profile=extension.extended_inverse_qft_profile(bits)
        assert tuple(tape)==actual['gates'];assert phase==actual['global_phase_correction_pi']
        assert norm==actual['norm_error_upper']==F(profile['norm_error_upper'])<F(1,64)
        assert omission==F(profile['omission_norm_upper'])
        assert kept==profile['controlled_phases_retained'] and lost==profile['controlled_phases_omitted']
        assert kept+lost==bits*(bits-1)//2
        assert lost==(max(0,bits-32)*max(0,bits-31))//2
        closed_omission=sum((F((bits-d)*355,113*(1<<d)) for d in range(32,bits)),F(0))
        assert omission==closed_omission
        counts=Counter(g[0] for g in tape);assert dict(counts)==profile['gate_counts']
        assert profile['T_count']==counts['T']+counts['Tdg']
        assert all(g[0] in ('H','CX','S','Sdg','T','Tdg') for g in tape)
        assert all(0<=wire<bits for gate in tape for wire in gate[1:])
        assert profile['chronological_tape_sha256']==hashlib.sha256(json.dumps(tuple(tape),separators=(',',':')).encode()).hexdigest()
        if bits<=32:
            # Compare the original emitted template at each shared width.
            import qft_synthesis
            before=qft_synthesis.inverse_qft_tape(bits)
            assert before['gates']==actual['gates'] and before['norm_error_upper']==norm
        records.append({'phase_bits':bits,'T_count':profile['T_count'],'norm_upper_fraction':str(norm),
            'omission_norm_upper_fraction':str(omission),'retained_controlled_phases':kept,
            'omitted_controlled_phases':lost,'tape_sha256':profile['chronological_tape_sha256']})
    rejection_checks=0
    for bad in (0,65,-1,True,1.5,F(3,2),'32'):
        try:extension.extended_inverse_qft_tape(bad)
        except ValueError:rejection_checks+=1
        else:raise AssertionError(('Invalid phase width accepted',bad))
    plans=[]
    for prior in (F(1),F(1,10),F(1,10**4)):
        for power in (3,6,9,12,15,18):
            epsilon=F(1,10**power);setup,profile=extension.compiled_extended_plan(epsilon,prior_upper=prior)
            if setup.zero_shots:
                assert prior/2<=epsilon;assert profile is None;continue
            M=setup.M;variance=prior*(1-prior) if prior<=F(1,2) else F(1,4)
            assert setup.variance_upper==variance
            assert setup.amplitude_error_upper+setup.postprocess_allowance<=epsilon
            if setup.phase_register_bits>1:
                # M/2 cannot meet the same stipulated variance-envelope criterion.
                assert 2*old.PI_UPPER*old.sqrt_upper_fraction(variance)/(M//2)+old.PI_UPPER**2/(M//2)**2+old.POSTPROCESS_ALLOWANCE>epsilon
            norm=F(profile['norm_error_upper']);rounded=setup.qft_norm_error_allowance
            assert norm<=rounded<norm+F(1,1<<80)
            success=8/F(355,113)**2-rounded
            assert setup.good_probability_lower==success>F(1,2)
            fail=exact_binomial(setup.repetitions,success)
            assert fail==setup.binomial_failure_upper<=setup.alpha
            if setup.repetitions>1:assert exact_binomial(setup.repetitions-2,success)>setup.alpha
            assert setup.total_controlled_grover_calls==setup.repetitions*(M-1)
            plans.append({'prior_upper_fraction':str(prior),'epsilon_fraction':str(epsilon),
                'phase_bits':setup.phase_register_bits,'repetitions':setup.repetitions,
                'success_lower_fraction':str(success),'binomial_failure_upper_fraction':str(fail)})
    family_plans=[]
    for alpha in (F(1,2700),F(1,5400)):
        setup,profile=extension.compiled_extended_plan(F(1,10**9),alpha,F(1,10))
        failure=exact_binomial(setup.repetitions,setup.good_probability_lower)
        assert failure==setup.binomial_failure_upper<=alpha
        assert setup.repetitions%2==1
        family_plans.append({'alpha_fraction':str(alpha),'repetitions':setup.repetitions,
            'binomial_failure_upper_fraction':str(failure)})
    assert 27*F(1,2700)==F(1,100)
    assert 27*2*F(1,5400)==F(1,100)
    remaining=F(1,200)
    for normaliser in (F(1,100),F(1),F(10000)):
        each=remaining/(2*normaliser)
        assert normaliser*(each+each)==remaining
    for bad_epsilon in (F(0),-F(1,10),F(1,10**61),F(1,2),F(1,10**30)):
        try:extension.extended_plan(bad_epsilon)
        except ValueError:rejection_checks+=1
        else:raise AssertionError(('Unsupported probability tolerance accepted',bad_epsilon))
    for kwargs in ({'alpha':0},{'alpha':1},{'prior_upper':-1},{'prior_upper':F(101,100)},
                   {'norm_error':-1},{'norm_error':F(1,63)}):
        try:extension.extended_plan(F(1,1000),**kwargs)
        except ValueError:rejection_checks+=1
        else:raise AssertionError(('Invalid plan parameter accepted',kwargs))
    reflection_records=[]
    for core in range(2,17):
        actual=extension.reflection_metrics(core)
        conjunctions=2*(core-1)
        expected=Counter({'X':2*core,'CZ':2,'Z':1,
            'H':2*conjunctions,'CX':6*conjunctions,'T':4*conjunctions,'Tdg':3*conjunctions})
        assert Counter(actual['gate_counts'])==expected
        assert actual['T']==7*conjunctions
        assert actual['explicit_T_depth_serial_upper']==4*conjunctions
        assert actual['conservative_5layer_T_depth_serial_upper']==5*conjunctions
        assert actual['expanded_primitive_operations']==sum(expected.values())
        assert actual['logical_qubits']==2*core
        assert actual['global_minus_Z_charged'] and actual['controlled_good_CZ_charged']
        reflection_records.append({'reflection_core_qubits':core,'counts':dict(expected),'T':actual['T']})
    # Unbalanced T/TD and S/SD counts expose incorrect use of twice A's
    # primitive multiset. The literal adjoint must swap gate types.
    example={'primitive_counts':{'H':11,'CX':13,'T':7,'Tdg':3,'S':2,'Sdg':1,'X':5},
        'reflection_core_qubits':5,'logical_qubits_upper':19,'uniform_H':6,
        'explicit_T_depth_serial_upper':10,'conservative_5layer_T_depth_serial_upper':10}
    profile_guard_checks=0
    invalid_profiles=[]
    for illegal in ('CCX','MZ','TD','arbitrary_rotation'):
        invalid_profiles.append({**example,'primitive_counts':{**example['primitive_counts'],illegal:1}})
    for invalid in (True,F(1,2),-1):
        invalid_profiles.append({**example,'primitive_counts':{**example['primitive_counts'],'T':invalid}})
    for key,invalid in (('reflection_core_qubits',True),('reflection_core_qubits',1),
                        ('logical_qubits_upper',4),('uniform_H',12),
                        ('explicit_T_depth_serial_upper',-1),('conservative_5layer_T_depth_serial_upper',9)):
        invalid_profiles.append({**example,key:invalid})
    for profile in invalid_profiles:
        try:extension.validate_oracle_profile(profile)
        except ValueError:profile_guard_checks+=1
        else:raise AssertionError(('Invalid financial-oracle profile accepted',profile))
    assert extension.validate_oracle_profile(example)
    setup,qft=extension.compiled_extended_plan(F(1,1000))
    charge=extension.charge_full_oracle(setup,example,qft)
    adjoint_swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    ap=Counter(example['primitive_counts']);inverse=Counter()
    for gate,count in ap.items():inverse[adjoint_swap.get(gate,gate)]+=count
    reflection=Counter(extension.reflection_metrics(5)['gate_counts']);cq=ap+inverse+reflection
    fourier=Counter(qft['gate_counts']);shots=setup.repetitions;calls=setup.M-1
    expected=ap+Counter({gate:calls*count for gate,count in cq.items()})+fourier;expected['H']+=setup.phase_register_bits
    for gate,count in expected.items():assert charge['metrics'][f'gate_{gate}']==shots*count,('Complete charge type',gate)
    wires=max(19,10)+setup.phase_register_bits
    assert charge['peak_logical_qubits']==wires
    assert charge['metrics']['T']==shots*(10+calls*(20+reflection['T']+reflection['Tdg'])+qft['T_count'])
    assert charge['metrics']['measurements']==shots*wires
    assert charge['metrics']['reset_X']==shots*wires
    assert charge['metrics']['expanded_primitive_operations']==shots*(sum(expected.values())+2*wires)
    for key,value in charge['per_repetition_metrics'].items():
        assert charge['metrics'][key]==shots*value,('Per-shot/total consistency',key)
    return {'status':'independent extended estimator audit passed','fourier_widths':records,'active_plans':plans,
        'small_fourier_matrix_diagnostic':numerical_fourier_diagnostic(extension),
        'family_confidence_and_split_normalisation_checks':family_plans,
        'reflection_profiles':reflection_records,'parameter_rejection_checks':rejection_checks,
        'oracle_profile_rejection_checks':profile_guard_checks,
        'complete_charge_test':{'shots':shots,'controlled_Q_calls_per_shot':calls,'logical_qubits':wires,
            'unbalanced_adjoint_type_counts_checked':True,'all_wire_observations_and_conditional_resets_charged':True},
        'source_sha256':digest(root/'integral_pricing_budget.py'),
        'input_sha256':{'results/qft_rotation_tapes.json':digest(tapes_path),'results/qft_rotation_certificates.json':digest(certificate_path)},
        'scope':'Emitted resource, norm and finite-confidence accounting; no physical QPE execution or practical quantum advantage claim'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default='/private/tmp/qbp-integral-value-step')
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_estimator_extension_audit.json'))
    args=parser.parse_args();result=evaluate(args.root);result['audit_source_sha256']=digest(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'widths':len(result['fourier_widths']),
        'active_plans':len(result['active_plans']),'parameter_rejections':result['parameter_rejection_checks']},indent=2))

if __name__=='__main__':main()
