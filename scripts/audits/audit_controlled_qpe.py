"""Independent controlled-reflection counts and corrected-instrument QPE audit.

Production comparisons replay retained additive gate multisets, not enormous
state vectors. The complete small lookup/QPE tapes are separate numerical
diagnostics on explicitly retained measurement records.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import random
import sys
import time


def audit(root):
    sys.path.insert(0,str(root))
    import architecture as arch
    import controlled_estimation_resources as cr
    import value_oracle as vo
    import quantum_gates as qg
    import qft_synthesis as qs
    import qpe_estimator as qe
    import numpy as np
    started=time.perf_counter()
    def counts(gates):
        return Counter(g[0] for g in gates)
    def expanded(gates):
        return counts(qg.expand_toffoli(gates))
    reflection_basis=0
    for n in range(2,10):
        own=qg.zero_phase(tuple(range(n)),tuple(range(n+1,2*n)),control=n)
        actual=cr.controlled_zero_reflection(n)
        assert own==list(actual)
        assert not any(g[0]=='X' and g[1]==n for g in actual)
        for value in range(1<<(n+1)):
            state=[(value>>j)&1 for j in range(n+1)]+[0]*(n-1)
            phase=1
            for g in actual:
                if g[0]=='X':state[g[1]]^=1
                elif g[0]=='CCX':state[g[3]]^=state[g[1]]&state[g[2]]
                elif g[0]=='CZ':phase*=(-1 if state[g[1]]&state[g[2]] else 1)
                else:raise ValueError(g)
            assert sum(state[j]<<j for j in range(n+1))==value
            assert not any(state[n+1:])
            expected=-1 if value==(1<<n) else 1
            assert phase==expected
            reflection_basis+=1
    retained=json.loads((root/'results'/'architecture_resource_candidates.json').read_text())['configurations']
    certificate_rows=json.loads((root/'results'/'finite_range_certificates.json').read_text())['configurations']
    def certificate_key(row):
        return tuple(row[k]for k in ('case_id','method','q','exp_address_bits',
                                    'log_address_bits','probability_address_bits'))
    certificate_map={certificate_key(r):r for r in certificate_rows}
    preliminary=qe.plan(F(1,100))
    ft_profile=qs.inverse_qft_profile(preliminary.phase_register_bits)
    setup=qe.plan(F(1,100),norm_error=ft_profile['norm_error_upper'])
    complete_scalar_checks=0
    production_checks=[]
    for row in retained:
        q=row['q'];body=vo.prepare_existing_body(row['resources'],q)
        maximum=certificate_map[certificate_key(row)]['max_integer']
        certified_bits=max(0,(maximum-1).bit_length())
        assert maximum<=1<<certified_bits
        flags=[('unit',vo.profile_from_prepared_range(body,q)),
               ('range',vo.profile_from_prepared_range(body,certified_bits)),
               ('clipped',vo.profile_from_prepared_clipped_range(body,max(0,q-2))),
               ('bucket',vo.profile_from_prepared_body(body,q//2))]
        for kind,flag in flags:
            profiles=cr.emitted_profiles(flag)
            n=flag['input_qubits']+flag['comparison_qubits']+1
            own_ref=qg.zero_phase(tuple(range(n)),tuple(range(n+1,2*n)),control=n)
            ac=Counter(flag['preparation_expanded_gate_counts'])
            expected=Counter({g:2*v for g,v in ac.items()})
            expected.update(expanded(own_ref));expected.update({'CZ':1,'Z':1})
            actual=Counter({k[5:]:v for k,v in profiles['CQ']['metrics'].items() if k.startswith('gate_')})
            assert expected==actual
            assert profiles['CQ']['metrics']['T']==expected['T']+expected['Tdg']
            assert profiles['CQ']['metrics']['T']==flag['grover_iterate_T_count']+14
            assert profiles['CQ']['metrics']['measurements']==expected['MZ']
            assert profiles['CQ']['metrics']['feedback_CZ']==expected['CZ_IF']
            assert profiles['CQ']['metrics']['reset_X']==expected['X_IF']
            assert profiles['CQ']['metrics']['expanded_primitive_operations']==sum(expected.values())
            complete=cr.complete_scalar(setup,flag,ft_profile)
            R=setup.repetitions;M=setup.M;m=setup.phase_register_bits
            A_T=ac['T']+ac['Tdg'];CQ_T=expected['T']+expected['Tdg']
            assert complete['T']==R*(A_T+(M-1)*CQ_T+ft_profile['T_count'])
            internal=R*(ac['MZ']+(M-1)*expected['MZ'])
            assert complete['metrics']['measurements']==internal+R*m
            assert complete['metrics']['all_measurements_with_reinitialisation_upper']==internal+R*(m+n)
            assert complete['metrics']['reinitialisation_conditional_X_upper']==R*(n+m)
            primitive=R*(sum(ac.values())+(M-1)*sum(expected.values())
                         +ft_profile['metrics']['expanded_primitive_operations']+m)
            assert complete['metrics']['expanded_primitive_operations']==primitive
            assert complete['metrics']['complete_primitive_operations_with_readout_and_reset_upper']==primitive+R*(m+n)+R*(n+m)
            complete_scalar_checks+=1
            production_checks.append({'candidate_id':row['candidate_id'],
                'architecture':row['architecture'],'flag_kind':kind,
                'comparison_bits':flag['comparison_qubits'],
                'active_financial_qubits':n,'controlled_T':profiles['CQ']['metrics']['T']})
    # Small literal lookup instruments, with their measurement records retained.
    def primitive_replay(operations,state,outcomes):
        state=state.copy();indices=np.arange(state.size)
        for op,wires,condition in operations:
            if op in ('CZ_IF','X_IF') and not outcomes[condition]:continue
            if op in ('X','X_IF','CX','CCX'):
                target=wires[-1];mask=np.ones(indices.size,bool)
                for w in wires[:-1]:mask &= ((indices>>w)&1).astype(bool)
                state=state[indices^(mask.astype(np.int64)<<target)]
            elif op=='H':
                low=indices[(indices&(1<<wires[0]))==0];high=low|(1<<wires[0])
                a,b=state[low].copy(),state[high].copy()
                state[low]=(a+b)/math.sqrt(2);state[high]=(a-b)/math.sqrt(2)
            elif op in ('T','TD','Tdg','S','SD','Sdg','Z','CZ','CZ_IF'):
                active=np.ones(indices.size,bool)
                for w in wires:active &= ((indices>>w)&1).astype(bool)
                phase={'T':np.exp(1j*math.pi/4),'TD':np.exp(-1j*math.pi/4),
                       'Tdg':np.exp(-1j*math.pi/4),'S':1j,'SD':-1j,
                       'Sdg':-1j,'Z':-1,'CZ':-1,'CZ_IF':-1}[op]
                state[active]*=phase
            elif op=='MZ':
                state[((indices>>wires[0])&1)!=outcomes[condition]]=0
                state*=math.sqrt(2)
            else:raise ValueError(op)
        return state
    instrument_records=[]
    for mode in arch.MODES:
        for m in (2,3):
            phase=tuple(range(6,6+m));n_total=6+m
            lookup=arch.lookup_plan((0,1),(2,),(1,0,0,0),mode,(3,4,5))
            A=[('H',0),('H',1),('LOOKUP',lookup)]
            gates=list(A)+[('H',w)for w in phase]
            for j,w in enumerate(phase):
                gates+=qg.grover(A,(0,1,2),2,(3,4),control=w)*(1<<j)
            ft=qs.inverse_qft_tape(m)
            gates+=[(g[0],*(phase[m-1-w] for w in g[1:]))for g in ft['gates']]
            number_records=sum(g[1].summary()['measurement_count'] for g in gates if g[0]=='LOOKUP')
            variants=([tuple([0]*number_records),tuple([1]*number_records),
                      tuple(j%2 for j in range(number_records))]
                      if number_records else [()])
            if number_records:
                rng=random.Random(58371+m)
                variants += [tuple(rng.randrange(2) for _ in range(number_records))for _ in range(3)]
            for branch in variants:
                bits=iter(branch);state=np.zeros(1<<n_total,complex);state[0]=1
                for g in gates:
                    if g[0]=='LOOKUP':
                        plan=g[1];outcomes={j:next(bits)for j in range(plan.summary()['measurement_count'])}
                        state=primitive_replay(plan.operations(),state,outcomes)
                    else:
                        state=primitive_replay([(g[0],tuple(g[1:]),None)],state,{})
                observed=[0.0]*(1<<m)
                for x,z in enumerate(state):
                    y=sum(((x>>w)&1)<<j for j,w in enumerate(phase))
                    observed[y]+=abs(z)**2
                target=qe.ideal_amplitude_distribution(.25,1<<m)
                tv=sum(abs(a-b)for a,b in zip(observed,target))/2
                assert tv<=float(ft['norm_error_upper'])+2e-11
                assert abs(np.vdot(state,state).real-1)<2e-11
                private=sum(abs(z)**2 for x,z in enumerate(state)if x&((1<<3)|(1<<4)|(1<<5)))
                assert private<2e-22
                instrument_records.append({'architecture':mode,'phase_bits':m,
                    'record_bits':number_records,'record':''.join(map(str,branch)),
                    'phase_law_total_variation':tv,'private_work_probability':private,
                    'scope':'literal expanded lookup primitives and specified normalised Kraus branches; all records exhaustive only in unitary modes'})
    return {'schema':'qbp-independent-controlled-qpe-audit-v1','status':'passed',
        'reflection_complete_basis_checks':reflection_basis,
        'production_controlled_multiset_checks':len(production_checks),
        'complete_scalar_resource_and_reset_checks':complete_scalar_checks,
        'retained_financial_candidate_rows':len(retained),
        'production_checks':production_checks,'instrument_records':instrument_records,
        'instrument_record_diagnostic_count':len(instrument_records),
        'measurement_scope':'six prescribed full record strings per measured case, not all exponentially many records; algebraic scalar-Kraus proof establishes the full channel',
        'financial_scope':'production emitted count multisets cover all 819 inherited profiles and four representative flag variants; small literal channel replay uses an explicitly synthetic one-bit lookup, not a production financial-state vector',
        'elapsed_seconds':time.perf_counter()-started,
        'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()for name in
            ('controlled_estimation_resources.py','value_oracle.py','quantum_gates.py','architecture.py',
             'qpe_estimator.py','qft_synthesis.py','results/qft_rotation_tapes.json','results/qft_rotation_certificates.json',
             'results/architecture_resource_candidates.json','results/finite_range_certificates.json')}}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    record=audit(args.root)
    output=args.output or args.root/'results'/'controlled_qpe_audit.json'
    output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k]for k in ('status','reflection_complete_basis_checks',
        'production_controlled_multiset_checks','instrument_record_diagnostic_count','elapsed_seconds')}))
