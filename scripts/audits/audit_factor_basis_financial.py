"""Independent financial/prior composition audit for seven exact factor bases.

The quantum producer and its frontend budget/prior routines are not imported.
All compositions below are rederived with exact Fraction arithmetic from
retained certified finite moments, scalar arithmetic, global norms and raw
financial transfers. Tape/ordinate and primitive certificates are linked to
separate independently retained audits; this record does not replace them.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()


def lo(iv):
    return F(iv['lower_exact_fraction'])


def hi(iv):
    return F(iv['upper_exact_fraction'])


def audit_informed(prior, previous, prep, certificate, financial):
    B=F(certificate['configuration']['monetary_signal_normaliser'])
    eta=F(certificate['probability_error_upper'])
    moment=prep['continuous_second_moment']
    remainder=F(prep['payoff_norms']['M2_upper'])*hi(moment['continuous_U2_interval'])/2
    c1=prep['coefficient_certificate']['c1']
    continuous_lower=lo(c1)-remainder
    continuous_upper=hi(c1)+remainder
    beta=F(financial['financial_transfer_upper'])
    ml,mh=lo(financial['mass_interval']),hi(financial['mass_interval'])
    assert 0<ml<=mh<=1
    finite_lower=max(F(0),continuous_lower-beta)/mh
    finite_upper=min(B,(continuous_upper+beta)/ml,F(previous['ideal_mean_upper']))
    assert 0<=finite_lower<=finite_upper<=B
    pl=max(F(0),finite_lower/B-eta)
    pu=min(F(1),F(previous['probability_upper']),finite_upper/B+eta)
    assert 0<=pl<=pu<=1
    expected={'probability_lower':pl,'probability_upper':pu,'finite_ideal_mean_lower':finite_lower,
        'finite_ideal_mean_upper':finite_upper,'ideal_mean_upper':finite_upper,
        'previous_probability_upper':F(previous['probability_upper']),
        'previous_ideal_mean_upper':F(previous['ideal_mean_upper']),
        'continuous_R2_mean_absolute_upper':remainder,
        'continuous_R1_mean_lower':continuous_lower,'continuous_R1_mean_upper':continuous_upper,
        'financial_transfer_upper':beta,'arithmetic_probability_enlargement':eta}
    for key,value in expected.items():
        assert F(prior[key])==value,key
    assert prior['native_positive'] is True and prior['continuous_R2_nonnegative'] is False
    assert prior['continuous_c1_interval']==c1
    assert prior['external_mass_interval']==financial['mass_interval']
    assert prior['continuous_moment_record_sha256']==digest(moment)
    return len(expected)


def aux_bias(kind,prep,cert):
    norm=prep['payoff_norms'];mu=prep['finite_moments']['centred_U_moments']
    m2,m3,m4,m6=(hi(mu[str(k)]) for k in (2,3,4,6))
    M2=F(norm['M2_upper']);M3=F(norm['C_r_upper']['3']);M4=F(norm['C_r_upper']['4']);M6=F(norm['C_r_upper']['6'])
    if kind=='boundary_r2':
        h=F(1,256)
        return h*h*(M3*m3/12+M4*m4/48)
    if kind=='stable_gauss_r2':
        constants=cert['configuration']['gauss_constants']
        node=sum((F(constants['w_'+suffix]['ideal_upper'])*F(constants['z_'+suffix]['error_upper']) for suffix in ('minus','plus')),F(0))
        weight=sum((F(constants['w_'+suffix]['error_upper']) for suffix in ('minus','plus')),F(0))
        return M6*m6/14400+M3*m3*node/2+M2*m2*weight/2
    assert kind in ('raw_direct_r2','positive_r1')
    return F(0)


def audit_budget(budget, financial, prep, certificate, kind):
    mass=financial['mass_interval'];ml,mh=lo(mass),hi(mass);mr=(mh-ml)/2
    assert budget['mass_interval']==mass
    coefficient=prep['coefficient_certificate']['c0' if kind=='positive_r1' else 'coefficient']
    cr=(hi(coefficient)-lo(coefficient))/2
    cfg=certificate['configuration'];B=F(cfg['monetary_signal_normaliser'])
    decode=F(cfg['monetary_probability_decode_factor'])
    arithmetic=decode*B*F(certificate['probability_error_upper'])
    auxiliary=aux_bias(kind,prep,certificate)
    beta=F(financial['financial_transfer_upper'])
    expected={'coefficient_radius_upper':cr,'external_factor_midpoint':(ml+mh)/2,
        'external_factor_upper':mh,'external_factor_radius':mr,
        'arithmetic_before_mass_upper':arithmetic,'financial_transfer_upper':beta,
        'monetary_probability_decode_factor':decode,
        'total_deterministic_undiscounted_price_error_upper':cr+beta+mh*(auxiliary+arithmetic)+mr*B}
    for key,value in expected.items():
        assert F(budget[key])==value,key
    assert F(budget['auxiliary_components']['bias_upper'])==auxiliary
    assert budget['coefficient_interval']==coefficient
    return len(expected)+1


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--step',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    resources=args.step/'resources';rp=resources/'results'
    certificate_path=args.step/'certificate/results/factor_basis_certificates.json'
    math_audit_path=args.step/'review/factor_basis_certificate_audit.json'
    resource_audit_path=rp/'factor_basis_quantum_independent_audit.json'
    data=json.loads((rp/'factor_basis_quantum_resources.json').read_text())
    certs=json.loads(certificate_path.read_text())
    math_audit=json.loads(math_audit_path.read_text())
    resource_audit=json.loads(resource_audit_path.read_text())
    assert math_audit['status']=='passed' and math_audit['certificate_sha256']==sha(certificate_path)
    assert resource_audit['status']=='passed' and resource_audit['resource_file_sha256']==sha(rp/'factor_basis_quantum_resources.json')
    assert data['certificate_sha256']==sha(certificate_path)
    protocol=json.loads((rp/'factor_basis_quantum_protocol.json').read_text())
    assert data['protocol_sha256']==sha(rp/'factor_basis_quantum_protocol.json')
    assert protocol['classical_price_answer_excluded_from_inputs'] is True
    assert data['alpha_per_price']=='1/3600' and data['selected_four_price_union_failure_upper']=='1/900'
    assert protocol['confidence_scope'].startswith('Alternative resource plans are deterministic hypothetical plans')
    by_cert={(r['original_model_id'],r['candidate_id']):r for r in certs['records']}
    prep_by_id={p['model']['model_id']:p for p in data['preparations']}
    assert len(prep_by_id)==7
    facts=[];priorchecks=budgetchecks=flagchecks=0
    for prep in data['preparations']:
        key=prep['original_model_id'],prep['candidate_id'];analytic=by_cert[key]
        assert prep['model']==analytic['rotated_model']
        direct=prep['R1_direct_global_financial_transfer'];ac=analytic['direct_global_transfer']
        for k in ('coordinates','linear_transport_upper','quadratic_transport_upper','tail_upper','financial_transfer_upper','mass_interval','rho1_mass_upper','rho2_mass_upper','other_coordinate_mass_upper'):
            assert direct[k]==ac[k],(key,k)
        selected=min((prep['R1_older_financial_transfer'],direct),key=lambda x:F(x['financial_transfer_upper']))
        chosen=prep['R1_selected_financial_transfer']
        for k,v in selected.items():
            assert chosen[k]==v,(key,k)
        assert chosen['mass_interval']==prep['R2_financial_transfer']['mass_interval']
        actual_bounds=prep['domain']['actual_coordinate_bounds'];cover=prep['domain']['coordinate_bounds']
        assert F(cover[0])<=F(actual_bounds[0])<F(actual_bounds[1])<=F(cover[1])
        assert F(prep['domain']['U_star'])==F(prep['domain']['T_star'])-1>=0
        assert prep['finite_moments']['model_id']==prep['model']['model_id']
        assert prep['domain']['vertices_evaluated']==2**prep['model']['outer_dimension']
        assert F(prep['domain']['T_star'])==max([F(1)]+[F(v['T_upper']) for v in prep['domain']['vertex_upper_enclosures']])
        identity=prep_by_id[prep['original_model_id']]
        # Exact Gaussian arguments are already proved in the all-58 audit.
        # Independent rounding encloses the same coefficients; do not mistake
        # a differing rational enclosure for a changed true coefficient.
        for coefficient in ('c0','c1','coefficient'):
            a,b=prep['coefficient_certificate'][coefficient],identity['coefficient_certificate'][coefficient]
            assert max(lo(a),lo(b))<=min(hi(a),hi(b)),(key,coefficient)
        a,b=prep['continuous_second_moment']['continuous_U2_interval'],identity['continuous_second_moment']['continuous_U2_interval']
        assert max(lo(a),lo(b))<=min(hi(a),hi(b))
        facts.append({'original_model_id':key[0],'candidate_id':key[1],
            'finite_U2_upper':str(hi(prep['finite_moments']['centred_U_moments']['2'])),
            'U_star':prep['domain']['U_star'],
            'selected_financial_transfer':chosen['financial_transfer_upper'],
            'finite_record_differs_from_identity':prep['finite_moments']['raw_T_moments']!=identity['finite_moments']['raw_T_moments'],
            'continuous_invariant_enclosures_consistent':True})
    for row in data['records']:
        prep=prep_by_id[row['model_id']];certificate=row['base_graph_certificate'];kind=row['kind'];cfg=certificate['configuration']
        U=F(prep['domain']['U_star']);M1=F(prep['payoff_norms']['M1_upper']);M2=F(prep['payoff_norms']['M2_upper'])
        assert F(cfg['U_star'])==U and F(cfg['retained_coordinate_bound'])==max(abs(F(x)) for x in prep['domain']['coordinate_bounds'])
        B1=min(F(1),M1*U);B2=M2*U*U/2
        assert F(cfg['positive_R1_normaliser'])==B1 and F(cfg['B2'])==B2
        B=F(cfg['monetary_signal_normaliser']);assert B==(B1 if kind=='positive_r1' else B2)
        financial=prep['R1_selected_financial_transfer'] if kind=='positive_r1' else prep['R2_financial_transfer']
        budget=row['complete_deterministic_budget']
        if budget is not None:
            budgetchecks+=audit_budget(budget,financial,prep,certificate,kind)
        if 'baseline_failure' in row:
            assert budget is None or F(budget['total_deterministic_undiscounted_price_error_upper'])>=F(1,100)
            continue
        moments=prep['finite_moments']['centred_U_moments'];eta=F(certificate['probability_error_upper']);previous=row['correction_prior']
        positive=cfg['output_encoding']=='positive'
        if kind=='positive_r1':
            mean=min(hi(prep['finite_moments']['correlated_underlying_gap_interval']),F(1),M1*hi(moments['1']))
            pu=min(F(1),mean/B+eta)
            assert F(previous['ideal_mean_upper'])==mean and F(previous['probability_upper'])==pu
            assert F(previous['underlying_gap_upper'])==hi(prep['finite_moments']['correlated_underlying_gap_interval'])
            assert F(previous['rounding_enlargement'])==eta
            priorchecks+=4
            priorchecks+=audit_informed(row['coefficient_informed_prior'],previous,prep,certificate,financial)
            oldfinancial=prep['R1_older_financial_transfer']
            budgetchecks+=audit_budget(row['older_financial_deterministic_budget'],oldfinancial,prep,certificate,kind)
            priorchecks+=audit_informed(row['older_coefficient_informed_prior'],previous,prep,certificate,oldfinancial)
        else:
            pu=min(F(1),hi(moments['2'])/(U*U)+(1 if positive else 2)*eta)
            assert F(previous['probability_upper'])==pu
            assert F(previous['ideal_envelope_probability_upper'])==hi(moments['2'])/(U*U)
            assert F(previous['rounding_enlargement'])==(1 if positive else 2)*eta
            priorchecks+=3
        assert previous['native_positive']==positive
        for strategy in row['strategies']:
            information=strategy['prior_information']
            chosen_budget=row['older_financial_deterministic_budget'] if information.startswith('older_financial') else budget
            if information in ('previous_proved_prior','older_financial_previous_prior'):
                prior=previous
            elif information=='coefficient_informed_direct_minimum':
                prior=row['coefficient_informed_prior']
            else:
                assert information=='older_financial_coefficient_informed_prior'
                prior=row['older_coefficient_informed_prior']
            out=strategy['result']
            if out['kind']=='rejected':
                assert F(chosen_budget['total_deterministic_undiscounted_price_error_upper'])>=F(1,100)
                continue
            C=hi(chosen_budget['mass_interval'])*B
            remaining=F(1,100)-F(chosen_budget['total_deterministic_undiscounted_price_error_upper'])
            if out['kind']=='deterministic_correction_prior':
                width=F(prior['probability_upper'])-F(prior.get('probability_lower','0'))
                radius=C*width/(2 if positive else 1)
                assert radius<=remaining
                assert F(out['total_price_error_upper'])==F(chosen_budget['total_deterministic_undiscounted_price_error_upper'])+radius
                continue
            assert out['kind']=='complete_quantum_resource_plan'
            assert out['correction_prior']==prior
            assert F(out['monetary_signal_factor_upper'])==C
            assert F(out['deterministic_error_upper'])==F(chosen_budget['total_deterministic_undiscounted_price_error_upper'])
            for flag in out['flags']:
                expected_prior=F(prior['probability_upper']) if positive or out['encoding']=='split' else F(1)
                assert F(flag['prior_upper'])==expected_prior
                assert F(flag['probability_tolerance'])==remaining/(C if positive else 2*C)
                assert F(flag['plan']['alpha'])==(F(1,7200) if out['encoding']=='split' else F(1,3600))
                flagchecks+=1
    pairs=[]
    for comparison in data['paired_comparison']:
        original=comparison['original_model_id']
        viable=[r for r in data['records'] if r['original_model_id']==original and r['selected_strategy']]
        best=min(viable,key=lambda r:(r['selected_strategy']['result']['T_count'],r['factor_basis_candidate_id'],r['kind'],r['candidate_index']))
        assert best['selected_strategy']==comparison['selected_strategy']
        assert best['factor_basis_candidate_id']==comparison['selected_candidate_id']
        identity_rows=[r for r in viable if r['factor_basis_candidate_id']=='identity']
        identity=min(identity_rows,key=lambda r:r['selected_strategy']['result']['T_count'])
        iresult=identity['selected_strategy']['result'];result=best['selected_strategy']['result']
        if result['kind']=='complete_quantum_resource_plan':
            phases=[(f['plan']['M'],f['plan']['repetitions']) for f in result['flags']]
            originalphases=[(f['plan']['M'],f['plan']['repetitions']) for f in iresult['flags']]
            assert phases==originalphases
            for a,b in zip(iresult['flags'],result['flags']):
                assert a['fourier_profile']==b['fourier_profile']
        else:
            phases=[]
        pairs.append({'original_model_id':original,'selected_candidate_id':comparison['selected_candidate_id'],
                      'selected_method':comparison['selected_method'],'phase_and_repetition_profile_unchanged':True,
                      'selected_phase_profiles':phases,'identity_T_count':comparison['identity_T_count'],
                      'selected_T_count':comparison['selected_T_count']})
    output={'schema':'qbp-independent-factor-basis-financial-composition-audit-v1','status':'passed',
        'prepared_bases_checked':len(facts),'financial_budget_scalar_checks':budgetchecks,
        'implemented_prior_scalar_checks':priorchecks,'flag_target_and_prior_checks':flagchecks,
        'preparations':facts,'paired_cost_interpretation':pairs,
        'source_sha256':sha(__file__),'resource_sha256':sha(rp/'factor_basis_quantum_resources.json'),
        'quantum_protocol_sha256':sha(rp/'factor_basis_quantum_protocol.json'),
        'certificate_sha256':sha(certificate_path),'mathematical_audit_sha256':sha(math_audit_path),
        'resource_audit_sha256':sha(resource_audit_path),
        'scope':'Independent exact composition of new finite moments and domains, raw transfer/cap normalisers, auxiliary monetary errors, coefficient-informed implemented probabilities and selected-estimator monetary targets. Continuous Gaussian invariance and financial transfer link to the all-58 mathematical audit; scalar/tape reconstruction links to its own independent resource audit. No new financial price, physical-device execution or originality certification.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:output[k] for k in ('status','prepared_bases_checked','financial_budget_scalar_checks','implemented_prior_scalar_checks','flag_target_and_prior_checks')}))

if __name__=='__main__':
    main()
