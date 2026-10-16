"""Bounded independent review of fresh financial and full QPE charge algebra.

No production helpers are imported or tables rebuilt. Emitted block
certification and score-function enclosures have separate full audits.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import comb, isqrt
from functools import lru_cache
import hashlib, json, time
from flint import arb, ctx

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE / 'quantum'
RESULTS = ROOT / 'results'
REPO = Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
EPS, ALPHA, PI, POST = F(1,100), F(1,3600), F(355,113), F(1,10**60)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()

def pair(record):
    return F(record['lower_exact_fraction']), F(record['upper_exact_fraction'])

def A(value):
    value=F(value)
    return arb(value.numerator)/value.denominator

def ball_bounds(value):
    def end(x):
        m,e=map(int,x.man_exp())
        return F(m*(1<<e)) if e>=0 else F(m,1<<-e)
    return end(value.lower()),end(value.upper())

def ceil(value):
    value = F(value)
    return -((-value.numerator)//value.denominator)

def square_root_upper(value):
    scale = 1 << 100
    k = ceil(value * scale * scale)
    r = isqrt(k)
    return F(r + (r*r < k), scale)

@lru_cache(None)
def binomial(n, success):
    return sum((F(comb(n,j)) * (1-success)**j * success**(n-j)
                for j in range((n+1)//2, n+1)), F(0))

def reflection(core):
    k = core - 1
    return {'CZ':2, 'X':2*core, 'H':4*k, 'CX':12*k,
            'Tdg':6*k, 'T':8*k, 'Z':1}

def financial_budget(row, preparation, financial):
    certificate = row['base_graph_certificate']
    cfg = certificate['configuration']
    b = row['complete_deterministic_budget']
    assert b is not None
    ml, mh = pair(financial['mass_interval'])
    mr = (mh-ml)/2
    base = preparation['coefficient_certificate']['c0' if row['kind']=='positive_r1' else 'coefficient']
    cr = (pair(base)[1]-pair(base)[0])/2
    B = F(cfg['monetary_signal_normaliser'])
    decode = F(cfg['monetary_probability_decode_factor'])
    eta = F(certificate['probability_error_upper'])
    aux = b['auxiliary_components']
    moment = {int(k):pair(v)[1] for k,v in preparation['finite_moments']['centred_U_moments'].items()}
    norm = preparation['payoff_norms']
    if row['kind']=='boundary_r2':
        expected_aux = F(1,1<<16) * (F(norm['C_r_upper']['3'])*moment[3]/12 + F(norm['C_r_upper']['4'])*moment[4]/48)
        assert F(aux['bias_upper']) == expected_aux
    elif row['kind']=='stable_gauss_r2':
        constants = aux['constants']
        node = sum((F(constants['w_'+s]['ideal_upper'])*F(constants['z_'+s]['error_upper']) for s in ('minus','plus')),F(0))
        weight = sum((F(constants['w_'+s]['error_upper']) for s in ('minus','plus')),F(0))
        remainder = F(norm['C_r_upper']['6'])*moment[6]/14400
        fixed = F(norm['C_r_upper']['3'])*moment[3]*node/2 + F(norm['M2_upper'])*moment[2]*weight/2
        assert F(aux['Gauss_remainder_upper']) == remainder
        assert F(aux['fixed_parameter_bias_upper']) == fixed
        assert F(aux['bias_upper']) == remainder + fixed
    else:
        assert F(aux['bias_upper']) == 0
    arithmetic = decode * B * eta
    expected = cr + F(financial['financial_transfer_upper']) + mh*(F(aux['bias_upper'])+arithmetic) + mr*B
    assert b['coefficient_interval'] == base and b['mass_interval'] == financial['mass_interval']
    assert F(b['coefficient_radius_upper']) == cr and F(b['external_factor_midpoint']) == (ml+mh)/2
    assert F(b['external_factor_upper']) == mh and F(b['external_factor_radius']) == mr
    assert F(b['arithmetic_before_mass_upper']) == arithmetic
    assert F(b['financial_transfer_upper']) == F(financial['financial_transfer_upper'])
    assert F(b['total_deterministic_undiscounted_price_error_upper']) == expected
    return expected

def previous_prior(row, preparation):
    cfg = row['base_graph_certificate']['configuration']
    eta = F(row['base_graph_certificate']['probability_error_upper'])
    B, U = F(cfg['monetary_signal_normaliser']), F(cfg['U_star'])
    moment = {int(k):pair(v)[1] for k,v in preparation['finite_moments']['centred_U_moments'].items()}
    prior = row['correction_prior']
    if row['kind']=='positive_r1':
        gap = pair(preparation['finite_moments']['correlated_underlying_gap_interval'])[1]
        mean = min(gap,F(1),F(cfg['M1'])*moment[1])
        assert F(prior['ideal_mean_upper']) == mean
        assert F(prior['probability_upper']) == min(F(1),mean/B+eta)
        assert prior['native_positive']
    else:
        positive = cfg['output_encoding']=='positive'
        assert prior['native_positive']==positive
        assert F(prior['probability_upper']) == min(F(1),moment[2]/U**2+(1 if positive else 2)*eta)

def informed_prior(record, row, preparation, financial):
    cfg = row['base_graph_certificate']['configuration']
    B, eta = F(cfg['monetary_signal_normaliser']), F(row['base_graph_certificate']['probability_error_upper'])
    c1lo,c1hi = pair(preparation['coefficient_certificate']['c1'])
    remainder = F(preparation['payoff_norms']['M2_upper'])*pair(preparation['continuous_second_moment']['continuous_U2_interval'])[1]/2
    beta = F(financial['financial_transfer_upper'])
    ml,mh = pair(financial['mass_interval'])
    lo,hi = max(F(0),c1lo-remainder-beta)/mh, min(B,(c1hi+remainder+beta)/ml,F(row['correction_prior']['ideal_mean_upper']))
    assert lo<=hi
    assert F(record['finite_ideal_mean_lower'])==lo and F(record['finite_ideal_mean_upper'])==hi
    assert F(record['probability_lower'])==max(F(0),lo/B-eta)
    assert F(record['probability_upper'])==min(F(1),F(row['correction_prior']['probability_upper']),hi/B+eta)
    assert record['continuous_moment_record_sha256']==digest(preparation['continuous_second_moment'])

def emitted_profile(profile):
    forward = Counter()
    for b in profile['blocks']:
        forward.update(b['primitive_counts'])
        assert b['logical_T_count']==b['primitive_counts'].get('T',0)+b['primitive_counts'].get('Tdg',0)
    assert forward==Counter(profile['primitive_counts_forward_body'])
    complete = Counter({g:2*n for g,n in forward.items()})
    complete.update(profile['comparator_primitive_counts'])
    complete['H']+=profile['uniform_H']
    assert complete==Counter(profile['primitive_counts'])
    assert profile['logical_T_count']==complete['T']+complete['Tdg']
    assert profile['lookup_invocations_forward_and_cleanup']==2*sum(b['lookup_invocations'] for b in profile['blocks'])
    assert profile['reflection_core_qubits']==profile['uniform_H']+1
    assert profile['active_qubits_after_clean_prep']==profile['reflection_core_qubits']+1
    assert profile['source_sha256']==sha(ROOT/'boundary_value_oracle.py')
    assert profile['retained_emitter_sha256']==sha(ROOT/'integral_value_oracle.py')

def complete_charge(flag, profile):
    plan, qft, charge = flag['plan'], flag['fourier_profile'], flag['complete_charge']
    m,M,R = plan['phase_register_bits'],plan['M'],plan['repetitions']
    assert M==1<<m and R>=1 and R%2==1
    epsilon,alpha,prior = map(F,(plan['epsilon'],plan['alpha'],plan['prior_upper']))
    assert epsilon==F(flag['probability_tolerance']) and prior==F(flag['prior_upper'])
    variance = prior*(1-prior) if prior<=F(1,2) else F(1,4)
    assert F(plan['variance_upper'])==variance
    error = 2*PI*square_root_upper(variance)/M + PI**2/M**2
    assert F(plan['amplitude_error_upper'])==error and F(plan['postprocess_allowance'])==POST
    assert error+POST<=epsilon
    assert M==2 or 2*PI*square_root_upper(variance)/(M//2)+PI**2/(M//2)**2+POST>epsilon
    norm = F(qft['norm_error_upper'])
    rounded = F(ceil(norm*(1<<80)),1<<80)
    assert F(plan['qft_norm_error_allowance'])==rounded and norm<=rounded<=F(1,64)
    success = 8/PI**2-rounded
    assert F(plan['good_probability_lower'])==success and success>F(1,2)
    fail = binomial(R,success)
    assert F(plan['binomial_failure_upper'])==fail<=alpha
    assert R==1 or binomial(R-2,success)>alpha
    assert qft['phase_register_bits']==m
    assert qft['T_count']==qft['gate_counts'].get('T',0)+qft['gate_counts'].get('Tdg',0)
    assert flag['financial_profile_sha256']==digest(profile)
    ap = Counter(profile['primitive_counts'])
    inverse = Counter()
    swap = {'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    for g,n in ap.items():inverse[swap.get(g,g)]+=n
    ref = reflection(profile['reflection_core_qubits'])
    assert Counter(charge['controlled_reflection']['gate_counts'])==Counter(ref)
    qp=ap+inverse;qp.update(ref)
    at=ap['T']+ap['Tdg'];rt=ref['T']+ref['Tdg'];qt=2*at+rt
    assert charge['complete_oracle_T']==at and charge['controlled_Q_T']==qt
    assert charge['controlled_reflection']['global_minus_Z_charged'] and charge['controlled_reflection']['controlled_good_CZ_charged']
    assert charge['controlled_reflection']['T']==rt
    peak=max(profile['logical_qubits_upper'],2*profile['reflection_core_qubits'])+m+qft.get('extra_ancillas',0)
    assert charge['peak_logical_qubits']==peak
    primitive=Counter(ap)
    primitive.update({g:(M-1)*n for g,n in qp.items()})
    primitive.update(qft['gate_counts']);primitive['H']+=m
    per=charge['per_repetition_metrics']
    for g,n in primitive.items():assert per['gate_'+g]==n
    assert per['T']==at+(M-1)*qt+qft['T_count']
    assert per['expanded_primitive_operations']==sum(primitive.values())+2*peak
    assert per['uniform_H']==profile['uniform_H']*(1+2*(M-1))+m
    assert per['measurements']==peak and per['reset_X']==peak
    for depth,refdepth in [('explicit_T_depth_serial_upper',8*(profile['reflection_core_qubits']-1)),('conservative_5layer_T_depth_serial_upper',10*(profile['reflection_core_qubits']-1))]:
        assert per[depth]==profile[depth]+(M-1)*(2*profile[depth]+refdepth)+qft['metrics'][depth]
    assert all(charge['metrics'][k]==R*v for k,v in per.items())
    assert charge['repetitions']==R and charge['controlled_Q_calls']==R*(M-1)
    assert charge['phase_register_measurements']==R*m and charge['complete_reinitialisations']==R

def main():
    start=time.perf_counter()
    ctx.prec=1024
    exact_mass=(A(6)/A(2).sqrt()).erf()
    protocol_path=RESULTS/'four_survivor_quantum_protocol.json'
    full=RESULTS/'four_survivor_quantum_resources.json'
    protocol=json.loads(protocol_path.read_text());data=json.loads(full.read_text())
    assert sha(protocol_path)==data['protocol_sha256']
    assert sha(ROOT/'four_survivor_quantum_resources.py')==protocol['source_sha256']==data['source_sha256']
    assert data['dependency_source_sha256']==protocol['dependency_source_sha256']
    for f,h in protocol['dependency_source_sha256'].items():assert sha(ROOT/f)==h
    assert data['input_sha256']==protocol['input_sha256']
    for f,h in protocol['input_sha256'].items():assert sha(RESULTS/f)==h and sha(REPO/'results'/f)==h
    models=json.loads((REPO/'results/dispersion_separation_protocol.json').read_text())['models']
    installed={m['model_id']:m for m in models}
    assert all(m==installed[m['model_id']] for m in protocol['models'])
    classical=json.loads((STAGE/'classical/results/four_case_classical_protocol.json').read_text())
    assert [m['model_id'] for m in protocol['models']]==classical['case_ids']
    assert protocol['classical_price_answer_excluded_from_inputs']
    assert F(data['monetary_radius'])==EPS and F(data['alpha_per_price'])==ALPHA
    assert F(data['simultaneous_four_price_union_failure_upper'])==4*ALPHA==F(1,900)
    assert sha(RESULTS/data['scalar_ordinate_archive'])==data['scalar_ordinate_archive_sha256']
    prep={p['model']['model_id']:p for p in data['preparations']}
    counts={'synthetic_models':len(prep),'direct_R1_coordinates':0,'financial_budgets':0,
            'priors':0,'emitted_profiles':0,'full_QPE_flag_charges':0,'strategies':0,'failures':0}
    for p in prep.values():
        model=p['model'];direct=p['R1_direct_global_financial_transfer'];norm=p['payoff_norms']
        assert model==installed[model['model_id']]
        w=list(map(F,model['weights']));R=[list(map(F,row)) for row in model['loadings']]
        m1,m2=F(norm['M1_upper']),F(norm['M2_upper'])
        linear=quadratic=F(0)
        for j,d in enumerate(direct['coordinates']):
            Aj=max(abs(row[j]) for row in R);r=sum(wi*row[j] for wi,row in zip(w,R))
            d1=m1*(Aj+abs(r));d2=(m1+m2)*(Aj*Aj+r*r)
            assert F(d['full_asset_loading_max'])==Aj and F(d['geometric_loading'])==r
            assert F(d['D1_upper'])==d1 and F(d['D2_upper'])==d2
            l=F(direct['other_coordinate_mass_upper'])*F(direct['rho1_mass_upper'])*d1
            q=F(direct['other_coordinate_mass_upper'])*F(direct['rho2_mass_upper'])*d2
            assert F(d['linear_error_upper'])==l and F(d['quadratic_error_upper'])==q
            linear+=l;quadratic+=q;counts['direct_R1_coordinates']+=1
        assert F(direct['linear_transport_upper'])==linear and F(direct['quadratic_transport_upper'])==quadratic
        assert F(direct['financial_transfer_upper'])==linear+quadratic+F(direct['tail_upper'])
        ml,mh=pair(direct['mass_interval'])
        mass_low,mass_high=ball_bounds(exact_mass**model['outer_dimension'])
        assert ml<=mass_low<=mass_high<=mh
        assert F(direct['tail_upper'])>=ball_bounds(1-exact_mass**model['outer_dimension'])[1]
        assert F(direct['other_coordinate_mass_upper'])>=ball_bounds(exact_mass**(model['outer_dimension']-1))[1]
        old=p['R1_older_financial_transfer'];selected=p['R1_selected_financial_transfer']
        assert F(selected['financial_transfer_upper'])==min(F(old['financial_transfer_upper']),F(direct['financial_transfer_upper']))
    for row in data['records']:
        p=prep[row['model_id']]
        assert row['scalar_bits'] in (12,14)
        assert set(row['node_address_bits'].values())=={row['scalar_bits']}
        assert row['scalar_bits']==12 or (row['kind']=='boundary_r2' and row['candidate_label']=='all14_declared_boundary_repair')
        assert row['base_graph_certificate']['configuration']['q']==28
        assert row['base_graph_certificate']['configuration']['b_y']==10
        assert F(row['base_graph_certificate']['configuration']['L'])==6
        financial=p['R1_selected_financial_transfer'] if row['kind']=='positive_r1' else p['R2_financial_transfer']
        if row['complete_deterministic_budget'] is not None:
            det=financial_budget(row,p,financial);counts['financial_budgets']+=1
        if 'baseline_failure' in row:
            assert row['scalar_bits']==12 and (row['complete_deterministic_budget'] is None or det>=EPS)
            counts['failures']+=1;continue
        previous_prior(row,p);counts['priors']+=1
        if row['kind']=='positive_r1':
            informed_prior(row['coefficient_informed_prior'],row,p,financial)
            oldrow=dict(row,complete_deterministic_budget=row['older_financial_deterministic_budget'])
            financial_budget(oldrow,p,p['R1_older_financial_transfer']);counts['financial_budgets']+=1
            informed_prior(row['older_coefficient_informed_prior'],row,p,p['R1_older_financial_transfer']);counts['priors']+=2
        for profiles in row['profiles'].values():
            for profile in profiles.values():emitted_profile(profile);counts['emitted_profiles']+=1
        for strategy in row['strategies']:
            result=strategy['result'];counts['strategies']+=1
            old=strategy['prior_information'].startswith('older_financial')
            budget=row['older_financial_deterministic_budget'] if old else row['complete_deterministic_budget']
            det=F(budget['total_deterministic_undiscounted_price_error_upper'])
            if result['kind']=='rejected':assert det>=EPS;continue
            if result['kind']=='deterministic_correction_prior':
                assert result['T_count']==0 and result['quantum_shots']==0
                assert F(result['total_price_error_upper'])==det+F(result['certified_correction_radius'])<=EPS
                continue
            cfg=row['base_graph_certificate']['configuration']
            C=F(budget['external_factor_upper'])*F(cfg['monetary_signal_normaliser'])
            info=strategy['prior_information']
            expected_prior=(row['coefficient_informed_prior'] if info=='coefficient_informed_direct_minimum'
                else row['older_coefficient_informed_prior'] if info=='older_financial_coefficient_informed_prior'
                else row['correction_prior'])
            assert result['correction_prior']==expected_prior
            positive=result['correction_prior']['native_positive']
            decode=1 if positive or strategy['strategy']=='split' else 2
            target=(EPS-det)/(C*(1 if positive else 2))
            alpha=ALPHA/(2 if strategy['strategy']=='split' and not positive else 1)
            total_T=total_R=0;stat=failure=F(0)
            expected_labels=['centred'] if positive or strategy['strategy']=='centred' else ['positive_part','negative_part']
            assert [flag['flag'] for flag in result['flags']]==expected_labels
            for flag in result['flags']:
                assert F(flag['probability_tolerance'])==target and F(flag['plan']['alpha'])==alpha
                expected_upper=F(1) if not positive and strategy['strategy']=='centred' else F(expected_prior['probability_upper'])
                assert F(flag['prior_upper'])==expected_upper
                profile=row['profiles'][strategy['architecture']][flag['flag']]
                complete_charge(flag,profile);counts['full_QPE_flag_charges']+=1
                total_T+=flag['complete_charge']['metrics']['T'];total_R+=flag['plan']['repetitions']
                stat+=F(flag['plan']['amplitude_error_upper'])+F(flag['plan']['postprocess_allowance'])
                failure+=F(flag['plan']['binomial_failure_upper'])
            stat*=C*decode
            assert result['T_count']==total_T and result['quantum_shots']==total_R
            assert F(result['monetary_signal_factor_upper'])==C
            assert F(result['deterministic_error_upper'])==det
            assert F(result['statistical_price_error_upper'])==stat
            assert F(result['total_price_error_upper'])==det+stat<=EPS
            assert F(result['failure_probability_upper'])==failure<=ALPHA
        valid=[s for s in row['strategies'] if s['result']['kind']!='rejected']
        assert row['selected_strategy']==min(valid,key=lambda s:(s['result']['T_count'],s['architecture'],s['strategy'],s['prior_information']))
    assert counts['emitted_profiles']==data['profiles_count'] and counts['full_QPE_flag_charges']==data['confidence_plan_count']
    for selected in data['selected_prices']:
        for kind,choice in selected['methods'].items():
            eligible=[r for r in data['records'] if r['model_id']==selected['model_id'] and r['kind']==kind and r['selected_strategy']]
            winner=min(eligible,key=lambda r:(r['selected_strategy']['result']['T_count'],r['candidate_index']))
            assert choice['candidate_index']==winner['candidate_index'] and choice['selected']==winner['selected_strategy']
        winner=min(selected['methods'],key=lambda k:selected['methods'][k]['selected']['result']['T_count'])
        assert selected['quantum_resource_winner']==winner
    summary=json.loads((RESULTS/'four_survivor_quantum_resource_summary.json').read_text())
    assert summary['full_result_sha256']==sha(full)
    record={'schema':'qbp-four-case-independent-quantum-algebra-audit-v1','status':'passed',
            'date':'2026-10-09','auditor_sha256':sha(__file__),'result_sha256':sha(full),
            'protocol_sha256':sha(protocol_path),'producer_sha256':data['source_sha256'],
            'counts':counts,'quantum_data_origin_verified':'Exact same four installed synthetic stress-model dictionaries as the classical protocol; no classical price answer or LSEG input.',
            'compiler_distinction':'The body cleanup re-emits self-inverse logical operations, so its expanded gate multiset matches the forward body. The complete A adjoint within controlled Q separately swaps T/Tdg and S/Sdg in expanded profiles.',
            'limits':'This bounded audit independently checks financial and prior algebra, full block/profile arithmetic, actual QPE phase/error/binomial linkage, reflection and reset charges and input/source hashes. It does not rebuild or revalidate all scalar tables, execute quantum hardware, prove classical superiority or establish original theory.',
            'licensed_data_access_or_deletion':False,'seconds_diagnostic':time.perf_counter()-start}
    (STAGE/'licence_review/final_quantum_algebra_audit.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Passed',counts)

if __name__=='__main__':
    main()
