"""Independent complete monetary, confidence and reset audit of tuned plans."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from math import comb,isqrt
import argparse
import gzip
import hashlib
import json
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def low(iv):return F(iv['lower_exact_fraction'])
def high(iv):return F(iv['upper_exact_fraction'])
def ball(x):
    from flint import arb
    x=F(x);return arb(x.numerator)/arb(x.denominator)


def first_quantile(path):
    with gzip.open(path,'rt') as stream:
        for line in stream:
            item=json.loads(line);certificate=item['certificate']
            if certificate['schema']=='qbp-integral-quantile-v1':
                assert hashlib.sha256(json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()).hexdigest()==item['sha256']
                return certificate
    raise AssertionError('No retained quantile certificate')


def exact_failure(shots,success):
    return sum((F(comb(shots,k))*(1-success)**k*success**(shots-k)
                for k in range((shots+1)//2,shots+1)),F(0))


def fresh_finite_moments(models,quantile,moment_records):
    from flint import arb,ctx
    ctx.prec=416;ys=[F(n['actual_coordinate']) for n in quantile['nodes']];records={}
    for model_id,model in models.items():
        if model['outer_dimension']!=1:continue
        weights=list(map(F,model['weights']));mu=list(map(F,model['log_drifts']));R=[F(row[0]) for row in model['loadings']]
        mg=sum(w*v for w,v in zip(weights,mu));rg=sum(w*v for w,v in zip(weights,R));lam=F(model['lambda'])
        d=[v-mg for v in mu];s=[v-rg for v in R]
        totals={j:arb(0) for j in (1,2,3)};underlying_gap=arb(0)
        for y in ys:
            ratio=sum((ball(w)*ball(di+si*y).exp() for w,di,si in zip(weights,d,s)),arb(0))
            U=ratio-1;assert U.lower()>=0
            for j in totals:totals[j]+=U**j
            underlying_gap+=100*ball(mg+rg*y+lam*lam/2).exp()*U
        for j in totals:
            totals[j]/=len(ys);stored=moment_records[model_id]['centred_U_moments'][str(j)]
            assert ball(low(stored)).lower()<=totals[j].lower()
            assert ball(high(stored)).upper()>=totals[j].upper()
        underlying_gap/=len(ys)
        records[model_id]={'moments':totals,'underlying_gap':underlying_gap,
            'finite_labels':len(ys),'scope':'Independent 416-bit direct sums on the unchanged actual outer grid'}
    return records


def reconstruct_auxiliary(H,lam,U,m2,m3,joint,hermite):
    h=F(1,256);E=hermite[1]/lam+hermite[2]/lam**2;F2=hermite[2]/lam+hermite[3]/lam**2
    D3,D4=joint[3,lam],joint[4,lam]
    beta=min(H*D3*m3*h/4,H*h*h*(D3*m3/12+D4*U*m3/48))
    strike=min(H*H*E*m2*h/(4*100*lam),H**3*m2*h*h*(F2/lam**2+E/lam)/(48*100**2))
    return beta+strike


def reflection(core):
    ccx=2*(core-1)
    return Counter({'X':2*core,'CZ':2,'Z':1,'H':2*ccx,'CX':6*ccx,'T':4*ccx,'Tdg':3*ccx})


def audit_charge(flag,profile):
    setup=flag['plan'];fourier=flag['fourier_profile'];charge=flag['complete_charge']
    m=int(setup['phase_register_bits']);M=int(setup['M']);shots=int(setup['repetitions'])
    assert 1<=m<=64 and M==1<<m and shots>=1 and shots%2==1
    assert fourier['phase_register_bits']==m
    norm=F(fourier['norm_error_upper']);rounded=F(setup['qft_norm_error_allowance'])
    scale=1<<80;ceiling=(norm*scale).numerator//(norm*scale).denominator
    if F(ceiling,scale)<norm:ceiling+=1
    assert rounded==F(ceiling,scale)
    success=8/F(355,113)**2-rounded
    assert F(setup['good_probability_lower'])==success>F(1,2)
    failure=exact_failure(shots,success)
    assert F(setup['binomial_failure_upper'])==failure<=F(setup['alpha'])
    if shots>1:assert exact_failure(shots-2,success)>F(setup['alpha'])
    prior=F(setup['prior_upper']);var=prior*(1-prior) if prior<=F(1,2) else F(1,4)
    assert F(setup['variance_upper'])==var
    # The stored root enclosure is assessed through its amplitude expression;
    # this independent inequality uses the exact square to avoid a float root.
    error=F(setup['amplitude_error_upper']);rest=error-F(355,113)**2/(M*M)
    assert rest>=0 and rest*rest>=4*F(355,113)**2*var/(M*M)
    scale_root=1<<100;scaled=var*scale_root*scale_root
    integer=(scaled.numerator+scaled.denominator-1)//scaled.denominator
    root=isqrt(integer)
    if root*root<integer:root+=1
    exact_root=F(root,scale_root)
    expected_error=2*F(355,113)*exact_root/M+F(355,113)**2/(M*M)
    assert error==expected_error
    allowance=F(setup['postprocess_allowance'])
    assert error+allowance<=F(setup['epsilon'])==F(flag['probability_tolerance'])
    if m>1:assert 2*F(355,113)*exact_root/(M//2)+F(355,113)**2/(M//2)**2+allowance>F(setup['epsilon'])
    assert F(setup['prior_upper'])==F(flag['prior_upper'])
    ap=Counter(profile['primitive_counts']);inverse=Counter()
    swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    for gate,count in ap.items():inverse[swap.get(gate,gate)]+=count
    rc=reflection(profile['reflection_core_qubits']);cq=ap+inverse+rc
    fp=Counter(fourier['gate_counts']);per=ap+Counter({g:(M-1)*n for g,n in cq.items()})+fp;per['H']+=m
    qubits=max(profile['logical_qubits_upper'],2*profile['reflection_core_qubits'])+m
    assert charge['peak_logical_qubits']==qubits
    assert charge['controlled_Q_calls']==shots*(M-1)
    assert charge['complete_oracle_T']==ap['T']+ap['Tdg']
    assert charge['controlled_Q_T']==cq['T']+cq['Tdg']
    for gate,count in per.items():assert charge['metrics'][f'gate_{gate}']==shots*count
    T=per['T']+per['Tdg']
    assert charge['metrics']['T']==shots*T
    assert charge['metrics']['measurements']==shots*qubits
    assert charge['metrics']['reset_X']==shots*qubits
    assert charge['metrics']['expanded_primitive_operations']==shots*(sum(per.values())+2*qubits)
    assert charge['phase_register_measurements']==shots*m
    assert charge['complete_reinitialisations']==shots
    for key,value in charge['per_repetition_metrics'].items():assert charge['metrics'][key]==shots*value
    explicit_A=profile['explicit_T_depth_serial_upper'];five_A=profile['conservative_5layer_T_depth_serial_upper']
    explicit_Q=2*explicit_A+8*(profile['reflection_core_qubits']-1)
    five_Q=2*five_A+10*(profile['reflection_core_qubits']-1)
    assert charge['metrics']['explicit_T_depth_serial_upper']==shots*(explicit_A+(M-1)*explicit_Q+fourier['T_count'])
    assert charge['metrics']['conservative_5layer_T_depth_serial_upper']==shots*(five_A+(M-1)*five_Q+fourier['T_count'])
    return shots*T,failure,error+F(setup['postprocess_allowance'])


def evaluate(root,path):
    root=Path(root);path=Path(path);data=json.loads(gzip.open(path,'rt').read())
    assert data['source_sha256']==sha(root/'integral_precision_selection.py')
    for name,digest in data['dependency_source_hashes'].items():assert sha(root/name)==digest
    for name,digest in data['input_sha256'].items():assert sha(root/'results'/name)==digest
    models={m['model_id']:m for m in json.loads((root/'results/generalisation_protocol.json').read_text())['models']}
    moments={r['model_id']:r for r in json.loads((root/'results/integral_finite_grid_moments.json').read_text())['records']}
    finances={(r['model_id'],F(r['cap'])):r for r in json.loads((root/'results/integral_financial_grid_assessment.json').read_text())['records']}
    norms=json.loads((root/'results/gaussian_derivative_norms.json').read_text())
    joint={(r['derivative_order'],F(r['lambda'])):F(r['norm_upper']) for r in norms['records'] if r['passed']}
    hermite={r['unit_derivative_order']:F(r['norm_upper']) for r in norms['unit_norms'] if r['passed']}
    archive=root/'results'/data['scalar_ordinate_archive'];assert data['scalar_ordinate_archive_sha256']==sha(archive)
    quantile=first_quantile(archive);flo,fhi=low(quantile['mass_interval']),high(quantile['mass_interval'])
    fresh=fresh_finite_moments(models,quantile,moments)
    epsilon=F(1,100);alpha=F(1,2700);counts=Counter();records=[]
    for index,candidate in enumerate(data['records']):
        model_id=candidate['model_id'];kind=candidate['kind'];H=F(candidate['cap']);lam=F(models[model_id]['lambda'])
        cert=candidate['base_graph_certificate'];cfg=cert['configuration'];budget=candidate['complete_deterministic_budget']
        B=F(cfg['monetary_signal_normaliser']);U=F(cfg['U_star']);eta=F(cert['probability_error_upper'])
        if cert['status']!='passed':
            assert all(s['result']['kind']=='rejected' for s in candidate['strategies']);counts['guard_rejected_candidates']+=1;continue
        m1,m2,m3=[high(moments[model_id]['centred_U_moments'][str(j)]) for j in (1,2,3)]
        finance=finances[model_id,H]['kinds']['R1' if kind=='positive_r1' else 'R2']
        coefficient_radius=(high(finance['coefficient_interval'])-low(finance['coefficient_interval']))/2
        arithmetic=B*eta*(1 if kind=='positive_r1' else 2)
        aux=reconstruct_auxiliary(H,lam,U,m2,m3,joint,hermite) if kind=='integral_r2' else F(0)
        transfer=F(finance['financial_transfer_including_tail_and_quantile_rounding_upper'])
        det=transfer+coefficient_radius+fhi*(arithmetic+aux)+(fhi-flo)*B/2
        assert F(budget['arithmetic_before_factor_upper'])==arithmetic
        assert F(budget['auxiliary_before_factor_upper'])==aux
        assert F(budget['coefficient_radius_upper'])==coefficient_radius
        assert F(budget['financial_transfer_including_tail_and_node_rounding_upper'])==transfer
        assert F(budget['external_factor_upper'])==fhi
        assert F(budget['total_deterministic_undiscounted_price_error_upper'])==det
        C=fhi*B;remaining=epsilon-det
        if kind=='positive_r1':
            prior_row=candidate['positive_r1_prior'];gap=prior_row['correlated_underlying_gap_interval']
            assert ball(low(gap)).lower()<=fresh[model_id]['underlying_gap'].lower()
            assert ball(high(gap)).upper()>=fresh[model_id]['underlying_gap'].upper()
            mean=min(H,H*joint[1,lam]*m1,high(gap));prior=min(F(1),mean/B+eta)
            assert F(prior_row['positive_R1_mean_upper'])==mean
            assert F(prior_row['probability_prior_upper'])==prior
            prior_radius=C*prior/2
        else:
            prior=min(F(1),m2/(U*U)+2*eta);prior_radius=C*prior
        for strategy in candidate['strategies']:
            result=strategy['result'];encoding=strategy['strategy'];architecture=strategy['architecture'];counts[result['kind']]+=1
            if det>=epsilon:
                assert result['kind']=='rejected';continue
            if prior_radius<=remaining:
                assert result['kind']=='deterministic_correction_prior'
                assert result['T_count']==result['quantum_shots']==0
                assert F(result['certified_correction_radius'])==prior_radius
                assert det+prior_radius<=epsilon
                continue
            assert result['kind']=='complete_quantum_resource_plan'
            labels=['centred'] if kind=='positive_r1' or encoding=='centred' else ['positive_part','negative_part']
            target=remaining/C if kind=='positive_r1' else remaining/(2*C)
            assert [flag['flag'] for flag in result['flags']]==labels
            totalT=0;failure=F(0);probability_error=F(0);shots=0
            for flag in result['flags']:
                expected_prior=prior if kind=='positive_r1' or encoding=='split' else F(1)
                expected_alpha=alpha/2 if encoding=='split' else alpha
                assert F(flag['prior_upper'])==expected_prior
                assert F(flag['probability_tolerance'])==target
                assert F(flag['plan']['alpha'])==expected_alpha
                profile=candidate['profiles'][architecture][flag['flag']]
                profile_hash=hashlib.sha256(json.dumps(profile,sort_keys=True,separators=(',',':')).encode()).hexdigest()
                assert profile_hash==flag['financial_profile_sha256']
                nt,fail,error=audit_charge(flag,profile);totalT+=nt;failure+=fail;probability_error+=error;shots+=flag['plan']['repetitions']
            money=C*probability_error*(2 if kind!='positive_r1' and encoding=='centred' else 1)
            assert F(result['statistical_price_error_upper'])==money
            assert F(result['total_price_error_upper'])==det+money<=epsilon
            assert F(result['failure_probability_upper'])==failure<=alpha
            assert result['T_count']==totalT and result['quantum_shots']==shots
        eligible=[s for s in candidate['strategies'] if s['result']['kind']!='rejected']
        expected=min(eligible,key=lambda s:(s['result']['T_count'],s['architecture'],s['strategy'])) if eligible else None
        assert candidate['selected_strategy']==expected
        records.append({'candidate_index':index,'model_id':model_id,'cap_fraction':str(H),'kind':kind,
            'deterministic_error_upper_fraction':str(det),'selected_T_count':expected['result']['T_count'] if expected else None})
    for price in data['selected_prices']:
        for kind,choice in price['methods'].items():
            candidates=[(i,r) for i,r in enumerate(data['records']) if r['model_id']==price['model_id'] and F(r['cap'])==F(price['cap']) and r['kind']==kind and r['selected_strategy'] is not None]
            if not candidates:assert choice is None;continue
            selected=min(candidates,key=lambda pair:pair[1]['selected_strategy']['result']['T_count'])
            assert choice['candidate_index']==selected[0] and choice['selected']==selected[1]['selected_strategy']
        assert 'zero quantum estimation calls' in price['complete_policy']
    classical={}
    hierarchy=json.loads((root/'results/gaussian_norm_pricing.json').read_text())['executed_hierarchy']['configurations']
    for row in hierarchy:
        item=row.get('selected')
        if not item or not item.get('certified'):continue
        iv=item['interval'];radius=(high(iv)-low(iv))/2
        assert radius<=epsilon
        classical[row['case_id']]={'source':'retained certified Taylor interval','radius_upper_fraction':str(radius)}
    continuous=json.loads((root/'results/continuous_quadrature.json').read_text())['configurations']
    for row in continuous:
        iv=row['undiscounted_price_interval'];radius=(high(iv)-low(iv))/2
        assert radius<=epsilon
        case=f"{row['model_id']}_cap{row['cap']}"
        if case not in classical or radius<F(classical[case]['radius_upper_fraction']):
            classical[case]={'source':'retained independently validated continuous interval','radius_upper_fraction':str(radius)}
    expected_cases={f"{m['model_id']}_cap{H}" for m in models.values() for H in (1,10,25)}
    assert set(classical)==expected_cases
    assert 27*alpha==F(1,100)
    return {'status':'independent precision monetary and complete estimator audit passed','candidate_count':len(data['records']),
        'strategy_kind_counts':dict(counts),'records':records,'fresh_finite_grid_models':len(fresh),
        'fresh_finite_label_evaluations':sum(r['finite_labels'] for r in fresh.values()),
        'retained_classical_zero_quantum_policy_intervals':classical,
        'comparison_sha256':sha(path),'scalar_ordinate_archive_sha256':sha(archive),
        'source_sha256':data['source_sha256'],'scope':'Exact monetary/confidence/count arithmetic and independently recomputed finite moment/prior premises; no physical quantum execution or advantage assertion'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--comparison',default=None)
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_precision_budget_audit.json'))
    args=parser.parse_args();result=evaluate(args.root,args.comparison or Path(args.root)/'results/integral_precision_comparison.json.gz');result['audit_source_sha256']=sha(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','candidate_count','strategy_kind_counts','fresh_finite_label_evaluations')}))

if __name__=='__main__':main()
