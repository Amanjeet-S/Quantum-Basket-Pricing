"""Own-formula reconstruction of finite budgets, tapes, confidence and charges.

The producer and its resource/confidence APIs are not imported. This audit
checks retained synthetic evidence and complete ideal-logical resource plans.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from math import comb
from functools import lru_cache
import gzip,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def lo(iv):return F(iv['lower_exact_fraction'])
def hi(iv):return F(iv['upper_exact_fraction'])
def clean(c):return {k:v for k,v in c.items() if v}
def failure(n,s):return sum((F(comb(n,j))*(1-s)**j*s**(n-j) for j in range((n+1)//2,n+1)),F(0))
PI=F(355,113); EPS=F(1,100); ALPHA=F(1,3600)

def fourier(bits,tapes,bounds):
    seq=[];error=F(0);phase=F(0);retained=omitted=0
    for j in range(bits//2):
        k=bits-1-j;seq.extend([('CX',j,k),('CX',k,j),('CX',j,k)])
    for j in range(bits-1,-1,-1):
        for k in range(bits-1,j,-1):
            d=k-j
            if d>31:error+=PI/F(1<<d);omitted+=1;continue
            rm=tapes['rotations'][f'-1:{d+1}'];rp=tapes['rotations'][f'1:{d+1}']
            minus=rm['chronological_gates'];plus=rp['chronological_gates']
            if minus or plus:
                seq.extend((g,k) for g in minus);seq.extend((g,j) for g in minus)
                seq.append(('CX',k,j));seq.extend((g,j) for g in plus);seq.append(('CX',k,j))
            error+=2*bounds[f'-1:{d+1}']+bounds[f'1:{d+1}'];retained+=1
            phase+=-F(1,1<<(d+2))+2*F(rm['phase_pi'])+F(rp['phase_pi'])
        seq.append(('H',j))
    return seq,error,phase,retained,omitted

def reflection(core):
    # Reconstruct the conjunction ladder and the literal exact CCX expansion.
    ctrl=core;work=list(range(core+1,2*core));body=[]
    body.append(('CCX',0,1,work[0]))
    for j in range(2,core):body.append(('CCX',work[j-2],j,work[j-1]))
    logical=[('CZ',ctrl,core-1)]+[('X',j) for j in range(core)]+body
    logical+=[('CZ',work[core-2],ctrl)]+body[::-1]+[('X',j) for j in range(core)]+[('Z',ctrl)]
    seq=[]
    for g in logical:
        if g[0]!='CCX':seq.append(g);continue
        _,a,b,t=g
        seq.extend([('H',t),('CX',b,t),('Tdg',t),('CX',a,t),('T',t),
            ('CX',b,t),('Tdg',t),('CX',a,t),('T',b),('T',t),
            ('H',t),('CX',a,b),('T',a),('Tdg',b),('CX',a,b)])
    return seq,2*(core-1)

def verify_profile(profile,known_tables):
    forward=Counter()
    for b in profile['blocks']:
        c=Counter(b['primitive_counts']);assert b['logical_T_count']==c['T']+c['Tdg']
        assert not set(c)-{'X','H','CX','CZ','Z','S','Sdg','T','Tdg'}
        if b['table_sha256'] is not None:assert b['table_sha256'] in known_tables
        forward.update(c)
    assert clean(forward)==clean(Counter(profile['primitive_counts_forward_body']))
    complete=Counter({g:2*n for g,n in forward.items()})
    complete.update(profile['comparator_primitive_counts']);complete['H']+=profile['uniform_H']
    assert clean(complete)==clean(Counter(profile['primitive_counts']))
    assert complete['T']+complete['Tdg']==profile['logical_T_count']
    assert profile['explicit_T_depth_serial_upper']==2*sum(b['explicit_T_depth_serial_upper'] for b in profile['blocks'])+profile_comparator_depth(profile,'explicit')
    assert profile['conservative_5layer_T_depth_serial_upper']==2*sum(b['conservative_5layer_T_depth_serial_upper'] for b in profile['blocks'])+profile_comparator_depth(profile,'five')
    return complete

@lru_cache(None)
def comparator(width):
    # Construct the subtraction ripple from Boolean operations independently.
    a=list(range(width));b=list(range(width,2*width));flag=2*width
    next_wire=flag+1;carry=next_wire;next_wire+=1;body=[]
    for ai,bi in zip(a,b):
        si=next_wire;next_wire+=1
        body.extend(('CX',x,si) for x in (ai,bi,carry))
        na=next_wire;next_wire+=1;body.extend([('X',na),('CX',ai,na)])
        nc=next_wire;next_wire+=1
        body.extend([('CCX',na,bi,nc),('CCX',na,carry,nc),('CCX',bi,carry,nc)])
        carry=nc
    logical=body+[('CX',carry,flag)]+body[::-1];sequence=[]
    five=[0]*next_wire
    for g in logical:
        level=max(five[w] for w in g[1:])+(5 if g[0]=='CCX' else 0)
        for w in g[1:]:five[w]=level
        if g[0]!='CCX':sequence.append(g);continue
        _,aa,bb,t=g
        sequence.extend([('H',t),('CX',bb,t),('Tdg',t),('CX',aa,t),('T',t),
            ('CX',bb,t),('Tdg',t),('CX',aa,t),('T',bb),('T',t),
            ('H',t),('CX',aa,bb),('T',aa),('Tdg',bb),('CX',aa,bb)])
    explicit=[0]*next_wire
    for g in sequence:
        level=max(explicit[w] for w in g[1:])+(g[0] in ('T','Tdg'))
        for w in g[1:]:explicit[w]=level
    return Counter(g[0] for g in sequence),max(explicit),max(five)

def profile_comparator_depth(profile,mode):
    counts,explicit,five=comparator(profile['comparison_bits']+1)
    assert clean(counts)==clean(Counter(profile['comparator_primitive_counts']))
    return explicit if mode=='explicit' else five

def verify_plan(flag,profile,tapes,bounds):
    plan=flag['plan'];R=plan['repetitions'];M=plan['M'];m=plan['phase_register_bits']
    assert M==1<<m and R>0 and R%2==1
    epsilon=F(flag['probability_tolerance']);prior=F(flag['prior_upper']);alpha=F(plan['alpha'])
    assert F(plan['epsilon'])==epsilon and F(plan['prior_upper'])==prior
    variance=prior*(1-prior) if prior<=F(1,2) else F(1,4)
    assert F(plan['variance_upper'])==variance and prior>2*epsilon
    seq,error,phase,retained,omitted=fourier(m,tapes,bounds);counts=Counter(g[0] for g in seq)
    qft=flag['fourier_profile'];assert digest(profile)==flag['financial_profile_sha256']
    assert qft['chronological_tape_sha256']==hashlib.sha256(json.dumps(seq,separators=(',',':')).encode()).hexdigest()
    assert error==F(qft['norm_error_upper']) and phase==F(qft['global_phase_correction_pi'])
    assert clean(counts)==clean(Counter(qft['gate_counts']))
    assert qft['T_count']==counts['T']+counts['Tdg']
    assert (retained,omitted)==(qft['controlled_phases_retained'],qft['controlled_phases_omitted'])
    scale=1<<80;scaled=error*scale;rounded=F(-(-scaled.numerator//scaled.denominator),scale)
    assert F(plan['qft_norm_error_allowance'])==rounded and rounded<=F(1,64)
    success=8/PI**2-rounded;assert success==F(plan['good_probability_lower'])
    fail=failure(R,success);assert fail==F(plan['binomial_failure_upper']) and fail<=alpha
    if R>1:assert failure(R-2,success)>alpha
    ae=F(plan['amplitude_error_upper']);root=(ae-PI**2/M**2)*M/(2*PI)
    assert root>=0 and root*root>=variance
    assert ae+F(plan['postprocess_allowance'])<=epsilon
    if M>2:assert 2*PI*root/(M//2)+PI**2/(M//2)**2+F(plan['postprocess_allowance'])>epsilon
    ap=Counter(profile['primitive_counts']);adj=Counter()
    swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
    for op,n in ap.items():adj[swap.get(op,op)]+=n
    refl,ccx=reflection(profile['reflection_core_qubits']);rc=Counter(g[0] for g in refl)
    qp=ap+adj+rc;per=ap+Counter({op:(M-1)*n for op,n in qp.items()})+counts;per['H']+=m
    charge=flag['complete_charge'];ref=charge['controlled_reflection']
    assert clean(rc)==clean(Counter(ref['gate_counts']))
    assert ref['chronological_tape_sha256']==hashlib.sha256(json.dumps(refl,separators=(',',':')).encode()).hexdigest()
    assert ref['T']==14*(profile['reflection_core_qubits']-1)
    peak=max(profile['logical_qubits_upper'],2*profile['reflection_core_qubits'])+m
    assert charge['peak_logical_qubits']==peak
    expected={'T':per['T']+per['Tdg'],'measurements':peak,'feedback_CZ':0,'reset_X':peak,
        'expanded_primitive_operations':sum(per.values())+2*peak,
        'uniform_H':profile['uniform_H']+(M-1)*2*profile['uniform_H']+m,
        'explicit_T_depth_serial_upper':profile['explicit_T_depth_serial_upper']+(M-1)*(2*profile['explicit_T_depth_serial_upper']+4*ccx)+qft['T_count'],
        'conservative_5layer_T_depth_serial_upper':profile['conservative_5layer_T_depth_serial_upper']+(M-1)*(2*profile['conservative_5layer_T_depth_serial_upper']+5*ccx)+qft['T_count']}
    expected.update({f'gate_{op}':n for op,n in per.items()})
    for op in ('H','S','Sdg','T','Tdg','CX'):expected.setdefault(f'gate_{op}',0)
    assert expected==charge['per_repetition_metrics']
    assert {k:R*v for k,v in expected.items()}==charge['metrics']
    assert charge['controlled_Q_calls']==R*(M-1)
    assert charge['complete_reinitialisations']==R and charge['phase_register_measurements']==R*m
    assert charge['complete_oracle_T']==ap['T']+ap['Tdg'] and charge['controlled_Q_T']==qp['T']+qp['Tdg']
    return charge['metrics']['T']

def main():
    path=RESULTS/'four_survivor_quantum_resources.json';data=json.loads(path.read_text())
    protocol=json.loads((RESULTS/data['protocol_file']).read_text())
    assert sha(RESULTS/data['protocol_file'])==data['protocol_sha256']
    assert sha(ROOT/'four_survivor_quantum_resources.py')==data['source_sha256']==protocol['source_sha256']
    for n,h in data['dependency_source_sha256'].items():assert sha(ROOT/n)==h
    for n,h in data['input_sha256'].items():assert sha(RESULTS/n)==h
    scalar_path=RESULTS/data['scalar_ordinate_archive'];assert sha(scalar_path)==data['scalar_ordinate_archive_sha256']
    known_tables=set();cert_count=0;node_count=0;reconstructible={}
    widths={profile['integer_word_width'] for row in data['records'] for family in row['profiles'].values() for profile in family.values()}
    with gzip.open(scalar_path,'rt') as stream:
        for line in stream:
            row=json.loads(line);c=row['certificate'];assert digest(c)==row['certificate_sha256'];cert_count+=1
            if c.get('schema')=='qbp-boundary-scalar-v1':
                values=[n['integer'] for n in c['nodes']];width=row['graph_word_width'];packed=[]
                for j,v in enumerate(values):
                    delta=values[j+1]-v if j+1<len(values) else 0
                    packed.append(v+(abs(delta)<<width)+(int(delta<0)<<(2*width)))
                maximum=F(0)
                for n in c['nodes']:
                    y=F(n['integer'],1<<row['q']);err=max(abs(y-F(n['lower'])),abs(y-F(n['upper'])))
                    assert err<=F(n['error_upper']);maximum=max(maximum,err);node_count+=1
                assert maximum<=F(c['nodal_error_upper'])
                assert F(c['output_error_upper'])==F(c['nodal_error_upper'])+F(c['linear_interpolation_error_upper'])+F(c['interpolation_floor_error_upper'])
                assert F(c['linear_interpolation_error_upper'])==F(c['second_derivative_upper'])*(F(c['span'])/(1<<c['table_bits']))**2/8
                # Width affects only the packed value/difference/sign fields.
                # The exact q28 ordinates and their certificates are unchanged.
                for word_width in widths:
                    words=[]
                    for j,v in enumerate(values):
                        delta=values[j+1]-v if j+1<len(values) else 0
                        assert max(v,abs(delta))<1<<word_width
                        words.append(v+(abs(delta)<<word_width)+(int(delta<0)<<(2*word_width)))
                    word_sha=hashlib.sha256(json.dumps(words,separators=(',',':')).encode()).hexdigest()
                    reconstructible[word_sha]={'certificate_sha256':row['certificate_sha256'],
                        'graph_word_width':word_width,'q':row['q'],'packed_table':words,'packed_table_sha256':word_sha}
            else:packed=[n['offset_integer'] for n in c['nodes']];node_count+=len(packed)
            assert packed==row['packed_table'] and digest(packed)==row['packed_table_sha256']
            # The production emitter hashes packed arrays without sorted-key JSON.
            known_tables.add(hashlib.sha256(json.dumps(packed,separators=(',',':')).encode()).hexdigest())
    wanted={b['table_sha256'] for row in data['records'] for family in row['profiles'].values()
            for profile in family.values() for b in profile['blocks'] if b['table_sha256'] is not None}
    missing=wanted-known_tables
    assert missing<=set(reconstructible),'A profile table cannot be recovered from retained ordinates'
    supplemental={'schema':'qbp-width-specific-packed-table-reconstruction-v1',
        'status':'passed','resource_file_sha256':sha(path),'original_scalar_archive_sha256':sha(scalar_path),
        'tables':[reconstructible[h] for h in sorted(missing)],
        'scope':'Exact integer repacking of retained synthetic scalar ordinates into the actual profile word width. Original certificates, archive and execution records remain unchanged.'}
    supplemental_path=RESULTS/'four_survivor_quantum_table_width_reconstruction.json.gz'
    with supplemental_path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as stream:
            stream.write(json.dumps(supplemental,separators=(',',':')).encode())
    known_tables.update(missing)
    preparations={p['model']['model_id']:p for p in data['preparations']};direct_count=0
    for p in preparations.values():
        model=p['model'];weights=list(map(F,model['weights']));loadings=[list(map(F,r)) for r in model['loadings']]
        rg=[sum(w*r[j] for w,r in zip(weights,loadings)) for j in range(model['outer_dimension'])]
        direct=p['R1_direct_global_financial_transfer'];norm=p['payoff_norms'];M1=F(norm['M1_upper']);M2=F(norm['M2_upper'])
        linear=quadratic=F(0)
        for j,row in enumerate(direct['coordinates']):
            A=max(abs(r[j]) for r in loadings);D1=M1*(A+abs(rg[j]));D2=(M1+M2)*(A*A+rg[j]*rg[j])
            assert F(row['full_asset_loading_max'])==A and F(row['geometric_loading'])==rg[j]
            assert F(row['D1_upper'])==D1 and F(row['D2_upper'])==D2
            l=F(direct['other_coordinate_mass_upper'])*F(direct['rho1_mass_upper'])*D1
            q=F(direct['other_coordinate_mass_upper'])*F(direct['rho2_mass_upper'])*D2
            assert l==F(row['linear_error_upper']) and q==F(row['quadratic_error_upper'])
            linear+=l;quadratic+=q;direct_count+=1
        assert linear==F(direct['linear_transport_upper']) and quadratic==F(direct['quadratic_transport_upper'])
        assert F(direct['financial_transfer_upper'])==linear+quadratic+F(direct['tail_upper'])
        assert F(direct['tail_upper'])>=1-hi(direct['mass_interval'])
        best=min(F(p['R1_older_financial_transfer']['financial_transfer_upper']),F(direct['financial_transfer_upper']))
        assert best==F(p['R1_selected_financial_transfer']['financial_transfer_upper'])
    tapes=json.loads((RESULTS/'qft_rotation_tapes.json').read_text());tc=json.loads((RESULTS/'qft_rotation_certificates.json').read_text())
    assert sha(RESULTS/'qft_rotation_tapes.json')==tc['tapes_sha256']
    bounds={r['key']:F(r['operator_norm_upper_rational']) for r in tc['certificates']}
    profiles=plans=budgets=complete=0;rejections=[]
    for row in data['records']:
        p=preparations[row['model_id']];certificate=row['base_graph_certificate'];budget=row['complete_deterministic_budget']
        assert certificate['precision_search_certificate_sha256']==digest({k:v for k,v in certificate.items() if k!='precision_search_certificate_sha256'})
        if budget:
            mass=budget['mass_interval'];base=p['coefficient_certificate']['c0' if row['kind']=='positive_r1' else 'coefficient']
            cr=(hi(base)-lo(base))/2;mh=hi(mass);mr=(hi(mass)-lo(mass))/2
            cfg=certificate['configuration'];B=F(cfg['monetary_signal_normaliser']);decode=F(cfg['monetary_probability_decode_factor'])
            arithmetic=decode*B*F(certificate['probability_error_upper'])
            for bud in [budget]+([row['older_financial_deterministic_budget']] if 'older_financial_deterministic_budget' in row else []):
                total=cr+F(bud['financial_transfer_upper'])+mh*(F(bud['auxiliary_components']['bias_upper'])+arithmetic)+mr*B
                assert total==F(bud['total_deterministic_undiscounted_price_error_upper'])
                assert cr==F(bud['coefficient_radius_upper']) and arithmetic==F(bud['arithmetic_before_mass_upper']);budgets+=1
        if 'baseline_failure' in row:
            assert budget is None or F(budget['total_deterministic_undiscounted_price_error_upper'])>=EPS
            rejections.append({'model_id':row['model_id'],'kind':row['kind'],'label':row['candidate_label']});continue
        for family in row['profiles'].values():
            for profile in family.values():verify_profile(profile,known_tables);profiles+=1
        for strategy in row['strategies']:
            result=strategy['result']
            if result['kind']=='rejected':continue
            if result['kind']=='deterministic_correction_prior':assert F(result['total_price_error_upper'])<=EPS;continue
            expected_total=F(0);T=0;failure_total=F(0)
            for flag in result['flags']:
                profile=row['profiles'][strategy['architecture']][flag['flag']]
                T+=verify_plan(flag,profile,tapes,bounds);plans+=1
                expected_total+=F(flag['plan']['amplitude_error_upper'])+F(flag['plan']['postprocess_allowance'])
                failure_total+=F(flag['plan']['binomial_failure_upper'])
                assert F(flag['plan']['alpha'])==(ALPHA/2 if result['encoding']=='split' else ALPHA)
            decoder=F(2) if result['encoding']=='centred' and not result['correction_prior']['native_positive'] else F(1)
            statistical=F(result['monetary_signal_factor_upper'])*decoder*expected_total
            assert statistical==F(result['statistical_price_error_upper'])
            assert F(result['deterministic_error_upper'])+statistical==F(result['total_price_error_upper'])<=EPS
            assert failure_total==F(result['failure_probability_upper'])<=ALPHA
            assert T==result['T_count'];complete+=1
        choices=[s for s in row['strategies'] if s['result']['kind']!='rejected']
        selected=min(choices,key=lambda x:(x['result']['T_count'],x['architecture'],x['strategy'],x['prior_information']))
        assert selected==row['selected_strategy']
    assert len(data['records'])==17 and len(data['selected_prices'])==4
    assert profiles==data['profiles_count']==80 and plans==data['confidence_plan_count']
    assert len(rejections)==1 and rejections[0]['model_id']=='n8_lambda0.2_scale3' and rejections[0]['kind']=='boundary_r2'
    for case in data['selected_prices']:
        assert len(case['methods'])==4 and all(case['methods'].values())
        values={k:v['selected']['result']['T_count'] for k,v in case['methods'].items()}
        assert case['quantum_resource_winner']==min(values,key=values.get)
    result={'schema':'qbp-independent-four-survivor-resources-audit-v1','status':'passed',
        'resource_file_sha256':sha(path),'producer_sha256':sha(ROOT/'four_survivor_quantum_resources.py'),
        'audit_source_sha256':sha(__file__),'protocol_sha256':data['protocol_sha256'],
        'checked_profiles':profiles,'checked_confidence_plans':plans,'checked_complete_strategies':complete,
        'checked_deterministic_budgets':budgets,'checked_direct_R1_coordinate_bounds':direct_count,
        'checked_scalar_certificates':cert_count,'checked_scalar_and_quantile_nodes':node_count,
        'width_specific_table_reconstruction_file':supplemental_path.name,
        'width_specific_table_reconstruction_sha256':sha(supplemental_path),
        'reconstructed_width_specific_packed_tables':len(missing),
        'retained_rejections':rejections,
        'independence':'Producer, compiler resource and estimator/confidence APIs are not imported. Exact rational formulas, literal retained rotation strings, separate CCX reflection expansion and table-word reconstruction are used.',
        'scope':'Complete ideal-logical counts and finite confidence, synthetic financial inputs and retained arithmetic certificates. No physical hardware execution, wall-clock speedup, novelty or global optimality claim.'}
    (RESULTS/'four_survivor_quantum_independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
