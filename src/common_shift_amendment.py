"""Sequential common-shift challenge using frozen financial preparations.

All circuit/scalar/profile artefacts are hash references. Registration precedes
every amended budget and cost. Previous sufficient certificates and complete
plans are retained, including failures and zero-shot analytical decisions.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import argparse, copy, gzip, hashlib, importlib.util, json, sys, time
sys.dont_write_bytecode = True

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(path):
    path=Path(path)
    return json.load(gzip.open(path,'rt')) if path.suffix=='.gz' else json.loads(path.read_text())
def write(path,value):
    path=Path(path)
    if path.exists(): raise FileExistsError('Preserve previous evidence: '+str(path))
    if path.suffix=='.gz':
        with path.open('wb') as handle:
            with gzip.GzipFile(fileobj=handle,mode='wb',mtime=0) as stream:
                encoder=json.JSONEncoder(sort_keys=True,separators=(',',':'),default=str)
                for block in encoder.iterencode(value): stream.write(block.encode())
                stream.write(b'\n')
    else: path.write_text(json.dumps(value,indent=2,default=str)+'\n')

RESOURCE_PATHS=(
 'raw_payoff_tight/raw_payoff_resources.json.gz',
 'raw_payoff_smoothing/raw_payoff_smoothing_resources.json',
 'expanded_raw_payoff/expanded_raw_payoff_resources.json.gz',
 'asymmetric_raw_payoff/asymmetric_raw_payoff_resources.json.gz',
 'strengthened_r1_resources.json','shared_threshold_resources.json',
 'factor_basis_quantum_resources.json')
ARCHIVE_PATHS=(
 'raw_payoff_tight/raw_payoff_scalar_ordinates.json.gz',
 'expanded_raw_payoff/expanded_raw_payoff_scalar_ordinates.json.gz',
 'asymmetric_raw_payoff/asymmetric_raw_payoff_scalar_ordinates.json.gz',
 'strengthened_r1_scalar_ordinates.json.gz','shared_threshold_scalar_ordinates.json.gz')
DEPENDENCIES=(
 'boundary_value_oracle.py','boundary_value_certificate.py',
 'integral_value_oracle.py','integral_value_certificate.py',
 'correction_information_priors.py','boundary_precision_selection.py',
 'integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py',
 'checks/audit_boundary_precision.py','check_expanded_precision_resources.py')
INPUTS=('gaussian_derivative_norms.json','qft_rotation_tapes.json','qft_rotation_certificates.json')

def bindings(args):
    repo=Path(args.repo)
    return {'producer_source_sha256':sha(__file__),
      'checked_common_shift_proof_sha256':sha(args.proof),
      'resources_sha256':{p:sha(repo/'results'/p) for p in RESOURCE_PATHS},
      'scalar_archives_sha256':{p:sha(repo/'results'/p) for p in ARCHIVE_PATHS},
      'dependencies_sha256':{p:sha(repo/p) for p in DEPENDENCIES},
      'inputs_sha256':{p:sha(repo/'results'/p) for p in INPUTS}}

def load_sources(repo):
    sources={p:read(repo/'results'/p) for p in RESOURCE_PATHS}
    pool=[]
    old=sources[RESOURCE_PATHS[0]]['records']
    old_amended=sources[RESOURCE_PATHS[1]]['records']
    for i,(base,amended) in enumerate(zip(old,old_amended)):
        pool.append(('raw_payoff',RESOURCE_PATHS[0],i,base,RESOURCE_PATHS[1],i,amended))
    for p in (RESOURCE_PATHS[2],RESOURCE_PATHS[3]):
        for i,(base,amended) in enumerate(zip(sources[p]['records'],sources[p]['amended_records'])):
            pool.append(('raw_payoff',p,i,base,p,i,amended))
    assert len(pool)==348
    for i,base in enumerate(sources['strengthened_r1_resources.json']['records']):
        pool.append(('numerical_conditional_price','strengthened_r1_resources.json',i,base,None,None,None))
    for i,base in enumerate(sources['shared_threshold_resources.json']['rows']):
        pool.append(('shared_threshold','shared_threshold_resources.json',i,base,None,None,None))
    return sources,pool

def record_binding(entry):
    method,path,index,base,amended_path,amended_index,amended=entry
    b={'method':method,'model_id':base['original_model_id'],'candidate_id':base['candidate_id'],
       'source_path':path,'source_row_index':index,'source_row_sha256':digest(base),
       'configuration_sha256':digest(base['graph_configuration']),
       'profiles_sha256':digest(base['profiles']),
       'retained_prior_sha256':digest(base.get('coefficient_informed_prior')),
       'retained_strategies_sha256':digest(base['strategies'])}
    if method=='raw_payoff':
        bg=base.get('PG_exp_bits',base.get('payoff_exp_bits'))
        ba=base.get('PA_exp_bits',base.get('payoff_exp_bits'))
        b.update(normal_bits=base['normal_bits'],PG_exp_bits=bg,PA_exp_bits=ba,
         graph_sha256=digest(base['graph_nodes']),amended_source_path=amended_path,
         amended_source_row_index=amended_index,amended_source_row_sha256=digest(amended),
         retained_amended_strategies_sha256=digest(amended['strategies']))
    elif method=='numerical_conditional_price':
        b.update(PA_bits=base['PA_bits'],PG_bits=base['PG_bits'],
           retained_body_certificate_sha256=digest(base['base_graph_certificate']))
    else: b['graph_sha256']=digest(base['graph_nodes'])
    return b

def register(args):
    repo,out=Path(args.repo),Path(args.output);out.mkdir(parents=True,exist_ok=True)
    sources,pool=load_sources(repo)
    known=sources['asymmetric_raw_payoff/asymmetric_raw_payoff_resources.json.gz']['selected_comparisons']
    protocol={'schema':'qbp-common-shift-protocol-v1',
      'registration_utc':datetime.now(timezone.utc).isoformat(),
      'research_record_date':'2026-10-11',**bindings(args),
      'registered_before_amended_budget_prior_or_cost_outcomes':True,
      'known_predecessors':{'disclosure':'All previous raw348, numerical64 and shared4 outcomes and the prospective common-shift derivation are known. This is a disclosed sequential certificate challenge, not a blind initial registration.',
        'latest_raw_comparisons':known,
        'numerical_comparisons':sources['strengthened_r1_resources.json']['selected_comparisons']},
      'fixed_pool':list(map(record_binding,pool)),
      'configurations':{'raw_payoff':348,'numerical_conditional_price':64,'shared_threshold':4},
      'architectures':['exact7t','unitary8t'],
      'common_shift_coefficient':'c0=nu0/lambda; Delta*=retained logT ideal upper=U_star; C=min(c0,nu1*Delta*/lambda^2). True gap is log(1+U)<=U_star.',
      'raw_formula':'eta=min(previous_pointwise,previous_smoothing,1,100(e_exp_G+e_exp_A)+101E_s+2D_N+C*(e_m+e_C+e_O)+c0*e_Delta). E_s=2L*abs(lambda_hat-lambda)+2^-q*indicator(drop>0). The shared log perturbation acts on h_Delta, whose derivative is a difference of two values in[0,101].',
      'numerical_formula':'For all retained B1=H=1 graphs: eta=min(previous_graph,1,e_PA_intrinsic+e_PG_intrinsic+C*e_m+c0*(e_Delta+e_(BA-BG))). Both decoder scales are exact1. Assert unclipped actual mG/mA offsets are within their scalar domains; all prior guards remain necessary.',
      'threshold_formula':'eta=min(previous_sandwich,1,beta_S+2D_N+2c0E_s+C*(e_m+e_C+e_O+e_strike_log)+c0*e_Delta). Full strike-log body error includes strike-ratio arithmetic. Normal rounding occurs only through D_N.',
      'prior_and_budget':'Retain coefficient interval, outer mass and continuous financial transfer exactly. Enlarge the same price-independent finite ideal R1 mean interval by each new eta. Rebuild canonical QPE allowance, exact amplification, Fourier and complete clean-circuit/reflection/reset charges from exact frozen profiles.',
      'retention':'Select minimum over every previous valid architecture/pointwise/smoothing plan and every amended plan. Earlier exclusions are retained. Every production failure and zero-shot decision is retained. Do not alter pool, graph, arithmetic/scalar precision, compiler or estimator after outcomes.',
      'missing_profile_policy':'Twelve numerical PA6 rows were previously budget-rejected and have no frozen profiles. Compute their amended budgets too. If newly admitted, report an uncosted admission and fail the complete-comparison status; do not silently exclude, infer free cost or emit a new graph.',
      'target':'Monetary radius1/100; alpha1/3600 for each price; same original risk-neutral model and analytic coefficient. No realised price supplied.',
      'scope':'Sufficient certificates and implemented finite-pool canonical-QPE logical resource caps. No new quantum primitive, theoretical priority, global optimum, physical runtime or practical quantum advantage.',
      'data_origin':'Exact rational synthetic models and certified scalar tables; no licensed/private observations or purchased data.'}
    write(out/'common_shift_protocol.json',protocol)
    print('Registered raw348/numerical64/shared4 before amended budgets or costs',flush=True)

def scalar_index(repo):
    index={}
    for p in ARCHIVE_PATHS:
        for a in read(repo/'results'/p):
            c=a['certificate']
            if c.get('name') is None:continue
            key=(c['name'],F(c['lo']),F(c['span']),tuple(sorted(c['params'].items())),c['table_bits'],c['q'])
            if key in index:assert index[key]==c
            index[key]=c
    return index

def coefficient(x,q):
    x=F(x)
    if x.denominator&(x.denominator-1)==0:return x.numerator,x.denominator.bit_length()-1
    t=x*(1<<q)+F(1,2);return t.numerator//t.denominator,q

def constants(row,nu):
    cfg=row['graph_configuration'];lam=F(cfg['lambda'])
    c0=nu[0]/lam;delta=F(cfg['U_star']);C=min(c0,nu[1]*delta/lam**2)
    return lam,c0,delta,C

def raw_mapping(row,archive,nu):
    graph={n['name']:n for n in row['graph_nodes']};sem=row['graph_semantic_nodes']
    cfg=row['graph_configuration'];cert=row['probability_certificate']
    body={n['node']:n for n in cert['body_arithmetic_certificate']['nodes']}
    assert cert['status']=='passed' and cert['body_arithmetic_certificate']['status']=='passed'
    err=lambda n:F(body[n]['absolute_output_error_upper'])
    lower=graph[sem['lower_endpoint']];assert lower['kind']=='add'
    inner=graph[lower['args'][1]];outer=graph[lower['args'][0]]
    assert inner['kind']=='scale' and inner['args']==[sem['normal_plus_L']]
    assert outer['kind']=='add' and outer['args'][0]==sem['mG_offset']
    cn=outer['args'][1];on=sem['common_offset'];assert graph[cn]['kind']==graph[on]['kind']=='constant'
    lam,c0,delta,C=constants(row,nu);q=cfg['q'];L=F(cfg['inner_normal_L'])
    assert F(inner['meta'])==lam
    I,d=coefficient(lam,q);Es=2*L*abs(F(I,1<<d)-lam)+(F(1,1<<q) if d else 0)
    exp=[]
    for side in ('PG','PA'):
        n=graph[sem[side+'_exp']];m=n['meta'];assert m['name']=='exp'
        cc=archive[m['name'],F(m['lo']),F(m['span']),tuple(sorted(m['parameters'].items())),n['bits'],q]
        exp.append(F(cc['output_error_upper']))
    Em=err(sem['mG_offset'])+err(cn)+err(on);edelta=err(sem['logT']);D=F(cert['inner_normal_CDF_error_upper'])
    arithmetic=100*sum(exp)+101*Es+C*Em+c0*edelta
    return {'c0':str(c0),'Delta_star':str(delta),'common_shift_coefficient':str(C),
      'Em':str(Em),'e_Delta':str(edelta),'E_s':str(Es),'D_N':str(D),
      'exp_errors':list(map(str,exp)),'common_shift_error':str(C*Em),
      'gap_error':str(c0*edelta),'arithmetic_upper':str(arithmetic),
      'auxiliary_bias_upper':str(2*D),'new_probability_error_upper':str(min(F(1),arithmetic+2*D))}

def numeric_graphs(repo,sources):
    sys.path.insert(0,str(repo))
    import boundary_value_oracle as o
    result={}
    factor=sources['factor_basis_quantum_resources.json']
    for r in sources['strengthened_r1_resources.json']['records']:
        key=r['original_model_id'],r['candidate_id']
        if key in result:continue
        prep=next(x for x in factor['preparations'] if (x['original_model_id'],x['candidate_id'])==key)
        g=o.build_graph(prep['model'],F(1),q=28,b_y=10,table_bits=12,outer_bound=6,
         U_star=F(prep['domain']['U_star']),retained_bound=max(abs(F(v)) for v in prep['domain']['coordinate_bounds']),kind='positive_r1')
        result[key]=g
    return result

def numerical_mapping(row,graph,archive,nu):
    cert=row['base_graph_certificate'];sem=row['graph_semantic_nodes'];cfg=row['graph_configuration']
    assert cert['status']=='passed' and not cert['failures']
    assert graph.semantic_nodes==sem and list(graph.nodes)==[r['node'] for r in cert['nodes']]
    assert F(cfg['monetary_signal_normaliser'])==F(cfg['cap'])==1
    body={n['node']:n for n in cert['nodes']};err=lambda n:F(body[n]['absolute_output_error_upper'])
    lam,c0,delta,C=constants(row,nu);mA=graph.nodes[sem['mA_offset']]
    assert mA.kind=='add' and mA.args[1]==sem['logT']
    head=graph.nodes[mA.args[0]];assert head.kind=='add' and head.args[0]==sem['mG_offset']
    cn=head.args[1];assert graph.nodes[cn].kind=='constant'
    exp=[]
    for side in ('PA','PG'):
        node=graph.nodes[sem[side]];assert node.kind=='function' and node.meta.name=='cap'
        raw=sem['mA_offset'] if side=='PA' else sem['mG_offset']
        clo,chi=map(F,body[raw]['actual_bounds'])
        assert 0<=clo<=chi<=node.meta.span
        assert graph.nodes[node.args[0]].kind=='clamp' and graph.nodes[node.args[0]].args==(raw,)
        cc=archive[node.meta.name,node.meta.lo,node.meta.span,
         tuple(sorted((k,str(v)) for k,v in node.meta.parameters)),row[side+'_bits'],cfg['q']]
        exp.append(F(cc['output_error_upper']))
    output=graph.nodes[sem['centred_probability']]
    while output.kind in ('clamp','clamp_signed'): output=graph.nodes[output.args[0]]
    assert output.kind=='signed_sub'
    for arg,side in zip(output.args,('PA','PG')):
        n=graph.nodes[arg];assert n.kind=='scale' and n.meta==1 and n.args==(sem[side],)
    Em=err(sem['mG_offset']);edelta=err(sem['logT'])+err(cn)
    ar=sum(exp)+C*Em+c0*edelta
    return {'c0':str(c0),'Delta_star':str(delta),'common_shift_coefficient':str(C),
      'Em':str(Em),'e_Delta_including_offset':str(edelta),'intrinsic_cap_errors':list(map(str,exp)),
      'decoder_scale_rounding_upper':'0','inactive_input_clamps_checked':True,
      'common_shift_error':str(C*Em),'gap_error':str(c0*edelta),
      'arithmetic_upper':str(ar),'auxiliary_bias_upper':'0','new_probability_error_upper':str(min(F(1),ar))}

def threshold_mapping(row,nu):
    graph={n['name']:n for n in row['graph_nodes']};sem=row['graph_semantic_nodes'];cfg=row['graph_configuration'];cert=row['probability_certificate']
    assert cert['status']=='passed' and cert['body_arithmetic_certificate']['status']=='passed'
    body={n['node']:n for n in cert['body_arithmetic_certificate']['nodes']};err=lambda n:F(body[n]['absolute_output_error_upper'])
    lower=graph[sem['lower_endpoint']];upper=graph[sem['upper_endpoint']]
    assert lower['kind']=='add' and lower['args'][0]==sem['mG_offset']
    assert upper['kind']=='add' and upper['args']==[sem['lower_endpoint'],sem['logT']]
    cn=lower['args'][1];threshold=graph[sem['threshold']];assert threshold['kind']=='sub'
    thresholdadd=graph[threshold['args'][0]];assert thresholdadd['kind']=='add' and thresholdadd['args'][0]==sem['strike_log']
    on=thresholdadd['args'][1];assert graph[cn]['kind']==graph[on]['kind']=='constant'
    inner=graph[threshold['args'][1]];assert inner['kind']=='scale' and inner['args']==[sem['normal_plus_L']]
    lam,c0,delta,C=constants(row,nu);assert F(inner['meta'])==lam
    I,d=coefficient(lam,cfg['q']);Es=2*F(cfg['inner_normal_L'])*abs(F(I,1<<d)-lam)+(F(1,1<<cfg['q']) if d else 0)
    common=err(sem['mG_offset'])+err(cn)+err(on)+err(sem['strike_log']);edelta=err(sem['logT'])
    D=F(cert['inner_normal_CDF_error_upper']);strike=F(cert['strike_midpoint_bias_upper'])
    ar=2*c0*Es+C*common+c0*edelta
    return {'c0':str(c0),'Delta_star':str(delta),'common_shift_coefficient':str(C),
      'common_error_including_strike':str(common),'e_Delta':str(edelta),'E_s':str(Es),
      'D_N':str(D),'strike_midpoint_bias_upper':str(strike),
      'common_shift_error':str(C*common),'gap_error':str(c0*edelta),
      'arithmetic_upper':str(ar),'auxiliary_bias_upper':str(2*D+strike),
      'new_probability_error_upper':str(min(F(1),ar+2*D+strike))}

def make_amended(entry,term,frontend):
    binding=record_binding(entry);method,path,index,base,ap,ai,oldamended=entry
    oldcert=base['base_graph_certificate'] if method=='numerical_conditional_price' else base['probability_certificate']
    previous=[F(oldcert['probability_error_upper'])]
    if oldamended:previous.append(F(oldamended['probability_certificate']['probability_error_upper']))
    eta=min(*previous,F(term['new_probability_error_upper']))
    budget=copy.deepcopy(base['complete_deterministic_budget']);A=F(budget['external_factor_upper'])
    det=F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*eta+F(budget['external_factor_radius'])
    selected_arithmetic=max(F(0),eta-F(term['auxiliary_bias_upper']))
    budget.update(arithmetic_before_mass_upper=str(selected_arithmetic),
      auxiliary_components={'bias_upper':term['auxiliary_bias_upper']},
      total_deterministic_undiscounted_price_error_upper=str(det),certificate_amendment='common_shift_coupled_conditional')
    oldprior=base.get('coefficient_informed_prior')
    assert oldprior is not None
    lo,hi=F(oldprior['finite_ideal_mean_lower']),F(oldprior['finite_ideal_mean_upper'])
    prior={'probability_lower':str(max(F(0),lo-eta)),'probability_upper':str(min(F(1),hi+eta)),
      'native_positive':True,'finite_ideal_mean_lower':str(lo),'finite_ideal_mean_upper':str(hi),
      'previous_mean_information_sha256':oldprior['previous_mean_information_sha256'],
      'retained_prior_sha256':digest(oldprior),
      'information_scope':'Same price-independent finite ideal R1 mean interval enlarged by minimum prior and checked common-shift errors.'}
    cert={'status':'passed','configuration':base['graph_configuration'],
      'probability_error_upper':str(eta),'monetary_value_error_upper':str(eta),
      'probability_actual_bounds':oldcert['probability_actual_bounds'],
      'retained_body_certificate_sha256':digest(oldcert),'common_shift_derivation':term,
      'previous_probability_error_bounds':list(map(str,previous)),
      'scope':'Uniform conditional mean certificate for exact unchanged preparation; previous integer range/domain/body guards retained.'}
    row={**binding,'probability_certificate':cert,'complete_deterministic_budget':budget,
      'coefficient_informed_prior':prior,'strategies':[]}
    if det>=F(1,100):row['rejection']={'reason':'Amended deterministic monetary tolerance exhausted','deterministic_error_upper':str(det)}
    elif not base['profiles']:row['rejection']={'reason':'Newly admitted row has no frozen emitted profiles','uncosted_admission':True}
    else:
        for mode in ('exact7t','unitary8t'):
            try: result=frontend.informed_complete_plan('positive',budget,cert,prior,base['profiles'][mode])
            except (ArithmeticError,AssertionError,ValueError) as e:
                result={'kind':'production_failure','reason':type(e).__name__,'message':str(e)}
            row['strategies'].append({'architecture':mode,'result':result})
    prior_strategies=[('previous_pointwise',s) for s in base['strategies']]
    if oldamended:prior_strategies.extend(('previous_smoothing',s) for s in oldamended['strategies'])
    current=[('common_shift',s) for s in row['strategies']]
    valid=[{'certificate_policy':tag,**s} for tag,s in prior_strategies+current
      if s['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')]
    if valid:row['selected_strategy']=min(valid,key=lambda s:(s['result']['T_count'],s['architecture'],s['certificate_policy']))
    row['minimum_previous_probability_error_upper']=str(min(previous))
    row['common_shift_selected_probability_error_upper']=str(eta)
    return row

def summary(rows):
    comparisons=[]
    for mid in dict.fromkeys(r['model_id'] for r in rows):
        methods={}
        for method in ('raw_payoff','numerical_conditional_price','shared_threshold'):
            group=[r for r in rows if r['model_id']==mid and r['method']==method]
            valid=[r for r in group if r.get('selected_strategy')]
            if not valid:raise AssertionError('No complete strategy: '+mid+' '+method)
            win=min(valid,key=lambda r:(r['selected_strategy']['result']['T_count'],
             r.get('normal_bits',0),r.get('PG_exp_bits',r.get('PG_bits',0)),r.get('PA_exp_bits',r.get('PA_bits',0)),r['selected_strategy']['architecture']))
            methods[method]={'selected_record_sha256':digest(win),
              'source_path':win['source_path'],'source_row_index':win['source_row_index'],
              'source_row_sha256':win['source_row_sha256'],'profiles_sha256':win['profiles_sha256'],
              'selected_strategy':win['selected_strategy'],'T_count':win['selected_strategy']['result']['T_count'],
              'amended_probability_error_upper':win['common_shift_selected_probability_error_upper'],
              'configurations':len(group),'admitted_complete_strategies':len(valid),
              'rejected_configurations':len(group)-len(valid),
              'configuration':{k:win[k] for k in ('normal_bits','PG_exp_bits','PA_exp_bits','PG_bits','PA_bits') if k in win}}
        rivals={k:methods[k]['T_count'] for k in ('raw_payoff','numerical_conditional_price')}
        strongest=min(rivals,key=lambda k:(rivals[k],k));T=rivals[strongest];E=methods['shared_threshold']['T_count']
        comparisons.append({'model_id':mid,'methods':methods,'strongest_rival':strongest,
         'strongest_rival_T_count':T,'shared_threshold_T_count':E,
         'saving_against_strongest_exact':str(F(T-E,T)) if T else None,
         'saving_against_strongest_per_cent_diagnostic':float(100*F(T-E,T)) if T else None})
    return comparisons

def run(args):
    repo,out=Path(args.repo),Path(args.output);protocol=read(out/'common_shift_protocol.json')
    frozen=bindings(args);assert all(protocol[k]==v for k,v in frozen.items())
    sources,pool=load_sources(repo);assert list(map(record_binding,pool))==protocol['fixed_pool']
    archive=scalar_index(repo);norms=read(repo/'results/gaussian_derivative_norms.json')
    nu={r:F(next(x['norm_upper'] for x in norms['unit_norms'] if x['unit_derivative_order']==r and x['passed'])) for r in (0,1)}
    graphs=numeric_graphs(repo,sources);sys.path.insert(0,str(repo));import correction_information_priors as frontend
    rows=[];start=time.monotonic()
    for entry in pool:
        method,path,index,base,ap,ai,oldamended=entry
        if method=='raw_payoff':term=raw_mapping(base,archive,nu)
        elif method=='numerical_conditional_price':term=numerical_mapping(base,graphs[base['original_model_id'],base['candidate_id']],archive,nu)
        else:term=threshold_mapping(base,nu)
        rows.append(make_amended(entry,term,frontend))
    comparisons=summary(rows)
    assert all(bindings(args)[k]==v for k,v in frozen.items())
    uncosted=[r for r in rows if r.get('rejection',{}).get('uncosted_admission')]
    outcome={'schema':'qbp-common-shift-resources-v1','status':'failed complete-comparison scope: uncosted admissions' if uncosted else 'completed retaining every amended outcome',
      'protocol_sha256':sha(out/'common_shift_protocol.json'),**frozen,
      'configurations':protocol['configurations'],'records':rows,'selected_comparisons':comparisons,
      'uncosted_admissions':len(uncosted),'generation_seconds_diagnostic':time.monotonic()-start,
      'scope':protocol['scope'],'data_origin':protocol['data_origin']}
    write(out/'common_shift_resources.json.gz',outcome)
    for r in comparisons:print(r['model_id'],{k:v['T_count'] for k,v in r['methods'].items()},r['saving_against_strongest_per_cent_diagnostic'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run'))
    p.add_argument('--repo',required=True);p.add_argument('--proof',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();register(a) if a.mode=='register' else run(a)
