"""Matched finite precision search and fully charged residual pricing plans.

All three financial methods receive the same finite outer law, monetary
accuracy and confidence. Local table search is a declared greedy finite pool,
not a global optimum. The actual circuits and complete estimator caps are
computed only from certified candidate programs. All trials are retained.
"""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict
import gzip,hashlib,json,sys,time
import integral_value_oracle as oracle
import integral_value_certificate as certificate
import integral_auxiliary_quadrature as auxiliary
import integral_pricing_budget as estimator
ROOT=Path(__file__).resolve().parent


def load(name):return json.loads((ROOT/'results'/name).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def upper(iv):return F(iv['upper_exact_fraction'])
def lower(iv):return F(iv['lower_exact_fraction'])


def build(model,cap,kind,domain,bits=None,encoding='centred'):
    graph=oracle.build_graph(model,cap,q=28,b_y=10,beta_bits=8,strike_bits=8,table_bits=12,
        outer_bound=6,U_star=F(domain['U_star']),D2=certificate.joint_norm(2,F(model['lambda'])),
        D1=certificate.joint_norm(1,F(model['lambda'])),kind=kind,
        retained_bound=max(abs(F(x)) for x in domain['coordinate_bounds']),output_encoding=encoding)
    for name,b in (bits or {}).items():
        if name not in graph.nodes or graph.nodes[name].kind!='function' or not 6<=b<=12:raise ValueError('Invalid local scalar precision')
        graph.nodes[name].bits=b
    graph.configuration['function_address_bits']={n:v.bits for n,v in graph.nodes.items() if v.kind=='function'}
    return graph


def bias(model,cap,kind,transfer,cert,factor,moments):
    old=certificate.complete_bias_budget(model,cap,kind,transfer,cert,factor,8,8,moments)
    if kind=='integral_r2':
        new=auxiliary.auxiliary_bias(model,cap,8,8,F(cert['configuration']['U_star']),moments[2],moments[3])
        old_aux=F(old['auxiliary_before_factor_upper']);new_aux=F(new['total_auxiliary_bias_upper'])
        old['superseded_linear_auxiliary_budget']=old['auxiliary_components']
        old['total_deterministic_undiscounted_price_error_upper']=str(F(old['total_deterministic_undiscounted_price_error_upper'])+F(old['external_factor_upper'])*(new_aux-old_aux))
        old['auxiliary_before_factor_upper']=str(new_aux);old['auxiliary_components']=new
    return old


def zero_cost(reason,price_radius):
    return {'kind':'deterministic_correction_prior','T_count':0,'quantum_shots':0,
        'certified_correction_radius':str(price_radius),'reason':reason}


def complete_plan(kind,encoding,budget,cert,moments,positive_prior,profiles):
    det=F(budget['total_deterministic_undiscounted_price_error_upper']);eps=F(1,100);alpha=F(1,2700)
    if det>=eps:return {'kind':'rejected','reason':'Full deterministic monetary envelope exhausts tolerance','deterministic_error_upper':str(det)}
    C=F(budget['external_factor_upper'])*F(cert['configuration']['monetary_signal_normaliser']);remaining=eps-det
    if C<=0:return zero_cost('Identically zero retained financial correction',F(0))
    eta=F(cert['probability_error_upper']);U=F(cert['configuration']['U_star'])
    rare=certificate.rare_flag_prior(moments[2],U,eta) if kind!='positive_r1' else None
    prior=F(positive_prior['probability_prior_upper']) if kind=='positive_r1' else F(rare['positive_prior_upper'])
    # Sharing this same correction prior prevents a weak centred baseline.
    prior_radius=C*prior/(2 if kind=='positive_r1' else 1)
    if prior_radius<=remaining:return zero_cost('Equally supplied finite correction interval meets remaining monetary tolerance',prior_radius)
    plans=[]
    if kind=='positive_r1':
        labels=['centred'];target=remaining/C;priors=[prior];alphas=[alpha]
    elif encoding=='centred':
        labels=['centred'];target=remaining/(2*C);priors=[F(1)];alphas=[alpha]
    else:
        labels=['positive_part','negative_part'];target=remaining/(2*C);priors=[prior,prior];alphas=[alpha/2,alpha/2]
    charges=[]
    for label,u,a in zip(labels,priors,alphas):
        setup,qft=estimator.compiled_extended_plan(target,a,u)
        profile=profiles[label]
        charge=estimator.charge_full_oracle(setup,profile,qft)
        plans.append({'flag':label,'probability_tolerance':str(target),'prior_upper':str(u),'plan':asdict(setup),
            'financial_profile_sha256':hashlib.sha256(json.dumps(profile,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            'fourier_profile':qft,'complete_charge':charge})
        charges.append(charge)
    money_stat=C*sum((p['plan']['amplitude_error_upper']+p['plan']['postprocess_allowance'] for p in plans),F(0))
    if kind!='positive_r1' and encoding=='centred':money_stat*=2
    assert det+money_stat<=eps
    failure=sum((p['plan']['binomial_failure_upper'] for p in plans),F(0));assert failure<=alpha
    return {'kind':'complete_quantum_resource_plan','encoding':encoding,'T_count':sum(p.get('metrics',{}).get('T',0) for p in charges),
        'quantum_shots':sum(p['plan']['repetitions'] for p in plans),'deterministic_error_upper':str(det),'statistical_price_error_upper':str(money_stat),
        'total_price_error_upper':str(det+money_stat),'failure_probability_upper':str(failure),
        'momentary_signal_factor_upper':str(C),'rare_flag_prior':rare,'flags':plans,
        'logical_qubits_upper':max(p.get('peak_logical_qubits',0) for p in charges),
        'execution_scope':'Constructed ideal-logical resource/confidence plan; no physical QPE price experiment was executed.'}


def main():
    start=time.perf_counter();data=load('generalisation_protocol.json');models={m['model_id']:m for m in data['models']}
    domains={x['model_id']:x for x in load('integral_finite_input_domains.json')['domains']}
    finances=load('integral_financial_grid_assessment.json')['records']
    moment_records={x['model_id']:x for x in load('integral_finite_grid_moments.json')['records']}
    mass=certificate.quantile_certificate(F(6),10,28)['mass_interval'];factor=(lower(mass),upper(mass))
    block_cache={};records=[];trials=[];scalar_seen=set();scalar_digest_by_header={};scalar_archive_path=ROOT/'results/integral_precision_scalar_archive.jsonl.gz'
    # Mathematical ordinate certificates can be cached independently of word width.
    original_scalar=certificate.scalar_certificate;scalar_cache={}
    def memo_scalar(name,lo,span,params,q,bits,output_offset=0):
        key=(name,str(lo),str(span),tuple(sorted((str(k),str(v)) for k,v in params.items())),q,bits,str(output_offset))
        if key not in scalar_cache:scalar_cache[key]=original_scalar(name,lo,span,params,q,bits,output_offset)
        return scalar_cache[key]
    certificate.scalar_certificate=memo_scalar
    with scalar_archive_path.open('wb') as raw:
        stream=gzip.GzipFile(fileobj=raw,mode='wb',mtime=0)
        for financial in finances:
            model=models[financial['model_id']];cap=financial['cap'];mid=model['model_id'];domain=domains[mid];tables_cache={}
            moments={r:upper(moment_records[mid]['centred_U_moments'][str(r)]) for r in (1,2,3)}
            for kind in ['integral_r2','direct_r2','positive_r1']:
                transfer=F(financial['kinds']['R1' if kind=='positive_r1' else 'R2']['financial_transfer_including_tail_and_quantile_rounding_upper'])
                local_programs={};evaluation_cache={}
                def evaluate(bits):
                    graph=build(model,cap,kind,domain,bits)
                    key=tuple(sorted(graph.configuration['function_address_bits'].items()))
                    if key in evaluation_cache:return evaluation_cache[key]
                    tables,ordinates=oracle.make_tables(graph,tables_cache,True)
                    cert=certificate.certify_graph(graph,tables,ordinates)
                    b=bias(model,cap,kind,transfer,cert,factor,moments) if cert['status']=='passed' else None
                    for ordinate in ordinates.values():
                        header=(ordinate['name'],ordinate['lo'],ordinate['span'],tuple(sorted((str(k),str(v)) for k,v in ordinate['params'].items())),ordinate['q'],ordinate['table_bits'],ordinate['output_offset'],ordinate['precision_bits']) if ordinate['schema']=='qbp-integral-scalar-v1' else ('quantile',ordinate['L'],ordinate['b_y'],ordinate['q'],ordinate['precision_bits'])
                        if header not in scalar_digest_by_header:
                            encoded=json.dumps(ordinate,sort_keys=True,separators=(',',':')).encode();digest=hashlib.sha256(encoded).hexdigest()
                            scalar_digest_by_header[header]=digest
                            stream.write(json.dumps({'sha256':digest,'certificate':ordinate},sort_keys=True,separators=(',',':')).encode()+b'\n');scalar_seen.add(digest)
                    cert['precision_search_certificate_sha256']=hashlib.sha256(json.dumps(cert,sort_keys=True,separators=(',',':')).encode()).hexdigest()
                    evaluation_cache[key]=(graph,tables,cert,b)
                    return evaluation_cache[key]
                initial=evaluate({});initial_bits=initial[0].configuration['function_address_bits']
                def retain(item,label):
                    g,t,c,b=item;key=tuple(sorted(g.configuration['function_address_bits'].items()))
                    if key not in local_programs:local_programs[key]={'item':item,'allocation_limits':[]}
                    local_programs[key]['allocation_limits'].append(label)
                retain(initial,'unchanged_all12_baseline')
                for limit in [F(1,1000),F(1,400),F(1,200)]:
                    current=initial;bits=dict(initial_bits)
                    if current[2]['status']!='passed' or F(current[3]['total_deterministic_undiscounted_price_error_upper'])>limit:
                        trials.append({'model_id':mid,'cap':cap,'kind':kind,'allocation_limit':str(limit),'status':'starting_program_outside_this_allocation','reason':'The declared greedy policy skips starting allocation violations; feasibility of coarser programs is not inferred'});continue
                    step=0
                    while True:
                        choices=[];old_bias=F(current[3]['total_deterministic_undiscounted_price_error_upper'])
                        for name,b in bits.items():
                            if b<=6:continue
                            proposed=dict(bits);proposed[name]=b-1
                            item=evaluate(proposed);det=F(item[3]['total_deterministic_undiscounted_price_error_upper']) if item[3] else None
                            passed=item[2]['status']=='passed' and det<=limit
                            score=F(1<<(b-1))/max(F(1,1<<28),det-old_bias) if passed else None
                            trial_index=len(trials)
                            trials.append({'model_id':mid,'cap':cap,'kind':kind,'allocation_limit':str(limit),'node':name,
                                'from_bits':b,'to_bits':b-1,'deterministic_error_upper':str(det) if det is not None else None,'accepted_by_certificate':passed,'chosen':False,'step':step,'current_width_map':dict(bits),'proposed_width_map':proposed,'graph_certificate_sha256':item[2]['precision_search_certificate_sha256'],'graph_failures':item[2]['failures'],'greedy_score':str(score) if score is not None else None})
                            if passed:
                                # Declared cell-saving/error proxy guides a greedy finite search.
                                choices.append((score,name,item,proposed,trial_index))
                        if not choices:break
                        choices.sort(key=lambda x:(-x[0],x[1]));_,name,current,bits,chosen_index=choices[0]
                        trials[chosen_index]['chosen']=True;step+=1
                    retain(current,str(limit))
                for candidate in local_programs.values():
                    graph,tables,cert,budget=candidate['item']
                    if budget is None:
                        records.append({'model_id':mid,'cap':cap,'kind':kind,'allocation_limits':candidate['allocation_limits'],'base_graph_certificate':cert,'selected_strategy':None,'reason':'Rejected numerical graph retained without a monetary budget'});continue
                    bits=graph.configuration['function_address_bits']
                    prior=certificate.positive_r1_prior(model,cap,F(6),10,28,F(domain['U_star']),F(cert['probability_error_upper']),moment_records[mid]) if kind=='positive_r1' else None
                    profiles={};strategies=[]
                    for architecture in ('exact7t','unitary8t'):
                        encodings=['centred'] if kind=='positive_r1' else ['centred','positive_part','negative_part']
                        profiles[architecture]={}
                        for encoding in encodings:
                            part=graph if encoding=='centred' else build(model,cap,kind,domain,bits,encoding)
                            pt=tables if encoding=='centred' else oracle.make_tables(part,tables_cache)
                            profiles[architecture][encoding]=oracle.counted_composition(part,pt,block_cache,architecture)
                        for strategy in (['centred'] if kind=='positive_r1' else ['centred','split']):
                            plan=complete_plan(kind,strategy,budget,cert,moments,prior,profiles[architecture]);strategies.append({'architecture':architecture,'strategy':strategy,'result':plan})
                    record={'model_id':mid,'cap':cap,'kind':kind,'allocation_limits':candidate['allocation_limits'],'node_address_bits':bits,
                        'base_graph_certificate':cert,'complete_deterministic_budget':budget,'positive_r1_prior':prior,'profiles':profiles,'strategies':strategies}
                    eligible=[x for x in strategies if x['result']['kind']!='rejected']
                    record['selected_strategy']=min(eligible,key=lambda x:(x['result']['T_count'],x['architecture'],x['strategy'])) if eligible else None
                    records.append(record)
                print('Certified and costed',mid,cap,kind,len(local_programs),'programs',flush=True)
                tables_cache.clear();scalar_cache.clear()
        stream.close()
    selected=[]
    for financial in finances:
        row={'model_id':financial['model_id'],'cap':financial['cap'],'methods':{}}
        for kind in ['integral_r2','direct_r2','positive_r1']:
            candidates=[x for x in records if x['model_id']==row['model_id'] and x['cap']==row['cap'] and x['kind']==kind and x['selected_strategy'] is not None]
            winner=min(candidates,key=lambda x:x['selected_strategy']['result']['T_count']) if candidates else None
            row['methods'][kind]={'candidate_index':records.index(winner),'selected':winner['selected_strategy']} if winner else None
        row['quantum_only_method']=min((k for k in row['methods'] if row['methods'][k]),key=lambda k:row['methods'][k]['selected']['result']['T_count'])
        row['complete_policy']='Classical price certificate already attains tolerance; zero quantum estimation calls.'
        selected.append(row)
    family=[]
    for model in data['models']:
        for cap in [1,10,25]:
            family.append({'model_id':model['model_id'],'cap':cap,'status':'certified continuous-transfer and resource plans retained' if model['outer_dimension']==1 else 'fixed-width higher-dimensional financial grid not certified to monetary tolerance; not ranked as an inferior algorithm','best_classical_policy':'Validated Taylor or continuous quadrature certificate; zero quantum estimation calls'})
    output={'schema':'qbp-independent-precision-comparison-v1','status':'complete certified finite-pool comparison','monetary_radius':'1/100','alpha_per_selected_price':'1/2700',
        'records':records,'search_trials':trials,'selected_prices':selected,'all27_scope':family,'scalar_ordinate_archive':str(scalar_archive_path.name),
        'scalar_ordinate_archive_sha256':sha(scalar_archive_path),'unique_scalar_certificates':len(scalar_seen),
        'source_sha256':sha(Path(__file__)),'dependency_source_hashes':{n:sha(ROOT/n) for n in ['integral_value_oracle.py','integral_value_certificate.py','integral_auxiliary_quadrature.py','integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py','certified_scalar.py','quantum_gates.py','value_oracle.py','architecture.py','reversible.py']},
        'input_sha256':{n:sha(ROOT/'results'/n) for n in ['generalisation_protocol.json','gaussian_derivative_norms.json','integral_finite_input_domains.json','integral_finite_grid_moments.json','integral_financial_grid_assessment.json','gaussian_norm_pricing.json','continuous_quadrature.json','taylor_pricing.json','qft_rotation_tapes.json','qft_rotation_certificates.json']},
        'elapsed_seconds_diagnostic':time.perf_counter()-start,'interpretation':'Greedy independently certified scalar address widths under three deterministic error allocations plus unchanged all12 baseline; no global precision optimum, physical quantum execution, superiority to classical pricing or theoretical priority claim.'}
    path=ROOT/'results/integral_precision_comparison.json.gz'
    with path.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(json.dumps(output,default=str,separators=(',',':')).encode())
    for row in selected:print(row['model_id'],row['cap'],{k:v['selected']['result']['T_count'] for k,v in row['methods'].items()},flush=True)
    print(json.dumps({'status':'passed','programs':len(records),'trials':len(trials),'seconds':time.perf_counter()-start}),flush=True)
if __name__=='__main__':main()
