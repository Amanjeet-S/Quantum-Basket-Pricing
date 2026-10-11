"""Frozen finite-pool comparison of four coherent curvature preparations.

All methods receive the same payoff-specific derivative norms, financial
transfer certificate, conditional coefficients, exact finite moments and
confidence policy. Greedy address allocation is a reproducible search policy,
not a proof of global optimality. Every trial, failure and emitted charge is
retained. Native positive outputs use their one-sided monetary decoder.
"""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict
import gzip,hashlib,json,time
import boundary_value_oracle as o
import boundary_value_certificate as c
import boundary_auxiliary_certificate as a
import integral_value_certificate as old
import integral_pricing_budget as e
ROOT=Path(__file__).resolve().parent
KINDS=('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1')
ALPHA=F(1,3600);EPS=F(1,100)
def load(n):return json.loads((ROOT/'results'/n).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(iv):return F(iv['upper_exact_fraction'])
def lower(iv):return F(iv['lower_exact_fraction'])
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def casekey(x):return x['model_id'],x['payoff'],x['cap']
def build(model,cap,kind,domain,bits=None,encoding='centred'):
    g=o.build_graph(model,cap,q=28,b_y=10,beta_bits=8,table_bits=12,outer_bound=6,
        U_star=F(domain['U_star']),kind=kind,retained_bound=max(abs(F(x)) for x in domain['coordinate_bounds']),output_encoding=encoding)
    for n,b in (bits or {}).items():
        if n not in g.nodes or g.nodes[n].kind!='function' or not 6<=b<=12:raise ValueError('Invalid scalar width')
        g.nodes[n].bits=b
    g.configuration['function_address_bits']={n:v.bits for n,v in g.nodes.items() if v.kind=='function'}
    return g

def budget(model,cap,kind,transfer,coef,cert,moments):
    mass=transfer['mass_interval'];ml,mh=lower(mass),upper(mass);mr=(mh-ml)/2
    base=coef['c0'] if kind=='positive_r1' else coef['coefficient'];cr=(upper(base)-lower(base))/2
    B=F(cert['configuration']['monetary_signal_normaliser']);decode=F(cert['configuration']['monetary_probability_decode_factor'])
    aux=a.auxiliary_bias(model['lambda'],cap,F(cert['configuration']['U_star']),moments[2],moments[3],kind=kind,MU4=moments[4],MU6=moments[6])
    arith=decode*B*F(cert['probability_error_upper']);grid=F(transfer['financial_transfer_upper'])
    total=cr+grid+mh*(F(aux['bias_upper'])+arith)+mr*B
    return {'coefficient_interval':base,'coefficient_radius_upper':str(cr),'mass_interval':mass,'external_factor_midpoint':str((ml+mh)/2),
        'external_factor_upper':str(mh),'external_factor_radius':str(mr),'auxiliary_components':aux,
        'arithmetic_before_mass_upper':str(arith),'financial_transfer_upper':str(grid),'monetary_probability_decode_factor':str(decode),
        'total_deterministic_undiscounted_price_error_upper':str(total)}

def prior(kind,cap,cert,moments,gap):
    cfg=cert['configuration'];B=F(cfg['monetary_signal_normaliser']);eta=F(cert['probability_error_upper']);U=F(cfg['U_star'])
    positive=cfg['output_encoding']=='positive'
    if kind=='positive_r1':
        bounds=[gap]
        if cap is not None:bounds.extend([F(cap),F(cfg['M1'])*moments[1]])
        mean=min(bounds);prob=min(F(1),mean/B+eta)
        return {'probability_upper':str(prob),'ideal_mean_upper':str(mean),'underlying_gap_upper':str(gap),'native_positive':True,'rounding_enlargement':str(eta)}
    prob=min(F(1),moments[2]/U**2+(1 if positive else 2)*eta)
    return {'probability_upper':str(prob),'ideal_envelope_probability_upper':str(moments[2]/U**2),'native_positive':positive,'rounding_enlargement':str((1 if positive else 2)*eta)}

def complete(encoding,bud,cert,pr,profiles):
    det=F(bud['total_deterministic_undiscounted_price_error_upper']);remaining=EPS-det
    if remaining<=0:return {'kind':'rejected','reason':'Deterministic monetary envelope exhausts tolerance','deterministic_error_upper':str(det)}
    C=F(bud['external_factor_upper'])*F(cert['configuration']['monetary_signal_normaliser']);p=F(pr['probability_upper']);positive=pr['native_positive']
    radius=C*p/(2 if positive else 1)
    if radius<=remaining:return {'kind':'deterministic_correction_prior','T_count':0,'quantum_shots':0,'certified_correction_radius':str(radius),
        'total_price_error_upper':str(det+radius),'reason':'Shared correction interval resolves the monetary tolerance; analytical decision, no quantum gain.'}
    if positive:labels=['centred'];target=remaining/C;priors=[p];alphas=[ALPHA];decode=1
    elif encoding=='centred':labels=['centred'];target=remaining/(2*C);priors=[F(1)];alphas=[ALPHA];decode=2
    else:labels=['positive_part','negative_part'];target=remaining/(2*C);priors=[p,p];alphas=[ALPHA/2,ALPHA/2];decode=1
    flags=[]
    for label,up,alpha in zip(labels,priors,alphas):
        plan,qft=e.compiled_extended_plan(target,alpha,up);profile=profiles[label];charge=e.charge_full_oracle(plan,profile,qft)
        flags.append({'flag':label,'probability_tolerance':str(target),'prior_upper':str(up),'plan':asdict(plan),
            'financial_profile_sha256':digest(profile),'fourier_profile':qft,'complete_charge':charge})
    stat=C*decode*sum((x['plan']['amplitude_error_upper']+x['plan']['postprocess_allowance'] for x in flags),F(0))
    failure=sum((x['plan']['binomial_failure_upper'] for x in flags),F(0))
    assert det+stat<=EPS and failure<=ALPHA
    return {'kind':'complete_quantum_resource_plan','encoding':encoding,'T_count':sum(x['complete_charge'].get('metrics',{}).get('T',0) for x in flags),
        'quantum_shots':sum(x['plan']['repetitions'] for x in flags),'deterministic_error_upper':str(det),'statistical_price_error_upper':str(stat),
        'total_price_error_upper':str(det+stat),'failure_probability_upper':str(failure),'monetary_signal_factor_upper':str(C),'correction_prior':pr,
        'flags':flags,'logical_qubits_upper':max(x['complete_charge'].get('peak_logical_qubits',0) for x in flags),
        'execution_scope':'Constructed ideal-logical resource/confidence plan; no physical QPE experiment.'}

def main():
    start=time.perf_counter();names=['boundary_value_oracle.py','boundary_value_certificate.py','boundary_derivative_certificate.py','boundary_auxiliary_certificate.py','integral_value_oracle.py','integral_value_certificate.py','integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py','certified_scalar.py','quantum_gates.py','value_oracle.py','architecture.py','reversible.py']
    source_hashes={n:sha(ROOT/n) for n in names};protocol=load('generalisation_protocol.json');models={m['model_id']:m for m in protocol['models']}
    dom={x['model_id']:x for x in load('integral_finite_input_domains.json')['domains']}
    grid=load('boundary_financial_grid_transfer_summary.json');trans={(*casekey(x),x['order']):x for x in grid['records']}
    oldmom={x['model_id']:x for x in load('integral_finite_grid_moments.json')['records']};newmom={x['model_id']:x for x in load('boundary_finite_moments.json')['records']}
    block_cache={};records=[];trials=[];seen=set();scalar_cache={};scalar_original=c.scalar_certificate
    def cached_scalar(name,lo,span,params,q,bits,output_offset=0):
        key=(name,str(lo),str(span),tuple(sorted((str(k),str(v)) for k,v in params.items())),q,bits,str(output_offset))
        if key not in scalar_cache:scalar_cache[key]=scalar_original(name,lo,span,params,q,bits,output_offset)
        return scalar_cache[key]
    c.scalar_certificate=cached_scalar
    archive=ROOT/'results/boundary_precision_scalar_archive.jsonl.gz'
    with archive.open('wb') as raw:
        stream=gzip.GzipFile(fileobj=raw,mode='wb',mtime=0)
        for coef in grid['coefficients']:
            mid,payoff,cap=casekey(coef);model=models[mid];domain=dom[mid];tables_cache={}
            moments={1:upper(oldmom[mid]['centred_U_moments']['1'])}
            moments.update({r:upper(newmom[mid]['centred_U_moments'][str(r)]) for r in (2,3,4,6)})
            gaprec=old.positive_r1_prior(model,cap if cap is not None else 1,F(6),10,28,F(domain['U_star']),F(0),oldmom[mid]);gap=upper(gaprec['correlated_underlying_gap_interval'])
            for kind in KINDS:
                transfer=trans[(mid,payoff,cap,1 if kind=='positive_r1' else 2)];local={};cache={}
                def evaluate(bits):
                    g=build(model,cap,kind,domain,bits);key=tuple(sorted(g.configuration['function_address_bits'].items()))
                    if key in cache:return cache[key]
                    tables,ords=o.make_tables(g,tables_cache,True);cert=c.certify_graph(g,tables,ords)
                    bud=budget(model,cap,kind,transfer,coef,cert,moments) if cert['status']=='passed' else None
                    for ordinate in ords.values():
                        header={k:v for k,v in ordinate.items() if k not in ('nodes',)}
                        hd=digest(header)
                        if hd not in seen:
                            blob=json.dumps(ordinate,sort_keys=True,separators=(',',':')).encode();ds=hashlib.sha256(blob).hexdigest();seen.add(hd)
                            stream.write(json.dumps({'sha256':ds,'certificate':ordinate},sort_keys=True,separators=(',',':')).encode()+b'\n')
                    cert['precision_search_certificate_sha256']=digest(cert);cache[key]=(g,tables,cert,bud);return cache[key]
                initial=evaluate({});initial_bits=initial[0].configuration['function_address_bits']
                def retain(item,label):
                    key=tuple(sorted(item[0].configuration['function_address_bits'].items()))
                    if key not in local:local[key]={'item':item,'allocation_limits':[]}
                    local[key]['allocation_limits'].append(label)
                retain(initial,'unchanged_all12_baseline')
                for limit in (F(1,1000),F(1,400),F(1,200)):
                    current=initial;bits=dict(initial_bits)
                    if current[3] is None or F(current[3]['total_deterministic_undiscounted_price_error_upper'])>limit:
                        trials.append({'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limit':str(limit),'status':'starting_program_outside_this_allocation'});continue
                    step=0
                    while True:
                        choices=[];previous=F(current[3]['total_deterministic_undiscounted_price_error_upper'])
                        for n,b in bits.items():
                            if b<=6:continue
                            proposed=dict(bits);proposed[n]=b-1;item=evaluate(proposed);det=F(item[3]['total_deterministic_undiscounted_price_error_upper']) if item[3] else None
                            passed=item[2]['status']=='passed' and det<=limit;score=F(1<<(b-1))/max(F(1,1<<28),det-previous) if passed else None
                            ix=len(trials);trials.append({'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limit':str(limit),'node':n,'from_bits':b,'to_bits':b-1,
                                'deterministic_error_upper':str(det) if det is not None else None,'accepted_by_certificate':passed,'chosen':False,'step':step,
                                'current_width_map':dict(bits),'proposed_width_map':proposed,'graph_certificate_sha256':item[2]['precision_search_certificate_sha256'],
                                'graph_failures':item[2]['failures'],'greedy_score':str(score) if score is not None else None})
                            if passed:choices.append((score,n,item,proposed,ix))
                        if not choices:break
                        choices.sort(key=lambda x:(-x[0],x[1]));_,n,current,bits,ix=choices[0];trials[ix]['chosen']=True;step+=1
                    retain(current,str(limit))
                for candidate in local.values():
                    g,tables,cert,bud=candidate['item'];bits=g.configuration['function_address_bits']
                    rec={'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limits':candidate['allocation_limits'],'node_address_bits':bits,
                        'base_graph_certificate':cert,'complete_deterministic_budget':bud,'profiles':{},'strategies':[],'selected_strategy':None}
                    if bud is not None:
                        pr=prior(kind,cap,cert,moments,gap);rec['correction_prior']=pr;positive=pr['native_positive']
                        for architecture in ('exact7t','unitary8t'):
                            profiles={}
                            for encoding in (('centred',) if positive else ('centred','positive_part','negative_part')):
                                pg=g if encoding=='centred' else build(model,cap,kind,domain,bits,encoding);pt=tables if encoding=='centred' else o.make_tables(pg,tables_cache)
                                profiles[encoding]=o.counted_composition(pg,pt,block_cache,architecture)
                            rec['profiles'][architecture]=profiles
                            for strategy in (('positive',) if positive else ('centred','split')):
                                rec['strategies'].append({'architecture':architecture,'strategy':strategy,'result':complete(strategy,bud,cert,pr,profiles)})
                        valid=[x for x in rec['strategies'] if x['result']['kind']!='rejected']
                        rec['selected_strategy']=min(valid,key=lambda x:(x['result']['T_count'],x['architecture'],x['strategy'])) if valid else None
                    records.append(rec)
                print('Costed',mid,payoff,cap,kind,len(local),'programs',flush=True)
                tables_cache.clear();scalar_cache.clear()
        stream.close()
    selected=[]
    for coef in grid['coefficients']:
        mid,payoff,cap=casekey(coef);row={'model_id':mid,'payoff':payoff,'cap':cap,'methods':{}}
        for kind in KINDS:
            valid=[(ix,x) for ix,x in enumerate(records) if casekey(x)==casekey(coef) and x['kind']==kind and x['selected_strategy']]
            winner=min(valid,key=lambda pair:(pair[1]['selected_strategy']['result']['T_count'],pair[0])) if valid else None
            row['methods'][kind]={'candidate_index':winner[0],'selected':winner[1]['selected_strategy']} if winner else None
        row['quantum_resource_winner']=min((k for k in KINDS if row['methods'][k]),key=lambda k:row['methods'][k]['selected']['result']['T_count'])
        selected.append(row)
    assert source_hashes=={n:sha(ROOT/n) for n in names},'A frozen dependency changed during execution'
    inputs=['boundary_study_protocol.json','generalisation_protocol.json','boundary_derivative_norms.json','boundary_finite_moments.json','boundary_financial_grid_transfer_summary.json','integral_finite_input_domains.json','integral_finite_grid_moments.json','qft_rotation_tapes.json','qft_rotation_certificates.json']
    result={'schema':'qbp-boundary-precision-comparison-v1','status':'complete finite-pool comparison','monetary_radius':str(EPS),'alpha_per_selected_price':str(ALPHA),
        'family_failure_probability_union_upper':'1/100','records':records,'search_trials':trials,'selected_prices':selected,'all36_scope':load('boundary_study_protocol.json')['all_financial_cases'],
        'unique_counted_blocks':len(block_cache),'unique_scalar_certificates':len(seen),'scalar_ordinate_archive':archive.name,'scalar_ordinate_archive_sha256':sha(archive),
        'source_sha256':sha(__file__),'dependency_source_hashes':source_hashes,'input_sha256':{n:sha(ROOT/'results'/n) for n in inputs},'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'interpretation':'Shared sharper norms and price-independent financial transfer, same classical law, four methods, both coherent architectures, native positive or signed decoders. Finite greedy pool, no global optimum, quantum execution, classical superiority or theoretical priority claim.'}
    out=ROOT/'results/boundary_precision_comparison.json.gz'
    with out.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(json.dumps(result,default=str,separators=(',',':')).encode())
    summary={k:v for k,v in result.items() if k not in ('records','search_trials')};summary['archive_sha256']=sha(out);summary['candidate_programs']=len(records);summary['search_trial_count']=len(trials)
    (ROOT/'results/boundary_precision_summary.json').write_text(json.dumps(summary,default=str,indent=2)+'\n')
    for row in selected:print(row['model_id'],row['payoff'],row['cap'],{k:v['selected']['result']['T_count'] if v else None for k,v in row['methods'].items()},flush=True)
    print('Completed',len(records),'programs',len(trials),'trials',time.perf_counter()-start,'seconds',flush=True)
if __name__=='__main__':main()
