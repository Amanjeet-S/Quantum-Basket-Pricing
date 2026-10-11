"""Matched complete pricing plans after factorised Gaussian financial transfer.

The frozen higher-dimensional financial law is unchanged. All methods share
coefficients, moments, payoff norms and the same exact rounded input law.
Both the original integer emitter and a proved exact simplification pass are
admitted. Precision allocation is a declared finite greedy pool, not a global
optimum. No financial price answer is used in the construction.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,time
import boundary_precision_selection as retained
import boundary_value_oracle as o
import boundary_value_certificate as c
from boundary_grid_certificate import interval,serial
ROOT=Path(__file__).resolve().parent
KINDS=retained.KINDS

def load(name):return json.loads((ROOT/'results'/name).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def upper(iv):return F(iv['upper_exact_fraction'])
def lower(iv):return F(iv['lower_exact_fraction'])
def key(row):return row['model_id'],row['payoff'],None if row['cap'] is None else int(F(row['cap']))
def digest(value):return retained.digest(value)

def main():
    import boundary_integer_optimisation as optimise
    start=time.perf_counter()
    deps=['boundary_precision_selection.py','boundary_integer_optimisation.py','higherdim_preparation.py','tensor_transfer.py','tensor_first_transfer.py','boundary_value_oracle.py','boundary_value_certificate.py','boundary_derivative_certificate.py','boundary_auxiliary_certificate.py','boundary_grid_certificate.py','integral_value_oracle.py','integral_value_certificate.py','integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py','certified_scalar.py','quantum_gates.py','value_oracle.py','architecture.py','reversible.py','gaussian_derivative_norm.py','taylor_pricing.py','continuous_quadrature.py','generalisation_models.py']
    inputs=['higherdim_protocol.json','higherdim_preparation.json','tensor_transfer.json','tensor_first_transfer.json','generalisation_protocol.json','boundary_derivative_norms.json','gaussian_derivative_norms.json','qft_rotation_tapes.json','qft_rotation_certificates.json']
    hashes={n:sha(ROOT/n) for n in deps};input_hashes={n:sha(ROOT/'results'/n) for n in inputs}
    models={m['model_id']:m for m in load('generalisation_protocol.json')['models']}
    preparation={r['model_id']:r for r in load('higherdim_preparation.json')['records']}
    cases=[case for row in preparation.values() for case in row['cases']]
    mass_one=preparation[cases[0]['model_id']]['moments']
    from higherdim_preparation import get_mgf
    mass=get_mgf(6,10,28).node_certificate['mass_interval'];masses={p:serial(interval(lower(mass),upper(mass))**p) for p in (3,7)}
    transfer={}
    for filename in ('tensor_transfer.json','tensor_first_transfer.json'):
        for row in load(filename)['records']:
            if row['model_id'] in preparation:
                record=dict(row);record['mass_interval']=masses[models[row['model_id']]['outer_dimension']]
                transfer[(*key(row),row['order'])]=record
    records=[];trials=[];block_cache={};seen=set();scalar_cache={};original=c.scalar_certificate
    def cached_scalar(name,lo,span,params,q,bits,output_offset=0):
        signature=(name,str(lo),str(span),tuple(sorted((str(k),str(v)) for k,v in params.items())),q,bits,str(output_offset))
        if signature not in scalar_cache:scalar_cache[signature]=original(name,lo,span,params,q,bits,output_offset)
        return scalar_cache[signature]
    c.scalar_certificate=cached_scalar
    archive=ROOT/'results/higherdim_precision_scalar_archive.jsonl.gz'
    with archive.open('wb') as raw:
        stream=gzip.GzipFile(fileobj=raw,mode='wb',mtime=0)
        for case in cases:
            mid,payoff,cap=key(case);model=models[mid];prep=preparation[mid];domain=prep['domain'];coef=case['coefficients']
            moments={r:upper(prep['moments']['centred_U_moments'][str(r)]) for r in (1,2,3,4,6)}
            gap=upper(prep['moments']['correlated_underlying_gap_interval']);tables_cache={}
            for kind in KINDS:
                financial=transfer[(mid,payoff,cap,1 if kind=='positive_r1' else 2)];local={};cache={}
                def evaluate(bits):
                    graph=retained.build(model,cap,kind,domain,bits);signature=tuple(sorted(graph.configuration['function_address_bits'].items()))
                    if signature in cache:return cache[signature]
                    tables,ordinates=o.make_tables(graph,tables_cache,True);cert=c.certify_graph(graph,tables,ordinates)
                    bud=retained.budget(model,cap,kind,financial,coef,cert,moments) if cert['status']=='passed' else None
                    for ordinate in ordinates.values():
                        header={k:v for k,v in ordinate.items() if k!='nodes'};hd=digest(header)
                        if hd not in seen:
                            blob=json.dumps(ordinate,sort_keys=True,separators=(',',':')).encode();seen.add(hd)
                            stream.write(json.dumps({'sha256':hashlib.sha256(blob).hexdigest(),'certificate':ordinate},sort_keys=True,separators=(',',':')).encode()+b'\n')
                    cert['precision_search_certificate_sha256']=digest(cert);cache[signature]=(graph,tables,cert,bud);return cache[signature]
                initial=evaluate({});initial_bits=initial[0].configuration['function_address_bits']
                def retain(item,label):
                    signature=tuple(sorted(item[0].configuration['function_address_bits'].items()))
                    if signature not in local:local[signature]={'item':item,'allocation_limits':[]}
                    local[signature]['allocation_limits'].append(label)
                retain(initial,'unchanged_all12_baseline')
                for limit in (F(1,1000),F(1,400),F(1,200)):
                    current=initial;bits=dict(initial_bits)
                    if current[3] is None or F(current[3]['total_deterministic_undiscounted_price_error_upper'])>limit:
                        trials.append({'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limit':str(limit),'status':'starting_program_outside_this_allocation'});continue
                    step=0
                    while True:
                        choices=[];old_error=F(current[3]['total_deterministic_undiscounted_price_error_upper'])
                        for name,bits_old in bits.items():
                            if bits_old<=6:continue
                            proposed=dict(bits);proposed[name]=bits_old-1;item=evaluate(proposed)
                            error=F(item[3]['total_deterministic_undiscounted_price_error_upper']) if item[3] else None
                            passed=item[2]['status']=='passed' and error<=limit;score=F(1<<(bits_old-1))/max(F(1,1<<28),error-old_error) if passed else None
                            ix=len(trials);trials.append({'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limit':str(limit),'node':name,'from_bits':bits_old,'to_bits':bits_old-1,
                                'deterministic_error_upper':str(error) if error is not None else None,'accepted_by_certificate':passed,'chosen':False,'step':step,
                                'current_width_map':dict(bits),'proposed_width_map':proposed,'graph_certificate_sha256':item[2]['precision_search_certificate_sha256'],
                                'graph_failures':item[2]['failures'],'greedy_score':str(score) if score is not None else None})
                            if passed:choices.append((score,name,item,proposed,ix))
                        if not choices:break
                        choices.sort(key=lambda x:(-x[0],x[1]));_,name,current,bits,ix=choices[0];trials[ix]['chosen']=True;step+=1
                    retain(current,str(limit))
                for candidate in local.values():
                    graph,tables,cert,bud=candidate['item'];bits=graph.configuration['function_address_bits']
                    record={'model_id':mid,'payoff':payoff,'cap':cap,'kind':kind,'allocation_limits':candidate['allocation_limits'],'node_address_bits':bits,
                        'base_graph_certificate':cert,'complete_deterministic_budget':bud,'profiles':{},'strategies':[],'selected_strategy':None}
                    if bud is not None:
                        pr=retained.prior(kind,cap,cert,moments,gap);record['correction_prior']=pr;positive=pr['native_positive']
                        for implementation in ('original','exact_integer_simplification'):
                            record['profiles'][implementation]={}
                            for architecture in ('exact7t','unitary8t'):
                                profiles={}
                                for encoding in (('centred',) if positive else ('centred','positive_part','negative_part')):
                                    pg=graph if encoding=='centred' else retained.build(model,cap,kind,domain,bits,encoding)
                                    pt=tables if encoding=='centred' else o.make_tables(pg,tables_cache)
                                    if implementation=='exact_integer_simplification':pg,proof=optimise.optimise_graph(pg,pt)
                                    else:proof={'applied':False,'scope':'Unmodified emitted integer graph retained as baseline.'}
                                    profile=o.counted_composition(pg,pt,block_cache,architecture)
                                    profile['integer_optimisation']=proof;profile['implementation']=implementation
                                    profile['integer_optimisation_source_sha256']=hashes['boundary_integer_optimisation.py']
                                    profiles[encoding]=profile
                                record['profiles'][implementation][architecture]=profiles
                                for strategy in (('positive',) if positive else ('centred','split')):
                                    record['strategies'].append({'implementation':implementation,'architecture':architecture,'strategy':strategy,
                                        'result':retained.complete(strategy,bud,cert,pr,profiles)})
                        accepted=[x for x in record['strategies'] if x['result']['kind']!='rejected']
                        record['selected_strategy']=min(accepted,key=lambda x:(x['result']['T_count'],x['implementation'],x['architecture'],x['strategy'])) if accepted else None
                    records.append(record)
                print('Costed',mid,payoff,cap,kind,len(local),'programs',flush=True)
                tables_cache.clear();scalar_cache.clear()
        stream.close()
    selected=[]
    for case in cases:
        mid,payoff,cap=key(case);row={'model_id':mid,'payoff':payoff,'cap':cap,'methods':{}}
        for kind in KINDS:
            eligible=[(ix,r) for ix,r in enumerate(records) if key(r)==key(case) and r['kind']==kind and r['selected_strategy']]
            winner=min(eligible,key=lambda pair:(pair[1]['selected_strategy']['result']['T_count'],pair[0])) if eligible else None
            row['methods'][kind]={'candidate_index':winner[0],'selected':winner[1]['selected_strategy']} if winner else None
        available=[kind for kind in KINDS if row['methods'][kind]]
        row['quantum_resource_winner']=min(available,key=lambda k:row['methods'][k]['selected']['result']['T_count']) if available else None
        selected.append(row)
    assert hashes=={n:sha(ROOT/n) for n in deps},'An executed dependency changed during costing'
    assert input_hashes=={n:sha(ROOT/'results'/n) for n in inputs},'An executed input changed during costing'
    result={'schema':'qbp-higherdim-complete-precision-comparison-v1','status':'completed with all rejections retained','monetary_radius':'1/100','alpha_per_selected_price':'1/3600',
        'union_over_original36_failure_upper':'1/100','records':records,'search_trials':trials,'selected_prices':selected,'unique_counted_blocks':len(block_cache),
        'unique_scalar_certificates':len(seen),'scalar_ordinate_archive':archive.name,'scalar_ordinate_archive_sha256':sha(archive),
        'source_sha256':sha(__file__),'dependency_source_hashes':hashes,'input_sha256':input_hashes,'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'interpretation':'New complete p3/p7 financial plans with positive transport-moment transfer, equally informed four methods, two coherent architectures and proved identical-output integer simplification. Finite greedy pool, no globally optimal algorithm, device execution, classical superiority or priority claim.'}
    out=ROOT/'results/higherdim_precision_comparison.json.gz'
    with out.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as z:z.write(json.dumps(result,default=str,separators=(',',':')).encode())
    summary={k:v for k,v in result.items() if k not in ('records','search_trials')};summary['archive_sha256']=sha(out);summary['candidate_programs']=len(records);summary['search_trial_count']=len(trials)
    (ROOT/'results/higherdim_precision_summary.json').write_text(json.dumps(summary,default=str,indent=2)+'\n')
    for row in selected:print(row['model_id'],row['payoff'],row['cap'],{k:v['selected']['result']['T_count'] if v else None for k,v in row['methods'].items()},flush=True)
    print('Completed',len(records),'programs',len(trials),'decisions',time.perf_counter()-start,'seconds',flush=True)
if __name__=='__main__':main()
