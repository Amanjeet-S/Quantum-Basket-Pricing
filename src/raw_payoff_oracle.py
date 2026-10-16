"""Registered, fully emitted common-normal capped-payoff R1 opponent.

Established monotone clamp decomposition avoids useless wide-range exponentials.
This is a finite financial implementation, not new quantum theory.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime, timezone
import argparse, copy, gzip, hashlib, json, sys, time
sys.dont_write_bytecode=True
ROOT=RESULTS=FACTOR=THRESHOLD=NUMERICAL=None
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def write(p,v):
    if p.exists():raise FileExistsError('Preserve existing evidence')
    p.write_text(json.dumps(v,indent=2,default=str)+'\n')
def configure(repo,factor,threshold,numerical,output):
    global ROOT,RESULTS,FACTOR,THRESHOLD,NUMERICAL,o,c,informed,shared
    ROOT=Path(repo).resolve();RESULTS=Path(output).resolve()
    FACTOR=Path(factor).resolve();THRESHOLD=Path(threshold).resolve();NUMERICAL=Path(numerical).resolve()
    RESULTS.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(ROOT))
    import boundary_value_oracle as o
    import boundary_value_certificate as c
    import correction_information_priors as informed
    import shared_threshold_oracle as shared
DEPENDENCIES=('boundary_value_oracle.py','boundary_value_certificate.py','integral_value_oracle.py',
 'integral_value_certificate.py','boundary_derivative_certificate.py','value_oracle.py','reversible.py',
 'architecture.py','quantum_gates.py','boundary_precision_selection.py','correction_information_priors.py',
 'integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py','certified_scalar.py','taylor_pricing.py',
 'continuous_quadrature.py','shared_threshold_oracle.py')
INPUTS=('gaussian_derivative_norms.json','boundary_derivative_norms.json','qft_rotation_tapes.json','qft_rotation_certificates.json')
def signatures():
    return {'source_sha256':sha(__file__),'factor_resources_sha256':sha(FACTOR),
      'threshold_resources_sha256':sha(THRESHOLD),'numerical_resources_sha256':sha(NUMERICAL),
      'dependencies_sha256':{n:sha(ROOT/n) for n in DEPENDENCIES},
      'inputs_sha256':{n:sha(ROOT/'results'/n) for n in INPUTS}}
def register():
    threshold=json.loads(THRESHOLD.read_text())
    write(RESULTS/'raw_payoff_protocol.json',{'schema':'qbp-raw-payoff-protocol-v1',
      'registered_before_fresh_costs':True,'registration_utc':datetime.now(timezone.utc).isoformat(),
      'research_record_date':'2026-10-10',**signatures(),
      'pairs':[[r['original_model_id'],r['candidate_id']] for r in threshold['rows']],
      'inner_normal_bits':[8,10,12],'narrow_exp_address_bits':[4,6,8,10],
      'architectures':['exact7t','unitary8t'],'configurations_per_pair':12,'configurations_total':48,
      'q':28,'outer_L':'6','outer_quantile_bits':10,'outer_scalar_table_bits':12,
      'inner_normal_L':'6','payoff_exp_span':'1/64','K':'100','H':'1',
      'representation':'Positive capped payoff difference on exactly the same normal label. Each payoff is clamp(100 positive_sub(exp(clamp(positive_sub(shifted_log,common_offset),1/64)),1),1).',
      'error_rule':'Graph absolute errors plus 2 times Gaussian Kolmogorov discrepancy; exact monotone clamps implement kinks and restrict exp to[0,1/64].',
      'selection':'Minimum complete canonical QPE T cap over the frozen pool and both architectures, retaining all failures.',
      'confidence':'Same alpha1/3600 and undiscounted monetary radius1/100.',
      'known_predecessors':'Shared threshold and strengthened numerical R1 outcomes are known. No fresh raw-payoff emitted cost was inspected.',
      'no_repairs':True,'financial_price_answers_excluded':True,
      'data_origin':'Synthetic exact rational no-arbitrage models only; no licensed observations.'})

def build_graph(model,domain,normal_bits,exp_bits):
    old=o.build_graph(model,F(1),q=28,b_y=10,table_bits=12,outer_bound=6,
       U_star=F(domain['U_star']),retained_bound=max(abs(F(v)) for v in domain['coordinate_bounds']),kind='positive_r1')
    mG0=old.semantic_nodes['mG_offset'];log0=old.semantic_nodes['logT'];BG=-old.nodes[old.semantic_nodes['PG']].meta.lo
    needed=set()
    def visit(n):
        if n in old.inputs or n in needed:return
        needed.add(n)
        for a in old.nodes[n].args:visit(a)
    visit(mG0);visit(log0)
    g=o.IntegralGraph(28,12);remap={n:g.inp(n,b) for n,b in old.inputs.items()}
    for name,node in old.nodes.items():
        if name not in needed:continue
        new=g.op(node.kind,*(remap[a] for a in node.args),meta=node.meta,bits=node.bits,
                 output_bits=0,bounds=old.bounds[name]);remap[name]=new
        if node.kind=='constant':g.constant_names[node.meta]=new
    mG,logT=remap[mG0],remap[log0];lam=F(model['lambda']);L=F(6)
    J=o.old.power_span(max(BG,lam*L)+F(1,64)+F(16,1<<g.q))
    label=g.inp('inner_normal_label',normal_bits)
    nplusL=g.op('quantile',label,meta=L,bits=normal_bits,bounds=(0,2*L))
    lower=g.add(g.add(mG,g.const(J-BG)),g.scale(nplusL,lam))
    upper=g.add(lower,logT);offset=g.const(J+lam*L)
    def payoff(endpoint):
        positive=g.positive_sub(endpoint,offset)
        bounded=g.clamp(positive,F(1,64))
        exponential=g.scalar('exp',bounded,0,F(1,64),output_bounds=(1,F(2)))
        g.nodes[exponential].bits=exp_bits
        value=g.clamp(g.scale(g.positive_sub(exponential,g.const(1)),100),1)
        return value,exponential,positive,bounded
    PG,eG,xG,bG=payoff(lower);PA,eA,xA,bA=payoff(upper)
    output=g.clamp(g.positive_sub(PA,PG),1)
    g.semantic_nodes={'mG_offset':mG,'logT':logT,'lower_endpoint':lower,'upper_endpoint':upper,
       'normal_plus_L':nplusL,'PG':PG,'PA':PA,'PG_exp':eG,'PA_exp':eA,
       'PG_positive_log':xG,'PA_positive_log':xA,'PG_bounded_log':bG,'PA_bounded_log':bA,
       'common_offset':offset,'probability':output}
    g.configuration=copy.deepcopy(old.configuration)
    g.configuration.update({'pricing_representation':'common_normal_raw_capped_payoff_difference',
       'monetary_signal_normaliser':'1','positive_R1_normaliser':'1','output_encoding':'positive',
       'monetary_probability_decode_factor':1,'normal_bits':normal_bits,'payoff_exp_bits':exp_bits,
       'inner_normal_L':'6','payoff_exp_span':'1/64','common_endpoint_offset':str(J+lam*L)})
    return g.finish(output)

def probability_certificate(g,body,ordinates,nu0):
    if body['status']!='passed':return {'status':'rejected','reason':'Arithmetic domain/range guard failure','body_arithmetic_certificate':body}
    normal=ordinates[g.semantic_nodes['normal_plus_L']];mass=normal['mass_interval'];a=F(mass['lower_exact_fraction'])
    eN=F(normal['maximum_node_error_upper']);N=1<<g.configuration['normal_bits']
    D=(1-a)/2+a/F(2*N)+nu0*eN
    arithmetic=F(body['probability_error_upper']);eta=min(F(1),arithmetic+2*D)
    return {'schema':'qbp-raw-payoff-probability-certificate-v1','status':'passed',
      'configuration':g.configuration,'body_arithmetic_certificate':body,'nodes':body['nodes'],
      'probability_error_upper':str(eta),'monetary_value_error_upper':str(eta),
      'probability_actual_bounds':body['probability_actual_bounds'],
      'arithmetic_probability_error_upper':str(arithmetic),'inner_normal_CDF_error_upper':str(D),
      'inner_normal_payoff_difference_bias_upper':str(2*D),'normal_quantile_rounding_error_upper':str(eN),
      'normal_mass_interval':mass,'decoded_observable':'Nonnegative difference of two capped payoffs on the same normal draw.',
      'positivity':'Represented upper endpoint equals represented lower plus nonnegative represented logT. Identical monotone payoff map and positive subtraction ensure a native [0,1] output.',
      'scope':'Every outer realised coordinate; normalized inner normal approximation charged by total variation<=2. Only original outer box mass is decoded externally.'}

def run():
    protocol=json.loads((RESULTS/'raw_payoff_protocol.json').read_text())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Frozen inputs changed')
    factor=json.loads(FACTOR.read_text());threshold=json.loads(THRESHOLD.read_text());numeric=json.loads(NUMERICAL.read_text())
    norms=json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())
    nu0=F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==0 and r['passed']))
    rows=[];comparisons=[];table_cache={};block_cache={};archive=[];seen=set();start=time.perf_counter()
    for original,candidate in protocol['pairs']:
        prep=next(p for p in factor['preparations'] if p['original_model_id']==original and p['candidate_id']==candidate)
        old=next(r for r in factor['records'] if r['original_model_id']==original and r['factor_basis_candidate_id']==candidate and r['kind']=='positive_r1')
        for bn,be in product(protocol['inner_normal_bits'],protocol['narrow_exp_address_bits']):
            row={'original_model_id':original,'candidate_id':candidate,'normal_bits':bn,'payoff_exp_bits':be,'profiles':{},'strategies':[]}
            try:
                g=build_graph(prep['model'],prep['domain'],bn,be);tables,ordinates=o.make_tables(g,table_cache,True)
                body=c.certify_graph(g,tables,ordinates);cert=probability_certificate(g,body,ordinates,nu0)
                row.update(model=prep['model'],graph_configuration=g.configuration,graph_semantic_nodes=g.semantic_nodes,
                    graph_nodes=[{'name':n,'kind':v.kind,'args':v.args,'meta':v.meta.as_dict() if isinstance(v.meta,o.ScalarSpec) else str(v.meta),
                        'bits':v.bits,'output_bits':v.output_bits,'ideal_bounds':list(map(str,g.bounds[n]))} for n,v in g.nodes.items()],
                    probability_certificate=cert)
                for n,cc in ordinates.items():
                    key=(digest(cc),g.width)
                    if key not in seen:
                        seen.add(key);tab=tables[o.table_key(n,g.nodes[n])]
                        archive.append({'certificate_sha256':key[0],'certificate':cc,'packed_table':tab,
                           'packed_table_sha256':digest(tab),'graph_word_width':g.width})
                if cert['status']!='passed':row['rejection']={'reason':cert['reason']}
                else:
                    eta=F(cert['probability_error_upper']);bud=copy.deepcopy(old['complete_deterministic_budget'])
                    mass=F(bud['external_factor_upper']);mr=F(bud['external_factor_radius']);ec=F(bud['coefficient_radius_upper']);beta=F(bud['financial_transfer_upper'])
                    det=ec+beta+mass*eta+mr
                    bud.update(arithmetic_before_mass_upper=cert['arithmetic_probability_error_upper'],
                      auxiliary_components={'inner_normal_payoff_difference_bias_upper':str(2*F(cert['inner_normal_CDF_error_upper'])),'bias_upper':str(2*F(cert['inner_normal_CDF_error_upper']))},
                      total_deterministic_undiscounted_price_error_upper=str(det))
                    prev=old['coefficient_informed_prior'];lo,hi=F(prev['finite_ideal_mean_lower']),F(prev['finite_ideal_mean_upper'])
                    prior={'probability_lower':str(max(F(0),lo-eta)),'probability_upper':str(min(F(1),hi+eta)),
                      'native_positive':True,'finite_ideal_mean_lower':str(lo),'finite_ideal_mean_upper':str(hi),
                      'previous_mean_information_sha256':digest(prev),
                      'information_scope':'Identical price-independent finite ideal R1 mean interval, enlarged by full raw-payoff conditional error.'}
                    row.update(complete_deterministic_budget=bud,coefficient_informed_prior=prior)
                    if det>=F(1,100):row['rejection']={'reason':'Deterministic monetary tolerance exhausted','deterministic_error_upper':str(det)}
                    else:
                        for mode in protocol['architectures']:
                            profile=o.counted_composition(g,tables,block_cache,mode);row['profiles'][mode]={'centred':profile}
                            plan=informed.informed_complete_plan('positive',bud,cert,prior,row['profiles'][mode])
                            row['strategies'].append({'architecture':mode,'result':plan})
                        row['selected_strategy']=min(row['strategies'],key=lambda s:(s['result']['T_count'],s['architecture']))
            except (ArithmeticError,ValueError,AssertionError) as error:
                row['rejection']={'reason':'Explicit production failure','type':type(error).__name__,'message':str(error)}
            rows.append(row)
            print(original,bn,be,'T',row.get('selected_strategy',{}).get('result',{}).get('T_count'),'reject',row.get('rejection',{}).get('reason'),flush=True)
        admitted=[r for r in rows if r['original_model_id']==original and r.get('selected_strategy')]
        tr=next(r for r in threshold['rows'] if r['original_model_id']==original)
        nr=next(r for r in numeric['selected_comparisons'] if r['original_model_id']==original)
        comparison={'original_model_id':original,'candidate_id':candidate,'admitted_configurations':len(admitted),
          'rejected_configurations':12-len(admitted),'fixed_threshold_selected_strategy':tr['selected_strategy'],
          'strengthened_numeric_selected_strategy':nr['selected_numerical_strategy']}
        if admitted:
            winner=min(admitted,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['payoff_exp_bits']))
            raw=winner['selected_strategy']['result']['T_count'];thresholdT=tr['selected_strategy']['result']['T_count'];numericT=nr['selected_numerical_strategy']['result']['T_count']
            comparison.update(selected_normal_bits=winner['normal_bits'],selected_payoff_exp_bits=winner['payoff_exp_bits'],
              selected_raw_payoff_strategy=winner['selected_strategy'],
              threshold_relative_T_saving_against_raw_exact=str(F(raw-thresholdT,raw)),
              raw_relative_T_saving_against_numeric_exact=str(F(numericT-raw,numericT)))
        comparisons.append(comparison)
    path=RESULTS/'raw_payoff_scalar_ordinates.json.gz'
    if path.exists():raise FileExistsError('Preserve archive')
    with path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as stream:stream.write(json.dumps(archive,separators=(',',':'),default=str).encode())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Frozen inputs changed during execution')
    write(RESULTS/'raw_payoff_resources.json',{'schema':'qbp-raw-payoff-resources-v1','status':'Completed retaining every failure',
       'protocol_sha256':sha(RESULTS/'raw_payoff_protocol.json'),**signatures(),'records':rows,'selected_comparisons':comparisons,
       'scalar_archive_file':path.name,'scalar_archive_sha256':sha(path),'unique_scalar_certificates_and_wordwidths':len(archive),
       'emitted_block_profiles':len(block_cache),'generation_seconds_diagnostic':time.perf_counter()-start,
       'scope':'Complete clean coherent financial preparations and canonical QPE resources with identical monetary accuracy and confidence. No hardware, quantum price or practical advantage.'})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run'));p.add_argument('--repo',required=True)
    p.add_argument('--factor-resources',required=True);p.add_argument('--threshold-resources',required=True)
    p.add_argument('--numerical-resources',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    configure(a.repo,a.factor_resources,a.threshold_resources,a.numerical_resources,a.output)
    register() if a.mode=='register' else run()
