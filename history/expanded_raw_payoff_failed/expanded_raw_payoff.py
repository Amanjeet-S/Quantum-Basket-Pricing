"""Separately registered widening of the checked raw-payoff precision pool.

Every previous outcome is disclosed as known. Only new integer exp addresses
11 through 14 are emitted. Old evidence is referenced by hash and remains intact.
Explicit paths support producer replay without repository mutation.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
import argparse,copy,gzip,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p):
 p=Path(p)
 return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def write(p,v,compact=False):
 p=Path(p)
 if p.exists():raise FileExistsError('Preserve previously retained evidence: '+str(p))
 data=json.dumps(v,separators=(',',':') if compact else None,indent=None if compact else 2,default=str)+'\n'
 if p.suffix=='.gz':
  with p.open('wb') as raw:
   with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as stream:stream.write(data.encode())
 else:p.write_text(data)
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def paths(args):
 repo=Path(args.repo).resolve()
 return {'repo':repo,'factor':repo/'results/factor_basis_quantum_resources.json',
 'threshold':repo/'results/shared_threshold_resources.json','numerical':repo/'results/strengthened_r1_resources.json',
 'old_raw':Path(args.old_raw).resolve(),'old_smoothing':Path(args.old_smoothing).resolve(),
 'proof':Path(args.proof).resolve(),'output':Path(args.output).resolve()}
def signatures(args,raw):
 p=paths(args);base=raw.signatures()
 return {'source_sha256':sha(__file__),'raw_builder_source_sha256':sha(p['repo']/'raw_payoff_oracle_tight.py'),
 'smoothing_source_sha256':sha(p['repo']/'raw_payoff_smoothing_amendment.py'),
 'old_raw_resource_sha256':sha(p['old_raw']/'raw_payoff_resources.json.gz'),
 'old_raw_protocol_sha256':sha(p['old_raw']/'raw_payoff_protocol.json'),
 'old_raw_archive_sha256':sha(p['old_raw']/'raw_payoff_scalar_ordinates.json.gz'),
 'old_smoothing_resource_sha256':sha(p['old_smoothing']/'raw_payoff_smoothing_resources.json'),
 'old_smoothing_protocol_sha256':sha(p['old_smoothing']/'raw_payoff_smoothing_protocol.json'),
 'old_smoothing_audit_sha256':sha(p['old_smoothing']/'raw_payoff_smoothing_audit.json'),
 'checked_smoothing_proof_sha256':sha(p['proof']),
 **{k:v for k,v in base.items() if k!='source_sha256'}}
def configured(args):
 p=paths(args);p['output'].mkdir(parents=True,exist_ok=True)
 raw=module('registered_existing_raw_builder',p['repo']/'raw_payoff_oracle_tight.py')
 raw.configure(p['repo'],p['factor'],p['threshold'],p['numerical'],p['output'])
 smooth=module('checked_existing_smoothing_mapping',p['repo']/'raw_payoff_smoothing_amendment.py')
 return p,raw,smooth
def register(args):
 p,raw,smooth=configured(args);threshold=read(p['threshold']);old=read(p['old_raw']/'raw_payoff_resources.json.gz');amended=read(p['old_smoothing']/'raw_payoff_smoothing_resources.json')
 assert len(old['records'])==len(amended['records'])==48
 audit=read(p['old_smoothing']/'raw_payoff_smoothing_audit.json');assert audit['status']=='passed' and audit['resources_sha256']==sha(p['old_smoothing']/'raw_payoff_smoothing_resources.json')
 write(p['output']/'expanded_raw_payoff_protocol.json',{'schema':'qbp-expanded-raw-payoff-protocol-v1',
 'registration_utc':datetime.now(timezone.utc).isoformat(),'research_record_date':'2026-10-10',
 **signatures(args,raw),'registered_before_new_scalar_and_profile_costs':True,
 'known_previous_outcomes':{'corrected_raw_configurations':48,'smoothing_amended_configurations':48,
 'old_exp_widths':[4,6,8,10],'old_winner_summaries':amended['selected_comparisons'],
 'disclosure':'All corrected raw-payoff and smoothing-first outcomes in the old 48-configuration pool, the fixed threshold and strengthened numerical outcomes are already known. This is a disclosed sequential extension, not an initial blind registration.'},
 'pairs':[[r['original_model_id'],r['candidate_id']] for r in threshold['rows']],
 'inner_normal_bits':[8,10,12],'new_narrow_exp_address_bits':[11,12,13,14],
 'combined_narrow_exp_address_bits':[4,6,8,10,11,12,13,14],
 'architectures':['exact7t','unitary8t'],'new_configurations':48,'combined_configurations':96,
 'q':28,'outer_quantile_bits':10,'outer_scalar_table_bits':12,'outer_L':'6','inner_normal_L':'6',
 'payoff_exp_span':'1/64','payoff_exp_output_range':['1','64/63'],'K':'100','H':'1',
 'representation':'Exactly the checked corrected monotone clamp graph, with only the narrow exponential address changed to each new registered integer width.',
 'error_rule':'Minimum of original pointwise conditional error and separately proved smoothing-first conditional error, eta=min(pointwise,min(1,100(e_exp_G+e_exp_A)+202E_s+2D_N+(nu0/lambda)(2(e_m+e_C+e_O)+e_Delta))). Intrinsic exponential errors exclude argument perturbations. Inner normal rounding enters only D_N in the smoothing rule.',
 'scalar_policy':'Fresh outward scalar certificate and packed table for each newly constructed DAG; all numeric ordinates retained. No repaired profile, sampled cost surrogate or financial price answer.',
 'full_profile_policy':'Fresh counted composition from literal existing primitive emitters, both architectures, then freshly rebuilt price-independent prior and canonical QPE/full charge plan.',
 'financial_and_confidence_policy':'Unchanged factor-basis models, exact smooth outer R1 transfer, coefficients and finite ideal mean interval; alpha1/3600 and undiscounted monetary radius1/100.',
 'selection':'Minimum complete canonical QPE T cap among exactly 96 old-plus-new configurations and both architectures; tie-break by normal bits, exp bits then architecture. Compare threshold to minimum of selected raw and strengthened numerical rival. Retain negative or zero savings.',
 'failure_policy':'Retain every graph/scalar/arithmetic/profile/QPE production failure and every exhausted deterministic budget. Do not silently omit failed configurations.',
 'boundary_policy':'A selected exp width14 is reported as an upper-boundary winner. No adaptive widening or global optimum claim is authorised by this protocol.',
 'data_origin':'Synthetic exact rational no-arbitrage models only; no licensed observations, hardware outcomes or quantum price estimates.'})
 print('Registered extension before all new outcomes',flush=True)
def make_amended(base,smooth,archive,nu0,frontend,old):
 binding={'model_id':base['original_model_id'],'candidate_id':base['candidate_id'],'normal_bits':base['normal_bits'],'exp_bits':base['payoff_exp_bits'],'graph_sha256':digest(base['graph_nodes']),'configuration_sha256':digest(base['graph_configuration']),'profile_sha256':digest(base['profiles'])}
 if base.get('rejection') and not base.get('complete_deterministic_budget'):
  return {**binding,'pointwise_row_sha256':digest(base),'rejection':copy.deepcopy(base['rejection']),'strategies':[]}
 term=smooth.mapping(base,archive,nu0);eta=F(term['selected_probability_error_upper']);bud=copy.deepcopy(base['complete_deterministic_budget']);previous=old['coefficient_informed_prior'];A=F(bud['external_factor_upper']);D=F(term['inner_normal_CDF_error_upper'])
 det=F(bud['coefficient_radius_upper'])+F(bud['financial_transfer_upper'])+A*eta+F(bud['external_factor_radius'])
 bud.update(arithmetic_before_mass_upper=str(min(F(base['probability_certificate']['arithmetic_probability_error_upper']),F(term['smoothing_arithmetic_error_upper']))),total_deterministic_undiscounted_price_error_upper=str(det),auxiliary_components={'bias_upper':str(2*D),'inner_normal_payoff_difference_bias_upper':str(2*D)},certificate_amendment='smoothing_first_uniform_conditional')
 prior={'probability_lower':str(max(F(0),F(previous['finite_ideal_mean_lower'])-eta)),'probability_upper':str(min(F(1),F(previous['finite_ideal_mean_upper'])+eta)),'native_positive':True,'finite_ideal_mean_lower':previous['finite_ideal_mean_lower'],'finite_ideal_mean_upper':previous['finite_ideal_mean_upper'],'previous_mean_information_sha256':digest(previous),'information_scope':'Unchanged price-independent finite ideal R1 mean interval, enlarged by minimum of checked pointwise and smoothing-first conditional error.'}
 cert={'status':'passed','configuration':base['graph_configuration'],'probability_error_upper':str(eta),'monetary_value_error_upper':str(eta),'probability_actual_bounds':base['probability_certificate']['probability_actual_bounds'],'retained_pointwise_certificate_sha256':digest(base['probability_certificate']),'smoothing_derivation':term,'scope':'Uniform conditional inner-law expectation error; original pointwise body certificate remains retained.'}
 row={**binding,'pointwise_row_sha256':digest(base),'probability_certificate':cert,'complete_deterministic_budget':bud,'coefficient_informed_prior':prior,'strategies':[]}
 if det>=F(1,100):row['rejection']={'reason':'Amended deterministic monetary tolerance exhausted','deterministic_error_upper':str(det)}
 else:
  for mode in ('exact7t','unitary8t'):
   if mode not in base['profiles']:
    result={'kind':'production_failure','reason':'Missing literal profile','message':'Retained profile production failure'}
   else:
    try:result=frontend.informed_complete_plan('positive',bud,cert,prior,base['profiles'][mode])
    except (ArithmeticError,ValueError,AssertionError) as e:result={'kind':'production_failure','reason':type(e).__name__,'message':str(e)}
   row['strategies'].append({'architecture':mode,'result':result})
  valid=[r for r in row['strategies'] if r['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')]
  if valid:row['selected_strategy']=min(valid,key=lambda s:(s['result']['T_count'],s['architecture']))
  else:row['rejection']={'reason':'No valid amended architecture plan'}
 return row
def combined_summary(old,new,threshold,numerical):
 allrows=[dict(r,evidence_generation='retained48',retained_record_sha256=digest(r)) for r in old['records']]+[dict(r,evidence_generation='new48',retained_record_sha256=digest(r)) for r in new]
 assert len(allrows)==96 and len({(r['model_id'],r['normal_bits'],r['exp_bits']) for r in allrows})==96
 comparisons=[]
 for mid in dict.fromkeys(r['model_id'] for r in allrows):
  cohort=[r for r in allrows if r['model_id']==mid];valid=[r for r in cohort if r.get('selected_strategy')]
  num=next(r['selected_numerical_strategy'] for r in numerical['selected_comparisons'] if r['original_model_id']==mid);thr=next(r['selected_strategy'] for r in threshold['rows'] if r['original_model_id']==mid)
  cr={'model_id':mid,'admitted_configurations':len(valid),'rejected_configurations':len(cohort)-len(valid),'numerical_T':num['result']['T_count'],'shared_threshold_T':thr['result']['T_count'],'numerical_selected_strategy':num,'threshold_selected_strategy':thr}
  if valid:
   win=min(valid,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['exp_bits'],r['selected_strategy']['architecture']))
   raw=win['selected_strategy']['result']['T_count'];strong=min(raw,num['result']['T_count']);event=thr['result']['T_count']
   cr.update(selected_normal_bits=win['normal_bits'],selected_exp_bits=win['exp_bits'],selected_evidence_generation=win['evidence_generation'],selected_record_sha256=win['retained_record_sha256'],selected_raw_strategy=win['selected_strategy'],raw_T=raw,strongest_registered_rival_T=strong,strongest_registered_rival='raw_payoff' if raw<num['result']['T_count'] else 'numerical_conditional_price',event_saving_against_raw_exact=str(F(raw-event,raw)) if raw else None,event_saving_against_strongest_exact=str(F(strong-event,strong)) if strong else None,event_saving_against_strongest_per_cent_diagnostic=float(100*F(strong-event,strong)) if strong else None,exp_upper_boundary=win['exp_bits']==14)
  comparisons.append(cr)
 # Bind exact old/new records rather than duplicate all old profile and plan evidence.
 bindings=[{'model_id':r['model_id'],'candidate_id':r['candidate_id'],'normal_bits':r['normal_bits'],'exp_bits':r['exp_bits'],'evidence_generation':r['evidence_generation'],'record_sha256':r['retained_record_sha256'],'T_count':r.get('selected_strategy',{}).get('result',{}).get('T_count'),'admitted':bool(r.get('selected_strategy'))} for r in allrows]
 return bindings,comparisons
def run(args):
 p,raw,smooth=configured(args);protocol_path=Path(args.protocol).resolve() if args.protocol else p['output']/'expanded_raw_payoff_protocol.json';protocol=read(protocol_path);frozen=signatures(args,raw)
 for k,v in frozen.items():assert protocol[k]==v,(k,'Frozen input changed')
 assert protocol['new_narrow_exp_address_bits']==[11,12,13,14] and protocol['inner_normal_bits']==[8,10,12]
 factor=read(p['factor']);threshold=read(p['threshold']);numeric=read(p['numerical']);oldamended=read(p['old_smoothing']/'raw_payoff_smoothing_resources.json')
 norms=read(p['repo']/'results/gaussian_derivative_norms.json');nu0=F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==0 and r['passed']))
 records=[];amended=[];cache={};block_cache={};archive=[];seen=set();start=time.monotonic()
 for original,candidate in protocol['pairs']:
  prep=next(r for r in factor['preparations'] if r['original_model_id']==original and r['candidate_id']==candidate)
  old=next(r for r in factor['records'] if r['original_model_id']==original and r['factor_basis_candidate_id']==candidate and r['kind']=='positive_r1')
  for bn,be in product(protocol['inner_normal_bits'],protocol['new_narrow_exp_address_bits']):
   row={'original_model_id':original,'candidate_id':candidate,'normal_bits':bn,'payoff_exp_bits':be,'profiles':{},'strategies':[]};smooth_index={}
   try:
    g=raw.build_graph(prep['model'],prep['domain'],bn,be);tables,ordinates=raw.o.make_tables(g,cache,True);body=raw.c.certify_graph(g,tables,ordinates);cert=raw.probability_certificate(g,body,ordinates,nu0)
    row.update(model=prep['model'],graph_configuration=g.configuration,graph_semantic_nodes=g.semantic_nodes,graph_nodes=[{'name':n,'kind':v.kind,'args':v.args,'meta':v.meta.as_dict() if isinstance(v.meta,raw.o.ScalarSpec) else str(v.meta),'bits':v.bits,'output_bits':v.output_bits,'ideal_bounds':list(map(str,g.bounds[n]))} for n,v in g.nodes.items()],probability_certificate=cert)
    for n,cc in ordinates.items():
     key=(digest(cc),g.width)
     if key not in seen:
      seen.add(key);table=tables[raw.o.table_key(n,g.nodes[n])]
      archive.append({'certificate_sha256':key[0],'certificate':cc,'packed_table':table,'packed_table_sha256':digest(table),'graph_word_width':g.width})
     if cc.get('name')=='exp' and F(cc['lo'])==0 and F(cc['span'])==F(1,64):
      smooth_index[('exp',F(cc['lo']),F(cc['span']),cc['table_bits'],cc['q'])]=cc
    if cert['status']!='passed':row['rejection']={'reason':cert['reason']}
    else:
     eta=F(cert['probability_error_upper']);bud=copy.deepcopy(old['complete_deterministic_budget']);A=F(bud['external_factor_upper']);D=F(cert['inner_normal_CDF_error_upper']);det=F(bud['coefficient_radius_upper'])+F(bud['financial_transfer_upper'])+A*eta+F(bud['external_factor_radius'])
     bud.update(arithmetic_before_mass_upper=cert['arithmetic_probability_error_upper'],auxiliary_components={'inner_normal_payoff_difference_bias_upper':str(2*D),'bias_upper':str(2*D)},total_deterministic_undiscounted_price_error_upper=str(det))
     prev=old['coefficient_informed_prior'];lo,hi=F(prev['finite_ideal_mean_lower']),F(prev['finite_ideal_mean_upper']);prior={'probability_lower':str(max(F(0),lo-eta)),'probability_upper':str(min(F(1),hi+eta)),'native_positive':True,'finite_ideal_mean_lower':str(lo),'finite_ideal_mean_upper':str(hi),'previous_mean_information_sha256':digest(prev),'information_scope':'Identical price-independent finite ideal R1 mean interval, enlarged by full raw-payoff pointwise conditional error.'}
     row.update(complete_deterministic_budget=bud,coefficient_informed_prior=prior)
     # Emit every arithmetic-admitted circuit, even if pointwise budget is exhausted,
     # because the registered smoothing certificate may admit the same circuit.
     for mode in protocol['architectures']:
      try:
       profile=raw.o.counted_composition(g,tables,block_cache,mode);row['profiles'][mode]={'centred':profile}
       if det>=F(1,100):result={'kind':'deterministic_budget_exhausted','deterministic_error_upper':str(det)}
       else:result=raw.informed.informed_complete_plan('positive',bud,cert,prior,row['profiles'][mode])
      except (ArithmeticError,ValueError,AssertionError) as e:result={'kind':'production_failure','reason':type(e).__name__,'message':str(e)}
      row['strategies'].append({'architecture':mode,'result':result})
     valid=[r for r in row['strategies'] if r['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')]
     if valid:row['selected_strategy']=min(valid,key=lambda s:(s['result']['T_count'],s['architecture']))
     else:row['rejection']={'reason':'Pointwise budget exhausted or no valid emitted architecture plan'}
   except (ArithmeticError,ValueError,AssertionError) as e:row['rejection']={'reason':'Explicit production failure','type':type(e).__name__,'message':str(e)}
   records.append(row)
   if row.get('graph_nodes'):
    try:amend=make_amended(row,smooth,smooth_index,nu0,raw.informed,old)
    except (ArithmeticError,ValueError,AssertionError) as e:amend={'model_id':original,'candidate_id':candidate,'normal_bits':bn,'exp_bits':be,'pointwise_row_sha256':digest(row),'strategies':[],'rejection':{'reason':'Explicit smoothing production failure','type':type(e).__name__,'message':str(e)}}
   else:amend={'model_id':original,'candidate_id':candidate,'normal_bits':bn,'exp_bits':be,'pointwise_row_sha256':digest(row),'strategies':[],'rejection':copy.deepcopy(row['rejection'])}
   amended.append(amend)
   print(original,bn,be,'pointwise_T',row.get('selected_strategy',{}).get('result',{}).get('T_count'),'amended_T',amend.get('selected_strategy',{}).get('result',{}).get('T_count'),'reject',amend.get('rejection'),flush=True)
 archive_path=p['output']/'expanded_raw_payoff_scalar_ordinates.json.gz';write(archive_path,archive,True)
 bindings,comparisons=combined_summary(oldamended,amended,threshold,numeric)
 for k,v in signatures(args,raw).items():assert protocol[k]==v,(k,'Frozen input changed during execution')
 write(p['output']/'expanded_raw_payoff_resources.json.gz',{'schema':'qbp-expanded-raw-payoff-resources-v1','status':'completed retaining every pointwise and smoothing outcome','protocol_sha256':sha(protocol_path),**frozen,'new_configurations':48,'combined_configurations':96,'records':records,'amended_records':amended,'combined_pool_bindings':bindings,'selected_comparisons':comparisons,'scalar_archive_file':archive_path.name,'scalar_archive_sha256':sha(archive_path),'unique_scalar_certificates_and_wordwidths':len(archive),'emitted_block_profiles':len(block_cache),'generation_seconds_diagnostic':time.monotonic()-start,'scope':'Disclosed sequential extension, explicit finite registered pool and sufficient error/canonical QPE logical upper caps. Old exact evidence remains referenced by hash. No global precision optimum, hardware quantum prices, physical advantage or originality claim.'},True)
 for c in comparisons:print('Combined96',c['model_id'],'winner',c.get('selected_normal_bits'),c.get('selected_exp_bits'),'raw_T',c.get('raw_T'),'event_saving_per_cent',c.get('event_saving_against_strongest_per_cent_diagnostic'),'boundary',c.get('exp_upper_boundary'),flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('mode',choices=('register','run','replay'));parser.add_argument('--repo',required=True);parser.add_argument('--old-raw',required=True);parser.add_argument('--old-smoothing',required=True);parser.add_argument('--proof',required=True);parser.add_argument('--output',required=True);parser.add_argument('--protocol');args=parser.parse_args()
 if args.mode=='register':
  if args.protocol:raise ValueError('Registration protocol path is output directory')
  register(args)
 else:
  if args.mode=='replay' and not args.protocol:raise ValueError('Replay requires exact previously registered --protocol')
  run(args)
