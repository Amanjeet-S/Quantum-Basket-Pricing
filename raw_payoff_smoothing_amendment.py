"""Proof-derived raw-payoff certificate amendment, preserving every circuit.

The registered finite pool and literal profiles are unchanged. Earlier
pointwise outcomes are declared as already known. The smooth conditional
price derivative controls common outer shifts; the small error in the inner
coefficient remains charged with the raw payoff slope.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,copy,gzip,hashlib,json,sys,time
sys.dont_write_bytecode=True
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
digest=lambda v:hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p):
 p=Path(p);return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def write(p,v):
 if p.exists():raise FileExistsError('Preserve previously retained evidence')
 p.write_text(json.dumps(v,indent=2,default=str)+'\n')
def coefficient_integer(x,q):
 x=F(x)
 if x.denominator&(x.denominator-1)==0:return x.numerator,x.denominator.bit_length()-1
 y=x*(1<<q)+F(1,2);return y.numerator//y.denominator,q

def mapping(row,archive,nu0):
 graph={n['name']:n for n in row['graph_nodes']};sem=row['graph_semantic_nodes'];q=row['graph_configuration']['q'];L=F(row['graph_configuration']['inner_normal_L']);lam=F(row['model']['lambda']);cert=row['probability_certificate'];body={n['node']:n for n in cert['body_arithmetic_certificate']['nodes']}
 lower=graph[sem['lower_endpoint']];assert lower['kind']=='add'
 inner=graph[lower['args'][1]];assert inner['kind']=='scale' and inner['args']==[sem['normal_plus_L']] and F(inner['meta'])==lam
 offset_add=graph[lower['args'][0]];assert offset_add['kind']=='add' and offset_add['args'][0]==sem['mG_offset']
 constant=graph[offset_add['args'][1]];assert constant['kind']=='constant'
 offset=graph[sem['common_offset']];assert offset['kind']=='constant'
 errors=lambda name:F(body[name]['absolute_output_error_upper'])
 em=errors(sem['mG_offset']);edelta=errors(sem['logT']);ec=errors(constant['name']);eo=errors(offset['name']);D=F(cert['inner_normal_CDF_error_upper'])
 integer,drop=coefficient_integer(lam,q);coefficient_error=abs(F(integer,1<<drop)-lam);floor=F(1,1<<q) if drop else F(0);Es=2*L*coefficient_error+floor
 exp_errors={}
 for side in ('PG','PA'):
  n=graph[sem[side+'_exp']];meta=n['meta'];assert n['kind']=='function' and meta['name']=='exp' and F(meta['lo'])==0 and F(meta['span'])==F(1,64)
  key=('exp',F(meta['lo']),F(meta['span']),n['bits'],q);ordinate=archive[key]
  exp_errors[side]=F(ordinate['output_error_upper'])
  assert F(meta['output_bounds'][0])==1 and F(meta['output_bounds'][1])==F(64,63)
 eexp=100*sum(exp_errors.values());scale_term=202*Es;outer_term=nu0/lam*(2*(em+ec+eo)+edelta);smooth_arith=eexp+scale_term+outer_term
 pointwise=F(cert['probability_error_upper']);smooth=min(F(1),smooth_arith+2*D);eta=min(pointwise,smooth)
 return {'outer_m_error_upper':str(em),'outer_logT_error_upper':str(edelta),'shared_shift_constant_error_upper':str(ec),'shared_offset_constant_error_upper':str(eo),'inner_scale_node':inner['name'],'outer_constant_node':constant['name'],'implemented_inner_coefficient':str(F(integer,1<<drop)),'inner_coefficient_error_upper':str(coefficient_error),'inner_scale_floor_allowance':str(floor),'inner_scale_error_upper':str(Es),'PG_intrinsic_exp_approximation_upper':str(exp_errors['PG']),'PA_intrinsic_exp_approximation_upper':str(exp_errors['PA']),'raw_exp_output_contribution_upper':str(eexp),'inner_scale_payoff_contribution_upper':str(scale_term),'smooth_outer_shift_contribution_upper':str(outer_term),'smoothing_arithmetic_error_upper':str(smooth_arith),'inner_normal_CDF_error_upper':str(D),'pointwise_probability_error_upper':str(pointwise),'smoothing_probability_error_upper':str(smooth),'selected_probability_error_upper':str(eta)}

def signatures(args):
 repo=Path(args.repo)
 return {'amendment_source_sha256':sha(__file__),'base_resources_sha256':sha(args.resources),'base_archive_sha256':sha(args.archive),'checked_smoothing_proof_sha256':sha(args.proof),'retained_independent_resource_review_sha256':sha(args.independent_review),'numerical_resources_sha256':sha(repo/'results/strengthened_r1_resources.json'),'threshold_resources_sha256':sha(repo/'results/shared_threshold_resources.json'),'mean_frontend_sha256':sha(repo/'correction_information_priors.py'),'confidence_frontend_sha256':sha(repo/'boundary_precision_selection.py'),'budget_charge_sha256':sha(repo/'integral_pricing_budget.py'),'gaussian_norms_sha256':sha(repo/'results/gaussian_derivative_norms.json')}

def register(args):
 resources=read(args.resources);assert resources['schema']=='qbp-raw-payoff-resources-tight-v2'
 independent=read(args.independent_review);assert independent['status']=='passed' and independent['resource_sha256']==sha(args.resources)
 output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
 write(output/'raw_payoff_smoothing_protocol.json',{'schema':'qbp-raw-payoff-smoothing-protocol-v1','registration_utc':datetime.now(timezone.utc).isoformat(),'research_record_date':'2026-10-10',**signatures(args),'registered_before_fresh_amended_costs':True,'pointwise_only_costs_already_known':True,'fixed_circuit_pool':[{'model_id':r['original_model_id'],'candidate_id':r['candidate_id'],'normal_bits':r['normal_bits'],'exp_bits':r['payoff_exp_bits'],'graph_sha256':digest(r['graph_nodes']),'configuration_sha256':digest(r['graph_configuration']),'profile_sha256':digest(r['profiles'])} for r in resources['records']],'configurations':48,'architectures':['exact7t','unitary8t'],'amendment_formula':'eta=min(pointwise_eta,min(1,100(e_exp_G+e_exp_A)+202 E_s+2D_N+(nu0/lambda)(2(e_m+e_C+e_O)+e_Delta)))','independent_coefficient_error_rule':'E_s=2L*abs(lambda_hat-lambda)+2^-q*indicator(drop>0), lambda_hat is the existing exact compiled fixed-point scaling coefficient. Normal quantile approximation appears only through the retained D_N, avoiding duplicate pointwise quantile perturbation in E_s.','intrinsic_exp_error_rule':'Use scalar certificate output_error_upper at PG_exp and PA_exp, excluding graph argument perturbations; both are charged.','outer_error_rule':'Read mG_offset and logT absolute_output_error_upper from unchanged body graph certificate. Read rounding errors of J-BG and J+lambda L constants separately from unchanged body graph certificate.','unchanged':['All financial models and finite laws','All 48 integer graphs, scalar tables, clean preparations and 96 emitted profiles','Existing smooth outer transfer and coefficient intervals','Same price-independent finite ideal mean intervals','Canonical QPE confidence alpha1/3600, monetary radius1/100 and architecture policy'],'selection':'Minimum complete canonical QPE T cap over exactly the retained 48 configurations and two architectures; retain all rejected outcomes. Compare against the minimum of this amended raw rival and the strengthened numerical rival. Negative or zero event saving is retained.','data_origin':'Synthetic no-arbitrage models and certified scalar data only; no licensed observations, physical QPE or price answers.'})

def run(args):
 repo=Path(args.repo);output=Path(args.output);protocol=read(output/'raw_payoff_smoothing_protocol.json');frozen=signatures(args)
 for k,v in frozen.items():assert protocol[k]==v
 sys.path.insert(0,str(repo));import correction_information_priors as frontend
 resources=read(args.resources);ords=read(args.archive);archive={}
 for r in ords:
  c=r['certificate']
  if c.get('name')=='exp' and F(c['lo'])==0 and F(c['span'])==F(1,64):
   key=('exp',F(c['lo']),F(c['span']),c['table_bits'],c['q'])
   if key in archive:assert archive[key]==c
   archive[key]=c
 norms=read(repo/'results/gaussian_derivative_norms.json');nu0=F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==0 and r['passed']))
 amended=[];start=time.monotonic()
 for base,binding in zip(resources['records'],protocol['fixed_circuit_pool']):
  assert digest(base['graph_nodes'])==binding['graph_sha256'] and digest(base['profiles'])==binding['profile_sha256'] and digest(base['graph_configuration'])==binding['configuration_sha256']
  term=mapping(base,archive,nu0);eta=F(term['selected_probability_error_upper']);bud=copy.deepcopy(base['complete_deterministic_budget']);old=base['coefficient_informed_prior'];A=F(bud['external_factor_upper']);det=F(bud['coefficient_radius_upper'])+F(bud['financial_transfer_upper'])+A*eta+F(bud['external_factor_radius']);D=F(term['inner_normal_CDF_error_upper']);ear=min(F(base['probability_certificate']['arithmetic_probability_error_upper']),F(term['smoothing_arithmetic_error_upper']))
  bud.update(arithmetic_before_mass_upper=str(ear),total_deterministic_undiscounted_price_error_upper=str(det),auxiliary_components={'bias_upper':str(2*D),'inner_normal_payoff_difference_bias_upper':str(2*D)},certificate_amendment='smoothing_first_uniform_conditional')
  prior={'probability_lower':str(max(F(0),F(old['finite_ideal_mean_lower'])-eta)),'probability_upper':str(min(F(1),F(old['finite_ideal_mean_upper'])+eta)),'native_positive':True,'finite_ideal_mean_lower':old['finite_ideal_mean_lower'],'finite_ideal_mean_upper':old['finite_ideal_mean_upper'],'previous_mean_information_sha256':old['previous_mean_information_sha256'],'information_scope':'Unchanged price-independent finite ideal R1 mean interval, enlarged by minimum of checked pointwise and smoothing-first conditional error bounds.'}
  cert={'status':'passed','configuration':base['graph_configuration'],'probability_error_upper':str(eta),'monetary_value_error_upper':str(eta),'probability_actual_bounds':base['probability_certificate']['probability_actual_bounds'],'retained_pointwise_certificate_sha256':digest(base['probability_certificate']),'smoothing_derivation':term,'scope':'Uniform conditional inner-law expectation error certificate; original pointwise body certificate remains unchanged and separately retained.'}
  row={**binding,'retained_pointwise_row_sha256':digest(base),'probability_certificate':cert,'complete_deterministic_budget':bud,'coefficient_informed_prior':prior,'strategies':[]}
  if det>=F(1,100):row['rejection']={'reason':'Amended deterministic monetary tolerance exhausted','deterministic_error_upper':str(det)}
  else:
   for mode in protocol['architectures']:
    try:result=frontend.informed_complete_plan('positive',bud,cert,prior,base['profiles'][mode])
    except (ArithmeticError,ValueError,AssertionError) as e:result={'kind':'production_failure','reason':type(e).__name__,'message':str(e)}
    row['strategies'].append({'architecture':mode,'result':result})
   valid=[s for s in row['strategies'] if s['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')]
   if valid:row['selected_strategy']=min(valid,key=lambda s:(s['result']['T_count'],s['architecture']))
   else:row['rejection']={'reason':'No valid amended architecture plan'}
  amended.append(row)
 numerical=read(repo/'results/strengthened_r1_resources.json');threshold=read(repo/'results/shared_threshold_resources.json');comparisons=[]
 for mid in dict.fromkeys(r['model_id'] for r in amended):
  group=[r for r in amended if r['model_id']==mid];valid=[r for r in group if r.get('selected_strategy')];nr=next(r for r in numerical['selected_comparisons'] if r['original_model_id']==mid);tr=next(r for r in threshold['rows'] if r['original_model_id']==mid);num=nr['selected_numerical_strategy']['result']['T_count'];thr=tr['selected_strategy']['result']['T_count']
  cr={'model_id':mid,'admitted_configurations':len(valid),'rejected_configurations':len(group)-len(valid),'numerical_T':num,'shared_threshold_T':thr,'numerical_selected_strategy':nr['selected_numerical_strategy'],'threshold_selected_strategy':tr['selected_strategy']}
  if valid:
   win=min(valid,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['exp_bits']));raw=win['selected_strategy']['result']['T_count'];strongest=min(num,raw)
   cr.update(selected_normal_bits=win['normal_bits'],selected_exp_bits=win['exp_bits'],selected_raw_strategy=win['selected_strategy'],raw_T=raw,strongest_registered_rival_T=strongest,strongest_registered_rival='raw_payoff' if raw<num else 'numerical_conditional_price',event_saving_against_raw_exact=str(F(raw-thr,raw)) if raw else None,event_saving_against_strongest_exact=str(F(strongest-thr,strongest)) if strongest else None,event_saving_against_strongest_per_cent_diagnostic=float(100*F(strongest-thr,strongest)) if strongest else None)
  comparisons.append(cr)
 for k,v in signatures(args).items():assert protocol[k]==v
 write(output/'raw_payoff_smoothing_resources.json',{'schema':'qbp-raw-payoff-smoothing-resources-v1','status':'completed retaining every amended outcome','protocol_sha256':sha(output/'raw_payoff_smoothing_protocol.json'),**frozen,'configurations':48,'records':amended,'selected_comparisons':comparisons,'generation_seconds_diagnostic':time.monotonic()-start,'scope':'Amended sufficient error bounds and canonical QPE upper caps only. The unchanged circuits and original plans are retained by hash reference. No global optimum, hardware, quantum prices, physical advantage or originality claim.'})
 for r in comparisons:print(r['model_id'],'raw_T',r.get('raw_T'),'event_saving_per_cent',r.get('event_saving_against_strongest_per_cent_diagnostic'),flush=True)

def audit(args):
 from collections import Counter
 from math import comb,isqrt
 import importlib.util
 repo=Path(args.repo);output=Path(args.output);protocol=read(output/'raw_payoff_smoothing_protocol.json');data=read(output/'raw_payoff_smoothing_resources.json');original=read(args.resources);ords=read(args.archive)
 assert data['protocol_sha256']==sha(output/'raw_payoff_smoothing_protocol.json')
 for k,v in signatures(args).items():assert protocol[k]==data[k]==v
 spec=importlib.util.spec_from_file_location('existing_independent_confidence',repo/'checks/audit_boundary_precision.py');independent=importlib.util.module_from_spec(spec);spec.loader.exec_module(independent)
 archive={}
 for a in ords:
  c=a['certificate']
  if c.get('name')=='exp' and F(c['lo'])==0 and F(c['span'])==F(1,64):
   q=c['q'];ulp=F(1,1<<q);N=1<<c['table_bits'];nodal=max(max(abs(F(n['integer'],1<<q)-F(n['lower'])),abs(F(n['integer'],1<<q)-F(n['upper']))) for n in c['nodes']);interp=F(c['second_derivative_upper'])*(F(c['span'])/N)**2/8
   assert F(c['nodal_error_upper'])==nodal and F(c['linear_interpolation_error_upper'])==interp and F(c['output_error_upper'])==nodal+interp+ulp
   key=('exp',F(c['lo']),F(c['span']),c['table_bits'],q)
   if key in archive:assert archive[key]==c
   archive[key]=c
 norms=read(repo/'results/gaussian_derivative_norms.json');nu0=F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==0 and r['passed']))
 tapes=read(repo/'results/qft_rotation_tapes.json')['rotations'];bounds={r['key']:F(r['operator_norm_upper_rational']) for r in read(repo/'results/qft_rotation_certificates.json')['certificates']}
 plans=deterministic=rejected=0
 assert len(data['records'])==48
 for r,old,binding in zip(data['records'],original['records'],protocol['fixed_circuit_pool']):
  assert all(r[k]==binding[k] for k in binding)
  assert r['retained_pointwise_row_sha256']==digest(old) and r['profile_sha256']==digest(old['profiles']) and r['graph_sha256']==digest(old['graph_nodes'])
  term=mapping(old,archive,nu0);assert term==r['probability_certificate']['smoothing_derivation'];eta=F(term['selected_probability_error_upper']);budget=r['complete_deterministic_budget'];oldbudget=old['complete_deterministic_budget']
  for k in ('coefficient_radius_upper','financial_transfer_upper','external_factor_upper','external_factor_radius','mass_interval','coefficient_interval','external_factor_midpoint'):assert budget[k]==oldbudget[k]
  A=F(budget['external_factor_upper']);det=F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*eta+F(budget['external_factor_radius']);assert det==F(budget['total_deterministic_undiscounted_price_error_upper'])
  pr=r['coefficient_informed_prior'];opr=old['coefficient_informed_prior'];assert F(pr['probability_lower'])==max(F(0),F(opr['finite_ideal_mean_lower'])-eta) and F(pr['probability_upper'])==min(F(1),F(opr['finite_ideal_mean_upper'])+eta) and pr['previous_mean_information_sha256']==opr['previous_mean_information_sha256']
  if r.get('rejection'):
   assert det>=F(1,100) or all(x['result']['kind'] not in ('complete_quantum_resource_plan','deterministic_correction_prior') for x in r['strategies']);rejected+=1;continue
  for strategy in r['strategies']:
   result=strategy['result'];profile=old['profiles'][strategy['architecture']]['centred'];upper=F(pr['probability_upper']);remaining=F(1,100)-det
   if result['kind']=='deterministic_correction_prior':
    radius=A*(upper-F(pr['probability_lower']))/2;assert result['T_count']==0 and F(result['total_price_error_upper'])==det+radius<=F(1,100);deterministic+=1;continue
   assert result['kind']=='complete_quantum_resource_plan' and len(result['flags'])==1
   flag=result['flags'][0];plan=flag['plan'];qft=flag['fourier_profile'];independent.check_plan(plan,remaining/A,F(1,3600),upper,qft,tapes,bounds);independent.charge_check(flag['complete_charge'],profile,qft,plan)
   assert flag['financial_profile_sha256']==digest(profile)
   stat=A*(F(plan['amplitude_error_upper'])+F(plan['postprocess_allowance']));assert F(result['statistical_price_error_upper'])==stat and F(result['total_price_error_upper'])==det+stat<=F(1,100)
   assert F(result['failure_probability_upper'])==F(plan['binomial_failure_upper'])<=F(1,3600) and result['T_count']==flag['complete_charge']['metrics']['T'];plans+=1
  valid=[x for x in r['strategies'] if x['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')];assert r['selected_strategy']==min(valid,key=lambda x:(x['result']['T_count'],x['architecture']))
 for comp in data['selected_comparisons']:
  valid=[r for r in data['records'] if r['model_id']==comp['model_id'] and r.get('selected_strategy')];win=min(valid,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['exp_bits']));assert win['selected_strategy']==comp['selected_raw_strategy'] and win['normal_bits']==comp['selected_normal_bits'] and win['exp_bits']==comp['selected_exp_bits']
  raw=win['selected_strategy']['result']['T_count'];num=comp['numerical_selected_strategy']['result']['T_count'];event=comp['threshold_selected_strategy']['result']['T_count'];strong=min(raw,num)
  assert comp['raw_T']==raw and comp['strongest_registered_rival_T']==strong and F(comp['event_saving_against_strongest_exact'])==F(strong-event,strong)
 write(output/'raw_payoff_smoothing_audit.json',{'schema':'qbp-raw-payoff-smoothing-audit-v1','status':'passed','source_sha256':sha(__file__),'resources_sha256':sha(output/'raw_payoff_smoothing_resources.json'),'protocol_sha256':sha(output/'raw_payoff_smoothing_protocol.json'),'prior_independent_resource_review_sha256':sha(args.independent_review),'existing_own_formula_confidence_auditor_sha256':sha(repo/'checks/audit_boundary_precision.py'),'fixed_circuit_bindings':48,'amended_QPE_and_full_charge_plans':plans,'deterministic_prior_plans':deterministic,'rejected_configurations':rejected,'scope':'Direct retained-field smoothing certificate recomputation, exact scalar error decomposition, unchanged graph/profile bindings, same independent own-formula QPE and complete charge replay. Mathematical justification is the separately checked supporting proof. Original pointwise plans and whole circuit replay remain retained.'})
 print('Amendment audit passed',plans,'QPE/full-charge plans',deterministic,'deterministic plans',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run','audit'));p.add_argument('--repo',required=True);p.add_argument('--raw-results',required=True);p.add_argument('--proof',required=True);p.add_argument('--independent-review',required=True);p.add_argument('--output',required=True);a=p.parse_args();a.resources=str(Path(a.raw_results)/'raw_payoff_resources.json.gz');a.archive=str(Path(a.raw_results)/'raw_payoff_scalar_ordinates.json.gz')
 register(a) if a.mode=='register' else run(a) if a.mode=='run' else audit(a)
