"""Disclosed all-pool deterministic fallback and worst-history upper caps.

This re-optimises a protected modern-controller guaranteed upper cap across
every already retained financial preparation and both architectures. It does
not simulate unknown finance amplitudes or equate a cap with realised cost.
"""
from pathlib import Path
from fractions import Fraction as F
from dataclasses import asdict
from datetime import datetime, timezone
import argparse,gzip,hashlib,importlib.util,json,sys
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p):
 p=Path(p);return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def write(p,v):
 p=Path(p)
 if p.exists():raise FileExistsError('Preserve previous panel: '+str(p))
 p.write_text(json.dumps(v,indent=2,default=str)+'\n')
def frozen(args):
 repo,out=Path(args.repo),Path(args.results)
 return {'source_sha256':sha(__file__),'common_shift_resources_sha256':sha(out/'common_shift_resources.json.gz'),
  'common_shift_protocol_sha256':sha(out/'common_shift_protocol.json'),
  'compiled_plan_and_charge_source_sha256':sha(repo/'integral_pricing_budget.py'),
  'existing_independent_formula_source_sha256':sha(repo/'check_expanded_precision_resources.py'),
  'controller_source_sha256':sha(repo/'chebyshev_controller.py'),
  'controller_resource_policy_sha256':sha(repo/'chebyshev_resources.py'),
  'qft_tapes_sha256':sha(repo/'results/qft_rotation_tapes.json'),
  'qft_certificates_sha256':sha(repo/'results/qft_rotation_certificates.json')}
def register(args):
 out=Path(args.results);d=read(out/'common_shift_resources.json.gz')
 bindings=[]
 for i,r in enumerate(d['records']):
  bindings.append({'record_index':i,'record_sha256':digest(r),'method':r['method'],'model_id':r['model_id'],
   'source_path':r['source_path'],'source_row_index':r['source_row_index'],'source_row_sha256':r['source_row_sha256'],
   'profiles_sha256':r['profiles_sha256'],'amended_strategies_sha256':digest(r['strategies']),
   'previous_minimum_strategy_sha256':digest(r.get('selected_strategy'))})
 write(out/'common_shift_poolwise_fallback_protocol.json',{'schema':'qbp-common-shift-poolwise-fallback-protocol-v1',
  'registration_utc':datetime.now(timezone.utc).isoformat(),**frozen(args),
  'registered_before_all_pool_half_alpha_fallback_outcomes':True,
  'known_outcomes':'All common-shift canonical-QPE results and the separately generated selected-preparation modern-controller diagnostics are already known. This is a disclosed sequential deterministic all-pool fallback panel, not a blind initial registration.',
  'record_bindings':bindings,'configurations':{'raw_payoff':348,'numerical_conditional_price':64,'shared_threshold':4},
  'architectures':['exact7t','unitary8t'],'complete_price_alpha':'1/3600','fallback_alpha':'1/7200',
  'fixed_fallback_rule':'For every admitted amended method/architecture, rebuild the canonical QPE fallback with precisely the same probability tolerance, upper prior and frozen clean financial profile, changing only alpha to1/7200. Retain exact Fourier and full circuit/reflection/cleanup/reset costs.',
  'fixed_worst_history_cap_rule':'T_worst=T_fallback+floor(T_fallback/10), using the protected controller pilot-T budget and fresh half-alpha fallback. This is a certified resource upper cap, not a claim that the cap is attained.',
  'selection':'Minimise the derived worst-history upper cap separately for each method/case over every retained admitted amended configuration and both architectures. Retain all previous and fresh rejections. Tie-break by normal bits, PG bits, PA bits, architecture and record index.',
  'structural_test':'Record whether M, Fourier profile and per-repetition charge stay unchanged, and repetitions25 become27, for every admitted architecture. Report any failure without adapting the pool or rule.',
  'independent_check':'Replay each fresh QFT/QPE error, odd binomial repetition minimum and complete charge with existing independent own-formula functions, setting their explicitly documented alpha to1/7200. The amendment producer and confidence frontend are excluded from this check.',
  'scope':'Pool-wise protected-controller guaranteed logical T upper caps for the frozen preparations. No financial means, history sampling, query-law novelty, globally optimal modern estimator, physical runtime or practical quantum advantage.'})
 print('Registered all416 settings/two architectures before808 fallback outcomes',flush=True)
def run(args):
 repo,out=Path(args.repo),Path(args.results);p=read(out/'common_shift_poolwise_fallback_protocol.json');d=read(out/'common_shift_resources.json.gz')
 assert all(p[k]==v for k,v in frozen(args).items())
 sys.path.insert(0,str(repo));from integral_pricing_budget import compiled_extended_plan,charge_full_oracle
 s=importlib.util.spec_from_file_location('independent_alpha_half',repo/'check_expanded_precision_resources.py');own=importlib.util.module_from_spec(s);s.loader.exec_module(own);own.ALPHA=F(1,7200)
 rotations=read(repo/'results/qft_rotation_tapes.json')['rotations'];bounds={r['key']:F(r['operator_norm_upper_rational']) for r in read(repo/'results/qft_rotation_certificates.json')['certificates']}
 sources={};rows=[];rejections=[];checks=0
 for row,binding in zip(d['records'],p['record_bindings']):
  assert digest(row)==binding['record_sha256']
  path=row['source_path']
  if path not in sources:sources[path]=read(repo/'results'/path)
  base=sources[path]['rows' if row['method']=='shared_threshold' else 'records'][row['source_row_index']]
  assert digest(base)==row['source_row_sha256'] and digest(base['profiles'])==row['profiles_sha256']
  if not row['strategies']:
   rejections.append({**binding,'rejection':row['rejection']});continue
  for strategy in row['strategies']:
   current=strategy['result']
   if current['kind']!='complete_quantum_resource_plan':
    rejections.append({**binding,'architecture':strategy['architecture'],'retained_non_QPE_outcome':current});continue
   flag=current['flags'][0];profile=base['profiles'][strategy['architecture']]['centred'];assert digest(profile)==flag['financial_profile_sha256']
   target,prior=F(flag['probability_tolerance']),F(flag['prior_upper'])
   setup,qft=compiled_extended_plan(target,F(1,7200),prior)
   charge=charge_full_oracle(setup,profile,qft);plan=asdict(setup)
   stat,failure=own.own_plan(plan,target,prior,qft,rotations,bounds);own.own_charge(charge,profile,qft,plan)
   factor=F(row['complete_deterministic_budget']['external_factor_upper']);det=F(row['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper'])
   assert det+factor*stat<=F(1,100) and failure<=F(1,7200)
   full,fallback=current['T_count'],charge['metrics']['T'];worst=fallback+fallback//10
   rec={**binding,'architecture':strategy['architecture'],'configuration':{k:row[k] for k in ('normal_bits','PG_exp_bits','PA_exp_bits','PG_bits','PA_bits') if k in row},
    'financial_profile_sha256':digest(profile),'probability_tolerance':str(target),'prior_upper':str(prior),
    'T_full_alpha':full,'T_half_alpha_fallback':fallback,'T_worst_history_upper':worst,
    'pilot_T_budget':fallback//10,'M_full':flag['plan']['M'],'M_half':plan['M'],
    'r_full':flag['plan']['repetitions'],'r_half':plan['repetitions'],
    'QFT_unchanged':digest(flag['fourier_profile'])==digest(qft),
    'per_repetition_charge_unchanged':charge['per_repetition_metrics']==flag['complete_charge']['per_repetition_metrics'],
    'full_plan_sha256':digest(flag['plan']),'half_plan_sha256':digest(plan),
    'half_plan':{k:str(v) if isinstance(v,F) else v for k,v in plan.items()},
    'half_fourier_profile_sha256':digest(qft),'half_full_charge_sha256':digest(charge),
    'half_failure_probability_upper':str(failure),'complete_monetary_error_upper':str(det+factor*stat),
    'frozen_profile_reused':True,'own_formula_replay_passed':True}
   rows.append(rec);checks+=1
 selected=[]
 for mid in dict.fromkeys(r['model_id'] for r in d['records']):
  methods={}
  for method in ('shared_threshold','raw_payoff','numerical_conditional_price'):
   group=[r for r in rows if r['model_id']==mid and r['method']==method]
   win=min(group,key=lambda r:(r['T_worst_history_upper'],r['configuration'].get('normal_bits',0),r['configuration'].get('PG_exp_bits',r['configuration'].get('PG_bits',0)),r['configuration'].get('PA_exp_bits',r['configuration'].get('PA_bits',0)),r['architecture'],r['record_index']))
   methods[method]={'selected_record_sha256':digest(win),'record_index':win['record_index'],'architecture':win['architecture'],'configuration':win['configuration'],'T_full_alpha':win['T_full_alpha'],'T_half_alpha_fallback':win['T_half_alpha_fallback'],'T_worst_history_upper':win['T_worst_history_upper'],'financial_profile_sha256':win['financial_profile_sha256']}
  rivals={m:methods[m]['T_worst_history_upper'] for m in ('raw_payoff','numerical_conditional_price')};strong=min(rivals,key=lambda m:(rivals[m],m));T=rivals[strong];E=methods['shared_threshold']['T_worst_history_upper']
  selected.append({'model_id':mid,'methods':methods,'strongest_rival':strong,'saving_against_strongest_exact':str(F(T-E,T)),'saving_against_strongest_per_cent_diagnostic':float(100*F(T-E,T))})
 assert all(p[k]==v for k,v in frozen(args).items())
 invariance={'all_M_unchanged':all(r['M_full']==r['M_half'] for r in rows),'all_QFT_unchanged':all(r['QFT_unchanged'] for r in rows),'all_per_repetition_charge_unchanged':all(r['per_repetition_charge_unchanged'] for r in rows),'all_repetitions25_to27':all((r['r_full'],r['r_half'])==(25,27) for r in rows),'all_half_cost_equals27_over25_full':all(F(r['T_half_alpha_fallback'])==F(27,25)*r['T_full_alpha'] for r in rows)}
 result={'schema':'qbp-common-shift-poolwise-fallback-v1','status':'passed own-formula full-pool replay','protocol_sha256':sha(out/'common_shift_poolwise_fallback_protocol.json'),**frozen(args),'records':rows,'rejections':rejections,'admitted_architecture_plans':checks,'rejected_configurations':len(rejections),'invariance':invariance,'selected_comparisons':selected,'scope':p['scope']}
 write(out/'common_shift_poolwise_fallback.json',result)
 print('Own-formula fallbacks',checks,'rejects',len(rejections),'invariance',invariance,flush=True)
 for x in selected:print(x['model_id'],{m:v['T_worst_history_upper'] for m,v in x['methods'].items()},x['saving_against_strongest_per_cent_diagnostic'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run'));p.add_argument('--repo',required=True);p.add_argument('--results',required=True);a=p.parse_args();register(a) if a.mode=='register' else run(a)
