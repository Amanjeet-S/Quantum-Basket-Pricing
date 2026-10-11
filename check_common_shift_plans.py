"""Own-formula complete-plan replay, excluding the amendment producer.

The existing independent resource library constructs QFT/QPE errors and full
clean financial-oracle/reflection charges. Primitive profiles are frozen hash
references, rather than freshly emitted or treated as free ideal queries.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, gzip, hashlib, importlib.util, json, sys
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def read(p):
 p=Path(p);return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def main(args):
 repo=Path(args.repo);out=Path(args.results);data=read(out/'common_shift_resources.json.gz');protocol=read(out/'common_shift_protocol.json')
 assert data['protocol_sha256']==sha(out/'common_shift_protocol.json')
 assert data['uncosted_admissions']==0 and data['status']=='completed retaining every amended outcome'
 for p,h in protocol['resources_sha256'].items():assert sha(repo/'results'/p)==h
 for p,h in protocol['scalar_archives_sha256'].items():assert sha(repo/'results'/p)==h
 for p,h in protocol['dependencies_sha256'].items():assert sha(repo/p)==h
 for p,h in protocol['inputs_sha256'].items():assert sha(repo/'results'/p)==h
 spec=importlib.util.spec_from_file_location('own_formula',repo/'check_expanded_precision_resources.py');own=importlib.util.module_from_spec(spec);spec.loader.exec_module(own)
 rotations=read(repo/'results/qft_rotation_tapes.json')['rotations'];bounds={r['key']:F(r['operator_norm_upper_rational']) for r in read(repo/'results/qft_rotation_certificates.json')['certificates']}
 sources={p:read(repo/'results'/p) for p in protocol['resources_sha256']}
 plans=zeros=rejects=0
 for row,binding in zip(data['records'],protocol['fixed_pool']):
  assert all(row[k]==v for k,v in binding.items())
  source=sources[row['source_path']];base=source['rows' if row['method']=='shared_threshold' else 'records'][row['source_row_index']]
  assert digest(base)==row['source_row_sha256'] and digest(base['profiles'])==row['profiles_sha256']
  cert,prior,budget=row['probability_certificate'],row['coefficient_informed_prior'],row['complete_deterministic_budget']
  old=base['complete_deterministic_budget'];oldprior=base['coefficient_informed_prior']
  for k in ('coefficient_interval','coefficient_radius_upper','mass_interval','external_factor_midpoint','external_factor_upper','external_factor_radius','financial_transfer_upper','monetary_probability_decode_factor'):assert budget[k]==old[k]
  assert F(cert['configuration']['monetary_signal_normaliser'])==F(cert['configuration']['cap'])==1 and F(budget['monetary_probability_decode_factor'])==1
  eta,A=F(cert['probability_error_upper']),F(budget['external_factor_upper'])
  expected=min(*(F(v) for v in cert['previous_probability_error_bounds']),F(cert['common_shift_derivation']['new_probability_error_upper']))
  assert eta==expected<=1 and eta>=0
  ar,bias=F(budget['arithmetic_before_mass_upper']),F(budget['auxiliary_components']['bias_upper']);assert ar+bias>=eta
  lo,hi=F(oldprior['finite_ideal_mean_lower']),F(oldprior['finite_ideal_mean_upper'])
  assert F(prior['finite_ideal_mean_lower'])==lo and F(prior['finite_ideal_mean_upper'])==hi
  assert F(prior['probability_lower'])==max(F(0),lo-eta) and F(prior['probability_upper'])==min(F(1),hi+eta)
  assert prior['retained_prior_sha256']==digest(oldprior) and prior['native_positive']
  det=F(old['coefficient_radius_upper'])+F(old['financial_transfer_upper'])+A*eta+F(old['external_factor_radius'])
  assert F(budget['total_deterministic_undiscounted_price_error_upper'])==det
  for s in row['strategies']:
   r=s['result']
   if r['kind']=='production_failure':raise AssertionError('Complete amendment plan failure')
   if r['kind']=='deterministic_correction_prior':
    radius=A*(F(prior['probability_upper'])-F(prior['probability_lower']))/2
    assert r['T_count']==0 and F(r['total_price_error_upper'])==det+radius<=F(1,100);zeros+=1;continue
   assert r['kind']=='complete_quantum_resource_plan' and len(r['flags'])==1 and det<F(1,100)
   flag=r['flags'][0];profile=base['profiles'][s['architecture']]['centred']
   assert flag['financial_profile_sha256']==digest(profile)
   target,upper=(F(1,100)-det)/A,F(prior['probability_upper'])
   assert F(flag['probability_tolerance'])==target and F(flag['prior_upper'])==upper
   statistical,failure=own.own_plan(flag['plan'],target,upper,flag['fourier_profile'],rotations,bounds)
   own.own_charge(flag['complete_charge'],profile,flag['fourier_profile'],flag['plan'])
   assert F(r['statistical_price_error_upper'])==A*statistical and F(r['total_price_error_upper'])==det+A*statistical<=F(1,100)
   assert F(r['failure_probability_upper'])==failure<=F(1,3600)
   assert r['T_count']==flag['complete_charge']['metrics']['T']
   assert r['quantum_shots']==flag['plan']['repetitions'] and r['logical_qubits_upper']==flag['complete_charge']['peak_logical_qubits'];plans+=1
  candidates=[{'certificate_policy':'previous_pointwise',**s} for s in base['strategies']]
  if row.get('amended_source_path'):
   source=sources[row['amended_source_path']]
   oldamended=source['amended_records' if row['amended_source_path']==row['source_path'] else 'records'][row['amended_source_row_index']]
   assert digest(oldamended)==row['amended_source_row_sha256']
   candidates.extend({'certificate_policy':'previous_smoothing',**s} for s in oldamended['strategies'])
  candidates.extend({'certificate_policy':'common_shift',**s} for s in row['strategies'])
  valid=[s for s in candidates if s['result']['kind'] in ('complete_quantum_resource_plan','deterministic_correction_prior')]
  if valid:assert row['selected_strategy']==min(valid,key=lambda s:(s['result']['T_count'],s['architecture'],s['certificate_policy']))
  else:
   assert row['rejection']['reason']=='Amended deterministic monetary tolerance exhausted' and det>=F(1,100);rejects+=1
 assert len(data['records'])==416
 for comparison in data['selected_comparisons']:
  for method,item in comparison['methods'].items():
   group=[r for r in data['records'] if r['model_id']==comparison['model_id'] and r['method']==method and r.get('selected_strategy')]
   win=min(group,key=lambda r:(r['selected_strategy']['result']['T_count'],r.get('normal_bits',0),r.get('PG_exp_bits',r.get('PG_bits',0)),r.get('PA_exp_bits',r.get('PA_bits',0)),r['selected_strategy']['architecture']))
   assert item['selected_record_sha256']==digest(win) and item['selected_strategy']==win['selected_strategy'] and item['T_count']==win['selected_strategy']['result']['T_count']
  rivals={m:comparison['methods'][m]['T_count'] for m in ('raw_payoff','numerical_conditional_price')};strong=min(rivals,key=lambda k:(rivals[k],k))
  T=rivals[strong];E=comparison['methods']['shared_threshold']['T_count'];assert comparison['strongest_rival']==strong and F(comparison['saving_against_strongest_exact'])==F(T-E,T)
 result={'schema':'qbp-common-shift-independent-plans-v1','status':'passed','source_sha256':sha(__file__),'resources_sha256':sha(out/'common_shift_resources.json.gz'),'protocol_sha256':sha(out/'common_shift_protocol.json'),'independent_library_sha256':sha(repo/'check_expanded_precision_resources.py'),'fixed_record_bindings':416,'amended_QPE_and_complete_charge_plans':plans,'deterministic_prior_plans':zeros,'rejected_configurations':rejects,'scope':'Own-formula QFT/QPE/full-charge replay and prior/budget/source/profile bindings. No producer or confidence/charging frontend imported. Mathematical certificate validity is separately reviewed.'}
 destination=Path(args.output)
 if destination.exists():raise FileExistsError('Preserve previous independent replay')
 destination.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--results',required=True);p.add_argument('--output',required=True);main(p.parse_args())
