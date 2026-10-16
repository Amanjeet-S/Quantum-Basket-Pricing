"""Read-only independent consistency review; does not rerun historical grids."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,hashlib,json
parser=argparse.ArgumentParser()
parser.add_argument('--root',type=Path,required=True)
parser.add_argument('--classical-protocol',type=Path)
parser.add_argument('--output',type=Path,default=Path(__file__).with_name('deterministic_evidence_review.json'))
args=parser.parse_args();S=args.root.resolve()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
a=load(S/'results/deterministic_pricing.json');b=load(S/'results/continuous_quadrature.json')
aa=load(S/'results/deterministic_pricing_audit.json');ba=load(S/'results/continuous_quadrature_audit.json')
pp=args.classical_protocol or S/'results/classical_pricing_protocol.json';p=load(pp);gp=S/'results/generalisation_protocol.json';g=load(gp)
assert a['protocol_sha256']==sha(pp) and b['protocol_sha256']==sha(gp)
assert len(p['all_candidate_profiles'])==len(a['tensor_inventory'])==273
assert a['enumerated_configurations']==len(a['configurations'])==47
assert a['total_enumerated_outputs']==sum(x['enumerated_states'] for x in a['configurations'])==330752
for n,h in a['source_hashes'].items():assert sha(S/n)==h,(n,h,sha(S/n))
assert b['source_sha256']==sha(S/'continuous_quadrature.py')
assert aa['result_sha256']==sha(S/'results/deterministic_pricing.json')
assert aa['continuous_result_sha256']==sha(S/'results/continuous_quadrature.json')
assert aa['audit_source_sha256']==sha(S/'checks/audit_deterministic_pricing.py')
assert ba['result_sha256']==sha(S/'results/continuous_quadrature.json')
assert ba['audit_source_sha256']==sha(S/'checks/audit_continuous_quadrature.py')
assert aa['complete_grid_configurations']==len(aa['checks'])==47
assert aa['independently_replayed_outputs']==sum(x['states'] for x in aa['checks'])==330752
assert len(b['configurations'])==len(b['independent_representation_checks'])==9
assert ba['scalar_cap_kernel_checks']==len(ba['checks'])==126
assert all(c['enclosures_overlap'] for c in ba['checks'])
assert len({(c['case_id'],c['y'],c['basket']) for c in ba['checks']})==126
models={m['model_id']:m for m in p['models']};rows={i['row']['candidate_id']:i for i in p['all_candidate_profiles']}
inv={i['candidate_id']:i for i in a['tensor_inventory']}
cont={r['case_id']:r for r in b['configurations']}
eligible=[]
for id,item in rows.items():
 row=item['row'];m=models[row['model_id']];ex=m['outer_dimension']*row['b_y']
 if row['method']=='direct':ex+=row['b_z']
 elif row['method']=='interval':ex+=row['b_u']
 assert int(inv[id]['tensor_states'])==1<<ex and inv[id]['grid_address_bits']==ex
 assert F(inv[id]['tensor_states_to_certified_MC_samples_exact_fraction'])==F(1<<ex,item['certified_samples'])
 expected=m['assets']==2 and row['method'] in ('full_global','full_cap') and (1<<ex)<=131072
 assert inv[id]['enumerated']==expected
 if expected:eligible.append(id)
assert set(eligible)=={r['candidate_id'] for r in a['configurations']}
ratios=[];radii=[];exact_checks=[]
for r in a['configurations']:
 row=r['row'];n=r['enumerated_states'];q=row['q'];s=int(r['integer_sum'])
 assert row==rows[r['candidate_id']]['row']
 assert 0<=s<=n*(1<<q)<(1<<64)
 mu=F(s,n*(1<<q));assert mu==F(r['normalised_finite_mean_exact_fraction'])
 point=F(r['geometric_correction']['undiscounted_value_decimal'])+row['cap']*mu
 beta=F(row['independent_deterministic_bias_upper_decimal'])
 assert beta<F(1,200)
 lo=F(r['undiscounted_price_interval_lower_exact_fraction']);hi=F(r['undiscounted_price_interval_upper_exact_fraction'])
 assert lo==point-beta and hi==point+beta
 assert F(r['undiscounted_price_interval_lower_decimal'])<=lo
 assert F(r['undiscounted_price_interval_upper_decimal'])>=hi
 assert r['confidence']==1 and r['statistical_error']==0
 ci=cont[r['case_id']]['undiscounted_price_interval']
 assert lo<=F(ci['lower_exact_fraction'])<=F(ci['upper_exact_fraction'])<=hi
 assert n<r['certified_MC_samples_same_prefix']
 work=r['online_algorithm_operation_components_upper'];assert sum(work.values())==r['online_algorithm_operation_upper_total']
 P=1+r['evaluator_profile']['nodes_by_kind']['function']+1;E=r['evaluator_profile']['scalar_lane_operation_upper_total']
 assert r['online_algorithm_operation_upper_total']==n*(E+3+2*P)-P-1
 ratios.append(F(n,r['certified_MC_samples_same_prefix']));radii.append(beta)
 exact_checks.append({'candidate_id':r['candidate_id'],'states':n,'exact_endpoint_and_sum_identity_checked':True,'continuous_interval_contained':True})
for winner in a['selected']:
 options=[r for r in a['configurations'] if r['case_id']==winner['case_id']]
 selected=min(options,key=lambda r:(r['online_algorithm_operation_upper_total'],r['setup']['table_entries_prepared'],r['row']['q'],r['method'],r['candidate_id']))
 assert selected['candidate_id']==winner['candidate_id']
 assert selected['evaluator_profile']['scalar_lane_operation_upper_total']==73
 assert selected['online_algorithm_operation_upper_total']==90*selected['enumerated_states']-8
widths=[];quad_checks=[]
for row in b['configurations']:
 assert row['precision_bits']==128 and row['normal_truncation']==8
 assert row['requested_absolute_integration_tolerance']=='1e-16'
 assert row['confidence']==1 and row['statistical_error']==0
 i=row['undiscounted_price_interval'];lo=F(i['lower_exact_fraction']);hi=F(i['upper_exact_fraction']);width=hi-lo
 assert width>0 and width<F(1,1000)
 assert F(i['lower_decimal'])<=lo and F(i['upper_decimal'])>=hi
 assert F(i['width_upper_decimal'])>=width
 assert F(i['midpoint_radius_upper_decimal'])>=(hi-lo)/2
 d=row['integration_diagnostics'];assert d['nonfinite_analytic_returns']>0
 assert d['nonfinite_analytic_returns']<=d['analytic_verification_requests']<=d['calls']
 ref=next(c['reference'] for c in b['independent_representation_checks'] if c['case_id']==row['case_id'])
 assert ref['precision_bits']==176 and ref['normal_truncation']==10
 assert ref['requested_absolute_integration_tolerance']=='1e-24'
 ri=ref['undiscounted_price_interval'];assert max(lo,F(ri['lower_exact_fraction']))<=min(hi,F(ri['upper_exact_fraction']))
 widths.append(width);quad_checks.append({'case_id':row['case_id'],'primary_bits':128,'primary_L':8,'reference_bits':176,'reference_L':10,'exact_interval_width_fraction':str(width),'analytic_branch_guard_rejections':d['nonfinite_analytic_returns']})
assert sum(r['integration_diagnostics']['calls'] for r in b['configurations'])==2112
for m in b['models']:
 w=list(map(F,m['weights']));rr=[F(r[0]) for r in m['loadings']];mu=list(map(F,m['log_drifts']));lam=F(m['lambda'])
 assert m['assets']==2 and m['outer_dimension']==1 and sum(w)==1 and min(w)>0 and lam>0
 assert all(mu[j]+(rr[j]**2+lam**2)/2==F(m['risk_free_rate'])-F(m['dividend_yield']) for j in range(2))
 assert F(m['maturity'])==1
 assert sum(w[j]*rr[j] for j in range(2))**2==F(m['geometric_outer_variance'])
with localcontext() as c:
 c.prec=100
 def dec(f):return str(Decimal(f.numerator)/Decimal(f.denominator))
 review={'status':'passed independent retained-evidence consistency review','date':'2026-10-06','scope':'Recomputes rational algebra, protocol inventory, selection and file provenance. Does not rerun historical exhaustive output evaluation or entire continuous price integration. Those executions remain separate hash-checked retained audits.',
 'checked_financial_prefixes':273,'checked_exact_price_records':47,'retained_full_interpreter_replay_outputs':330752,'checked_continuous_price_records':9,'retained_positive_kernel_checks':126,'total_primary_integrand_calls':2112,
 'finite_radius_max_exact_fraction':str(max(radii)),'finite_radius_max_decimal':dec(max(radii)),
 'finite_states_to_certified_MC_samples_min_exact_fraction':str(min(ratios)),'finite_states_to_certified_MC_samples_max_exact_fraction':str(max(ratios)),
 'primary_continuous_width_min_exact_fraction':str(min(widths)),'primary_continuous_width_max_exact_fraction':str(max(widths)),
 'primary_continuous_width_min_decimal':dec(min(widths)),'primary_continuous_width_max_decimal':dec(max(widths)),
 'exact_checks':exact_checks,'continuous_checks':quad_checks,
 'source_hashes':{str(p.relative_to(S)):sha(p) for p in [S/'deterministic_pricing.py',S/'continuous_quadrature.py',S/'proof_sources/deterministic_baseline_proof.tex',S/'checks/audit_deterministic_pricing.py',S/'checks/audit_continuous_quadrature.py']},
 'result_hashes':{str(p.relative_to(S)):sha(p) for p in [S/'results/deterministic_pricing.json',S/'results/deterministic_pricing_audit.json',S/'results/continuous_quadrature.json',S/'results/continuous_quadrature_audit.json']},
 'protocol_hashes':{'classical_pricing_protocol.json':sha(pp),'generalisation_protocol.json':sha(gp)},
 'review_source_sha256':sha(Path(__file__)),
 'limitations':'Both continuous representations and the scalar integral audit share Arb. Numerical-library correctness is a premise. No empirical sample, calibrated market validation, classical complexity lower bound, theoretical novelty or physical quantum runtime advantage follows.'}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(review,indent=2)+'\n')
print(json.dumps({k:v for k,v in review.items() if k not in ('source_hashes','result_hashes','protocol_hashes','exact_checks','continuous_checks')},indent=2))
