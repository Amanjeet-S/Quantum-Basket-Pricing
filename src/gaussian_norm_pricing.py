"""Matched finite basket certificates using joint Gaussian derivative norms.

Full degree paths reuse retained exact coefficient intervals; their geometry
moments are recomputed and charged. The new stopping hierarchy is executed
with its own fresh geometry bank and actually evaluated coefficients.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json,time
from flint import ctx
from continuous_quadrature import frac_arb,exact_endpoint,serialise_interval
import taylor_order_benchmark as policy
from taylor_pricing import model_parameters
from gaussian_derivative_norm import upper_norm
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
TOLERANCE=F(1,100)

def joint_radius(bank,lam,cap,k):
    bank.ensure(k+1)
    return max(F(0),exact_endpoint(cap*frac_arb(upper_norm(k+1,lam))*bank.centred[k+1]/factorial(k+1),True))

def main():
    ctx.prec=224;start=time.perf_counter()
    protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text())
    old=json.loads((ROOT/'results/taylor_pricing.json').read_text())
    previous={r['case_id']:r for r in old['configurations']}
    old_hierarchy=json.loads((ROOT/'results/taylor_order_benchmark.json').read_text())['family_runs'][2]
    old_choices={r['case_id']:r for r in old_hierarchy['configurations']}
    models={m['model_id']:m for m in protocol['models']};banks={};paths=[]
    path_start=time.perf_counter()
    for row in old['configurations']:
        model=models[row['model_id']];w,d,s,mg,rg,lam=model_parameters(model)
        key=(w,d,s,rg);bank=banks.setdefault(key,policy.Geometry(key));orders=[];clo=chi=F(0)
        for o in row['orders']:
            k=o['order'];c=o['coefficient_interval']
            clo+=F(c['lower_exact_fraction']);chi+=F(c['upper_exact_fraction'])
            rem=joint_radius(bank,lam,row['cap'],k)
            lo=max(F(0),clo-rem);hi=min(F(row['cap']),chi+rem);assert lo<=hi
            old_it=o['undiscounted_price_interval'];olo,ohi=F(old_it['lower_exact_fraction']),F(old_it['upper_exact_fraction'])
            assert max(lo,olo)<=min(hi,ohi)
            orders.append({'order':k,'remainder_upper':str(rem),'old_remainder_upper':o['remainder_upper_exact_fraction'],
               'interval':serialise_interval(lo,hi),'certified':(hi-lo)/2<=TOLERANCE,
               'old_interval_radius':str((ohi-olo)/2),'radius_ratio_to_old':str((hi-lo)/(ohi-olo)),
               'coefficient_interval_sum':{'lower':str(clo),'upper':str(chi)}})
        good=[o for o in orders if o['certified']];best=min(good or orders,key=lambda o:F(o['interval']['upper_exact_fraction'])-F(o['interval']['lower_exact_fraction']))
        paths.append({'case_id':row['case_id'],'assets':row['assets'],'cap':row['cap'],'lambda':str(lam),
           'orders':orders,'first_certifying_order':good[0]['order'] if good else None,
           'old_first_certifying_order':row['first_certifying_order'],'certified':bool(good),'best':best})
    retightening={'status':'recomputed moment geometry and reused retained certified coefficients',
       'kernel_evaluations_executed':0,'coefficient_interval_source_sha256':sha(ROOT/'results/taylor_pricing.json'),
       'geometry_exponentials_executed':sum(b.exponentials for b in banks.values()),
       'seconds_diagnostic':time.perf_counter()-path_start,'configurations':paths}
    hierarchy=[];banks={};original=policy.radius;run_start=time.perf_counter()
    try:
        policy.radius=joint_radius
        for model in protocol['models']:
            w,d,s,mg,rg,lam=model_parameters(model);key=(w,d,s,rg);bank=banks.setdefault(key,policy.Geometry(key))
            for cap in protocol['caps']:
                row=policy.evaluate(model,cap,'monotone_then_screen',bank)
                old_row=old_choices[row['case_id']]
                row['remainder_policy']='certified joint Hermite-Gaussian norm'
                row['old_hierarchy_order']=old_row['selected']['order'];row['old_hierarchy_kernel_calls']=old_row['kernel_calls']
                hierarchy.append(row)
    finally:policy.radius=original
    experiment={'execution_scope':'new hierarchy, own fresh banks, same frozen 27 contracts',
       'configurations':hierarchy,'kernel_calls':sum(r['kernel_calls'] for r in hierarchy),
       'geometry_exponentials':sum(b.exponentials for b in banks.values()),
       'certified_cases':sum(r['selected']['certified'] for r in hierarchy),
       'seconds_diagnostic':time.perf_counter()-run_start}
    ratios=[F(r['best']['radius_ratio_to_old']) for r in paths]
    result={'schema':'qbp-joint-gaussian-norm-pricing-v1','status':'completed',
      'precision_bits':224,'degree_limit':8,'radius_target':str(TOLERANCE),
      'retightened_paths':retightening,'executed_hierarchy':experiment,
      'summary':{'certified_paths':sum(r['certified'] for r in paths),
         'first_degree_reductions':sum(r['first_certifying_order'] is not None and r['old_first_certifying_order'] is not None and r['first_certifying_order']<r['old_first_certifying_order'] for r in paths),
         'minimum_best_radius_ratio':float(min(ratios)),'maximum_best_radius_ratio':float(max(ratios)),
         'old_hierarchy_kernel_calls':old_hierarchy['kernel_calls'],'new_hierarchy_kernel_calls':experiment['kernel_calls'],
         'old_hierarchy_geometry_terms':old_hierarchy['geometry_exponentials'],'new_hierarchy_geometry_terms':experiment['geometry_exponentials'],
         'executed_hierarchy_kernel_ratio':str(F(experiment['kernel_calls'],old_hierarchy['kernel_calls']))},
      'source_hashes':{n:sha(ROOT/n) for n in ['gaussian_norm_pricing.py','gaussian_derivative_norm.py','taylor_pricing.py','taylor_order_benchmark.py','continuous_quadrature.py']},
      'input_hashes':{n:sha(ROOT/n) for n in ['results/gaussian_derivative_norms.json','results/taylor_pricing.json','results/taylor_order_benchmark.json','results/generalisation_protocol.json']},
      'interpretation':'A sharper certified remainder for the same finite coefficients. Norm generation is additional shared setup and is not omitted from total-work interpretation. No new quantum baseline or physical advantage.'}
    (ROOT/'results/gaussian_norm_pricing.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'summary':result['summary'],'hierarchy':{k:v for k,v in experiment.items() if k!='configurations'}},indent=2))
if __name__=='__main__':main()
