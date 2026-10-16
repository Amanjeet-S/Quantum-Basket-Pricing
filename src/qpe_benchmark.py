"""Complete canonical-QPE price budgets on the immutable financial pool.

No quantum hardware is executed. The primary budgets use actual retained
Clifford+T inverse-QFT norms and exact binomial amplification. All single
flags and dyadic buckets receive identical proved financial support/moments.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import gzip,hashlib,json,time
import value_oracle as vo
import dyadic_plan as dp
from clipping_certificate import clipped_tail_bound
import qpe_estimator as qe
import qft_synthesis as qft
import controlled_estimation_resources as cr
ROOT=Path(__file__).resolve().parent
ALPHA=F(1,100)
ALLOCATIONS=('equal','cost_moment')
class UnsupportedScalarPlan(ValueError):pass

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(r):return tuple(r[k] for k in ('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits'))
def pack(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [pack(v) for v in x]
    return x

def scalar(e,p,alpha,flag):
    if p<=2*e:return {'status':'prior_zero_shot','epsilon':e,'prior':p,'alpha':alpha,'estimate':p/2,'T':0,'runs':0,'controlled_Q_calls':0,'metrics':{},'phase_register_bits':0,'peak_logical_qubits':0}
    if not dp.MIN_EPSILON<=e<F(1,2):raise UnsupportedScalarPlan('Probability tolerance outside declared domain')
    plan,profile=qe.compiled_plan(e,alpha,p)
    cost=cr.complete_scalar(plan,flag,profile)
    return {'status':'estimated','epsilon':e,'prior':p,'alpha':alpha,
        **cost,'plan':plan.record(),'qft_norm_error_upper':profile['norm_error_upper'],
        'A_T':vo.resource_api(flag)['A']['T'],'CQ_T':vo.resource_api(flag)['Q']['T']+14,
        'qft_T':profile['T_count'],'confidence_scope':'fresh independent complete QPE runs; certified approximate inverse QFT'}

def flag_candidate(kind,c,q,u,v,rho,maximum,e,body):
    level=F(1 << c,1 << q)
    if kind=='range' and maximum>1 << c:return {'kind':kind,'c':c,'status':'invalid_support'}
    tail=clipped_tail_bound(rho,u,v,level) if kind=='clipped' else {'upper':F(0),'case':'no_tail'}
    D=tail['upper']
    if D/2>=e:return {'kind':kind,'c':c,'status':'tail_exhausts_budget','tail_upper':D}
    mean=min(u,level,dp.sharp_grid_mean_upper(q,v,upper_integer=min(maximum,1 << c))) if maximum else F(0)
    prior=min(F(1),mean/level);eps=(e-D/2)/level
    flag=vo.profile_from_prepared_clipped_range(body,c) if kind=='clipped' else vo.profile_from_prepared_range(body,c)
    try:out=scalar(eps,prior,ALPHA,flag)
    except UnsupportedScalarPlan as err:return {'kind':kind,'c':c,'status':'unsupported_scalar_plan','reason':str(err),'tail_upper':D}
    out.update(kind=kind,c=c,level=level,tail_upper=D,tail_case=tail['case'],tail_midpoint=D/2,total_mean_error_upper=level*eps+D/2)
    assert out['total_mean_error_upper']<=e
    return out

def dyadic(allocation,q,u,v,maximum,e,body):
    s=dp.default_prefix(q,u,v,e,max_integer=maximum);buckets,tail=dp.partition(q,s,u,v,True,maximum)
    profiles={b.index:vo.profile_from_prepared_body(body,b.index) for b in buckets}
    resources={i:{'A':vo.resource_api(p)['A'],'Q':cr.emitted_profiles(p)['CQ']['metrics']} for i,p in profiles.items()}
    weights=dp.allocation_weights(buckets,resources,allocation,metric='T')
    errors=dp.capped_allocation(buckets,e-tail/2,weights,enforce_probability_floor=allocation!='equal')
    active=[i for i,(b,err) in enumerate(zip(buckets,errors)) if err<b.weight*b.prior/2];alpha=ALPHA/len(active) if active else F(0)
    records=[];totals={};T=runs=calls=peak=0
    for b,err in zip(buckets,errors):
        eps=err/b.weight
        if err>=b.weight*b.prior/2:
            out={'status':'prior_zero_shot','T':0,'runs':0,'controlled_Q_calls':0,'metrics':{},'phase_register_bits':0,'peak_logical_qubits':0,'estimate':b.prior/2}
        else:out=scalar(eps,b.prior,alpha,profiles[b.index])
        rec={**out,'index':b.index,'weight':b.weight,'prior':b.prior,'epsilon':eps,'alpha':alpha if out['runs'] else F(0),'mean_error_allowance':err};records.append(rec)
        T+=out['T'];runs+=out['runs'];calls+=out['controlled_Q_calls'];peak=max(peak,out['peak_logical_qubits'])
        for m,z in out['metrics'].items():totals[m]=totals.get(m,0)+z
    assert sum((r['alpha'] for r in records),F(0))<=ALPHA
    totalerror=tail/2+sum(errors,F(0));assert totalerror<=e
    return {'kind':'dyadic','allocation':allocation,'status':'complete','T':T,'runs':runs,'controlled_Q_calls':calls,'metrics':totals,'peak_logical_qubits':peak,'omitted_prefix_bits':s,'omitted_mean_upper':tail,'omitted_midpoint':tail/2,'buckets':records,'total_mean_error_upper':totalerror}

def main():
    started=time.perf_counter();rp=ROOT/'results/architecture_resource_candidates.json';fp=ROOT/'results/finite_range_certificates.json';raw=json.loads(rp.read_text());cert={key(r):r for r in json.loads(fp.read_text())['configurations']};rows=[]
    choose=lambda seq:min((r for r in seq if 'T' in r),key=lambda r:(r['T'],r.get('runs',0),r.get('c',99)))
    for index,row in enumerate(raw['configurations']):
        bound=cert[key(row)];q=row['q'];maximum=bound['max_integer'];u=F(bound['fair_mean_upper_fraction']);v=F(bound['fair_second_moment_upper_fraction']);rho=F(maximum,1 << q);e=F(1,200*row['cap']);body=vo.prepare_existing_body(row['resources'],q)
        plain=flag_candidate('range',q,q,u,v,rho,maximum,e,body)
        ranges=[flag_candidate('range',c,q,u,v,rho,maximum,e,body) for c in range(q+1)]
        clipped=[flag_candidate('clipped',c,q,u,v,rho,maximum,e,body) for c in range(q+1)]
        dyads=[]
        for allocation in ALLOCATIONS:
            try:dyads.append(dyadic(allocation,q,u,v,maximum,e,body))
            except UnsupportedScalarPlan as err:dyads.append({'kind':'dyadic','allocation':allocation,'status':'unsupported_scalar_plan','reason':str(err)})
        bestR,bestC,bestD=choose(ranges),choose(clipped),choose(dyads);bestS=choose([plain,bestR,bestC]);bestAll=choose([bestS,bestD])
        out={k:row[k] for k in ('case_id','model_id','method','candidate_id','architecture','cap','q','exp_address_bits','log_address_bits','probability_address_bits')}
        out.update(mean_upper=u,second_moment_upper=v,range=rho,max_integer=maximum,estimation_tolerance=e,alpha=ALPHA,deterministic_price_error_upper=row['independent_deterministic_bias_upper_decimal'],plain=plain,range_candidates=ranges,clipping_candidates=clipped,dyadic_candidates=dyads,best_range=bestR,best_clipped=bestC,best_dyadic=bestD,best_single_flag=bestS,best_all=bestAll)
        rows.append(pack(out))
        if (index+1)%100==0:print(json.dumps({'profiles_complete':index+1,'seconds':round(time.perf_counter()-started,3)}),flush=True)
    selected=[]
    for case in sorted({r['case_id'] for r in rows}):
        pool=[r for r in rows if r['case_id']==case]
        for family in ('plain','best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
            row=min(pool,key=lambda r:(r[family]['T'],r[family].get('runs',0),r['candidate_id'],r['architecture']));c=row[family]
            selected.append({'case_id':case,'family':family,'candidate_id':row['candidate_id'],'architecture':row['architecture'],'method':row['method'],'T':c['T'],'runs':c['runs'],'controlled_Q_calls':c['controlled_Q_calls'],'peak_logical_qubits':c['peak_logical_qubits'],'c':c.get('c'),'kind':c['kind'],'allocation':c.get('allocation')})
    ix={(x['case_id'],x['family']):x for x in selected};cases=sorted({x['case_id'] for x in selected});ratios={}
    for f in ('best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
        vals=[F(ix[(c,f)]['T'],ix[(c,'plain')]['T']) for c in cases]
        ratios[f]={'improvements':sum(x<1 for x in vals),'ties':sum(x==1 for x in vals),'losses':sum(x>1 for x in vals),'minimum_ratio':float(min(vals)),'maximum_ratio':float(max(vals))}
    result={'status':'fully charged canonical QPE comparison; no quantum hardware execution','scope':'all819 inherited financial/architecture profiles, same support/moments and monetary tolerances; per-price99% confidence; finite candidate minimum only','source_hashes':{n:sha(ROOT/n) for n in ['qpe_benchmark.py','controlled_estimation_resources.py','qpe_estimator.py','qft_synthesis.py','value_oracle.py','dyadic_plan.py','clipping_certificate.py']},'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [rp,fp,ROOT/'results/qft_rotation_tapes.json',ROOT/'results/qft_rotation_certificates.json']},'selection_protocol':{'alpha':'1/100','confidence':'exact binomial median using actual emitted inverse-QFT norm bound','range_and_clip_exponents':'all0..q, with support/tail/scalar-domain rejections retained','dyadic_allocations':ALLOCATIONS,'dyadic_prefix':'same quarter-budget omitted-prefix rule','primary_metric':'complete logical T including controlled reflections and certified inverse QFT','qubits':'conservative serial reusable-work bounds','controlled_sign':'Z on phase control included','wire_convention':'certified inverse QFT wire0MSB; controlled powers2**(m-1-j) on phasewirej'},'configurations':rows,'selected_cases':selected,'summary':{'profiles':len(rows),'selected_cases':len(selected),'ratios_to_QPE_plain':ratios},'construction_seconds_diagnostic':time.perf_counter()-started}
    payload=(json.dumps(pack(result),indent=2)+'\n').encode();archive=ROOT/'results/qpe_benchmark.json.gz';archive.write_bytes(gzip.compress(payload,compresslevel=9,mtime=0))
    summary={'status':result['status'],'summary':result['summary'],'selected_cases':selected,'archive_sha256':sha(archive),'uncompressed_sha256':hashlib.sha256(payload).hexdigest(),'uncompressed_bytes':len(payload),'compressed_bytes':archive.stat().st_size,'construction_seconds_diagnostic':result['construction_seconds_diagnostic']};(ROOT/'results/qpe_benchmark_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(result['summary'],indent=2),flush=True)
if __name__=='__main__':main()
