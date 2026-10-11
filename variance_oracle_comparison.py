"""Costed finite-output mean estimators with fair range and clipping opponents.

All methods receive the same actual-output support and proved moment bounds.
Primary T caps use the same coupled MIQAE scheduling DP, or the same product
cap for the three/five controller. Other resources use conservative product
caps. No quantum hardware or physical wall-clock advantage is claimed.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import hashlib,json,time,argparse
from functools import lru_cache
import value_oracle as vo
import dyadic_plan as dp
import clipping_certificate as clip
import full_batch_miqae as miq
import iterative_estimator as three
from certified_scalar import ceil_fraction
ROOT=Path(__file__).resolve().parent
ALPHA=F(1,100)
CONTROLLERS=('three_five','miqae')
ALLOCATIONS=('equal','cost_moment')

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def key(r):return tuple(r[k] for k in ('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits'))
def pack(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):pack(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [pack(v) for v in x]
    return x

def thin(p):return {'A':{'T':p['A']['T']},'Q':{'T':p['Q']['T']},'qubits':p['qubits']}
def charge(profile,n,k,T):
    names=set(profile['A'])|set(profile['Q'])
    result={m:n*profile['A'].get(m,0)+k*profile['Q'].get(m,0) for m in names}
    result['T']=T;result['final_flag_measurements']=n;result['complete_reinitialisations']=n
    return result

def scalar(controller,epsilon,prior,p):
    if prior<=2*epsilon:
        return {'status':'prior_zero_shot','epsilon':epsilon,'prior':prior,
                'estimate_without_shots':prior/2,'shots':0,'Q_calls':0,'T':0,
                'qubits':0,'charges':charge(p,0,0,0)}
    if not dp.MIN_EPSILON<=epsilon<F(1,2):raise ValueError('scalar tolerance outside supported range')
    s=(miq if controller=='miqae' else three).plan(epsilon,ALPHA,prior)
    n=s.total_shots_upper;k=s.total_grover_calls_upper
    T=dp.scalar_cost(s,controller,p['A']['T'],p['Q']['T'])
    return {'status':'estimated','epsilon':epsilon,'prior':prior,'shots':n,'Q_calls':k,'T':T,
            'qubits':p['qubits'],'charges':charge(p,n,k,T),
            'scalar_plan_input':{'epsilon':epsilon,'alpha':ALPHA,'prior':prior},
            'A_T':p['A']['T'],'Q_T':p['Q']['T']}

def candidate(controller,kind,c,q,u,v,rho,maximum,e,body):
    r=F(1 << c,1 << q)
    if kind=='range' and maximum>1 << c:return {'kind':kind,'c':c,'status':'invalid_support'}
    tail=clip.clipped_tail_bound(rho,u,v,r) if kind=='clipped' else {'upper':F(0),'case':'no_tail'}
    D=tail['upper']
    if D/2>=e:return {'kind':kind,'c':c,'status':'tail_exhausts_budget','tail_upper':D,'tail_case':tail['case']}
    bound=min(u,r,dp.sharp_grid_mean_upper(q,v,upper_integer=min(maximum,1 << c))) if maximum else F(0)
    prior=min(F(1),bound/r)
    profile=vo.resource_api(vo.profile_from_prepared_clipped_range(body,c) if kind=='clipped' else vo.profile_from_prepared_range(body,c))
    try:s=scalar(controller,(e-D/2)/r,prior,profile)
    except ValueError as error:return {'kind':kind,'c':c,'status':'unsupported_scalar_plan','reason':str(error),'tail_upper':D}
    s.update(kind=kind,c=c,range=r,tail_upper=D,tail_case=tail['case'],tail_midpoint=D/2,
             total_mean_error_upper=r*s['epsilon']+D/2,
             moment_bounds_shared_with_all_opponents=True)
    assert s['total_mean_error_upper']<=e
    return s

def dyadic(controller,allocation,q,u,v,e,maximum,body):
    profiles={b:vo.resource_api(vo.profile_from_prepared_body(body,b)) for b in range(q+1)}
    try:s=dp.build_plan(q,u,v,e,ALPHA,{b:thin(p) for b,p in profiles.items()},
        controller=controller,allocation=allocation,metric='T',sharp_priors=True,max_integer=maximum)
    except ValueError as error:return {'kind':'dyadic','allocation':allocation,'status':'unsupported_scalar_plan','reason':str(error)}
    T=s['aggregate_metrics_upper'].get('T',0)
    charges={};blackboxT=0;records=[]
    for b in s['buckets']:
        n,k=b['maximum_shots'],b['maximum_Q_calls'];i=b['index'];bt=b['metric_caps'].get('T',0)
        cost=charge(profiles[i],n,k,bt)
        for m,z in cost.items():charges[m]=charges.get(m,0)+z
        bp=vo.resource_api(vo.profile_from_prepared_body(body,i),blackbox=True)
        blackboxT+=n*bp['A']['T']+k*bp['Q']['T']
        records.append({z:b[z] for z in b if z!='scalar_plan'})
    assert charges.get('T',0)==T
    return {'kind':'dyadic','allocation':allocation,'status':s['status'],'T':T,
        'shots':s['total_shots_upper'],'Q_calls':s['total_Q_calls_upper'],
        'qubits':s['peak_qubits_upper'],'charges':charges,'blackbox_value_T_cap':blackboxT,
        'total_mean_error_upper':s['total_error_upper'],'mean_upper':s['mean_upper'],
        'second_moment_upper':v,'max_integer':maximum,
        'omitted_prefix_bits':s.get('omitted_prefix_bits'),
        'omitted_mean_interval':s.get('omitted_mean_interval'),
        'estimate_without_shots':s.get('estimate_without_shots'),
        'buckets':records}

def classical(r,e):
    rho=F(r['finite_lattice_range_fraction']);u=F(r['fair_mean_upper_fraction']);v=F(r['fair_second_moment_upper_fraction'])
    m=min(u,rho/2,v/rho) if rho else F(0);variance=m*(rho-m)
    # Bernstein with |X-EX|<=rho. Log allowance is strictly outward.
    with localcontext() as c:
        c.prec=100;ln=F(Decimal(200).ln()+Decimal('1e-80'))
    B=ceil_fraction((2*variance+2*rho*e/3)*ln/(e*e))
    H=ceil_fraction(rho*rho*ln/(2*e*e))
    return {'range':rho,'mean_upper':u,'second_moment_upper':v,'variance_upper':variance,
        'Bernstein_sample_cap':B,'Hoeffding_sample_cap':H,'sample_cap':min(B,H),
        'status':'prospective bound only; no new stream executed; ideal independent samples'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--start',type=int,default=0);parser.add_argument('--stop',type=int,default=819);parser.add_argument('--part',action='store_true');parser.add_argument('--merge',action='store_true');args=parser.parse_args()
    start=time.perf_counter();rp=ROOT/'results/architecture_resource_candidates.json';fp=ROOT/'results/finite_range_certificates.json'
    raw=json.loads(rp.read_text());ranges=json.loads(fp.read_text());cert={key(r):r for r in ranges['configurations']}
    rows=[]
    for index,row in enumerate(raw['configurations']):
        if args.merge or not args.start<=index<args.stop:continue
        r=cert[key(row)];q=row['q'];maximum=r['max_integer'];rho=F(maximum,1 << q)
        u=F(r['fair_mean_upper_fraction']);v=F(r['fair_second_moment_upper_fraction']);e=F(1,200*row['cap'])
        body=vo.prepare_existing_body(row['resources'],q)
        result={k:row[k] for k in ('case_id','model_id','method','candidate_id','architecture','cap','q','exp_address_bits','log_address_bits','probability_address_bits')}
        result.update(mean_upper=u,second_moment_upper=v,max_integer=maximum,range=rho,
                      normalised_estimation_budget=e,alpha=ALPHA,
                      deterministic_price_error_upper=row['independent_deterministic_bias_upper_decimal'],
                      classical_range_budget=classical(r,e),controllers={})
        for controller in CONTROLLERS:
            plain=candidate(controller,'range',q,q,u,v,rho,maximum,e,body)
            rng=[candidate(controller,'range',c,q,u,v,rho,maximum,e,body) for c in range(q+1)]
            clipped=[candidate(controller,'clipped',c,q,u,v,rho,maximum,e,body) for c in range(q+1)]
            dyads=[dyadic(controller,a,q,u,v,e,maximum,body) for a in ALLOCATIONS]
            valid=lambda seq:[x for x in seq if 'T' in x]
            choose=lambda seq:min(valid(seq),key=lambda x:(x['T'],x['shots'],x.get('c',q)))
            best_range=choose(rng);best_clip=choose(clipped);best_dyad=choose(dyads)
            best_single=choose([plain,best_range,best_clip]);best_all=choose([best_single,best_dyad])
            result['controllers'][controller]={'plain':plain,'range_candidates':rng,'clipping_candidates':clipped,
                'dyadic_candidates':dyads,'best_range':best_range,'best_clipped':best_clip,
                'best_dyadic':best_dyad,'best_single_flag':best_single,'best_all':best_all}
        rows.append(pack(result))
        if (index+1)%3==0:
            (ROOT/f'results/variance_oracle_checkpoint_{args.start}_{args.stop}.json').write_text(json.dumps({'completed':len(rows),'source_sha256':digest(Path(__file__)),'rows':rows},separators=(',',':'))+'\n')
        if (index+1)%3==0:miq.plan.cache_clear()
        if (index+1)%25==0:print(json.dumps({'profiles_complete':index+1,'elapsed_seconds':round(time.perf_counter()-start,2)}),flush=True)
    if args.part:
        part={'source_sha256':digest(Path(__file__)),'input_hashes':{str(p.relative_to(ROOT)):digest(p) for p in (rp,fp)},'start':args.start,'stop':args.stop,'rows':rows,'elapsed_seconds':time.perf_counter()-start}
        (ROOT/f'results/variance_oracle_part_{args.start}_{args.stop}.json').write_text(json.dumps(part,indent=2)+'\n');print(json.dumps({'part_complete':len(rows),'start':args.start,'stop':args.stop}),flush=True);return
    execution_segments=[]
    if args.merge:
        for lo,hi in ((0,273),(273,546),(546,819)):
            part=json.loads((ROOT/f'results/variance_oracle_part_{lo}_{hi}.json').read_text())
            assert part['source_sha256']==digest(Path(__file__))
            assert part['input_hashes']=={str(p.relative_to(ROOT)):digest(p) for p in (rp,fp)}
            assert len(part['rows'])==hi-lo;rows.extend(part['rows']);execution_segments.append({k:part[k] for k in ('start','stop','elapsed_seconds','source_sha256')})
        assert len(rows)==819
    selected=[]
    for case in sorted({r['case_id'] for r in rows}):
        group=[r for r in rows if r['case_id']==case]
        for name in ('plain','best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
            pool=[(r,c,r['controllers'][c][name]) for r in group for c in CONTROLLERS]
            r,c,s=min(pool,key=lambda z:(z[2]['T'],z[2]['shots'],z[0]['candidate_id'],z[1]))
            selected.append({'case_id':case,'family':name,'candidate_id':r['candidate_id'],
                'method':r['method'],'architecture':r['architecture'],'controller':c,
                'T':s['T'],'shots':s['shots'],'Q_calls':s['Q_calls'],'qubits':s['qubits'],
                'c':s.get('c'),'allocation':s.get('allocation'),'deterministic_price_error_upper':r['deterministic_price_error_upper']})
    selected_dict={(r['case_id'],r['family']):r for r in selected}
    ratios={}
    for family in ('best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
        vals=[F(selected_dict[(case,family)]['T'],selected_dict[(case,'plain')]['T']) for case in sorted({r['case_id'] for r in rows})]
        ratios[family]={'strict_improvements':sum(x<1 for x in vals),'ties':sum(x==1 for x in vals),
            'worse':sum(x>1 for x in vals),'minimum_ratio':float(min(vals)),'maximum_ratio':float(max(vals))}
    dyad_vs_single=[F(selected_dict[(case,'best_dyadic')]['T'],selected_dict[(case,'best_single_flag')]['T']) for case in sorted({r['case_id'] for r in rows})]
    out={'status':'complete certified resource-cap comparison; no quantum execution',
         'confidence_scope':'each selected price at 99 per cent under ideal fresh shots; not simultaneous 27-price confidence',
         'scope':'all 819 inherited emitted profiles; finite candidate minimum, not global optimality; same strengthened moments/support/controller/T scheduling caps for all families',
         'selection_protocol':{'controllers':CONTROLLERS,'dyadic_allocations':ALLOCATIONS,'dyadic_prefix':'default quarter-budget midpoint prefix only',
             'range_exponents':'all c=0,...,q subject to actual-output support','clip_exponents':'all c=0,...,q subject to tail and scalar constraints',
             'primary_resource':'logical T; MIQAE coupled admissible-chain DP; three/five product cap',
             'other_resources':'conservative maximum shots*A + maximum Q calls*Q; black-box value alternative uses this product cap; flag measurements and resets additional'},
         'source_hashes':{n:digest(ROOT/n) for n in ('variance_oracle_comparison.py','value_oracle.py','dyadic_plan.py','clipping_certificate.py','full_batch_miqae.py','iterative_estimator.py')},
         'input_hashes':{str(p.relative_to(ROOT)):digest(p) for p in (rp,fp)},
         'configurations':rows,'selected_cases':selected,
         'summary':{'profiles':len(rows),'selected_records':len(selected),'family_ratios_to_plain':ratios,
             'dyadic_vs_best_single_flag':{'wins':sum(x<1 for x in dyad_vs_single),'ties':sum(x==1 for x in dyad_vs_single),
                 'losses':sum(x>1 for x in dyad_vs_single),'minimum_ratio':float(min(dyad_vs_single)),'maximum_ratio':float(max(dyad_vs_single))}},
         'execution_segments':execution_segments,'elapsed_seconds':time.perf_counter()-start,
         'novelty_boundary':'dyadic reductions and moment-option clipping bounds are established; project-specific finite circuit/range/cost coupling, priority unresolved'}
    (ROOT/'results/variance_oracle_comparison.json').write_text(json.dumps(pack(out),indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2),flush=True)
if __name__=='__main__':main()
