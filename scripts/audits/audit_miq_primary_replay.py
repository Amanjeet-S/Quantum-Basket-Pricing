"""Independent complete scalar scheduling replay for retained representatives.

Previously verified scalar interval primitives and Boolean gate emitters are
retained premises. No MIQAE
controller, production cost DP or variance-oracle comparator is imported.
The multiplier bounds, shot table and three objective DPs are reconstructed.
"""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
import argparse,gzip,hashlib,json
from certified_scalar import sqrt_upper_fraction,pi_rational_interval,inverse_probability_interval


def floor_f(x):return x.numerator//x.denominator
def ceil_f(x):return -((-x.numerator)//x.denominator)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def replay(epsilon,alpha,prior,a_cost,q_cost):
    epsilon,alpha,prior=map(F,(epsilon,alpha,prior))
    if prior<=2*epsilon:return {'shots':0,'Q_calls':0,'T':0,'table_entries':0}
    guarded=min(F(1),prior+F(1,10**18))
    variance=guarded*(1-guarded) if guarded<=F(1,2) else F(1,4)
    gamma=2*sqrt_upper_fraction(variance)
    M=pi_rational_interval()[1]*gamma/(4*(epsilon-F(1,10**60)))
    maximum=floor_f(M)
    if maximum%2==0:maximum-=1
    hi=inverse_probability_interval(prior)[1]
    first=floor_f(1/hi)
    if first%2==0:first-=1
    first=max(1,first)
    radius=F(69,1000);table=[]
    with localcontext() as ctx:
        ctx.prec=200
        for m in range(first,maximum+1,2):
            if m!=first and m<3*first:continue
            argument=3*M/(alpha*m)
            log_upper=F((Decimal(argument.numerator)/Decimal(argument.denominator)).ln()+Decimal('1e-60'))
            assert 0<log_upper<100
            n=ceil_f(log_upper/(2*radius*radius))
            table.append((m,n))
    # A separately coded sparse predecessor scan, maximising shot, Q and T
    # objectives independently. Every predecessor must be <=current/3.
    objectives=((1,0),(0,1),(a_cost,q_cost));dps=[{} for _ in objectives]
    eligible=[None]*len(objectives);best=[0]*len(objectives);cursor=0
    for m,n in table:
        while cursor<len(table) and table[cursor][0]<=m//3:
            old=table[cursor][0]
            for j,dp in enumerate(dps):
                if old in dp and (eligible[j] is None or dp[old]>eligible[j]):eligible[j]=dp[old]
            cursor+=1
        for j,(ca,cq) in enumerate(objectives):
            local=n*(ca+cq*((m-1)//2))
            if m==first:value=local
            elif eligible[j] is not None:value=eligible[j]+local
            else:continue
            dps[j][m]=value;best[j]=max(best[j],value)
    return {'shots':best[0],'Q_calls':best[1],'T':best[2],'table_entries':len(table),
            'initial_multiplier':first,'maximum_odd_multiplier':maximum,
            'shot_table_sha256':hashlib.sha256(repr(tuple(table)).encode()).hexdigest()}


def representatives(data,source_root):
    rows=data.get('configurations',data.get('rows',[]));result=[]
    # Include the largest retained Grover budgets in each financial method
    # and scalar family. Selection itself uses no outcomes.
    for method in sorted({r['method'] for r in rows}):
        for family in ('plain','best_range','best_clipped'):
            pool=[(r,r['controllers']['miqae'][family]) for r in rows if r['method']==method]
            r,s=max(pool,key=lambda z:(z[1]['Q_calls'],z[1]['T'],z[0]['candidate_id']))
            if s['status']!='estimated':continue
            result.append({'case_id':r['case_id'],'method':method,'architecture':r['architecture'],
                'candidate_id':r['candidate_id'],'family':family,
                'input':s['scalar_plan_input'],'A_T':s['A_T'],'Q_T':s['Q_T'],
                'expected':{k:s[k] for k in ('shots','Q_calls','T')}})
    from audit_variance_oracle_comparison import independent_profile
    original=json.loads((source_root/'results/architecture_resource_candidates.json').read_text())['configurations']
    original={(r['candidate_id'],r['architecture']):r['resources'] for r in original}
    for method in sorted({r['method'] for r in rows}):
        pool=[(r,d,b) for r in rows if r['method']==method for d in r['controllers']['miqae']['dyadic_candidates'] if 'T' in d for b in d['buckets'] if b['status']=='estimated']
        if not pool:continue
        r,d,b=max(pool,key=lambda z:(z[2]['maximum_Q_calls'],z[2]['metric_caps'].get('T',0),z[0]['candidate_id']))
        A,Q,qubits=independent_profile(original[(r['candidate_id'],r['architecture'])],r['q'],'dyadic',b['index'])
        result.append({'case_id':r['case_id'],'method':method,'architecture':r['architecture'],
            'candidate_id':r['candidate_id'],'family':'dyadic_largest_bucket',
            'bucket':b['index'],'allocation':d['allocation'],
            'input':{'epsilon':b['probability_tolerance'],'alpha':b['alpha'],'prior':b['prior_upper']},
            'A_T':A['T'],'Q_T':Q['T'],
            'expected':{'shots':b['maximum_shots'],'Q_calls':b['maximum_Q_calls'],'T':b['metric_caps']['T']}})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args();container=args.input.read_bytes();raw=gzip.decompress(container) if args.input.suffix=='.gz' else container;data=json.loads(raw);rows=[]
    for r in representatives(data,args.input.parents[1]):
        inputs=r['input'];actual=replay(inputs['epsilon'],inputs['alpha'],inputs['prior'],r['A_T'],r['Q_T'])
        assert all(actual[k]==v for k,v in r['expected'].items()),(r,actual)
        rows.append({**r,'replay':actual})
    out={'status':'passed','input_sha256':hashlib.sha256(container).hexdigest(),'uncompressed_sha256':hashlib.sha256(raw).hexdigest(),
         'audit_source_sha256':digest(Path(__file__)),'representatives':rows,
         'scope':'independent complete finite shot-table and coupled-T recurrence; no production planner/controller imports; verified scalar interval primitives retained as premises'}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'representatives':len(rows),'total_table_entries':sum(r['replay']['table_entries'] for r in rows)},indent=2))
