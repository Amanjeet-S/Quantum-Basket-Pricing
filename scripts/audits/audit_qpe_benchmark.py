"""Independent result audit of complete QPE budgets and finite Taylor policy.

This script does not import either benchmark runner or its estimator planner.
It checks finite formulas, actual retained primitive rotation lists, emitted
predicate gates, frozen certificates, selection minima and classical counters.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import math
import sys
import time

PI=F(355,113);XI=F(1,10**60);ALPHA=F(1,100)

def ceil(x):return -((-x.numerator)//x.denominator)
def upper_sqrt(x,bits=100):
    x=F(x);scale=1<<bits;n=x*scale*scale;k=math.isqrt(n.numerator//n.denominator)
    return F(k+(k*k*n.denominator<n.numerator),scale)
def clip_root(x):
    x=F(x);a,b=math.isqrt(x.numerator),math.isqrt(x.denominator)
    if a*a==x.numerator and b*b==x.denominator:return F(a,b)
    return upper_sqrt(x)
def sharp_mean(q,v,low=1,high=None):
    high=(1<<q)if high is None else high;d=F(1,1<<q);lo,hi=low*d,high*d
    if v<lo*lo:return v/lo
    if v>=hi*hi:return hi
    z=v*(1<<(2*q));j=math.isqrt(z.numerator//z.denominator);x,y=j*d,(j+1)*d
    return (v+x*y)/(x+y)
def tail_bound(rho,u,v,r):
    u=min(u,rho);v=min(v,rho*u,rho*rho)
    if rho<=r or u==0 or v==0:return F(0)
    if r==0 and u*u<=v:return u
    if r==0 or(v>4*r*r and u*u>=v):return clip_root(v)-r
    if rho<2*r:return(rho-r)*v/(rho*rho)
    if v<=4*r*r and u>=v/(2*r):return v/(4*r)
    assert 2*r<=v/u<=rho
    return u-r*u*u/v

@lru_cache(None)
def binomial(n,s):
    return sum((F(math.comb(n,j))*(1-s)**j*s**(n-j)for j in range((n+1)//2,n+1)),F(0))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bounds(x):return F(x['lower_exact_fraction']),F(x['upper_exact_fraction'])
def key(r):return tuple(r[k]for k in('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits'))
def expand(raw):
    out=Counter()
    for gate,n in raw.items():
        if gate=='CCX':out.update({'H':2*n,'CX':6*n,'T':4*n,'Tdg':3*n})
        else:out[{'TD':'Tdg','SD':'Sdg'}.get(gate,gate)]+=n
    return out


def main_audit(root):
    sys.path.insert(0,str(root))
    import value_oracle as vo
    import quantum_gates as qg
    # The predicate compiler is the independently audited existing Boolean
    # graph emitter. Financial A/CQ aggregate helpers are not imported.
    @lru_cache(None)
    def predicate(kind,q,c):
        g=vo.ValueCircuit(store_gates=False,mode='exact7t');value=g.alloc(q+1)
        width=(c+1 if c<q else 0)if kind=='bucket' else c
        comparison=g.alloc(width);flag=g.alloc(1)[0]
        if kind=='bucket':vo.bucket_predicate_into(g,value,comparison,c,flag)
        elif kind=='clipped':vo.clipped_range_predicate_into(g,value,comparison,c,flag)
        else:vo.range_predicate_into(g,value,comparison,c,flag)
        metrics=g.resources()
        private=g.wire_count-(q+1+width+1)
        return expand(metrics['primitive_gate_counts']),width,private,metrics['explicit_7T_sequence_T_depth'],metrics['conservative_5layer_CCX_T_depth_upper']
    @lru_cache(None)
    def old_reflection_depth(n):
        gates=qg.zero_phase(tuple(range(n)),tuple(range(n,2*n-2)))
        return qg.resources(qg.expand_toffoli(gates),max(n,2*n-2))['T_depth_assigned_3layer_Toffoli_or_explicit_T_upper']
    tapes=json.loads((root/'results/qft_rotation_tapes.json').read_text())
    cert=json.loads((root/'results/qft_rotation_certificates.json').read_text())
    assert cert['tapes_sha256']==sha(root/'results/qft_rotation_tapes.json')
    rotation_norm={r['key']:F(r['operator_norm_upper_rational'])for r in cert['certificates']}
    @lru_cache(None)
    def fourier(m):
        out=Counter({'H':m,'CX':3*(m//2)});error=F(0)
        for distance in range(1,m):
            k=distance+1;neg=tapes['rotations'][f'-1:{k}']['chronological_gates']
            pos=tapes['rotations'][f'1:{k}']['chronological_gates'];copies=m-distance
            if neg or pos:
                out.update({g:copies*(2*Counter(neg)[g]+Counter(pos)[g])for g in set(neg)|set(pos)})
                out['CX']+=2*copies
            error+=copies*(2*rotation_norm[f'-1:{k}']+rotation_norm[f'1:{k}'])
        rounded=F(ceil(error*(1<<80)),1<<80)
        assert error<=rounded<=F(1,64)
        return out,error,rounded
    @lru_cache(None)
    def primitive_profiles(candidate_id,architecture,kind,q,c):
        row=raw_map[(candidate_id,architecture)];body=Counter()
        for block in row['resources']['blocks']:body.update(block['primitive_gate_counts'])
        raw_body=body.copy();body=expand(body);pred,width,private,p4,p5=predicate(kind,q,c)
        L=row['resources']['input_qubits'];n=L+width+1
        ac=Counter({g:2*z for g,z in body.items()});ac.update(pred);ac['H']+=L+width
        cq=Counter({g:2*z for g,z in ac.items()})
        # n original logical wires, one positive phase control, n-1 scratch.
        cq.update({'X':2*n,'H':4*(n-1),'CX':12*(n-1),'T':8*(n-1),
                   'Tdg':6*(n-1),'CZ':2,'Z':1})
        body_peak=row['resources']['logical_qubits_serial_upper_bound']-(q+2)
        peak=max(body_peak+width+1+private,2*n-2)
        a4=2*(raw_body['T']+raw_body['TD']+4*raw_body['CCX'])+p4
        a5=2*(raw_body['T']+raw_body['TD']+5*raw_body['CCX'])+p5
        return ac,cq,L,width,n,peak,a4,a5
    report={'schema':'qbp-independent-complete-qpe-benchmark-audit-v1','status':'passed',
            'scalar_records_checked':0,'scalar_zero_shots_checked':0,
            'single_flag_candidates_checked':0,'support_rejections_checked':0,
            'tail_rejections_checked':0,'unsupported_records':[],
            'dyadic_candidates_checked':0,'bucket_records_checked':0,
            'exact_median_confidence_pairs_checked':set(),'fourier_widths_checked':set()}
    started=time.perf_counter()
    with gzip.open(root/'results/qpe_benchmark.json.gz','rt')as f:data=json.load(f)
    summary=json.loads((root/'results/qpe_benchmark_summary.json').read_text())
    assert summary['archive_sha256']==sha(root/'results/qpe_benchmark.json.gz')
    assert summary['compressed_bytes']==(root/'results/qpe_benchmark.json.gz').stat().st_size
    h=hashlib.sha256();length=0
    with gzip.open(root/'results/qpe_benchmark.json.gz','rb')as f:
        while chunk:=f.read(1<<20):h.update(chunk);length+=len(chunk)
    assert summary['uncompressed_sha256']==h.hexdigest()and summary['uncompressed_bytes']==length
    for name,digest in data['source_hashes'].items():assert sha(root/name)==digest,(name,'source hash')
    for name,digest in data['input_hashes'].items():assert sha(root/name)==digest,(name,'input hash')
    raw=json.loads((root/'results/architecture_resource_candidates.json').read_text())['configurations']
    raw_map={(r['candidate_id'],r['architecture']):r for r in raw}
    frozen={key(r):r for r in json.loads((root/'results/finite_range_certificates.json').read_text())['configurations']}
    assert len(data['configurations'])==len(raw)==819
    def scalar(rec,row,kind,c,epsilon,prior,alpha):
        report['scalar_records_checked']+=1
        if rec['status']=='prior_zero_shot':
            assert prior<=2*epsilon;assert F(rec['estimate'])==prior/2
            assert rec['T']==rec['runs']==rec['controlled_Q_calls']==0
            assert not rec['metrics'];report['scalar_zero_shots_checked']+=1;return
        assert rec['status']=='estimated',rec['status'];assert prior>2*epsilon
        assert F(rec['epsilon'])==epsilon and F(rec['prior'])==prior and F(rec['alpha'])==alpha
        p=rec['plan'];M=p['M'];m=p['phase_register_bits'];R=p['repetitions']
        assert p['confidence_rule']=='exact_binomial'
        assert M==1<<m and rec['phase_register_bits']==m
        var=prior*(1-prior)if prior<=F(1,2)else F(1,4)
        err=2*PI*upper_sqrt(var)/M+PI*PI/(M*M)
        assert F(p['variance_upper'])==var and F(p['amplitude_error_upper'])==err
        assert err+XI<=epsilon
        assert 2*PI*upper_sqrt(var)/(M//2)+PI*PI/(M//2)**2+XI>epsilon
        assert F(p['epsilon'])==epsilon and F(p['prior_upper'])==prior and F(p['alpha'])==alpha
        assert F(p['postprocess_allowance'])==XI
        ft,actual_error,rounded=fourier(m);s=8/(PI*PI)-rounded
        assert F(rec['qft_norm_error_upper'])==actual_error
        assert F(p['qft_norm_error_allowance'])==rounded and F(p['good_probability_lower'])==s
        assert R%2==1 and R>=1 and binomial(R,s)<=alpha
        if R>1:assert binomial(R-2,s)>alpha
        assert F(p['binomial_failure_upper'])==binomial(R,s)
        if(m,str(alpha),R)not in report['exact_median_confidence_pairs_checked']:
            assert all(binomial(j,s)>alpha for j in range(1,R,2))
        report['exact_median_confidence_pairs_checked'].add((m,str(alpha),R))
        report['fourier_widths_checked'].add(m)
        ac,cq,L,width,n,peak,a4,a5=primitive_profiles(row['candidate_id'],row['architecture'],kind,row['q'],c)
        aT=ac['T']+ac['Tdg'];qT=cq['T']+cq['Tdg'];fT=ft['T']+ft['Tdg']
        expected=Counter({g:R*(ac[g]+(M-1)*cq[g]+ft[g])for g in set(ac)|set(cq)|set(ft)})
        expected['H']+=R*m
        metric=rec['metrics'];actual=Counter({k[5:]:z for k,z in metric.items()if k.startswith('gate_')})
        assert actual==expected
        assert rec['A_T']==aT and rec['CQ_T']==qT and rec['qft_T']==fT
        assert rec['T']==metric['T']==R*(aT+(M-1)*qT+fT)
        assert rec['controlled_Q_calls']==R*(M-1)and rec['runs']==R
        assert metric['expanded_primitive_operations']==sum(expected.values())
        assert metric['measurements']==expected['MZ']+R*m
        assert metric['feedback_CZ']==expected['CZ_IF']and metric['reset_X']==expected['X_IF']
        assert metric['uniform_H']==R*((L+width)*(1+2*(M-1))+m)
        assert metric['all_measurements_with_reinitialisation_upper']==expected['MZ']+R*(m+n)
        assert metric['reinitialisation_additional_Z_measurements_upper']==R*n
        assert metric['reinitialisation_conditional_X_upper']==R*(m+n)
        assert metric['complete_primitive_operations_with_readout_and_reset_upper']==sum(expected.values())+2*R*(m+n)
        assert rec['phase_register_measurements']==R*m and rec['complete_reinitialisations']==R
        assert rec['active_financial_qubits']==n
        assert rec['peak_logical_qubits']==peak+2+m
        assert metric['explicit_T_depth_serial_upper']==R*(a4+(M-1)*(2*a4+old_reflection_depth(n)+8)+fT)
        assert metric['conservative_5layer_T_depth_serial_upper']==R*(a5+(M-1)*(2*a5+10*(n-2)+10)+fT)
    def choose(records):return min((x for x in records if 'T'in x),key=lambda x:(x['T'],x.get('runs',0),x.get('c',99)))
    for row in data['configurations']:
        original=raw_map[(row['candidate_id'],row['architecture'])];bound=frozen[key(row)]
        q=row['q'];maximum=bound['max_integer'];u=F(bound['fair_mean_upper_fraction']);v=F(bound['fair_second_moment_upper_fraction'])
        rho=F(maximum,1<<q);e=F(1,200*row['cap'])
        assert F(row['mean_upper'])==u and F(row['second_moment_upper'])==v
        assert F(row['range'])==rho and row['max_integer']==maximum
        assert F(row['estimation_tolerance'])==e and F(row['alpha'])==ALPHA
        assert row['deterministic_price_error_upper']==original['independent_deterministic_bias_upper_decimal']
        assert F(row['deterministic_price_error_upper'])<=F(1,200)
        assert len(row['range_candidates'])==len(row['clipping_candidates'])==q+1
        for kind,records in(('range',row['range_candidates']),('clipped',row['clipping_candidates'])):
            for c,rec in enumerate(records):
                report['single_flag_candidates_checked']+=1
                assert rec['kind']==kind and rec['c']==c
                if kind=='range'and maximum>1<<c:
                    assert rec['status']=='invalid_support';report['support_rejections_checked']+=1;continue
                D=tail_bound(rho,u,v,F(1<<c,1<<q))if kind=='clipped'else F(0)
                assert F(rec['tail_upper'])==D
                if D/2>=e:
                    assert rec['status']=='tail_exhausts_budget';report['tail_rejections_checked']+=1;continue
                r=F(1<<c,1<<q);mean=min(u,r,sharp_mean(q,v,high=min(maximum,1<<c)))if maximum else F(0)
                prior=min(F(1),mean/r);epsilon=(e-D/2)/r
                if rec['status']=='unsupported_scalar_plan':
                    assert not F(1,10**6)<=epsilon<F(1,2)
                    report['unsupported_records'].append((row['candidate_id'],kind,c));continue
                scalar(rec,row,kind,c,epsilon,prior,ALPHA)
                assert F(rec['level'])==r and F(rec['tail_midpoint'])==D/2
                assert F(rec['total_mean_error_upper'])==r*epsilon+D/2==e
        assert row['plain']==row['range_candidates'][q]
        for d in row['dyadic_candidates']:
            assert d['status']=='complete';report['dyadic_candidates_checked']+=1
            s=d['omitted_prefix_bits'];tail=min(F(min((1<<s)-1,maximum),1<<q),u,upper_sqrt(v))
            assert F(d['omitted_mean_upper'])==tail and F(d['omitted_midpoint'])==tail/2
            assert tail/2<=e/4
            if s<q:
                nexttail=min(F(min((1<<(s+1))-1,maximum),1<<q),u,upper_sqrt(v))
                assert nexttail/2>e/4
            assert[b['index']for b in d['buckets']]==list(range(s,q+1))
            active=[];totalmetrics=Counter();totalT=totalR=totalQ=peak=0;error=tail/2
            for b in d['buckets']:
                j=b['index'];weight=F(1)if j==q else F(1<<(j+1),1<<q)
                highest=min((1<<(j+1))-1,maximum)
                prior=(min(v,u)if maximum==1<<q else F(0))if j==q else(F(0)if highest<1<<j else min(F(highest,1<<(j+1)),2*v/(weight*weight),u/weight,sharp_mean(q,v,1<<j,highest)/weight))
                assert F(b['weight'])==weight and F(b['prior'])==prior
                allowance=F(b['mean_error_allowance']);epsilon=allowance/weight
                assert F(b['epsilon'])==epsilon and 0<=allowance<=weight*prior/2
                error+=allowance
                if allowance<weight*prior/2:active.append(b)
            budget=e-tail/2;caps=[F(b['weight'])*F(b['prior'])/2 for b in d['buckets']]
            floors=[min(cap,F(b['weight'])*F(1,10**6))if d['allocation']!='equal'else F(0)
                    for b,cap in zip(d['buckets'],caps)]
            weights=[]
            for b in d['buckets']:
                if d['allocation']=='equal':weights.append(F(1));continue
                prior=F(b['prior']);w=F(b['weight'])
                if prior==0:weights.append(F(1,1<<100));continue
                ac,cq,*_=primitive_profiles(row['candidate_id'],row['architecture'],'bucket',q,b['index'])
                cost=max(1,ac['T']+ac['Tdg']+cq['T']+cq['Tdg'])
                weights.append(max(F(1,1<<100),upper_sqrt(cost*w*upper_sqrt(prior))))
            values=[F(b['mean_error_allowance'])for b in d['buckets']]
            assert all(floor<=value<=cap for floor,value,cap in zip(floors,values,caps))
            if sum(caps,F(0))<=budget:assert values==caps
            else:
                assert sum(values,F(0))==budget
                multipliers=[(value-floor)/weight for value,floor,weight,cap in zip(values,floors,weights,caps)if value<cap]
                assert multipliers and len(set(multipliers))==1
                lam=multipliers[0]
                assert all(cap-floor<=lam*weight for value,floor,weight,cap in zip(values,floors,weights,caps)if value==cap)
            for b in d['buckets']:
                j=b['index'];prior=F(b['prior']);epsilon=F(b['epsilon']);alpha=ALPHA/len(active)if b in active else F(0)
                scalar(b,row,'bucket',j,epsilon,prior,alpha)
                assert F(b['alpha'])==alpha
                totalT+=b['T'];totalR+=b['runs'];totalQ+=b['controlled_Q_calls'];peak=max(peak,b['peak_logical_qubits'])
                totalmetrics.update(b['metrics']);report['bucket_records_checked']+=1
            assert sum((F(b['alpha'])for b in d['buckets']),F(0))<=ALPHA
            assert error==F(d['total_mean_error_upper'])<=e
            assert d['T']==totalT and d['runs']==totalR and d['controlled_Q_calls']==totalQ
            assert d['peak_logical_qubits']==peak and Counter(d['metrics'])==totalmetrics
        assert row['best_range']==choose(row['range_candidates'])
        assert row['best_clipped']==choose(row['clipping_candidates'])
        assert row['best_dyadic']==choose(row['dyadic_candidates'])
        assert row['best_single_flag']==choose([row['plain'],row['best_range'],row['best_clipped']])
        assert row['best_all']==choose([row['best_single_flag'],row['best_dyadic']])
    selections={}
    for case in sorted({r['case_id']for r in data['configurations']}):
        pool=[r for r in data['configurations']if r['case_id']==case]
        for family in('plain','best_range','best_clipped','best_single_flag','best_dyadic','best_all'):
            winner=min(pool,key=lambda r:(r[family]['T'],r[family].get('runs',0),r['candidate_id'],r['architecture']))
            selections[(case,family)]=winner
    assert len(data['selected_cases'])==162
    for rec in data['selected_cases']:
        winner=selections[(rec['case_id'],rec['family'])];w=winner[rec['family']]
        for k in('candidate_id','architecture','method'):assert rec[k]==winner[k]
        for k in('T','runs','controlled_Q_calls','peak_logical_qubits','kind'):assert rec[k]==w[k]
        assert rec['c']==w.get('c')and rec['allocation']==w.get('allocation')
    assert summary['selected_cases']==data['selected_cases']
    assert summary['summary']==data['summary']
    for family,rec in data['summary']['ratios_to_QPE_plain'].items():
        ratios=[F(r[family]['T'],selections[(case,'plain')]['plain']['T'])for(case,f),r in selections.items()if f==family]
        assert rec['improvements']==sum(x<1 for x in ratios)
        assert rec['ties']==sum(x==1 for x in ratios)and rec['losses']==sum(x>1 for x in ratios)
        assert rec['minimum_ratio']==float(min(ratios))and rec['maximum_ratio']==float(max(ratios))
    report['financial_profiles_checked']=819;report['selected_minima_checked']=162
    report['exact_median_confidence_pairs_checked']=len(report['exact_median_confidence_pairs_checked'])
    report['fourier_widths_checked']=sorted(report['fourier_widths_checked'])
    report['source_hashes_verified']=len(data['source_hashes']);report['input_hashes_verified']=len(data['input_hashes'])
    report['elapsed_seconds']=time.perf_counter()-started
    report['scope']='independent finite formulas and primitive-count reconstruction; no benchmark or planner import, no physical quantum execution'
    report['archive_sha256']=summary['archive_sha256']
    report['auditor_sha256']=sha(Path(__file__))
    return report


def taylor_audit(root):
    sys.path.insert(0,str(root))
    from taylor_pricing import model_parameters
    from flint import arb,ctx
    from continuous_quadrature import exact_endpoint
    data=json.loads((root/'results/taylor_order_benchmark.json').read_text())
    protocol=json.loads((root/'results/generalisation_protocol.json').read_text())
    assert data['protocol_sha256']==sha(root/'results/generalisation_protocol.json')
    for name,digest in data['source_hashes'].items():assert sha(root/name)==digest
    models={m['model_id']:m for m in protocol['models']}
    old={r['case_id']:r for r in json.loads((root/'results/taylor_pricing.json').read_text())['configurations']}
    quad={r['case_id']:r for r in json.loads((root/'results/continuous_quadrature.json').read_text())['configurations']}
    full={};counters={};pricechecks=[];monotonechecks=[]
    def fresh_monotone(rec):
        model=models[rec['case_id'].rsplit('_cap',1)[0]]
        w,d,s,mg,rg,lam=model_parameters(model);H=rec['cap'];saved=ctx.prec
        ctx.prec=384
        try:
            def a(x):
                x=F(x);return arb(x.numerator)/arb(x.denominator)
            first=sum((a(wi)*(a(di)+a(sum(z*z for z in si))/2).exp()
                       for wi,di,si in zip(w,d,s)),arb(0))
            assert exact_endpoint(first,False)>1
            L=1/(2*arb.pi()).sqrt()/a(lam)
            cap_slope=F(100+H,H)
            if exact_endpoint(L,False)>=cap_slope:L=a(cap_slope)
            elif exact_endpoint(L,True)>cap_slope:
                # This branch would require an interval minimum; none of
                # the frozen models straddles the rational slope crossover.
                raise AssertionError('Independent slope crossover interval unresolved')
            gap=a(H)*L*first.log();gap_upper=exact_endpoint(gap,True)
            sigma=a(lam*lam+sum(z*z for z in rg)).sqrt();mu=a(mg)
            Phi=lambda x:(1+(x/arb(2).sqrt()).erf())/2
            def call(K):
                x=(mu-(a(K)/100).log())/sigma
                return 100*(mu+sigma*sigma/2).exp()*Phi(x+sigma)-a(K)*Phi(x)
            geo=call(100)-call(100+H)
            glo,ghi=exact_endpoint(geo,False),exact_endpoint(geo,True)
        finally:ctx.prec=saved
        recorded=F(rec['selected']['remainder_upper'])if rec['selected'].get('one_sided_monotonicity')else F(rec['one_sided_screen_upper_exact_fraction'])
        assert recorded>=gap_upper>=0
        selected_zero=rec['selected'].get('one_sided_monotonicity',False)
        simple=rec['selected']if selected_zero else rec.get('one_sided_preliminary_interval')
        if simple is not None:
            assert recorded<=F(1,50)
            lo,hi=bounds(simple['interval'])
            assert lo<=max(F(0),glo)
            assert hi>=min(F(H),ghi+gap_upper)
            assert simple['certified']==((hi-lo)/2<=F(1,100))
            if not selected_zero:assert not simple['certified']
        else:assert recorded>F(1,50)
        monotonechecks.append({'case_id':rec['case_id'],
            'selected_one_sided_control':selected_zero,'recorded_gap_upper':str(recorded),
            'independent_384bit_gap_upper':str(gap_upper),
            'independent_geometric_price_interval':{'lower':str(glo),'upper':str(ghi)},
            'scalar_preliminary_evaluated':simple is not None,
            'first_moment_and_control_formula_checked':True})
        return simple is not None
    def geometry_key(record):
        model=models[record['case_id'].rsplit('_cap',1)[0]]
        w,d,s,mg,rg,lam=model_parameters(model)
        return w,d,s,rg
    def check_record(rec,bank_degree):
        n=rec['assets'];limit=data['maximum_order'];key=geometry_key(rec)
        count=lambda k:math.comb(n+k,k)
        preliminary=False
        if rec['policy']=='monotone_then_screen':preliminary=fresh_monotone(rec)
        if rec['policy']=='monotone_then_screen'and rec['selected'].get('one_sided_monotonicity'):
            selected=rec['selected'];required=1;calls=1
            assert rec['screens']==[]and rec['evaluated_orders']==[selected]
            assert selected['certified']and selected['order']==0
        elif rec['policy']=='full_search':
            required=limit+1;calls=count(limit)
            assert[o['order']for o in rec['evaluated_orders']]==list(range(limit+1))
            good=[o for o in rec['evaluated_orders']if o['certified']]
            selected=min(good or rec['evaluated_orders'],key=lambda o:bounds(o['interval'])[1]-bounds(o['interval'])[0])
        else:
            screens=rec['screens'];assert[s['order']for s in screens]==list(range(screens[-1]['order']+1))
            evaluated=[];calls=0
            for s in screens:
                qualified=F(s['remainder_upper'])<=F(1,100)
                assert s['remainder_qualified']==qualified
                compute=qualified or s['order']==limit
                assert s['coefficient_evaluated']==compute
                if compute:evaluated.append(s['order']);calls+=count(s['order'])
                full_rem=full[rec['case_id']]['evaluated_orders'][s['order']]['remainder_upper']if rec['case_id']in full else None
                if full_rem is not None:assert F(s['remainder_upper'])==F(full_rem)
            assert[o['order']for o in rec['evaluated_orders']]==evaluated
            good=[o for o in rec['evaluated_orders']if o['certified']]
            if good:
                selected=good[0];assert screens[-1]['order']==selected['order']
            else:
                selected=rec['evaluated_orders'][-1];assert screens[-1]['order']==limit
            required=screens[-1]['order']+1
            if rec['policy']=='monotone_then_screen':calls+=int(preliminary)
        assert rec['selected']==selected and rec['kernel_calls']==calls
        old_degree=bank_degree.get(key,-1);next_degree=max(old_degree,required)
        new=count(next_degree)-(count(old_degree)if old_degree>=0 else 0)
        assert rec['new_geometry_exponentials']==new
        bank_degree[key]=next_degree
        for o in rec['evaluated_orders']:
            lo,hi=bounds(o['interval']);assert 0<=lo<=hi<=rec['cap']
            assert o['certified']==((hi-lo)/2<=F(1,100))
        return calls,new
    for family in data['family_runs']:
        bank={};calls=geometry=0
        for r in family['configurations']:
            c,g=check_record(r,bank);calls+=c;geometry+=g
            if family['policy']=='full_search':full[r['case_id']]=r
            else:
                lo,hi=bounds(r['selected']['interval']);olo,ohi=bounds(old[r['case_id']]['selected_undiscounted_price_interval'])
                assert max(lo,olo)<=min(hi,ohi)
                assert r['selected']['certified']==old[r['case_id']]['certified_within_001']
                flo,fhi=bounds(full[r['case_id']]['selected']['interval']);assert max(lo,flo)<=min(hi,fhi)
                contains=None
                if r['case_id']in quad:
                    qlo,qhi=bounds(quad[r['case_id']]['undiscounted_price_interval'])
                    assert lo<=qlo<=qhi<=hi;contains=True
                pricechecks.append({'case_id':r['case_id'],'old_Taylor_enclosure_overlaps':True,
                    'high_precision_quadrature_interval_contained':contains})
        assert family['kernel_calls']==calls and family['geometry_exponentials']==geometry
        assert family['certified_cases']==sum(r['selected']['certified']for r in family['configurations'])
        counters[family['policy']]={'kernel_calls':calls,'geometry_exponentials':geometry,'certified_cases':family['certified_cases']}
    for r in data['cold_standalone_presets']:check_record(r,{})
    assert all(x['certified_cases']==25 for x in counters.values())
    assert F(data['summary']['screened_to_full_kernel_ratio'])==F(counters['remainder_screen']['kernel_calls'],counters['full_search']['kernel_calls'])
    assert F(data['summary']['hierarchy_to_full_kernel_ratio'])==F(counters['monotone_then_screen']['kernel_calls'],counters['full_search']['kernel_calls'])
    assert data['summary']['monotone_screen_certificates']==sum(r['selected'].get('one_sided_monotonicity',False)for r in data['family_runs'][2]['configurations'])==3
    return {'status':'passed','family_contracts':27,'family_policy_runs':81,
            'cold_standalone_runs':len(data['cold_standalone_presets']),
            'counters':counters,'price_comparison_checks':pricechecks,
            'high_precision_quadrature_containment_checks':sum(r['high_precision_quadrature_interval_contained']is True for r in pricechecks),
            'monotone_control_checks':monotonechecks,
            'monotone_first_moment_and_control_precision_bits':384,
            'scope':'result and source-hash audit of finite order policies and exact combinatorial work; fresh 384-bit direct first moments and independent two-call Black-Scholes differences audit the hierarchy; existing validated Taylor coefficient/remainder theorem retained; no new kernel integration performed'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=main_audit(args.root);result['taylor_order_audit']=taylor_audit(args.root)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items()if k!='taylor_order_audit'},indent=2))
    print(json.dumps(result['taylor_order_audit']['counters'],indent=2))
