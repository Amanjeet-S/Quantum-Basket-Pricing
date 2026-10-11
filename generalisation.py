"""Reversible financial graphs and precision search for frozen synthetic models."""
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import encoding_resources as e
from resource_study import propagate, coupling_bound
import general_kernels as kernels

D=Decimal
ROOT=Path(__file__).resolve().parent


def weighted_sum(g,values,weights):
    terms=[g.scale(v,w) for v,w in zip(values,weights)]
    result=terms[0]
    for term in terms[1:]:result=g.add(result,term)
    return result


def affine(g,y,coefficients,drift,bias,z=None,lam=None):
    with localcontext() as ctx:
        ctx.prec=200
        constant=D(bias)+drift-5*sum(coefficients)
        if z is not None:constant-=5*lam
    result=g.const(constant)
    for yy,coefficient in zip(y,coefficients):
        if coefficient:
            term=g.scale(yy,abs(coefficient))
            result=g.add(result,term) if coefficient>0 else g.sub(result,term)
    if z is not None:result=g.add(result,g.scale(z,lam))
    return result


def build_graph(model,method,q,by,bz,bu,exp_bits,log_bits,prob_bits,cap):
    if q<max(exp_bits,log_bits,prob_bits,bu+1):raise ValueError('Insufficient fractional width')
    g=e.Graph(q);g.width=q+7
    g.kernel_bound_mode='cap_aware' if method=='full_cap' else 'global'
    w=list(map(D,model['weights']));r=[list(map(D,rr)) for rr in model['loadings']]
    mu=list(map(D,model['log_drifts']));lam=D(model['lambda']);h=D(cap)/100
    y=[g.op('quantile',g.inp(f'iy{j}',by),bits=by) for j in range(model['outer_dimension'])]
    with localcontext() as ctx:
        ctx.prec=200
        mg=sum(a*b for a,b in zip(w,mu))
        rg=[sum(w[i]*r[i][j] for i in range(model['assets'])) for j in range(len(y))]
    z=g.op('quantile',g.inp('iz',bz),bits=bz) if method=='direct' else None
    logs=[affine(g,y,rr,m,4,z,lam) for rr,m in zip(r,mu)]
    prices=[g.function('general_exp',m,exp_bits) for m in logs]
    ba=weighted_sum(g,prices,w)
    if method=='direct':
        geom=g.function('general_exp',affine(g,y,rg,mg,4,z,lam),exp_bits)
        res=g.clamp(g.sub(g.payoff(ba,h),g.payoff(geom,h)),h)
        g.output=g.clamp(g.scale(res,D(1)/h),1)
    else:
        ma=g.sub(g.function('general_log',g.sub(ba,g.const('.0625')),log_bits),g.const(2))
        geom=affine(g,y,rg,mg,2)
        if method=='interval':
            u=g.op('uniform',g.inp('iu',bu),bits=bu)
            strike=g.function(f'general_strike:{cap}',u,log_bits)
            probs=[g.function(f'general_cdf:{lam}',g.sub(m,strike),prob_bits) for m in (ma,geom)]
        elif method in ('full_global','full_cap'):
            probs=[g.function(f'general_cap:{lam}:{cap}',m,prob_bits) for m in (ma,geom)]
        else:raise ValueError(method)
        g.output=g.clamp(g.sub(*probs),1)
    return g


def make_tables(graph,cache,certify=False):
    audits=[]
    for node in graph.nodes.values():
        if node.kind not in ('quantile','function'):continue
        kind='quantile' if node.kind=='quantile' else node.meta
        key=(kind,node.bits,graph.q)
        if key in cache:continue
        if kind=='quantile':
            cache[key],residual=e.quantile_table(node.bits,graph.q,certify)
            audits.append({'bits':node.bits,'q':graph.q,'entries':len(cache[key]),'exhaustive_inverse_enclosure':certify,'largest_cdf_residual':residual})
        else:cache[key]=kernels.function_table(kind,node.bits,graph.q,graph.width)
    return audits


def distribution(model,method,cap,q,epsilon):
    with localcontext() as ctx:
        ctx.prec=70
        r=[list(map(D,rr)) for rr in model['loadings']];w=list(map(D,model['weights']))
        rg=[sum(w[i]*r[i][j] for i in range(len(w))) for j in range(len(r[0]))]
        financial_lip,_=kernels.probability_constants(D(model['lambda']),D(cap)/100,True)
        shared_log_price_lip=min(100+D(cap),D(cap)*financial_lip)
        # All methods approximate the same continuously inner-averaged price.
        # Use its strongest checked outer-coordinate bound for every opponent.
        lips=[shared_log_price_lip*(max(abs(rr[j]) for rr in r)+abs(rg[j])) for j in range(len(rg))]
        allowance=D(str(epsilon))*D('.4')-D('1e-12')
        bu=q-1 if method=='interval' else 0;choices=[]
        for by in range(6,19):
            outer=sum(min(D('1.5')*cap/D(1<<by),lip*coupling_bound(by)) for lip in lips)
            if method=='direct':
                for bz in range(6,19):
                    inner=min(D(cap)/D(1<<bz),2*D(model['lambda'])*(100+D(cap))*coupling_bound(bz))
                    if outer+inner<=allowance:choices.append((len(rg)*(1<<by)+(1<<bz),by,bz,bu,outer,inner))
            else:
                inner=D(cap)/D(1<<bu) if method=='interval' else D(0)
                if outer+inner<=allowance:choices.append((len(rg)*(1<<by),by,0,bu,outer,inner))
        if not choices:return None
        _,by,bz,bu,outer,inner=min(choices)
        return {'b_y':by,'b_z':bz,'b_u':bu,
            'outer_distribution_bound':math.nextafter(float(outer),math.inf),
            'inner_distribution_bound':math.nextafter(float(inner),math.inf),
            'distribution_price_bound':math.nextafter(float(outer+inner),math.inf)}


def proxy(model,method,q,p):
    n=model['assets'];exps=n+(method=='direct');logs=0 if method=='direct' else 1+(method=='interval')
    probs=0 if method=='direct' else 2
    return 28*(exps*(1<<p[0])+logs*(1<<p[1])+probs*(1<<p[2]))+28*(n+probs+logs)*(q+7)**2


def feasible_candidates(model,method,cap,epsilon):
    found=[];budget=epsilon*.1
    for q in range(18,31):
        dist=distribution(model,method,cap,q,epsilon)
        if dist is None:continue
        p=[4,0,0] if method=='direct' else [4,4,4]
        def assess(pp):
            g=build_graph(model,method,q,dist['b_y'],dist['b_z'],dist['b_u'],*pp,cap)
            cc=propagate(g,cap)
            return cc['arithmetic_price_error_upper'],cc
        for _ in range(12):
            try:current,cert=assess(p);break
            except AssertionError as exc:
                if exc.args[0][0] not in ('function address out of domain','register overflow certificate'):raise
                p[0]+=1
                if p[1]:p[1]+=1
        else:raise AssertionError('Unable to certify initial prefix')
        for _ in range(45):
            if current<=budget:break
            nexts=[]
            for i in range(3):
                if not p[i] or p[i]>=min(15,q-1):continue
                pp=p.copy();pp[i]+=1
                error,cc=assess(pp)
                extra=proxy(model,method,q,pp)-proxy(model,method,q,p)
                nexts.append(((current-error)/max(extra,1),error,cc,pp))
            if not nexts:break
            _,current,cert,p=max(nexts,key=lambda t:t[0])
        if current<=budget:
            found.append({'q':q,'precision':p,'distribution':dist,'certificate':cert,'proxy':proxy(model,method,q,p)})
    if not found:raise ValueError(('No feasible configuration',model['model_id'],method,cap))
    return sorted(found,key=lambda t:t['proxy'])


def select(model,method,cap,epsilon,tables,blocks,global_winner=None):
    feasible=feasible_candidates(model,method,cap,epsilon)
    shortlist=feasible[:2]
    if global_winner is not None:
        p=[global_winner[k] for k in ('exp_address_bits','log_address_bits','probability_address_bits')]
        q=global_winner['q'];dist=distribution(model,method,cap,q,epsilon)
        shortlist.append({'q':q,'precision':p,'distribution':dist,'proxy':proxy(model,method,q,p),'certificate':None})
    candidates=[]
    for row in shortlist:
        q=row['q'];p=row['precision'];dist=row['distribution']
        graph=build_graph(model,method,q,dist['b_y'],dist['b_z'],dist['b_u'],*p,cap)
        cert=propagate(graph,cap)
        assert cert['arithmetic_price_error_upper']<=epsilon*.1
        make_tables(graph,tables)
        resources=e.counted_composition(graph,tables,blocks)
        candidates.append((resources['logical_T_count_exact_7T_Toffoli'],q,p,dist,cert,resources))
    _,q,p,dist,cert,resources=min(candidates,key=lambda t:(t[0],t[1]))
    return {'model_id':model['model_id'],'case_id':f"{model['model_id']}_cap{cap}",
        'method':method,'cap':cap,'price_tolerance':epsilon,'q':q,
        'exp_address_bits':p[0],'log_address_bits':p[1],'probability_address_bits':p[2],
        **dist,**cert,'resources':resources,
        'bias_bound_including_correction':math.nextafter(dist['distribution_price_bound']+cert['arithmetic_price_error_upper']+1e-12,math.inf),
        'reserved_statistical_price_allowance':epsilon*.5,
        'feasible_q_configurations':len(feasible),
        'exact_counted_candidates':[{'T':t[0],'q':t[1],'function_address_bits':t[2]} for t in candidates],
        'selection':'best emitted T count among two proxy-shortlisted feasible candidates; cap-aware also includes the global winner; not globally optimal'}


def geometric_correction(model,cap):
    with localcontext() as ctx:
        ctx.prec=70
        w=list(map(D,model['weights']));mu=sum(a*D(b) for a,b in zip(w,model['log_drifts']))
        variance=D(model['lambda'])**2+D(model['geometric_outer_variance']);sd=variance.sqrt();h=D(cap)/100
        x0=mu/sd;x1=(mu-(1+h).ln())/sd
        val=100*((mu+variance/2).exp()*(kernels.clipped_cdf(x0+sd)-kernels.clipped_cdf(x1+sd))
                 -(kernels.clipped_cdf(x0)-kernels.clipped_cdf(x1))+h*kernels.clipped_cdf(x1))
        # In this family geometric arguments are all within the ±8 series.
        assert max(map(abs,[x0,x1,x0+sd,x1+sd]))<8
        return {'model_id':model['model_id'],'cap':cap,'undiscounted_value_decimal':str(val),'absolute_allowance':1e-12}


def main():
    p=argparse.ArgumentParser();p.add_argument('--model-ids',nargs='+');p.add_argument('--caps',nargs='+',type=int)
    p.add_argument('--output',type=Path,default=ROOT/'results'/'generalisation_resources.json')
    args=p.parse_args();protocol_path=ROOT/'results'/'generalisation_protocol.json'
    protocol=json.loads(protocol_path.read_text());models=protocol['models']
    if args.model_ids:models=[m for m in models if m['model_id'] in args.model_ids]
    caps=args.caps or protocol['caps'];epsilon=protocol['price_tolerances'][0]
    start=time.perf_counter();tables={};blocks={};rows=[];quantile_audits=[]
    for model in models:
        for cap in caps:
            winner=None
            for method in protocol['methods']:
                row=select(model,method,cap,epsilon,tables,blocks,winner if method=='full_cap' else None)
                if method=='full_global':winner=row
                rows.append(row)
                print('Completed',row['case_id'],method,row['resources']['logical_T_count_exact_7T_Toffoli'],flush=True)
    retained_tables={}
    for row in rows:
        model=next(m for m in models if m['model_id']==row['model_id'])
        graph=build_graph(model,row['method'],row['q'],row['b_y'],row['b_z'],row['b_u'],row['exp_address_bits'],row['log_address_bits'],row['probability_address_bits'],row['cap'])
        quantile_audits+=make_tables(graph,retained_tables,True)
    out={'status':'broader constructed evaluator comparison; estimator and theoretical novelty unverified',
        'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),'models':models,
        'configurations':rows,'data':'stipulated synthetic parameters only',
        'quantile_audits':quantile_audits,'table_manifest':[{'kernel':k[0],'address_bits':k[1],'q':k[2],
            'entries':len(v),'sha256':hashlib.sha256(json.dumps(v,separators=(',',':')).encode()).hexdigest()} for k,v in sorted(retained_tables.items())],
        'classical_workflow_seconds':time.perf_counter()-start,
        'geometric_corrections':[geometric_correction(m,c) for m in models for c in caps],
        'source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ['generalisation.py','generalisation_models.py','general_kernels.py','encoding_resources.py','resource_study.py','reversible.py','reflections.py']}}
    args.output.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
