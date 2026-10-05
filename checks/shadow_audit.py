"""Independent vector integer replay of the retained root graphs.

The graph and its QROM constants are read-only inputs.  This program does not
call Graph.evaluate for its main computations.  It checks every packed function
endpoint using a separately implemented 100-digit CDF enclosure, reconstructs
all integer operations using a split product that cannot overflow uint64, and
compares against positive integral and direct residual references.  A scalar
Graph.evaluate comparison checks that the independent evaluator implements the
retained graph.  The million-sample tests are diagnostics, not uniform proofs.
"""
from decimal import Decimal, localcontext
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import encoding_resources as er
from endpoint_cert import decimal_cdf, rounded_enclosure, EPS
from reference import ncdf, full_positive, model, quantiles

D=Decimal
ENCODING_SOURCE_SHA256=hashlib.sha256((ROOT/'encoding_resources.py').read_bytes()).hexdigest()

def vector_integer(graph,inputs,tables):
    q=graph.q; width=graph.width; mask=np.uint64((1<<width)-1)
    unitmask=np.uint64((1<<q)-1)
    vals={k:np.asarray(v,dtype=np.uint64) for k,v in inputs.items()}
    for name,node in graph.nodes.items():
        a=[vals[k] for k in node.args]
        if node.kind=='constant': v=np.uint64(node.meta)
        elif node.kind=='add': v=a[0]+a[1]
        elif node.kind in ('sub','positive_sub'): v=a[0]-a[1]
        elif node.kind=='shift': v=a[0]>>np.uint64(node.meta)
        elif node.kind=='mul':
            # (ab)>>q exactly, avoiding a full-width overflowing product.
            ah,al=a[0]>>np.uint64(q),a[0]&unitmask
            bh,bl=a[1]>>np.uint64(q),a[1]&unitmask
            v=(ah*bh<<np.uint64(q))+ah*bl+al*bh+(al*bl>>np.uint64(q))
        elif node.kind=='scale':
            integer,drop,_=node.meta
            if drop==0:
                v=a[0]*np.uint64(integer)
            elif drop==q:
                # The scale factor is a compiled classical integer, not an
                # input register. It need not share the signed safe range.
                ah,al=a[0]>>np.uint64(q),a[0]&unitmask
                bh,bl=np.uint64(integer)>>np.uint64(q),np.uint64(integer)&unitmask
                v=(ah*bh<<np.uint64(q))+ah*bl+al*bh+(al*bl>>np.uint64(q))
            else:raise ValueError(('unrecognised constant scale',node.meta))
        elif node.kind=='quantile':
            v=np.asarray(tables[('quantile',node.bits,q)],dtype=np.uint64)[a[0].astype(np.int64)]
        elif node.kind=='uniform':v=(2*a[0]+1)*np.uint64((1<<q)>>(node.bits+1))
        elif node.kind=='function':
            span_bits=2 if node.meta in ('exp','log') else 4
            f=q+span_bits-node.bits
            index=(a[0]>>np.uint64(f)).astype(np.int64)
            rem=a[0]&np.uint64((1<<f)-1)
            packed=tables[(node.meta,node.bits,q)]
            base=np.asarray([x&int(mask) for x in packed],dtype=np.uint64)[index]
            delta=np.asarray([x>>width for x in packed],dtype=np.uint64)[index]
            v=base+(delta*rem>>np.uint64(f))
        elif node.kind=='clamp':
            signed=a[0].astype(np.int64)
            signed=np.where(a[0]<(1<<(width-1)),signed,signed-(1<<width))
            v=np.clip(signed,0,node.meta).astype(np.uint64)
        else:raise ValueError(node.kind)
        vals[name]=v&mask
    return vals[graph.output]

def ideal(kind,iy,iz,iu,by,bz,bu,cap):
    y=quantiles(by)[iy]
    logs=np.stack([-.015-.1*y,-.03+.2*y],axis=-1)
    ba=100*np.exp(logs).mean(axis=-1); bg=100*np.exp(logs.mean(axis=-1))
    if kind=='direct':
        z=quantiles(bz)[iz]
        f=lambda b:np.clip(b*np.exp(.2*z)-100,0,cap)
        return np.maximum(0,f(ba)-f(bg))
    if kind=='interval':
        strike=100+cap*(iu+.5)/(1<<bu)
        return cap*(ncdf(np.log(ba/strike)/.2)-ncdf(np.log(bg/strike)/.2))
    # The positive layer-cake integral is independent of the retained full
    # formula's ordered CDF intervals and signed integer subtraction.
    return full_positive(ba,bg,.2,100,cap,192)

def certify_packed(key,packed):
    kind,bits,q=key; n=1<<bits; width=q+6; mask=(1<<width)-1
    endpoints=[x&mask for x in packed]+[(packed[-1]&mask)+(packed[-1]>>width)]
    assert len(endpoints)==n+1
    with localcontext() as ctx:
        ctx.prec=100
        for j,stored in enumerate(endpoints):
            if kind=='exp':x=D(-2)+D(4)*j/n;value,error=x.exp(),EPS
            elif kind=='log':x=D('.25')+D(4)*j/n;value,error=x.ln()+2,EPS
            elif kind=='cdf':x=D(-8)+D(16)*j/n;value,error=decimal_cdf(x)
            else:raise ValueError(kind)
            assert stored==rounded_enclosure(value,error,q),(kind,bits,q,j)
    assert all(x<=y for x,y in zip(endpoints,endpoints[1:]))
    return {'function':kind,'address_bits':bits,'fraction_bits':q,
            'enclosed_endpoint_count':n+1,'monotone_integer_endpoints':True,
            'independent_decimal_absolute_enclosure':str(EPS)}

def certify_quantile_table(key,stored):
    _,bits,q=key;n=1<<bits;scale=D(1<<q);unit=1/scale
    assert len(stored)==n
    assert all(x+y==10*(1<<q) for x,y in zip(stored,stored[::-1]))
    with localcontext() as ctx:
        ctx.prec=100
        for j in range(n//2,n):
            centre=D(stored[j])/scale-5
            target=(D(j)+D('.5'))/n
            lo,elo=decimal_cdf(centre-unit)
            hi,ehi=decimal_cdf(centre+unit)
            assert lo+elo<=target<=hi-ehi,('independent inverse enclosure',bits,q,j)
    return {'gaussian_address_bits':bits,'fraction_bits':q,'entries':n,
            'exhaustive_positive_half_and_exact_symmetric_negative_half':True,
            'coordinate_error_enclosure':float(unit),
            'independent_decimal_absolute_cdf_enclosure':str(EPS)}


def audit_case(cfg,cache,samples,rng):
    kind=cfg['method'];q=cfg['q'];cap=cfg['cap'];b=(cfg['b_y'],cfg['b_z'],cfg['b_u'])
    p=(cfg['exp_address_bits'],cfg['log_address_bits'],cfg['cdf_address_bits'])
    graph=er.build_graph(kind,q,*b,*p,cap);er.make_tables(graph,cache)
    coarse=er.build_graph(kind,q,6,6,6,*p,cap);er.make_tables(coarse,cache)
    if kind=='full':
        ci={'iy':np.arange(64,dtype=np.int64)}
    else:
        yy,zz=np.meshgrid(np.arange(64),np.arange(64),indexing='ij')
        ci={'iy':yy.ravel(),('iz' if kind=='direct' else 'iu'):zz.ravel()}
    def evaluate(g,inp,bb):
        result=vector_integer(g,inp,cache)
        exact=ideal(kind,inp['iy'],inp.get('iz',np.zeros_like(inp['iy'])),
                    inp.get('iu',np.zeros_like(inp['iy'])),*bb,cap)
        monetary=result.astype(float)*cap/(1<<q)
        diff=np.abs(monetary-exact)
        return result,{'states':len(exact),'maximum_absolute_monetary_error':float(diff.max()),
                      'mean_signed_monetary_error':float((monetary-exact).mean()),
                      'mean_ideal_residual':float(exact.mean()),'mean_fixed_residual':float(monetary.mean())}
    cres,cd=evaluate(coarse,ci,(6,6,6))
    assert cd['maximum_absolute_monetary_error']<=cfg['arithmetic_price_error_upper']
    # Root scalar evaluator is a comparison, not the main implementation.
    count=min(1024,len(cres));ind=np.linspace(0,len(cres)-1,count,dtype=int)
    for j in ind:
        assert int(cres[j])==coarse.evaluate({k:int(v[j]) for k,v in ci.items()},cache)
    # Actual retained grid: full conditioning exhaustively tests every Y.
    # Other methods combine systematic extreme/coarse pairs and uniform samples.
    if kind=='full':pi={'iy':np.arange(1<<b[0],dtype=np.int64)}
    else:
        inner='iz' if kind=='direct' else 'iu';ib=b[1] if kind=='direct' else b[2]
        yy=np.unique(np.r_[0,(1<<b[0])-1,((np.arange(64)+.5)*(1<<b[0])/64).astype(int)])
        zz=np.unique(np.r_[0,(1<<ib)-1,((np.arange(64)+.5)*(1<<ib)/64).astype(int)])
        ymesh,zmesh=np.meshgrid(yy,zz,indexing='ij')
        pi={'iy':np.r_[ymesh.ravel(),rng.integers(0,1<<b[0],samples)],
            inner:np.r_[zmesh.ravel(),rng.integers(0,1<<ib,samples)]}
    pres,pd=evaluate(graph,pi,b)
    assert pd['maximum_absolute_monetary_error']<=cfg['arithmetic_price_error_upper']
    count=min(1024,len(pres));ind=np.linspace(0,len(pres)-1,count,dtype=int)
    for j in ind:
        assert int(pres[j])==graph.evaluate({k:int(v[j]) for k,v in pi.items()},cache)
    return {'method':kind,'cap':cap,'price_tolerance':cfg['price_tolerance'],
            'retained_fraction_bits':q,'table_address_bits':p,'production_input_bits':b,
            'parent_uniform_arithmetic_price_bound':cfg['arithmetic_price_error_upper'],
            'coarse_finite_grid':cd,'retained_grid_diagnostic':pd,
            'retained_grid_exhaustive':kind=='full','scalar_graph_crosschecks':count+min(1024,len(cres)),
            'all_checks_passed':True}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--samples',type=int,default=1000000)
    ap.add_argument('--config',type=Path,default=ROOT/'results'/'encoding_resources.json')
    ap.add_argument('--output',type=Path,default=ROOT/'results'/'independent_encoding_audit.json')
    ap.add_argument('--certify-quantiles',action='store_true')
    args=ap.parse_args()
    start=time.perf_counter();config=args.config;config_bytes=config.read_bytes();data=json.loads(config_bytes)
    cache={};rng=np.random.default_rng(20261005);rows=[]
    for cfg in data['configurations']:
        rows.append(audit_case(cfg,cache,args.samples,rng))
        print('Passed integer replay',cfg['method'],cfg['price_tolerance'],flush=True)
    tables=[];qcerts=[]
    for key,packed in cache.items():
        if key[0]!='quantile':
            tables.append(certify_packed(key,packed))
            print('Passed independent endpoint enclosure',key,flush=True)
        elif args.certify_quantiles:
            qcerts.append(certify_quantile_table(key,packed))
            print('Passed independent quantile enclosure',key,flush=True)
    ref=json.loads((ROOT/'results'/'continuous_reference.json').read_text())
    assert ENCODING_SOURCE_SHA256==hashlib.sha256((ROOT/'encoding_resources.py').read_bytes()).hexdigest(),'retained encoding changed during audit'
    assert config_bytes==config.read_bytes(),'configuration changed during audit'
    out={'status':'independent retained integer replay and endpoint enclosure passed',
         'config_sha256':hashlib.sha256(config_bytes).hexdigest(),
         'root_encoding_sha256':ENCODING_SOURCE_SHA256,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seed':20261005,'random_samples_per_two_input_configuration':args.samples,
         'integer_replay':rows,'function_tables':tables,'quantile_tables':qcerts,
         'continuous_reference':ref['continuous_reference'],
         'finite_reference':ref['finite'],
         'limitations':['empirical replay does not prove a uniform arithmetic error bound',
             'retained parent analytical certificates still require a separate proof review',
             *([] if args.certify_quantiles else ['large inverse-normal QROM constants are not exhaustively independently enclosed here']),
             'continuous quadrature and variance refinements are numerical diagnostics, not directed enclosures',
             'no mean estimator or quantum speedup is tested'],
         'elapsed_seconds':time.perf_counter()-start}
    path=args.output;path.write_text(json.dumps(out,indent=2)+'\n');print(path)

if __name__=='__main__':main()
