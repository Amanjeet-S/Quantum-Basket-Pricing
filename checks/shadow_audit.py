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
from endpoint_cert import decimal_cdf, rounded_enclosure, capped_kernel, EPS
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
            if node.meta.startswith('capped_') or node.meta=='scaled_cdf':span_bits=1
            elif node.meta.startswith('strike_log_'):span_bits=0
            else:span_bits=2 if node.meta in ('exp','log') else 4
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
    if kind in ('interval','interval_kernel'):
        strike=100+cap*(iu+.5)/(1<<bu)
        return cap*(ncdf(np.log(ba/strike)/.2)-ncdf(np.log(bg/strike)/.2))
    # The positive layer-cake integral is independent of the retained full
    # formula's ordered CDF intervals and signed integer subtraction.
    return full_positive(ba,bg,.2,100,cap,192)

def certify_packed(key,packed):
    kind,bits,q=key; n=1<<bits; width=q+6; mask=(1<<width)-1
    endpoints=[x&mask for x in packed]+[(packed[-1]&mask)+(packed[-1]>>width)]
    assert len(endpoints)==n+1
    maximum_error=EPS
    with localcontext() as ctx:
        ctx.prec=100
        for j,stored in enumerate(endpoints):
            if kind=='exp':x=D(-2)+D(4)*j/n;value,error=x.exp(),EPS
            elif kind=='log':x=D('.25')+D(4)*j/n;value,error=x.ln()+2,EPS
            elif kind=='cdf':x=D(-8)+D(16)*j/n;value,error=decimal_cdf(x)
            elif kind.startswith('capped_'):
                x=D(-1)+D(2)*j/n;value,error=capped_kernel(x,D(kind.split('_',1)[1]))
            elif kind=='scaled_cdf':x=D(-1)+D(2)*j/n;value,error=decimal_cdf(5*x)
            elif kind.startswith('strike_log_'):
                x=D(j)/n;h=D(kind.split('_',2)[2])/100
                value,error=1+(1+h*x).ln(),EPS
            else:raise ValueError(kind)
            maximum_error=max(maximum_error,error)
            assert stored==rounded_enclosure(value,error,q),(kind,bits,q,j)
    assert all(x<=y for x,y in zip(endpoints,endpoints[1:]))
    return {'function':kind,'address_bits':bits,'fraction_bits':q,
            'enclosed_endpoint_count':n+1,'monotone_integer_endpoints':True,
            'independent_decimal_absolute_enclosure':str(maximum_error)}

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


def special_table_audit(key,packed):
    """Adversarial primitive input cells, independent of the Gaussian law."""
    kind,bits,q=key;n=1<<bits;unit=1<<q
    if kind.startswith('capped_') or kind=='scaled_cdf':span=2
    elif kind.startswith('strike_log_'):span=1
    else:return None
    boundaries=np.arange(n+1,dtype=np.int64)*(span*unit//n)
    address=np.unique(np.clip(np.r_[boundaries-1,boundaries,boundaries+1],0,span*unit-1))
    graph=er.Graph(q);name=graph.inp('x',q+2)
    graph.output=graph.function(kind,name,bits)
    values=vector_integer(graph,{'x':address},{key:packed})
    exact=[]
    with localcontext() as ctx:
        ctx.prec=100
        for a in address:
            if kind.startswith('capped_'):
                value,_=capped_kernel(D(int(a))/unit-1,D(kind.split('_',1)[1]))
            elif kind=='scaled_cdf':value,_=decimal_cdf(5*(D(int(a))/unit-1))
            else:
                h=D(kind.split('_',2)[2])/100
                value=1+(1+h*D(int(a))/unit).ln()
            exact.append(float(value))
    observed=float(np.max(abs(values.astype(float)/unit-np.asarray(exact))))
    curvature=(25/math.sqrt(2*math.pi*math.e) if not kind.startswith('strike_log_')
               else (float(kind.split('_',2)[2])/100)**2)
    primitive_bound=curvature*(span/n)**2/8+1.5/unit+1e-43
    assert observed<=primitive_bound+1e-14,(key,observed,primitive_bound)
    assert np.all(values[:-1]<=values[1:])
    for a,v in zip(address,values):
        assert int(v)==graph.evaluate({'x':int(a)},{key:packed})
    record={'function':kind,'address_bits':bits,'fraction_bits':q,
            'cell_boundary_and_adjacent_unit_states':len(address),
            'maximum_absolute_output_error_diagnostic':observed,
            'global_interpolation_rounding_bound':primitive_bound,
            'monotone_integer_interpolator':True,'scalar_graph_agreement':True}
    if kind.startswith('capped_'):
        cap=float(kind.split('_',1)[1]);h=cap/100
        m=address/unit-1
        quadratures=[]
        for order in (64,192):
            nodes,weights=np.polynomial.legendre.leggauss(order)
            strikes=1+h*(nodes+1)/2
            quadratures.append(ncdf((m[:,None]-np.log(strikes)[None,:])/.2)@weights/2)
        refinement=float(np.max(abs(quadratures[0]-quadratures[1])))
        formula_difference=float(np.max(abs(np.asarray(exact)-quadratures[1])))
        assert refinement<1e-10 and formula_difference<1e-10
        record['positive_quadrature']={'orders':[64,192],'states':len(address),
            'maximum_refinement_difference':refinement,
            'maximum_100_digit_formula_vs_positive_integral_difference':formula_difference,
            'qualification':'quadrature is an independent numerical diagnostic, not a directed error enclosure'}
    return record


def audit_new_uncomputation():
    """Complete financial tapes at deliberately coarse precisions."""
    from itertools import product
    records=[];tables={}
    for kind in ('full_kernel','interval_kernel'):
        for cap in (1,10,25):
            graph=er.build_graph(kind,8,2,2,2,3,3,4,cap)
            er.make_tables(graph,tables)
            circuit,inputs,comparison,flag=er.compile_small(graph,tables)
            actual=circuit.resources();composed=er.counted_composition(graph,tables,{})
            for label,key in [('X','x'),('CNOT','cnot'),('Toffoli','toffoli')]:
                assert composed['expanded_gate_counts'][label]==actual[key]
            assert composed['logical_qubits_serial_upper_bound']>=actual['wires']
            tested=0
            for indices in product(*[range(1<<b) for b in graph.inputs.values()]):
                mapping=dict(zip(graph.inputs,indices));threshold=graph.evaluate(mapping,tables)
                levels=sorted(set([0,min(threshold,255),max(0,threshold-1),255]))
                for level in levels:
                    state=circuit.run({**{inputs[k]:v for k,v in mapping.items()},comparison:level})
                    assert state[flag]==int(level<threshold),(kind,cap,mapping,level)
                    assert all(circuit.value(state,inputs[k])==v for k,v in mapping.items())
                    assert circuit.value(state,comparison)==level
                    reserved={flag,*comparison,*[w for reg in inputs.values() for w in reg]}
                    assert all(v==0 for i,v in enumerate(state) if i not in reserved)
                    tested+=1
            records.append({'method':kind,'cap':cap,'fraction_bits':8,
                'complete_basis_state_checks':tested,'inputs_preserved':True,
                'all_work_and_comparison_high_bit_reset':True,'exact_flag_predicate':True,
                'complete_vs_composed_X_CNOT_Toffoli_agreement':True,
                'qualification':'coarse gate correctness, not financial price accuracy'})
    return records


def audit_case(cfg,cache,samples,rng):
    kind=cfg['method'];q=cfg['q'];cap=cfg['cap'];b=(cfg['b_y'],cfg['b_z'],cfg['b_u'])
    p=(cfg['exp_address_bits'],cfg['log_address_bits'],cfg['cdf_address_bits'])
    graph=er.build_graph(kind,q,*b,*p,cap);er.make_tables(graph,cache)
    coarse=er.build_graph(kind,q,6,6,6,*p,cap);er.make_tables(coarse,cache)
    if kind in ('full','full_kernel'):
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
    if kind in ('full','full_kernel'):pi={'iy':np.arange(1<<b[0],dtype=np.int64)}
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
            'retained_grid_exhaustive':kind in ('full','full_kernel'),'scalar_graph_crosschecks':count+min(1024,len(cres)),
            'all_checks_passed':True}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--samples',type=int,default=1000000)
    ap.add_argument('--config',type=Path,default=ROOT/'results'/'conditioning_resources.json')
    ap.add_argument('--output',type=Path,default=ROOT/'results'/'conditioning_independent_audit.json')
    ap.add_argument('--certify-quantiles',action='store_true')
    ap.add_argument('--new-variant-gates',action='store_true')
    args=ap.parse_args()
    start=time.perf_counter();config=args.config;config_bytes=config.read_bytes();data=json.loads(config_bytes)
    cache={};rng=np.random.default_rng(20261005);rows=[]
    for cfg in data['configurations']:
        rows.append(audit_case(cfg,cache,args.samples,rng))
        print('Passed integer replay',cfg['method'],cfg['price_tolerance'],flush=True)
    tables=[];qcerts=[];special_checks=[]
    for key,packed in cache.items():
        if key[0]!='quantile':
            tables.append(certify_packed(key,packed))
            print('Passed independent endpoint enclosure',key,flush=True)
            if key[0].startswith(('capped_','strike_log_')) or key[0]=='scaled_cdf':
                special_checks.append(special_table_audit(key,packed))
        elif args.certify_quantiles:
            qcerts.append(certify_quantile_table(key,packed))
            print('Passed independent quantile enclosure',key,flush=True)
    gates=audit_new_uncomputation() if args.new_variant_gates else []
    ref=json.loads((ROOT/'results'/'continuous_reference.json').read_text())
    assert ENCODING_SOURCE_SHA256==hashlib.sha256((ROOT/'encoding_resources.py').read_bytes()).hexdigest(),'retained encoding changed during audit'
    assert config_bytes==config.read_bytes(),'configuration changed during audit'
    out={'status':'independent retained integer replay and endpoint enclosure passed',
         'config_sha256':hashlib.sha256(config_bytes).hexdigest(),
         'root_encoding_sha256':ENCODING_SOURCE_SHA256,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seed':20261005,'random_samples_per_two_input_configuration':args.samples,
         'integer_replay':rows,'function_tables':tables,'quantile_tables':qcerts,
         'special_primitive_boundary_and_quadrature_checks':special_checks,
         'new_variant_complete_gate_checks':gates,
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
