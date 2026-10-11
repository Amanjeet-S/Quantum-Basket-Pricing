"""Independent finite integer and financial audit of the broader constructor.

The retained graph/tables are inputs, not the main evaluator. A separate
uint64 interpreter uses split products. Endpoint enclosures use an independent
100-digit normal series and charged Mills tails. Gaussian production samples
and selected complete gate tapes are diagnostics; uniform error claims rely
on the separately checked premises and analytical certificate.
"""
from decimal import Decimal, localcontext
from itertools import product
from pathlib import Path
from statistics import NormalDist
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import generalisation as gen
import encoding_resources as er
from classical_reference_generalisation import (normal_cdf,normal_interval,log_baskets,
    conditional_positive,conditional_closed,exact_model_checks,geometric_positive)

D=Decimal;EPS=D('1e-50');ND=NormalDist()


def arctan(x):
    total=D(0);term=x;n=0
    while True:
        total+=term/D(2*n+1);term*=-x*x;n+=1
        if abs(term/D(2*n+1))<D('1e-96'):return total


with localcontext() as ctx:
    ctx.prec=100
    PI=16*arctan(D(1)/5)-4*arctan(D(1)/239)
    PHI0=1/(2*PI).sqrt()


def cdf_enclosure(x):
    """Explicit series remainder on [-8,8], explicit Mills tail outside."""
    x=D(x)
    with localcontext() as ctx:
        ctx.prec=100
        if abs(x)>=8:
            error=PHI0*(-x*x/2).exp()/abs(x)+EPS
            return (D(0) if x<0 else D(1)),error
        if not x:return D('.5'),D(0)
        term=x;total=term;n=0;peak=abs(term)
        while True:
            nxt=-term*x*x*D(2*n+1)/(2*D(n+1)*D(2*n+3));n+=1
            peak=max(peak,abs(nxt),abs(total))
            if n>32 and abs(nxt)<D('1e-80'):
                assert n<250 and peak<D('1e30')
                return D('.5')+PHI0*total,EPS+abs(PHI0*nxt)
            total+=nxt;term=nxt


def endpoint_enclosure(kind,x):
    x=D(x)
    with localcontext() as ctx:
        ctx.prec=100
        if kind=='general_exp':return (x-4).exp(),EPS
        if kind=='general_log':return (x+D('.0625')).ln()+4,EPS
        parts=kind.split(':')
        if parts[0]=='general_strike':return (1+D(parts[1])/100*x).ln(),EPS
        lam=D(parts[1]);m=x-2
        if parts[0]=='general_cdf':return cdf_enclosure(m/lam)
        h=D(parts[2])/100;a=(1+h).ln();t0=m/lam;t1=(m-a)/lam
        p0,e0=cdf_enclosure(t0);p1,e1=cdf_enclosure(t1)
        a0,ea0=cdf_enclosure(t0+lam);a1,ea1=cdf_enclosure(t1+lam)
        growth=(m+lam*lam/2).exp()
        value=(growth*(a0-a1)-(p0-p1)+h*p1)/h
        error=(growth*(ea0+ea1)+e0+e1+h*e1)/h+D('1e-44')
        # No projection or cumulative maximum in the independent reference.
        # The interval encloses the true continuous normal expectation.
        return value,error


def span(kind):
    if kind in ('general_exp','general_log'):return 8
    if kind.startswith('general_strike:'):return 1
    if kind.startswith(('general_cdf:','general_cap:')):return 4
    raise ValueError(kind)


def vector_integer(graph,inputs,tables,retain=False):
    q=graph.q;width=graph.width;mask=np.uint64((1<<width)-1);umask=np.uint64((1<<q)-1)
    vals={k:np.asarray(v,dtype=np.uint64) for k,v in inputs.items()}
    def mul(a,b):
        ah,al=a>>np.uint64(q),a&umask;bh,bl=b>>np.uint64(q),b&umask
        return (ah*bh<<np.uint64(q))+ah*bl+al*bh+(al*bl>>np.uint64(q))
    for name,node in graph.nodes.items():
        a=[vals[k] for k in node.args]
        if node.kind=='constant':v=np.uint64(node.meta)
        elif node.kind=='add':v=a[0]+a[1]
        elif node.kind in ('sub','positive_sub'):v=a[0]-a[1]
        elif node.kind=='shift':v=a[0]>>np.uint64(node.meta)
        elif node.kind=='mul':v=mul(*a)
        elif node.kind=='scale':
            integer,drop,_=node.meta
            if drop==q:v=mul(a[0],np.uint64(integer))
            elif drop==0:v=a[0]*np.uint64(integer)
            else:raise AssertionError(node.meta)
        elif node.kind=='quantile':v=np.asarray(tables[('quantile',node.bits,q)],np.uint64)[a[0].astype(np.int64)]
        elif node.kind=='uniform':v=(2*a[0]+1)*np.uint64((1<<q)>>(node.bits+1))
        elif node.kind=='function':
            f=q+int(math.log2(span(node.meta)))-node.bits
            index=(a[0]>>np.uint64(f)).astype(np.int64);rem=a[0]&np.uint64((1<<f)-1)
            assert np.all(index<(1<<node.bits)),('independent function address',node.meta,name,index.max())
            table=tables[(node.meta,node.bits,q)]
            base=np.asarray([p&int(mask) for p in table],np.uint64)[index]
            delta=np.asarray([p>>width for p in table],np.uint64)[index]
            # The actually retained slopes and remainders must fit the
            # independent interpreter's uint64 product, checked explicitly.
            assert int(delta.max())*((1<<f)-1)<1<<64
            v=base+(delta*rem>>np.uint64(f))
        elif node.kind=='clamp':
            signed=np.where(a[0]<(1<<(width-1)),a[0].astype(np.int64),a[0].astype(np.int64)-(1<<width))
            v=np.clip(signed,0,node.meta).astype(np.uint64)
        else:raise ValueError(node.kind)
        vals[name]=v&mask
    return vals if retain else vals[graph.output]


def quantiles(bits):
    return np.fromiter((ND.inv_cdf((j+.5)/(1<<bits)) for j in range(1<<bits)),float,count=1<<bits)


def ideal(model,method,inputs,by,bz,bu,cap):
    yy=quantiles(by);y=np.column_stack([yy[inputs[f'iy{j}']] for j in range(model['outer_dimension'])])
    ma,mg=log_baskets(model,y);lam=float(model['lambda'])
    if method=='direct':
        z=quantiles(bz)[inputs['iz']]
        return np.maximum(np.clip(100*np.expm1(ma+lam*z),0,cap)-np.clip(100*np.expm1(mg+lam*z),0,cap),0)
    if method=='interval':
        k=np.log1p(cap/100*(inputs['iu']+.5)/(1<<bu))
        return cap*normal_interval((mg-k)/lam,(ma-k)/lam)
    return conditional_positive(ma,mg,lam,cap,96)


def certify_table(key,packed,width):
    kind,bits,q=key;n=1<<bits;mask=(1<<width)-1;unit=D(1)/D(1<<q)
    ends=[p&mask for p in packed]+[(packed[-1]&mask)+(packed[-1]>>width)]
    assert len(ends)==n+1 and all(x<=y for x,y in zip(ends,ends[1:]))
    allowance=D('2e-12') if kind.startswith('general_cap:') else D('6.4e-16') if kind.startswith('general_cdf:') else D('1e-38')
    max_error=D(0);max_ref_error=D(0);not_unique=0
    with localcontext() as ctx:
        ctx.prec=100
        for j,stored in enumerate(ends):
            val,error=endpoint_enclosure(kind,D(span(kind))*j/n)
            observed=abs(D(stored)*unit-val)
            # The entire independently enclosed interval must lie within
            # the declared ordinate allowance, not merely intersect it.
            assert observed+error<=unit/2+allowance,(key,j,str(observed),str(error))
            max_error=max(max_error,observed+error);max_ref_error=max(max_ref_error,error)
            if error>EPS:not_unique+=1
    return {'kernel':kind,'address_bits':bits,'q':q,'endpoint_count':n+1,
            'absolute_ordinate_allowance':str(allowance),'maximum_stored_vs_independent_enclosure':str(max_error),
            'largest_independent_reference_enclosure':str(max_ref_error),
            'monotone_integer_endpoints':True,'qualification':'bounded rounding; no claim of unique exact rounding after charged tails or projection'}


def certify_quantile(key,stored):
    _,bits,q=key;n=1<<bits;unit=D(1)/D(1<<q)
    assert all(a+b==10*(1<<q) for a,b in zip(stored,stored[::-1]))
    with localcontext() as ctx:
        ctx.prec=100
        for j in range(n//2,n):
            centre=D(stored[j])*unit-5;target=(D(j)+D('.5'))/n
            lo,elo=cdf_enclosure(centre-unit);hi,ehi=cdf_enclosure(centre+unit)
            assert lo+elo<=target<=hi-ehi,(key,j)
    return {'bits':bits,'q':q,'entries':n,'exhaustive_positive_half_with_exact_negative_symmetry':True,
            'inverse_coordinate_error_upper':float(unit)}


def critical_states(kind,key,table,width):
    """All cell boundaries +/- one unit, with selected positive integrals."""
    _,bits,q=key;n=1<<bits;unit=1<<q;sp=span(kind)
    edges=np.arange(n+1,dtype=np.int64)*(sp*unit//n)
    x=np.unique(np.clip(np.r_[edges-1,edges,edges+1],0,sp*unit-1))
    g=er.Graph(q);g.width=width;inp=g.inp('x',q+3);g.output=g.function(kind,inp,bits)
    v=vector_integer(g,{'x':x},{key:table})
    assert np.all(v[:-1]<=v[1:]);record={'kernel':kind,'bits':bits,'q':q,'cell_boundary_and_adjacent_states':len(x),
        'integer_interpolator_monotone':True}
    # The financial cap kernels include worst tail transitions at ±8lambda,
    # near strike thresholds and zero, besides evenly spaced signed inputs.
    if kind.startswith('general_cap:'):
        _,sl,sh=kind.split(':');lam=float(sl);cap=int(sh);h=cap/100
        m=np.unique(np.clip(np.r_[np.linspace(-2,2,257),[-8*lam,8*lam,0,math.log1p(h),
            -8*lam-lam*lam,8*lam-lam*lam,math.log1p(h)-8*lam,math.log1p(h)+8*lam]],-2,2))
        nodes64,w64=np.polynomial.legendre.leggauss(64);nodes192,w192=np.polynomial.legendre.leggauss(192)
        a64=normal_cdf((m[:,None]-np.log1p(h*(nodes64+1)/2))/lam)@(w64/2)
        a192=normal_cdf((m[:,None]-np.log1p(h*(nodes192+1)/2))/lam)@(w192/2)
        refs=np.array([float(endpoint_enclosure(kind,D(str(z))+2)[0]) for z in m])
        assert np.max(abs(a64-a192))<1e-10
        assert np.max(abs(refs-a192))<3e-12
        record['independent_positive_integral']={'states':len(m),'orders':[64,192],
            'maximum_refinement':float(np.max(abs(a64-a192))),
            'maximum_decimal_formula_vs_integral':float(np.max(abs(refs-a192))),
            'qualification':'numerical diagnostic, not directed quadrature enclosure'}
    return record


def systematic_inputs(graph,samples,rng):
    names=list(graph.inputs);bits=list(graph.inputs.values());d=len(names)
    patterns=list(product((0,1),repeat=d))
    # Every outer Gaussian corner, extreme inner values, centre and coarse
    # axial slices. No tensor exhaustive claim for dimensions three/seven.
    rows=[[v*((1<<b)-1) for v,b in zip(pat,bits)] for pat in patterns]
    centre=[(1<<b)//2 for b in bits];rows.append(centre)
    for j,b in enumerate(bits):
        for index in np.unique(((np.arange(64)+.5)*(1<<b)/64).astype(int)):
            row=centre.copy();row[j]=int(index);rows.append(row)
    system=np.asarray(rows,np.int64)
    random=np.column_stack([rng.integers(0,1<<b,samples) for b in bits])
    both=np.vstack([system,random])
    return {name:both[:,j] for j,name in enumerate(names)},len(system)


def audit_case(model,row,cache,rng,samples):
    args=[row[k] for k in ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')]
    g=gen.build_graph(model,row['method'],*args,row['cap']);gen.make_tables(g,cache)
    inputs,systematic=systematic_inputs(g,samples,rng)
    if model['outer_dimension']==1 and row['method'].startswith('full'):
        # Retain exhaustive Y for the two-asset cases in addition to the
        # common prescribed sample design.
        inputs['iy0']=np.unique(np.r_[inputs['iy0'],np.arange(1<<row['b_y'])])
    result=vector_integer(g,inputs,cache)
    ref=ideal(model,row['method'],inputs,row['b_y'],row['b_z'],row['b_u'],row['cap'])
    monetary=result.astype(float)*row['cap']/(1<<row['q']);diff=abs(monetary-ref)
    assert diff.max()<=row['arithmetic_price_error_upper']+1e-12,(row['case_id'],row['method'],diff.max())
    ids=np.linspace(0,len(result)-1,min(256,len(result)),dtype=int)
    for j in ids:assert int(result[j])==g.evaluate({k:int(v[j]) for k,v in inputs.items()},cache)
    # Independent finite-grid price expectation for n2 full conditioning.
    finite=None
    if model['outer_dimension']==1 and row['method'].startswith('full'):
        idx={'iy0':np.arange(1<<row['b_y'])};rr=ideal(model,row['method'],idx,row['b_y'],0,0,row['cap'])
        finite={'outer_grid_states':len(rr),'positive_integral_residual_mean':float(rr.mean()),
                'with_continuous_geometric_correction':float(rr.mean()+geometric_positive(model,row['cap']))}
    return {'case_id':row['case_id'],'method':row['method'],'states':len(ref),'systematic_states':systematic,
            'maximum_absolute_monetary_error':float(diff.max()),
            'mean_signed_monetary_error':float((monetary-ref).mean()),
            'analytical_arithmetic_price_error_upper':row['arithmetic_price_error_upper'],
            'scalar_graph_comparisons':len(ids),'independent_integer_agreement':True,
            'n2_finite_expectation':finite,'qualification':'sampling diagnoses implementation; does not prove uniform financial error'}


def gate_replays(models):
    records=[]
    selected=[next(m for m in models if m['assets']==n and m['lambda']==lam)
              for n,lam in [(2,'0.05'),(4,'0.2'),(8,'0.4')]]
    for model in selected:
        for method in ('direct','interval','full_cap'):
            cap=25 if model['assets']==2 else 10 if model['assets']==4 else 1
            g=gen.build_graph(model,method,9,1,1,1,3,4,4,cap);tables={};gen.make_tables(g,tables)
            circuit,inputs,comparison,flag=er.compile_small(g,tables)
            actual=circuit.resources();composed=er.counted_composition(g,tables,{})
            for label,k in [('X','x'),('CNOT','cnot'),('Toffoli','toffoli')]:assert composed['expanded_gate_counts'][label]==actual[k]
            assert composed['logical_qubits_serial_upper_bound']>=actual['wires']
            names=list(g.inputs);mapping_set=list(product((0,1),repeat=len(names)))
            if len(mapping_set)>32:
                rr=np.random.default_rng(6347);keep=np.unique(np.r_[0,len(mapping_set)-1,rr.integers(0,len(mapping_set),30)])
                mapping_set=[mapping_set[j] for j in keep]
            count=0
            for mapping_values in mapping_set:
                mapping=dict(zip(names,mapping_values));threshold=g.evaluate(mapping,tables)
                for level in sorted(set([0,511,min(threshold,511),max(0,min(threshold-1,511))])):
                    state=circuit.run({**{inputs[k]:v for k,v in mapping.items()},comparison:level})
                    assert state[flag]==int(level<threshold)
                    assert all(circuit.value(state,inputs[k])==v for k,v in mapping.items())
                    assert circuit.value(state,comparison)==level
                    reserved={flag,*comparison,*[x for reg in inputs.values() for x in reg]}
                    assert all(v==0 for j,v in enumerate(state) if j not in reserved)
                    count+=1
            records.append({'model_id':model['model_id'],'method':method,'cap':cap,'fraction_bits':9,
                'basis_comparator_checks':count,'input_registers_preserved':True,'all_work_reset':True,
                'flag_predicate_correct':True,'gate_counts_equal_composed_blocks':True,
                'qualification':'selected coarse complete tapes, not production price accuracy'})
    return records


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,default=ROOT/'results/generalisation_resources.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/generalisation_independent_audit.json')
    p.add_argument('--samples',type=int,default=16384);p.add_argument('--certify-quantiles',action='store_true')
    p.add_argument('--skip-gates',action='store_true');args=p.parse_args();start=time.perf_counter()
    config=json.loads(args.config.read_text())
    sourcefiles=[ROOT/f for f in config['source_hashes']]
    before={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sourcefiles}
    own_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    models=config['models'];metrics=[exact_model_checks(m) for m in models]
    assert config['protocol_sha256']==hashlib.sha256((ROOT/'results/generalisation_protocol.json').read_bytes()).hexdigest()
    cache={};rng=np.random.default_rng(20261005);cases=[]
    for row in config['configurations']:
        model=next(m for m in models if m['model_id']==row['model_id'])
        cases.append(audit_case(model,row,cache,rng,args.samples));print('Audited',row['case_id'],row['method'],flush=True)
    endpoints=[];quantile_records=[];critical=[]
    for key,table in cache.items():
        if key[0]=='quantile':
            if args.certify_quantiles:quantile_records.append(certify_quantile(key,table))
        else:
            endpoints.append(certify_table(key,table,key[2]+7));critical.append(critical_states(key[0],key,table,key[2]+7))
    gates=[] if args.skip_gates else gate_replays(models)
    after={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sourcefiles};assert before==after,'Source changed during independent audit'
    assert own_sha==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'Audit source changed during execution'
    for f,expected in config.get('source_hashes',{}).items():
        path=ROOT/f
        if path.exists():assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,('configuration source hash mismatch',f)
    out={'status':'independent retained-graph diagnostics and constant enclosures passed',
         'configuration_sha256':hashlib.sha256(args.config.read_bytes()).hexdigest(),
         'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'classical_source_sha256':hashlib.sha256((Path(__file__).parent/'classical_reference_generalisation.py').read_bytes()).hexdigest(),
         'source_hashes_unchanged_during_audit':before,'model_checks':metrics,'configurations':cases,
         'endpoint_enclosures':endpoints,'quantile_enclosures':quantile_records,'primitive_boundary_diagnostics':critical,
         'selected_complete_gate_replays':gates,'summary':{'financial_configurations':len(cases),
             'finite_production_states':sum(c['states'] for c in cases),
             'scalar_graph_checks':sum(c['scalar_graph_comparisons'] for c in cases),
             'function_endpoints':sum(c['endpoint_count'] for c in endpoints),
             'quantile_entries':sum(c['entries'] for c in quantile_records),
             'primitive_cell_boundary_states':sum(c['cell_boundary_and_adjacent_states'] for c in critical),
             'complete_gate_basis_checks':sum(c['basis_comparator_checks'] for c in gates)},
         'qualifications':['No exhaustive production-gate replay for the high-dimensional cases',
             'No end-to-end quantum estimation or hardware confidence claim',
             'Ordinate rounding is bounded after Mills clipping, projection and cumulative maximum',
             'Financial interpolation certificates and gate correctness are distinct statements'],
         'elapsed_seconds':time.perf_counter()-start}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(args.output)


if __name__=='__main__':main()
