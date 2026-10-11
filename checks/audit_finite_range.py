"""Independent finite-support, convex-vertex and lattice-moment audit.

The range generator is not imported. Reference calculations use a different
log-sum-exp evaluation and 450 decimal digits; production values come from
the existing exact NumPy integer evaluator and are checked against Python
integer graph replay. Those point checks supplement the uniform proof.
"""
from decimal import Decimal,localcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import math
import sys
import time
sys.dont_write_bytecode=True
D=Decimal
REFERENCE_ERROR=D('1e-350')
SOURCE_CENTRE_ERROR=D('1e-200')
KEY_FIELDS=('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits')


def key(row):return tuple(row[k] for k in KEY_FIELDS)
def dec(x):
    x=F(x);return D(x.numerator)/D(x.denominator)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def phi_reference():
    # pi=6 asin(1/2), independently evaluated as a positive series whose
    # successive term ratio is below1/4. Its tail is below8*next_term.
    with localcontext() as ctx:
        ctx.prec=450
        total,term,k=D(0),D('.5'),0
        while True:
            total+=term
            nxt=term*D((2*k+1)**2)/D(8*(k+1)*(2*k+3))
            if 8*nxt<D('1e-420'):break
            term=nxt;k+=1
        return 1/(12*total).sqrt()


def geometry_reference(model,maximum):
    w=list(map(F,model['weights']));r=[list(map(F,row)) for row in model['loadings']]
    mu=list(map(F,model['log_drifts']));p=model['outer_dimension']
    assert len(w)==len(r)==len(mu)==model['assets'] and sum(w)==1 and all(x>0 for x in w)
    assert all(len(row)==p for row in r) and p<=7 and len(w)<=8
    assert F(model['strike'])==F(model['spot'])==100
    rg=[sum(ww*rr[j] for ww,rr in zip(w,r)) for j in range(p)]
    mg=sum(ww*mm for ww,mm in zip(w,mu));vertices=[]
    with localcontext() as ctx:
        ctx.prec=450
        for signs in product((-1,1),repeat=p):
            arguments=[mm-mg+maximum*sum((rr[j]-rg[j])*signs[j] for j in range(p)) for mm,rr in zip(mu,r)]
            assert sum(ww*z for ww,z in zip(w,arguments))==0
            shift=max(arguments)
            # Scaled log-sum-exp differs from the generator's raw exp sum.
            total=sum(dec(ww)*dec(z-shift).exp() for ww,z in zip(w,arguments))
            J=dec(shift)+total.ln();assert J>=-REFERENCE_ERROR
            vertices.append((signs,J))
        maximum_J=max(J for _,J in vertices)
        maximum_argument=max(abs(mm-mg)+maximum*sum(abs(rr[j]-rg[j]) for j in range(p)) for mm,rr in zip(mu,r))
        assert maximum_argument<=20
        mg_upper=mg+maximum*sum(map(abs,rg))
    return maximum_J,vertices,maximum_argument,mg_upper,(w,r,mu,rg,mg)


def gaussian_mean_reference(model,row,geometry):
    w,r,mu,rg,mg=geometry;lam=F(model['lambda']);cap=F(row['cap']);h=cap/F(model['strike'])
    with localcontext() as ctx:
        ctx.prec=450
        EA=D(model['spot'])*sum(dec(ww)*(dec(mm)+dec(sum(z*z for z in rr)+lam*lam)/2).exp() for ww,mm,rr in zip(w,mu,r))
        EG=D(model['spot'])*(dec(mg)+dec(sum(z*z for z in rg)+lam*lam)/2).exp()
        J0=sum(dec(ww)*(dec(mm-mg)+dec(sum((a-b)**2 for a,b in zip(rr,rg)))/2).exp() for ww,mm,rr in zip(w,mu,r)).ln()
        slope=min(phi_reference()/dec(lam),dec((1+h)/h))
        mean=min(D(1),(EA-EG)/dec(cap),slope*J0)
        numerical=(D.from_float(row['arithmetic_price_error_upper'])+D.from_float(row['distribution_price_bound']))/dec(cap)
        return min(D(1),mean+numerical)


def chord_inverse(q,v,maximum):
    """Independent rational bisection of the finite squared-value hull."""
    v=F(v);delta=F(1,1<<q)
    if maximum==0:return F(0)
    if v>=maximum*maximum*delta*delta:return maximum*delta
    low,high=0,maximum
    while high-low>1:
        middle=(low+high)//2
        if middle*middle*delta*delta<=v:low=middle
        else:high=middle
    # For every integer a, (a-low)(a-low-1)>=0, so the following
    # linear moment certificate holds pointwise throughout the lattice.
    x,y=low*delta,(low+1)*delta
    return (v+x*y)/(x+y)


def cdf_reference(x):
    """Unclipped positive-series Gaussian CDF, independent of the tables.

    integral_0^t exp(-z²/2)dz = exp(-t²/2) sum t^(2k+1)/(2k+1)!!.
    The first omitted term divided by1-next_ratio encloses the positive tail.
    This diagnostic reference uses180digits and a conservative1e-70 envelope.
    """
    if x==0:return D('.5')
    t=abs(x);assert t<=20
    square=t*t;term,total,k=t,t,0
    while True:
        nxt=term*square/D(2*k+3)
        next_ratio=square/D(2*k+5)
        if next_ratio<1 and nxt/(1-next_ratio)<D('1e-120'):break
        total+=nxt;term=nxt;k+=1
    integral=phi_reference()*(-square/2).exp()*total
    return D('.5')+integral if x>0 else D('.5')-integral


def ideal_reference(model,src,address,tables):
    q=src['q'];w=list(map(F,model['weights']));loadings=[list(map(F,r)) for r in model['loadings']]
    drifts=list(map(F,model['log_drifts']));lam=F(model['lambda']);cap=F(src['cap'])
    y=[F(tables[('quantile',src['b_y'],q)][address[f'iy{j}']],1<<q)-5 for j in range(model['outer_dimension'])]
    logs=[mu+sum(a*b for a,b in zip(rr,y)) for mu,rr in zip(drifts,loadings)]
    with localcontext() as ctx:
        ctx.prec=180
        mg=sum(dec(ww*ll) for ww,ll in zip(w,logs));ma=sum(dec(ww)*dec(ll).exp() for ww,ll in zip(w,logs)).ln()
        if src['method']=='direct':
            z=F(tables[('quantile',src['b_z'],q)][address['iz']],1<<q)-5
            A=D(model['spot'])*(ma+dec(lam*z)).exp();G=D(model['spot'])*(mg+dec(lam*z)).exp()
            payoff=lambda stock:min(max(stock-D(model['strike']),D(0)),dec(cap))
            return (payoff(A)-payoff(G))/dec(cap)
        h=dec(cap/F(model['strike']));l=dec(lam)
        if src['method']=='interval':
            u=F(2*address['iu']+1,1<<(src['b_u']+1));strike=(1+h*dec(u)).ln()
            return cdf_reference((ma-strike)/l)-cdf_reference((mg-strike)/l)
        edge=(1+h).ln()
        def capped(m):
            a,b=m/l,(m-edge)/l
            return ((m+l*l/2).exp()*(cdf_reference(a+l)-cdf_reference(b+l))
                    -cdf_reference(a)+cdf_reference(b)+h*cdf_reference(b))/h
        return capped(ma)-capped(mg)


def analytic_audit(source_root,data,certificate,moments,er):
    source={key(r):r for r in data['configurations']};moment={key(r):r for r in moments['configurations']}
    models={m['model_id']:m for m in data['models']};supports={};geometries={};records=[]
    def actual_support(bits,q):
        if (bits,q) not in supports:
            raw,_=er.quantile_table(bits,q,False);values=[int(x)-5*(1<<q) for x in raw]
            assert all(a+b==0 for a,b in zip(values,reversed(values)))
            assert all(a<=b for a,b in zip(values,values[1:]))
            maximum=F(values[-1],1<<q);assert maximum<F(19,4)
            supports[(bits,q)]=(maximum,hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest())
        return supports[(bits,q)]
    assert len(source)==len(moment)==len(certificate['configurations'])==273
    for row in certificate['configurations']:
        src=source[key(row)];mom=moment[key(row)];model=models[row['model_id']]
        maximum,sha=actual_support(src['b_y'],src['q'])
        assert F(row['quantile_support_maximum_fraction'])==maximum and row['quantile_table_sha256']==sha
        geomkey=(row['model_id'],src['b_y'],src['q'])
        if geomkey not in geometries:geometries[geomkey]=geometry_reference(model,maximum)
        J,vertices,max_arg,mg,geometry=geometries[geomkey]
        assert F(row['maximum_exponential_argument_exact'])==max_arg
        assert F(row['geometric_log_base_maximum_exact'])==mg
        assert row['vertex_count']==len(vertices)
        assert D(row['maximum_log_gap_upper_decimal'])>=J+REFERENCE_ERROR
        assert abs(D(row['maximum_log_gap_centre_decimal'])-J)<=SOURCE_CENTRE_ERROR
        win=tuple(row['maximising_vertex_signs'])
        assert any(signs==win and abs(z-J)<=SOURCE_CENTRE_ERROR for signs,z in vertices)
        with localcontext() as ctx:
            ctx.prec=450
            cap=F(row['cap']);h=cap/F(model['strike']);lam=F(model['lambda'])
            assert F(1,20)<=lam<=F(2,5) and F(1,100)<=h<=F(1,4)
            assert 18<=src['q']<=30
            slope=phi_reference()/dec(lam)
            if row['method'] in ('full_global','full_cap'):slope=min(slope,dec((1+h)/h))
            if row['method']=='direct':slope=dec((1+h)/h)
            assert D(row['normalised_slope_upper_decimal'])>=slope+REFERENCE_ERROR
            candidates={'range':D(1),'log_gap':slope*J}
            if row['method']=='direct':
                zmax,zsha=actual_support(src['b_z'],src['q'])
                assert F(row['inner_quantile_support_maximum_fraction'])==zmax
                assert row['inner_quantile_table_sha256']==zsha
                assert abs(mg+lam*zmax)<=20
                candidates['positive_AG_gap']=D(model['spot'])*(dec(mg+lam*zmax)).exp()*(J.exp()-1)/dec(cap)
            elif row['method'] in ('full_global','full_cap'):
                assert abs(mg+lam*lam/2)<=20
                candidates['positive_AG_gap']=D(model['spot'])*(dec(mg+lam*lam/2)).exp()*(J.exp()-1)/dec(cap)
            assert set(candidates)==set(row['ideal_grid_range_candidates_decimal'])
            for name,value in candidates.items():
                retained=D(row['ideal_grid_range_candidates_decimal'][name])
                assert retained==1 if name=='range' else retained>=value+REFERENCE_ERROR
            eta=F.from_float(src['arithmetic_price_error_upper'])/cap
            assert F(row['body_pointwise_error_upper_decimal'])==eta
            assert F(mom['body_pointwise_normalised_error_upper_decimal'])==eta
            assert D(row['finite_output_range_upper_decimal'])>=min(D(1),min(candidates.values())+dec(eta))
            upper=F(row['finite_output_range_upper_decimal']);scaled=upper*(1<<src['q'])
            maximum_integer=scaled.numerator//scaled.denominator
            assert row['max_integer']==maximum_integer and 0<=maximum_integer<=1<<src['q']
            rho=F(maximum_integer,1<<src['q'])
            assert F(row['finite_lattice_range_fraction'])==rho
            legacy=D(row['legacy_mean_upper_decimal']);refmean=gaussian_mean_reference(model,src,geometry)
            assert legacy>=refmean if legacy==1 else legacy>=refmean+REFERENCE_ERROR
            mean=min(F(legacy),rho)
            second=min(F(mom['finite_output_second_moment_upper_decimal']),rho*mean,rho*rho)
            mean=min(mean,chord_inverse(src['q'],second,maximum_integer))
            second=min(second,rho*mean)
            assert F(row['fair_mean_upper_fraction'])==mean
            assert F(row['fair_second_moment_upper_fraction'])==second
        records.append({'case_id':row['case_id'],'method':row['method'],'q':src['q'],
            'exp_address_bits':row['exp_address_bits'],'log_address_bits':row['log_address_bits'],
            'probability_address_bits':row['probability_address_bits'],
            'exact_support_checked':True,'vertices_checked':len(vertices),
            'maximum_log_gap_reference_decimal':str(J),'max_integer':maximum_integer,
            'finite_lattice_range_fraction':str(rho),'exact_lattice_moment_transfer_checked':True})
    return records,supports,source,models


def production_replay(data,certificate,source,models,gen,ce,np):
    # Select the strongest retained range in each asset-count/method cell.
    # This stresses the certificates that materially tighten the competitors.
    chosen={}
    for row in certificate['configurations']:
        group=(models[row['model_id']]['assets'],row['method'])
        score=F(row['finite_lattice_range_fraction'])
        if group not in chosen or score<chosen[group][0]:chosen[group]=(score,row)
    cache={};rng=np.random.default_rng(710619);records=[]
    for _,row in sorted(chosen.values(),key=lambda z:(models[z[1]['model_id']]['assets'],z[1]['method'])):
        src=source[key(row)];model=models[row['model_id']]
        graph=gen.build_graph(model,src['method'],*[src[k] for k in ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],src['cap'])
        gen.make_tables(graph,cache)
        engine=ce.FiniteEvaluator(graph,cache);names=list(graph.inputs);bits=list(graph.inputs.values())
        addresses=[[p*((1<<b)-1) for p,b in zip(pattern,bits)] for pattern in product((0,1),repeat=len(names))]
        centre=[1<<(b-1) for b in bits];addresses.append(centre)
        for j,b in enumerate(bits):
            for value in sorted({0,1,(1<<b)-2,(1<<b)-1,(1<<(b-1))-1,1<<(b-1),(1<<(b-1))+1}):
                a=centre.copy();a[j]=value;addresses.append(a)
        addresses.extend(np.column_stack([rng.integers(0,1<<b,size=4096,dtype=np.uint64) for b in bits]).tolist())
        matrix=np.asarray(addresses,dtype=np.uint64);inputs={name:matrix[:,j] for j,name in enumerate(names)}
        values=engine.evaluate(inputs);assert int(values.max(initial=0))<=row['max_integer']
        # Python unbounded integer replay independently checks the vectorised
        # evaluator at all vertices and a retained set of random points.
        selected=sorted(set(range(min(1<<len(names),len(addresses))))|set(rng.integers(0,len(addresses),64).tolist()))
        maximum_body_error=D(0);eta=F.from_float(src['arithmetic_price_error_upper'])/F(src['cap'])
        for j in selected:
            address={name:int(matrix[j,k]) for k,name in enumerate(names)}
            expected=graph.evaluate(address,cache)
            assert int(values[j])==expected
            ideal=ideal_reference(model,src,address,cache)
            with localcontext() as ctx:
                ctx.prec=180
                error=abs(dec(F(expected,1<<src['q']))-ideal)
                assert error+D('1e-70')<=dec(eta)
                maximum_body_error=max(maximum_body_error,error)
        records.append({'case_id':src['case_id'],'method':src['method'],'q':src['q'],
            'maximum_certified_integer':row['max_integer'],'maximum_observed_integer':int(values.max(initial=0)),
            'production_points':len(addresses),'Python_integer_cross_checks':len(selected),
            'latent_vertex_points':1<<len(names),'random_points':4096,
            'input_addresses_sha256':hashlib.sha256(matrix.tobytes()).hexdigest(),
            'finite_output_sha256':hashlib.sha256(values.tobytes()).hexdigest(),
            'independent_ideal_body_point_checks':len(selected),
            'maximum_independent_body_normalised_error_decimal':str(maximum_body_error),
            'retained_uniform_body_error_fraction':str(eta),
            'ideal_reference':'180digit direct exponentials or unclipped positive-series CDF;1e-70 numerical allowance',
            'table_manifest':engine.table_manifest,
            'scope':'production-precision finite-value diagnostics; not an empirical proof of the uniform range'})
        print('Production range replay',src['case_id'],src['method'],len(addresses),flush=True)
    return records


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source-root',type=Path)
    parser.add_argument('--skip-production',action='store_true')
    here=Path(__file__).resolve().parent;default=here if (here/'finite_range.py').exists() else here.parent
    output_root=here if (here/'value_oracle.py').exists() else default
    parser.add_argument('--output',type=Path,default=output_root/'results/finite_range_independent_audit.json')
    args=parser.parse_args();root=args.source_root or default;sys.path.insert(0,str(root))
    er=importlib.import_module('encoding_resources');gen=importlib.import_module('generalisation')
    ce=importlib.import_module('classical_engine');np=importlib.import_module('numpy')
    source_path=root/'results/architecture_resource_candidates.json';moment_path=root/'results/classical_variance_certificates.json'
    range_path=root/'results/finite_range_certificates.json'
    data=json.loads(source_path.read_text());moments=json.loads(moment_path.read_text());certificate=json.loads(range_path.read_text())
    assert certificate['configuration_sha256']==digest(source_path) and certificate['moment_sha256']==digest(moment_path)
    before={name:digest(root/name) for name in certificate['source_hashes']};assert before==certificate['source_hashes']
    started=time.perf_counter()
    with localcontext() as ctx:
        ctx.prec=450
        records,supports,source,models=analytic_audit(root,data,certificate,moments,er)
    production=[] if args.skip_production else production_replay(data,certificate,source,models,gen,ce,np)
    assert {name:digest(root/name) for name in before}==before
    out={'status':'independent exact support, convex vertices, outward ranges and lattice moments passed',
        'range_certificate_sha256':digest(range_path),'configuration_sha256':digest(source_path),'moment_sha256':digest(moment_path),
        'source_hashes':before,'audit_sha256':digest(Path(__file__)),
        'reference_decimal_precision':450,'reference_absolute_error_premise':str(REFERENCE_ERROR),
        'source_220digit_centre_comparison_allowance':str(SOURCE_CENTRE_ERROR),
        'normal_density_reference':'pi=6 asin(1/2), positive series with ratio below1/4 and tail below8 next_term',
        'configurations':records,'quantile_supports':[{'bits':b,'q':q,'support_maximum_fraction':str(v),'actual_table_sha256':sha} for (b,q),(v,sha) in sorted(supports.items())],
        'production_replays':production,'elapsed_seconds':time.perf_counter()-started,
        'summary':{'range_rows':len(records),'vertices_checked':sum(x['vertices_checked'] for x in records),
            'actual_quantile_supports':len(supports),'production_graphs':len(production),
            'production_points':sum(x['production_points'] for x in production)},
        'limits':['uniform result requires the retained body-error certificate and finite-moment certificate',
            'production point checks supplement the proof and do not establish an exhaustive support claim',
            'finite rounded quantile support does not bound the continuous Gaussian model',
            '450-digit exp/ln/sqrt are correctly rounded Decimal primitives; stated reference envelope is an explicit numerical premise']}
    args.output.parent.mkdir(exist_ok=True,parents=True);args.output.write_text(json.dumps(out,indent=2)+'\n');print(args.output)


if __name__=='__main__':main()
