"""Bounded degree-three continuous price test with the transport norm inputs.

Established Gaussian tilting and finite Taylor bounds. Full continuous ratio
moments replace hybrid transport moments: they are additional scalar queries,
not a price answer, but must not be claimed literally present in raw transport
moments. No earlier higher-degree search is rerun.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import hashlib,json,time,platform
import flint
from flint import arb,ctx
from taylor_pricing import model_parameters,exact_geometry_terms,dot
from boundary_derivative_certificate import payoff_norms
from boundary_grid_certificate import ball,serial,endpoints,interval,phi,coefficient_certificate
ROOT=Path(__file__).resolve().parent
PRECISION=512

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def geometry(model):
    w,d,s,mg,rg,lam=model_parameters(model);terms=exact_geometry_terms(w,d,s,rg,4)
    weighted=[];raw=[]
    for row in terms:
        out=[(alpha,ball(coefficient)*ball(exponent).exp(),shift) for alpha,coefficient,exponent,shift in row]
        weighted.append(out);raw.append(sum((a for _,a,_ in out),arb(0)))
    centred=[]
    for k in range(5):
        v=sum(((-1)**(k-j)*comb(k,j)*raw[j] for j in range(k+1)),arb(0));lo,hi=endpoints(v)
        if hi<0:raise ArithmeticError('A nonnegative continuous U moment is excluded')
        centred.append(interval(max(F(0),lo),hi))
    return {'weighted':weighted,'moments':centred,'terms':terms,'parameters':(w,d,s,mg,rg,lam)}

def boundary23(m,sigma,cap):
    z=m/sigma;f=phi(z);c2=100*f/sigma;c3=100*(-z/sigma-2)*f/sigma
    if cap is not None:
        h=ball(cap);ell=(1+h/100).log();z2=(m-ell)/sigma;f2=phi(z2)
        c2-=(100+h)*f2/sigma;c3-=(100+h)*(-z2/sigma-2)*f2/sigma
    return c2,c3

def evaluate(model,cap,g):
    start=time.perf_counter();w,d,s,mg,rg,lam=g['parameters'];sigma=ball(lam*lam+dot(rg,rg)).sqrt()
    coeff=coefficient_certificate(model,cap,'vanilla' if cap is None else 'capped',1)
    c0=interval(F(coeff['c0']['lower_exact_fraction']),F(coeff['c0']['upper_exact_fraction']))
    c1=interval(F(coeff['c1']['lower_exact_fraction']),F(coeff['c1']['upper_exact_fraction']))
    norms=payoff_norms(lam,cap);sums=[];kernel_calls=0;density_calls=0
    for row in g['weighted'][:4]:
        total2=total3=arb(0)
        for _,factor,shift in row:
            v2,v3=boundary23(ball(mg+shift),sigma,cap);total2+=factor*v2;total3+=factor*v3
            kernel_calls+=1;density_calls+=1 if cap is None else 2
        sums.append((total2,total3))
    cs=[c0,c1]
    for k in (2,3):cs.append(sum(((-1)**(k-j)*comb(k,j)*sums[j][k-2] for j in range(k+1)),arb(0))/factorial(k))
    partial=c0;orders=[];discount=ball(-F(model['risk_free_rate'])*F(model['maturity'])).exp()
    for k in (1,2,3):
        partial+=cs[k];M=F(norms['C_r_upper'][str(k+1)])
        radius=M*endpoints(g['moments'][k+1])[1]/factorial(k+1)
        lo,hi=endpoints(partial);lo=max(F(0),lo-(F(0) if cap is None and k==1 else radius));hi+=radius
        if cap is not None:hi=min(F(cap),hi)
        if lo>hi:raise ArithmeticError('Price enclosure empty')
        price=interval(lo,hi);disc=price*discount
        orders.append({'order':k,'coefficient_interval':serial(cs[k]),'partial_sum_interval':serial(partial),
          'remainder_M_upper':str(M),'continuous_ratio_moment_interval':serial(g['moments'][k+1]),
          'remainder_upper':str(radius),'continuous_price_interval':serial(price),
          'undiscounted_radius_upper':str((hi-lo)/2),'certifies_undiscounted_cent':(hi-lo)/2<=F(1,100),
          'discounted_price_interval':serial(disc),
          'required_kernel_terms':sum(len(row) for row in g['terms'][:k+1])})
    return {'model_id':model['model_id'],'payoff':'vanilla' if cap is None else 'capped','cap':None if cap is None else str(cap),
      'orders':orders,'first_successful_order':next((r['order'] for r in orders if r['certifies_undiscounted_cent']),None),
      'additional_C2_C3_kernel_evaluations':kernel_calls,'additional_Gaussian_density_evaluations':density_calls,
      'c0_c1_input_certificate':coeff,'computation_seconds_diagnostic':time.perf_counter()-start}

def main():
    ctx.prec=PRECISION;start=time.perf_counter();protocol=ROOT/'results/generalisation_protocol.json';models=json.loads(protocol.read_text())['models']
    cache={};records=[];geometry_evidence=[]
    for model in models:
        w,d,s,mg,rg,lam=model_parameters(model);key=(w,d,s,rg)
        fresh=key not in cache;begin=time.perf_counter()
        if fresh:cache[key]=geometry(model)
        g={**cache[key],'parameters':(w,d,s,mg,rg,lam)}
        if fresh:
            geometry_evidence.append({'assets':model['assets'],'continuous_U_moments':{str(k):serial(v) for k,v in enumerate(g['moments'])},
              'Gaussian_geometry_exponential_evaluations':sum(map(len,g['terms'])),'monomial_counts_by_degree':list(map(len,g['terms'])),
              'construction_seconds_diagnostic':time.perf_counter()-begin})
        for cap in (F(1),F(10),F(25),None):records.append(evaluate(model,cap,g))
    result={'schema':'qbp-degree3-same-information-classical-test-v1','status':'complete; unsuccessful degree caps retained',
      'precision_bits':PRECISION,'degree_cap':3,'remainder_degree_cap':4,'records':records,'geometry_evidence':geometry_evidence,
      'source_sha256':sha(__file__),'protocol_sha256':sha(protocol),
      'dependencies':{name:sha(ROOT/name) for name in ('taylor_pricing.py','boundary_derivative_certificate.py','boundary_grid_certificate.py')},
      'norm_input_sha256':sha(ROOT/'results/boundary_derivative_norms.json'),
      'generation_seconds_diagnostic':time.perf_counter()-start,'python_version':platform.python_version(),'python_flint_version':flint.__version__,
      'scope':'Same frozen model and global conditional boundary norms as transport. Degree<=4 continuous exponential geometry plus additional tilted Gaussian kernels. Hybrid transport moments are not literally continuous price coefficients. No financial grid, prior basket-price answer, old higher-degree search or quantum circuit is executed.'}
    out=ROOT/'results/residual_same_information_taylor.json';out.write_text(json.dumps(result,indent=2,default=str)+'\n')
    print('degree3 certified',sum(r['first_successful_order'] is not None for r in records),'of',len(records),'seconds',result['generation_seconds_diagnostic'],flush=True)
if __name__=='__main__':main()
