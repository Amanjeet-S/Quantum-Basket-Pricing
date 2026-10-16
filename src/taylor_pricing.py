"""Finite-order, outward-certified common-factor basket Taylor prices.

Frozen exact rational Gaussian models, no empirical data. The order cap is 8.
No Taylor-series convergence claim is made; lognormal moments can grow.
"""
from pathlib import Path
from fractions import Fraction
from math import comb, factorial
from functools import lru_cache
import hashlib,json,time,platform,sys
import flint
from flint import arb,ctx
from continuous_quadrature import frac_arb, exact_endpoint, serialise_interval
ROOT=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
MAX_ORDER=8
PRECISION_BITS=224


def compositions(total,n):
    if n==1:
        yield (total,);return
    for first in range(total+1):
        for rest in compositions(total-first,n-1):yield (first,)+rest


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def model_parameters(model):
    w=tuple(map(Fraction,model['weights']))
    r=tuple(tuple(map(Fraction,a)) for a in model['loadings'])
    mu=tuple(map(Fraction,model['log_drifts']));lam=Fraction(model['lambda'])
    assert model['spot']==100 and model['strike']==100 and Fraction(model['maturity'])==1
    n=len(w);p=len(r[0]);rate=Fraction(model['risk_free_rate'])-Fraction(model['dividend_yield'])
    assert n==model['assets'] and p==model['outer_dimension'] and sum(w)==1 and min(w)>0 and lam>0
    assert all(mu[i]+(dot(r[i],r[i])+lam*lam)/2==rate for i in range(n))
    mg=dot(w,mu);rg=tuple(sum(w[i]*r[i][j] for i in range(n)) for j in range(p))
    d=tuple(x-mg for x in mu);s=tuple(tuple(x-y for x,y in zip(a,rg)) for a in r)
    assert dot(w,d)==0 and all(sum(w[i]*s[i][j] for i in range(n))==0 for j in range(p))
    return w,d,s,mg,rg,lam


def signed_stirling(order):
    table=[[1]]
    for k in range(1,order+1):
        prior=table[-1];row=[0]*(k+1)
        for j in range(1,k+1):row[j]=(prior[j-1] if j-1<len(prior) else 0)-(k-1)*(prior[j] if j<len(prior) else 0)
        table.append(row)
    return table


def normal_pdf(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()

def normal_cdf(x):return (1+(x/arb(2).sqrt()).erf())/2


def hermites(x,maximum):
    values=[arb(1)]
    if maximum>0:values.append(x)
    for j in range(1,maximum):values.append(x*values[-1]-j*values[-2])
    return values


def kernel_derivatives(m,scale,cap,order):
    """P^(ell), ell=0..order, closed form, all in outward Arb arithmetic."""
    ell=(1+arb(cap)/100).log();x0=m/scale;x1=(m-ell)/scale
    z0=x0+scale;z1=x1+scale
    q0=normal_cdf(z0)-normal_cdf(z1)
    e=(m+scale*scale/2).exp()
    values=[100*(e*q0-(normal_cdf(x0)-normal_cdf(x1))+(arb(cap)/100)*normal_cdf(x1))]
    if order==0:return values
    h0=hermites(z0,max(0,order-2));h1=hermites(z1,max(0,order-2))
    f0=normal_pdf(z0);f1=normal_pdf(z1)
    q=[q0]
    for j in range(1,order):q.append(((-1)**(j-1))*(h0[j-1]*f0-h1[j-1]*f1)/(scale**j))
    for derivative in range(1,order+1):
        values.append(100*e*sum(comb(derivative-1,j)*q[j] for j in range(derivative)))
    return values


def Fourier_constants(order):
    b=[1/(2*arb.pi()).sqrt(),1/arb.pi()]
    for j in range(2,order):b.append((j-1)*b[j-2])
    return b[:order]


def terms_for_geometry(model,max_moment):
    """Each multinomial coefficient, Gaussian exponent and shift is exact."""
    w,d,s,mg,rg,lam=model_parameters(model)
    return exact_geometry_terms(w,d,s,rg,max_moment)


@lru_cache(maxsize=12)
def exact_geometry_terms(w,d,s,rg,max_moment):
    terms=[]
    for j in range(max_moment+1):
        row=[]
        for alpha in compositions(j,len(w)):
            multinomial=factorial(j)//__import__('math').prod(factorial(a) for a in alpha)
            weight=Fraction(multinomial)
            for wi,a in zip(w,alpha):weight*=wi**a
            exposure=tuple(sum(alpha[i]*s[i][v] for i in range(len(w))) for v in range(len(rg)))
            exponent=dot(alpha,d)+dot(exposure,exposure)/2
            shift=dot(rg,exposure)
            row.append((alpha,weight,exponent,shift))
        assert len(row)==comb(len(w)+j-1,j)
        terms.append(row)
    return terms



@lru_cache(maxsize=12)
def prepared_geometry(w,d,s,rg,max_order,precision_bits):
    # Precision is part of the key; cached balls are never reused at a higher
    # advertised precision. All current production calls use 224 bits.
    assert ctx.prec==precision_bits
    terms=exact_geometry_terms(w,d,s,rg,max_order+1)
    moments_raw=[];weighted=[]
    for row in terms:
        translated=[(a,frac_arb(weight)*frac_arb(exponent).exp(),shift) for a,weight,exponent,shift in row]
        weighted.append(translated);moments_raw.append(sum((value for _,value,_ in translated),arb(0)))
    centred=[arb(1)]
    for k in range(1,max_order+2):
        value=sum(((-1)**(k-j))*comb(k,j)*moments_raw[j] for j in range(k+1))
        assert exact_endpoint(value,True)>=0
        centred.append(value)
    return terms,weighted,centred


def evaluate_family(model,caps,max_order=MAX_ORDER):
    start=time.perf_counter();ctx.prec=PRECISION_BITS
    w,d,s,mg,rg,lam=model_parameters(model)
    cache_hits_before=prepared_geometry.cache_info().hits
    terms,weighted,centred=prepared_geometry(w,d,s,rg,max_order,PRECISION_BITS)
    cache_hit=prepared_geometry.cache_info().hits>cache_hits_before
    stirling=signed_stirling(max_order+1);b=Fourier_constants(max_order+1)
    scale=frac_arb(lam*lam+dot(rg,rg)).sqrt();mg_arb=frac_arb(mg)
    geometric_setup=time.perf_counter()-start
    rows=[];work_start=time.perf_counter();evaluations=0
    for cap in caps:
        cap_start=time.perf_counter();sums=[]
        for j in range(max_order+1):
            accumulated=[arb(0)]*(max_order+1)
            for _,factor,shift in weighted[j]:
                derivatives=kernel_derivatives(mg_arb+frac_arb(shift),scale,cap,max_order)
                evaluations+=1
                for ell in range(max_order+1):accumulated[ell]+=factor*derivatives[ell]
            sums.append(accumulated)
        coefficients=[kernel_derivatives(mg_arb,scale,cap,0)[0]]
        for k in range(1,max_order+1):
            coefficient=arb(0)
            for ell in range(1,k+1):
                tilted=sum(((-1)**(k-j))*comb(k,j)*sums[j][ell] for j in range(k+1))
                coefficient+=stirling[k][ell]*tilted
            coefficients.append(coefficient/factorial(k))
        partial=arb(0);orders=[]
        for k in range(max_order+1):
            partial+=coefficients[k];r=k+1
            c=sum(abs(stirling[r][ell])*b[ell-1]/(frac_arb(lam)**ell) for ell in range(1,r+1))
            radius=cap*c*centred[r]/factorial(r)
            radius_upper=max(Fraction(0),exact_endpoint(radius,True))
            lower=exact_endpoint(partial,False)-radius_upper;upper=exact_endpoint(partial,True)+radius_upper
            # Exact cap range can intersect any mathematically valid enclosure.
            lower=max(Fraction(0),lower);upper=min(Fraction(cap),upper)
            assert lower<=upper
            half=(upper-lower)/2
            discount=(-frac_arb(Fraction(model['risk_free_rate'])*Fraction(model['maturity']))).exp()
            orders.append({'order':k,'coefficient_ball':coefficients[k].str(50,more=True),
                'coefficient_interval':serialise_interval(exact_endpoint(coefficients[k],False),exact_endpoint(coefficients[k],True)),
                'partial_sum_ball':partial.str(50,more=True),'remainder_constant_ball':c.str(50,more=True),
                'centred_moment_ball':centred[r].str(50,more=True),
                'remainder_upper_exact_fraction':str(radius_upper),
                'undiscounted_price_interval':serialise_interval(lower,upper),
                'discounted_price_interval':serialise_interval(exact_endpoint(frac_arb(lower)*discount,False),exact_endpoint(frac_arb(upper)*discount,True)),
                'meets_undiscounted_radius_001':half<=Fraction(1,100)})
        successful=[o for o in orders if o['meets_undiscounted_radius_001']]
        selected=min(successful or orders,key=lambda o:Fraction(o['undiscounted_price_interval']['upper_exact_fraction'])-Fraction(o['undiscounted_price_interval']['lower_exact_fraction']))
        rows.append({'case_id':f"{model['model_id']}_cap{cap}",'model_id':model['model_id'],'assets':model['assets'],'cap':cap,
            'orders':orders,'selected_order':selected['order'],'certified_within_001':bool(successful),
            'first_certifying_order':successful[0]['order'] if successful else None,
            'first_certifying_order_terms':sum(len(t) for t in terms[:successful[0]['order']+1]) if successful else None,
            'selected_undiscounted_price_interval':selected['undiscounted_price_interval'],
            'work_seconds_diagnostic':time.perf_counter()-cap_start,
            'kernel_calls':sum(len(t) for t in terms[:max_order+1]),
            'selected_order_terms':sum(len(t) for t in terms[:selected['order']+1])})
        print(rows[-1]['case_id'],'order',selected['order'],'radius',selected['undiscounted_price_interval']['midpoint_radius_upper_decimal'],'successful',bool(successful),flush=True)
    return rows,{'model_id':model['model_id'],'geometry_setup_seconds_diagnostic':geometric_setup,
        'geometry_terms_including_remainder':sum(map(len,terms)),
        'geometry_cache_hit':cache_hit,
        'new_geometry_exponential_evaluations':0 if cache_hit else sum(map(len,terms)),
        'kernel_calls':evaluations,
        'work_seconds_diagnostic':time.perf_counter()-work_start,
        'moment_intervals':[serialise_interval(exact_endpoint(m,False),exact_endpoint(m,True)) for m in centred]}


def main():
    start=time.perf_counter();protocol_path=ROOT/'results/generalisation_protocol.json';protocol=json.loads(protocol_path.read_text())
    rows=[];work=[]
    for model in protocol['models']:
        output,diagnostic=evaluate_family(model,protocol['caps']);rows.extend(output);work.append(diagnostic)
    out={'status':'finite-order certified deterministic Gaussian basket approximation',
        'prespecified_maximum_order':MAX_ORDER,'precision_bits':PRECISION_BITS,
        'price_unit':'unspecified price currency unit','data':'All 27 frozen synthetic models; no empirical data.',
        'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_source_sha256':{'continuous_quadrature.py':hashlib.sha256((ROOT/'continuous_quadrature.py').read_bytes()).hexdigest()},
        'python_flint_version':flint.__version__,'python_version':platform.python_version(),
        'configurations':rows,'model_work_diagnostics':work,
        'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'successful_prices':sum(r['certified_within_001'] for r in rows),
        'unresolved_prices':sum(not r['certified_within_001'] for r in rows),
        'claim_limit':'Finite-order enclosure only; no convergent infinite series, new generic Taylor theory, market accuracy or quantum runtime comparison is claimed.'}
    (ROOT/'results/taylor_pricing.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['successful_prices','unresolved_prices','elapsed_seconds_diagnostic']},indent=2))

if __name__=='__main__':main()
