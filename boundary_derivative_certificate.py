"""Certified payoff-specific log-price Greek norms for call boundaries.

The falling log-price derivative operator is established mathematics. This
module certifies its concrete capped and vanilla Gaussian boundary kernels.
No new stochastic asset model or generic option-pricing priority is asserted.
Generation uses the existing python-flint runtime; retained API reads stdlib.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json,time
import gaussian_derivative_norm as polynomial
ROOT=Path(__file__).resolve().parent
PRECISION=512
ROOT_BITS=160


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ball(x):
    from flint import arb
    x=F(x);return arb(x.numerator)/arb(x.denominator)

def endpoints(x):
    from continuous_quadrature import exact_endpoint
    if not x.is_finite():raise ArithmeticError('Nonfinite payoff boundary enclosure')
    return exact_endpoint(x,False),exact_endpoint(x,True)

def interval(lo,hi):
    from flint import arb
    lo,hi=F(lo),F(hi)
    if lo>hi:raise ValueError('Reversed interval')
    return ball((lo+hi)/2)+arb(0,ball((hi-lo)/2).upper())

def serial(x):
    lo,hi=endpoints(x);return {'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)}

def polynomial_ball(p,x):
    from flint import arb
    value=arb(0)
    for coefficient in reversed(p):value=value*x+ball(coefficient)
    return value


def density(z):
    from flint import arb
    return (-z*z/2).exp()/(2*arb.pi()).sqrt()

def cdf(z):
    from flint import arb
    return (1+(z/arb(2).sqrt()).erf())/2


def boundary_polynomial(order,lam):
    lam=F(lam)
    if type(order) is not int or not 2<=order<=6 or lam<=0:raise ValueError('Supported boundary orders are2through6 and lambda is positive')
    p=(F(1),)
    for r in range(2,order):
        p=polynomial.add(polynomial.scale(polynomial.add(polynomial.derivative(p),polynomial.scale(polynomial.shift(p),-1)),1/lam),polynomial.scale(p,-r))
    return p


def polynomial_norm(order,lam,bits=ROOT_BITS):
    from flint import ctx
    ctx.prec=PRECISION;lam=F(lam);p=boundary_polynomial(order,lam)
    q=polynomial.add(polynomial.derivative(p),polynomial.scale(polynomial.shift(p),-1))
    roots,chain,bound,splits=polynomial.isolate(q,bits)
    lo_norm=hi_norm=F(0);rows=[]
    for lo,hi in roots:
        x=interval(lo,hi);value=abs(polynomial_ball(p,x)*density(x))
        sample=abs(polynomial_ball(p,ball((lo+hi)/2))*density(ball((lo+hi)/2)))
        lower=max(F(0),endpoints(sample)[0]);upper=endpoints(value)[1]
        lo_norm=max(lo_norm,lower);hi_norm=max(hi_norm,upper)
        rows.append({'lower':str(lo),'upper':str(hi),'variation_left':polynomial.variations(chain,lo),
            'variation_right':polynomial.variations(chain,hi),'norm_upper':str(upper),'sample_lower':str(lower)})
    if order==2:
        direct=density(ball(0));hi_norm=min(hi_norm,endpoints(direct)[1]);lo_norm=max(lo_norm,endpoints(direct)[0])
    assert 0<lo_norm<=hi_norm
    return {'order':order,'lambda':str(lam),'polynomial':[str(x) for x in p],
        'stationary_polynomial':[str(x) for x in q],'sturm_chain':[[str(x) for x in row] for row in chain],
        'cauchy_bound':str(bound),'root_isolation_bits':bits,'root_intervals':rows,
        'real_stationary_points':len(rows),'subdivision_steps':splits,'precision_bits':PRECISION,
        'norm_lower':str(lo_norm),'norm_upper':str(hi_norm),'passed':True}


def monotone_root(function,lo,hi,margin,bits=ROOT_BITS,domain=None):
    """Decreasing sign root, including outward uncertain model constants."""
    lo,hi,margin=map(F,(lo,hi,margin));target=F(1,1<<bits)
    assert margin>0 and endpoints(function(lo))[0]>0 and endpoints(function(hi))[1]<0
    steps=0;ambiguous=0
    while hi-lo>target:
        mid=(lo+hi)/2;fl,fu=endpoints(function(mid));steps+=1
        if fl>0:lo=mid
        elif fu<0:hi=mid
        else:
            ambiguous+=1
            radius=2*max(abs(fl),abs(fu))/margin+target/4
            nl,nh=max(lo,mid-radius),min(hi,mid+radius)
            if domain is not None:assert domain(nl,nh)
            # The positive monotonicity margin proves the root inclusion;
            # fresh endpoint signs also certify the retained bracket directly.
            assert endpoints(function(nl))[0]>0 and endpoints(function(nh))[1]<0
            if nh-nl>=hi-lo:raise ArithmeticError('Precision insufficient for the requested root bracket')
            lo,hi=nl,nh
        if steps>4*bits+2048:raise ArithmeticError('Root isolation did not reach its declared tolerance')
    if domain is not None:assert domain(lo,hi)
    return {'lower':str(lo),'upper':str(hi),'width':str(hi-lo),'steps':steps,'ambiguous_midpoint_refinements':ambiguous,
        'left_function_interval':serial(function(lo)),'right_function_interval':serial(function(hi)),
        'strict_decreasing_derivative_margin':str(margin),'root_isolation_bits':bits}


def capped_C1(m,lam,K,H):
    lam,K,H=map(ball,(lam,K,H));ell=(1+H/K).log()
    return K*(m+lam*lam/2).exp()*(cdf((m+lam*lam)/lam)-cdf((m-ell+lam*lam)/lam))


def capped_C2_z(z,lam,K,H):
    lam=F(lam);K,H=map(F,(K,H));c=1+H/K;a=ball(c).log()/ball(lam)
    return ball(K/lam)*(density(z)-ball(c)*density(z-a))


def capped_second_norm(lam,H,K=100):
    from flint import ctx
    ctx.prec=PRECISION;lam,H,K=map(F,(lam,H,K))
    if min(lam,H,K)<=0:raise ValueError('Positive volatility, cap and strike required')
    c=1+H/K;ell=ball(c).log();a=ell/ball(lam);alo,ahi=endpoints(a)
    def equation(z):
        z=ball(z);return (z/(z-a)).log()-a*z-ell+a*a/2
    left=F(-1)
    while endpoints(equation(left))[0]<=0:left*=2
    right=F(-1,2)
    while endpoints(equation(right))[1]>=0:right/=2
    negative=monotone_root(equation,left,right,alo,domain=lambda lo,hi:hi<0)
    delta=F(1,2);left=ahi+delta
    while endpoints(equation(left))[0]<=0:delta/=2;left=ahi+delta
    delta=F(1);right=ahi+delta
    while endpoints(equation(right))[1]>=0:delta*=2;right=ahi+delta
    positive=monotone_root(equation,left,right,alo,domain=lambda lo,hi:lo>ahi)
    peaks=[];lower=upper=F(0)
    for label,root in [('positive_maximum',negative),('negative_minimum',positive)]:
        lo,hi=F(root['lower']),F(root['upper']);value=capped_C2_z(interval(lo,hi),lam,K,H)
        sample=abs(capped_C2_z(ball((lo+hi)/2),lam,K,H));pl=max(F(0),endpoints(sample)[0]);pu=endpoints(abs(value))[1]
        assert endpoints(value)[0]>0 if label=='positive_maximum' else endpoints(value)[1]<0
        lower=max(lower,pl);upper=max(upper,pu)
        peaks.append({'label':label,'stationary_root':root,'value_interval':serial(value),'absolute_peak_lower':str(pl),'absolute_peak_upper':str(pu),
            'stationary_derivative_interval':serial(-interval(lo,hi)*density(interval(lo,hi))+ball(c)*(interval(lo,hi)-a)*density(interval(lo,hi)-a))})
    return {'norm_lower':str(lower),'norm_upper':str(upper),'c':str(c),'log_width_interval':serial(ell),
        'scaled_width_interval':serial(a),'stationary_peaks':peaks,'passed':True,
        'proof':'Exactly one strictly decreasing log-ratio root on each outside component; no stationary point between zero and a; both Gaussian tails vanish.'}


def capped_first_norm(lam,H,K=100):
    from flint import ctx
    ctx.prec=PRECISION;lam,H,K=map(F,(lam,H,K))
    if min(lam,H,K)<=0:raise ValueError('Positive volatility, cap and strike required')
    ell=ball(1+H/K).log();elllo,ellhi=endpoints(ell)
    dmin=density(ball((ellhi+lam*lam)/lam))/ball(lam)
    gamma=endpoints(ball(elllo/(lam*lam))*dmin)[0]
    assert gamma>0
    def score(m):
        mm=ball(m);p=capped_C1(mm,lam,K,H)
        assert p.lower()>0
        return (p+capped_C2_z(mm/ball(lam),lam,K,H))/p
    root=monotone_root(score,F(0),ellhi,gamma,domain=lambda lo,hi:lo>=0 and hi<=ellhi)
    lo,hi=F(root['lower']),F(root['upper']);peak=capped_C1(interval(lo,hi),lam,K,H)
    sample=capped_C1(ball((lo+hi)/2),lam,K,H)
    lower=max(F(0),endpoints(sample)[0]);upper=endpoints(peak)[1]
    assert 0<lower<=upper
    return {'norm_lower':str(lower),'norm_upper':str(upper),'log_width_interval':serial(ell),
        'stationary_root':root,'peak_interval':serial(peak),'uniform_log_concavity_margin_lower':str(gamma),
        'second_derivative_at_root_interval':serial(peak+capped_C2_z(interval(lo,hi)/ball(lam),lam,K,H)),
        'passed':True,'proof':'Truncated-normal integration by parts gives variance strictly below lambda squared, so the positive slope has one peak inside the log-strike interval.'}


@lru_cache(None)
def generic_norms():
    x=json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())
    if x['status']!='passed':raise ValueError('Unaccepted generic Gaussian derivative norms')
    return {(F(r['lambda']),r['derivative_order']):F(r['norm_upper']) for r in x['records'] if r['passed']}


@lru_cache(None)
def retained():
    data=json.loads((ROOT/'results/boundary_derivative_norms.json').read_text())
    if data['status']!='passed':raise ValueError('Unaccepted payoff boundary norms')
    return data


def payoff_norms(lam,cap=None,K=100):
    lam,K=map(F,(lam,K));data=retained();key='vanilla' if cap is None else 'capped';cap=None if cap is None else F(cap)
    row=next(r for r in data['records'] if r['payoff']==key and F(r['lambda'])==lam and F(r['K'])==K and (cap is None or F(r['cap'])==cap))
    return row


def main():
    from flint import ctx
    ctx.prec=PRECISION;start=time.perf_counter();records=[];poly=[];generic=generic_norms()
    hermite={r['unit_derivative_order']:F(r['norm_upper']) for r in json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())['unit_norms'] if r['passed']}
    for lam in (F(1,20),F(1,5),F(2,5)):
        norms={r:polynomial_norm(r,lam) for r in range(2,7)};poly.extend(norms.values())
        K=F(100);v2=ball(K/lam)*density(ball(0));v2lo,v2hi=endpoints(v2)
        vanilla={str(r):str(K/lam*F(norms[r]['norm_upper'])) for r in range(2,7)};vanilla['2']=str(v2hi)
        records.append({'payoff':'vanilla','K':str(K),'cap':None,'lambda':str(lam),'M1_upper':None,
            'M2_lower':str(v2lo),'M2_upper':str(v2hi),'C_r_upper':vanilla,
            'first_derivative_global_bounded':False,'second_remainder_nonnegative':True,
            'scalar_boundary_density_first_upper':str(K*F(next(u['norm_upper'] for u in json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())['unit_norms'] if u['unit_derivative_order']==1))/(lam*lam*v2hi)),
            'scalar_boundary_density_second_upper':str(K*F(next(u['norm_upper'] for u in json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())['unit_norms'] if u['unit_derivative_order']==2))/(lam**3*v2hi)),
            'passed':True,'scope':'Known nonnegative vanilla cash-Gamma boundary kernel; unbounded raw payoff with finite Gaussian moments, not a new risk law.'})
        for H in (F(1),F(10),F(25)):
            c1=capped_first_norm(lam,H,K);c2=capped_second_norm(lam,H,K)
            m1=min(F(c1['norm_upper']),H*generic[lam,1],K+H);m2=min(F(c2['norm_upper']),H*generic[lam,2],(2*K+H)*F(norms[2]['norm_upper'])/lam)
            assert F(c1['norm_lower'])<=m1 and F(c2['norm_lower'])<=m2
            bounds={str(r):str(min(H*generic[lam,r],(2*K+H)*F(norms[r]['norm_upper'])/lam)) for r in range(3,7)};bounds['2']=str(m2)
            rows={'payoff':'capped','K':str(K),'cap':str(H),'lambda':str(lam),'M1_lower':c1['norm_lower'],'M1_upper':str(m1),
                'M2_lower':c2['norm_lower'],'M2_upper':str(m2),'C_r_upper':bounds,'C1_root_certificate':c1,'C2_root_certificate':c2,
                'generic_H_D1_upper':str(H*generic[lam,1]),'generic_H_D2_upper':str(H*generic[lam,2]),
                'competent_first_derivative_upper':str(min(K+H,H*generic[lam,1])),
                'M1_to_generic_ratio_upper':str(m1/(H*generic[lam,1])),'M2_to_generic_ratio_upper':str(m2/(H*generic[lam,2])),
                'scalar_cap_first_upper':str(m1),'scalar_cap_second_upper':str(min(m1+m2,H*hermite[1]/lam**2)),
                'scalar_cap_prime_first_upper':str(min(m1+m2,H*hermite[1]/lam**2)),
                'scalar_cap_prime_second_upper':str(min(m1+3*m2+F(bounds['3']),H*hermite[2]/lam**3)),
                'scalar_boundary_h_first_upper':str(min(F(bounds['3'])+2*m2,(2*K+H)*hermite[1]/lam**2,H*(hermite[1]/lam**2+hermite[2]/lam**3))/(2*m2)),
                'scalar_boundary_h_second_upper':str(min(F(bounds['4'])+5*F(bounds['3'])+4*m2,(2*K+H)*hermite[2]/lam**3,H*(hermite[2]/lam**3+hermite[3]/lam**4))/(2*m2)),
                'passed':True,'scope':'Sharp first/second payoff-specific Gaussian norms; higher norms are certified minima of generic and boundary-triangle envelopes shared across all methods.'}
            records.append(rows);print('lambda',lam,'H',H,'M1 ratio',float(F(rows['M1_to_generic_ratio_upper'])),'M2 ratio',float(F(rows['M2_to_generic_ratio_upper'])),flush=True)
    result={'schema':'qbp-payoff-boundary-derivative-norm-v1','status':'passed','records':records,'boundary_polynomial_norms':poly,
        'root_isolation_bits':ROOT_BITS,'precision_bits':PRECISION,'source_sha256':sha(__file__),
        'input_sha256':{'results/gaussian_derivative_norms.json':sha(ROOT/'results/gaussian_derivative_norms.json')},
        'dependency_source_sha256':{name:sha(ROOT/name) for name in ['gaussian_derivative_norm.py','continuous_quadrature.py']},
        'generation_seconds_diagnostic':time.perf_counter()-start,
        'priority':'Falling log-price Greek operators and vanilla/CPWA option pricing are established; the retained project-specific root and error certificates imply no world novelty.'}
    (ROOT/'results/boundary_derivative_norms.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Certified',len(records),'payoff cases and',len(poly),'boundary polynomial norms',flush=True)

if __name__=='__main__':main()
