"""Continuous ordinary-call analytical screens in every frozen model.

No financial-grid, finite-grid moment or price-reference file is an input.
The screen is an established convex Taylor consequence, not a new theorem.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from datetime import datetime,timezone
import hashlib,json,time
from flint import arb,ctx
from continuous_quadrature import frac_arb,exact_endpoint,serialise_interval

ROOT=Path(__file__).resolve().parent
PRECISION=512
ENDPOINT_BITS=640

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def endpoints(a):
    def outward(x,up):
        x=F(x);den=1<<ENDPOINT_BITS;v=x*den;k=v.numerator//v.denominator
        return F(k+int(up and v.denominator!=1),den)
    return outward(exact_endpoint(a,False),False),outward(exact_endpoint(a,True),True)

def serial(a):return serialise_interval(*endpoints(a))

def geometry(model):
    w=list(map(F,model['weights']));rr=[list(map(F,r)) for r in model['loadings']];mu=list(map(F,model['log_drifts']));lam=F(model['lambda'])
    if min(w)<=0 or sum(w)!=1 or lam<=0:raise ValueError('Positive centred financial model required')
    p=len(rr[0]);dot=lambda a,b:sum(x*y for x,y in zip(a,b));mg=dot(w,mu)
    rg=[sum(wi*r[j] for wi,r in zip(w,rr)) for j in range(p)]
    d=[a-mg for a in mu];s=[[a-b for a,b in zip(r,rg)] for r in rr]
    if dot(w,d)!=0 or any(sum(wi*r[j] for wi,r in zip(w,s))!=0 for j in range(p)):raise ValueError('Centred geometry failed')
    for drift,row in zip(mu,rr):
        if drift!=F(model['risk_free_rate'])-F(model['dividend_yield'])-(lam*lam+dot(row,row))/2:raise ValueError('Drift does not describe the frozen risk-neutral GBM')
    return w,d,s,mg,rg,lam

class Arithmetic:
    def __init__(self):
        self.counts=Counter();self.counts['square_root_calls']=2
        self.sqrt2=arb(2).sqrt();self.sqrt2pi=(2*arb.pi()).sqrt()
    def exp(self,a):self.counts['exponential_calls']+=1;return a.exp()
    def cdf(self,a):self.counts['Gaussian_CDF_calls']+=1;return (1+(a/self.sqrt2).erf())/2
    def sqrt(self,a):self.counts['square_root_calls']+=1;return a.sqrt()


def screen(model,tolerance=F(1,100)):
    ctx.prec=PRECISION;start=time.perf_counter();w,d,s,mg,rg,lam=geometry(model);K=F(model['strike'])
    if K!=100 or F(model['spot'])!=100:raise ValueError('Frozen spot/strike100 is required')
    a=Arithmetic();dot=lambda x,y:sum(u*v for u,v in zip(x,y));first_terms=[];first=arb(0)
    for wi,di,si in zip(w,d,s):
        exponent=di+dot(si,si)/2;value=a.exp(frac_arb(exponent));term=frac_arb(wi)*value;first+=term
        first_terms.append({'weight':str(wi),'exponent':str(exponent),'exponential_interval':serial(value),'weighted_term_interval':serial(term)})
    second=arb(0);second_terms=[]
    for i in range(len(w)):
        for j in range(i,len(w)):
            vec=[x+y for x,y in zip(s[i],s[j])];exponent=d[i]+d[j]+dot(vec,vec)/2
            multiplicity=1 if i==j else 2;coefficient=multiplicity*w[i]*w[j]
            term=frac_arb(coefficient)*a.exp(frac_arb(exponent));second+=term
            second_terms.append({'i':i,'j':j,'multiplicity':multiplicity,'coefficient':str(coefficient),'exponent':str(exponent),'term_interval':serial(term)})
    U2=second-2*first+1;ulo,uhi=endpoints(U2)
    if uhi<0:raise ArithmeticError('Nonnegative continuous moment has a negative enclosure')
    variance=lam*lam+dot(rg,rg);sigma=a.sqrt(frac_arb(variance))
    base_growth=frac_arb(K)*a.exp(frac_arb(mg+variance/2));base_slope=base_growth*a.cdf(frac_arb(mg)/sigma+sigma)
    c0=base_slope-frac_arb(K)*a.cdf(frac_arb(mg)/sigma);c1=-base_slope;tilts=[]
    forward_exponent=F(model['risk_free_rate'])-F(model['dividend_yield'])
    common_forward=frac_arb(K)*a.exp(frac_arb(forward_exponent))
    for i,(wi,di,si) in enumerate(zip(w,d,s)):
        shifted=mg+dot(rg,si)
        # Exact risk-neutral forward equality cancels the two Gaussian-tilt
        # exponentials, rather than evaluating them separately.
        assert di+dot(si,si)/2+shifted+variance/2==forward_exponent
        term=frac_arb(wi)*common_forward*a.cdf(frac_arb(shifted)/sigma+sigma);c1+=term
        tilts.append({'i':i,'shifted_geometric_mean':str(shifted),'term_interval':serial(term)})
    coefficient=c0+c1;clo,chi=endpoints(coefficient)
    M2=frac_arb(K)/(frac_arb(lam)*a.sqrt2pi);mlo,mhi=endpoints(M2)
    moment_upper=max(F(0),uhi);correction=mhi*moment_upper/2
    undiscounted_lo=clo;undiscounted_hi=chi+correction
    discount=a.exp(frac_arb(-F(model['risk_free_rate'])*F(model['maturity'])));dl,du=endpoints(discount)
    products=[d*v for d in (dl,du) for v in (undiscounted_lo,undiscounted_hi)];lo,hi=min(products),max(products)
    point=(lo+hi)/2;radius=(hi-lo)/2
    a.counts['first_moment_terms']=len(first_terms);a.counts['symmetric_second_moment_terms']=len(second_terms)
    a.counts['Gaussian_tilt_terms']=len(tilts);a.counts['raw_ordered_second_moment_terms_replaced_by_symmetry']=len(w)**2
    return {'schema':'qbp-continuous-convex-call-screen-v1','model_id':model['model_id'],'assets':model['assets'],
        'outer_dimension':model['outer_dimension'],'lambda':str(lam),'effective_geometric_variance':str(variance),
        'continuous_T_mean_interval':serial(first),'continuous_T_square_mean_interval':serial(second),
        'continuous_U_square_mean_interval':serial(U2),'continuous_U_square_nonnegative_upper':str(moment_upper),
        'c0_interval':serial(c0),'c1_interval':serial(c1),'coefficient_interval':serial(coefficient),
        'global_M2_interval':serial(M2),'undiscounted_convex_correction_upper':str(correction),
        'undiscounted_price_interval':serialise_interval(undiscounted_lo,undiscounted_hi),
        'discount_interval':serial(discount),'price_interval':serialise_interval(lo,hi),'price_midpoint':str(point),
        'price_radius_upper':str(radius),'price_tolerance':str(tolerance),'stops_without_sampling':radius<=tolerance,
        'status':'certified tolerance stop' if radius<=tolerance else 'retained analytical interval exceeds tolerance',
        'first_moment_terms':first_terms,'second_moment_terms':second_terms,'coefficient_tilt_terms':tilts,
        'operation_counts':dict(a.counts),'runtime_seconds_diagnostic':time.perf_counter()-start,
        'precision_bits':PRECISION,'endpoint_lattice_bits':ENDPOINT_BITS,
        'scope':'Continuous Gaussian moments and geometric coefficients only. The established nonnegative second remainder gives this analytical interval, not a new generic theorem or empirical price.'}


def interval_from_record(record):
    lo,hi=F(record['lower_exact_fraction']),F(record['upper_exact_fraction']);mid=(lo+hi)/2;rad=(hi-lo)/2
    return frac_arb(mid)+arb(0,frac_arb(rad).upper())


def main():
    start=time.perf_counter();source=sha(__file__);path=ROOT/'results/generalisation_protocol.json';ph=sha(path)
    data=json.loads(path.read_text());rows=[screen(m) for m in data['models']]
    if len(rows)!=9 or sha(__file__)!=source or sha(path)!=ph:raise RuntimeError('The complete frozen family or source changed')
    result={'schema':'qbp-continuous-convex-call-screen-family-v1','status':'complete nine-model analytical screen',
        'created_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':source,'input_sha256':{'results/generalisation_protocol.json':ph},
        'dependency_sha256':{'continuous_quadrature.py':sha(ROOT/'continuous_quadrature.py')},
        'library_versions':{'python_flint':__import__('flint').__version__},'records':rows,'model_count':9,
        'tolerance_stop_count':sum(r['stops_without_sampling'] for r in rows),
        'total_runtime_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Established convex Taylor consequence: c0+c1 <= undiscounted V <= c0+c1 + K E(U^2)/(2 lambda sqrt(2 pi)), applied to all nine ordinary-call models. No finite financial grid, financial-grid certificate, finite-grid moment or continuous price answer is an input.'}
    (ROOT/'results/boundary_convex_screening.json').write_text(json.dumps(result,indent=2)+'\n')
    for r in rows:print(r['model_id'],'radius',float(F(r['price_radius_upper'])),'stops',r['stops_without_sampling'])
    return result

if __name__=='__main__':main()
