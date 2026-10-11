"""Independent formal-series, strike-integral and one-factor coefficient checks."""
import sys,json,time,hashlib
from pathlib import Path
from fractions import Fraction
from math import comb,factorial
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from flint import arb,acb,arb_series,acb_series,ctx
from taylor_pricing import kernel_derivatives,model_parameters,terms_for_geometry,Fourier_constants,signed_stirling
from continuous_quadrature import frac_arb,exact_endpoint,serialise_interval
ROOT=Path(__file__).resolve().parents[1]


def series_price(m,sigma,cap,maximum,complex_series=False):
    """Formal differentiation of the original price formula, independent AD."""
    ctx.cap=maximum+1
    constructor=acb_series if complex_series else arb_series
    x=constructor([m,1]);scale=constructor([sigma]);h=constructor([arb(cap)/100]);l=(1+h).log()
    def Phi(z):return (1+(z/arb(2).sqrt()).erf())/2
    z0=x/scale;z1=(x-l)/scale
    return 100*((x+scale*scale/2).exp()*(Phi(z0+scale)-Phi(z1+scale))
        -(Phi(z0)-Phi(z1))+h*Phi(z1))


def check_overlap(a,b):
    if isinstance(b,acb):assert b.imag.contains(0);b=b.real
    assert max(exact_endpoint(a,False),exact_endpoint(b,False))<=min(exact_endpoint(a,True),exact_endpoint(b,True))


def one_factor_tail(model,order,truncation):
    w,d,s,mg,rg,lam=model_parameters(model);assert len(rg)==1
    out=arb(0)
    # Direct binomial terms, not production compositions or moments.
    for j in range(order+1):
        for a in range(j+1):
            b=j-a;v=a*s[0][0]+b*s[1][0]
            factor=frac_arb(Fraction(comb(j,a))*w[0]**a*w[1]**b)
            factor*=frac_arb(a*d[0]+b*d[1]+v*v/2).exp()
            shift=frac_arb(v);root2=arb(2).sqrt()
            tail=((arb(truncation)+shift)/root2).erfc()/2+((arb(truncation)-shift)/root2).erfc()/2
            out+=comb(order,j)*factor*tail
    return out


def coefficient_integral(model,cap,order,truncation=12):
    w,d,s,mg,rg,lam=model_parameters(model);assert len(rg)==1
    mua=[frac_arb(Fraction(x)) for x in model['log_drifts']]
    ra=[frac_arb(Fraction(x[0])) for x in model['loadings']]
    wa=[frac_arb(x) for x in w];mg=frac_arb(mg);rg=frac_arb(rg[0]);la=frac_arb(lam)
    st=signed_stirling(order)[order];calls=0
    def integrand(y,analytic):
        nonlocal calls;calls+=1
        m=mg+rg*y
        derivative=series_price(m,la,cap,order,True)
        if order==0:out=derivative[0]
        else:
            T=sum(wa[i]*(mua[i]+ra[i]*y-m).exp() for i in range(2))
            falling=sum(st[j]*factorial(j)*derivative[j] for j in range(1,order+1))
            out=(T-1)**order*falling/factorial(order)
        return out*(-y*y/2).exp()/(2*arb.pi()).sqrt()
    result=acb.integral(integrand,-truncation,truncation,abs_tol=arb('1e-28'),rel_tol=arb('1e-28'),eval_limit=100000,depth_limit=40)
    assert result.is_finite() and result.imag.contains(0)
    if order==0:tail=cap*(arb(truncation)/arb(2).sqrt()).erfc()
    else:
        b=Fourier_constants(order);C=sum(abs(st[j])*b[j-1]/la**j for j in range(1,order+1))
        tail=cap*C*one_factor_tail(model,order,truncation)/factorial(order)
    lo=exact_endpoint(result.real,False)-exact_endpoint(tail,True)
    hi=exact_endpoint(result.real,True)+exact_endpoint(tail,True)
    return lo,hi,calls


def main():
    start=time.perf_counter();ctx.prec=288
    protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text())
    output=json.loads((ROOT/'results/taylor_pricing.json').read_text())
    rows={r['case_id']:r for r in output['configurations']}
    formal=[];strike=[];coeff=[]
    for sigma in [Fraction(1,20),Fraction(1,5),Fraction(2,5),Fraction(7,13)]:
        for m in [Fraction(-2),Fraction(-1,7),Fraction(0),Fraction(2,9),Fraction(3)]:
            for cap in [1,10,25]:
                ma=frac_arb(m);sa=frac_arb(sigma)
                derivatives=kernel_derivatives(ma,sa,cap,8);independent=series_price(ma,sa,cap,8)
                for ell in range(9):check_overlap(derivatives[ell],factorial(ell)*independent[ell])
                formal.append({'m':str(m),'sigma':str(sigma),'cap':cap,'derivatives_checked':9})
                if m in [Fraction(-1,7),Fraction(2,9)] and sigma in [Fraction(1,20),Fraction(2,5)]:
                    for ell in [0,1,2,5,8]:
                        def integrand(x,analytic):
                            z=(ma-x)/sa
                            if ell==0:base=(1+(z/arb(2).sqrt()).erf())/2
                            else:
                                # Literal polynomial recursion in the integrand.
                                hermite=[acb(1),z]
                                for v in range(1,ell-1):hermite.append(z*hermite[-1]-v*hermite[-2])
                                base=(-1)**(ell-1)*hermite[ell-1]*(-z*z/2).exp()/(2*arb.pi()).sqrt()/sa**ell
                            return 100*x.exp()*base
                        hi=(1+arb(cap)/100).log()
                        integral=acb.integral(integrand,0,hi,abs_tol=arb('1e-32'),rel_tol=arb('1e-32'),eval_limit=100000,depth_limit=40)
                        assert integral.is_finite();check_overlap(derivatives[ell],integral)
                        strike.append({'m':str(m),'sigma':str(sigma),'cap':cap,'derivative':ell,'positive_integrand':ell in [0,1]})
    for model in protocol['models']:
        if model['assets']!=2:continue
        for cap in protocol['caps']:
            for k in range(9):
                lo,hi,calls=coefficient_integral(model,cap,k)
                interval=rows[f"{model['model_id']}_cap{cap}"]['orders'][k]['coefficient_interval']
                assert max(lo,Fraction(interval['lower_exact_fraction']))<=min(hi,Fraction(interval['upper_exact_fraction']))
                coeff.append({'case_id':f"{model['model_id']}_cap{cap}",'order':k,'direct_gaussian_coefficient_interval':serialise_interval(lo,hi),'integrand_calls':calls})
        print('audited coefficients',model['model_id'],flush=True)
    polynomial=[];exposure_checks=[]
    for asset_count in [2,4,8]:
        model=next(m for m in protocol['models'] if m['assets']==asset_count)
        w,d,s,mg,rg,lam=model_parameters(model)
        expansion={(0,)*asset_count:Fraction(1)}
        production=terms_for_geometry(model,9)
        for j in range(10):
            actual={a:c for a,c,e,v in production[j]}
            assert expansion==actual
            assert sum(expansion.values())==1
            polynomial.append({'assets':asset_count,'degree':j,'multinomial_terms':len(actual),'recurrence_weights_match':True})
            indices=sorted(actual)
            for alpha in {indices[0],indices[-1],indices[len(indices)//2]}:
                # Independent pair-covariance expansion instead of the exposure norm.
                quadratic=sum(alpha[i]*alpha[v]*sum(s[i][z]*s[v][z] for z in range(len(rg))) for i in range(asset_count) for v in range(asset_count))
                expo=sum(alpha[i]*d[i] for i in range(asset_count))+quadratic/2
                shifted=sum(alpha[i]*sum(rg[z]*s[i][z] for z in range(len(rg))) for i in range(asset_count))
                match=next((e,v) for a,c,e,v in production[j] if a==alpha)
                assert match==(expo,shifted)
                exposure_checks.append({'assets':asset_count,'degree':j,'alpha':list(alpha),'pair_covariance_identity_matches':True})
            if j<9:
                multiplied={}
                for alpha,c in expansion.items():
                    for i,wi in enumerate(w):
                        beta=list(alpha);beta[i]+=1;beta=tuple(beta)
                        multiplied[beta]=multiplied.get(beta,Fraction(0))+c*wi
                expansion=multiplied
    old=json.loads((ROOT/'results/continuous_quadrature.json').read_text())
    whole=[]
    for row in old['configurations']:
        reference=row['undiscounted_price_interval']
        for order in rows[row['case_id']]['orders']:
            current=order['undiscounted_price_interval']
            assert Fraction(current['lower_exact_fraction'])<=Fraction(reference['lower_exact_fraction'])<=Fraction(reference['upper_exact_fraction'])<=Fraction(current['upper_exact_fraction'])
            whole.append({'case_id':row['case_id'],'order':order['order'],'retained_continuous_quadrature_interval_contained':True})
    out={'status':'passed','formal_series_checks':formal,'strike_integral_checks':strike,
        'one_factor_coefficient_quadrature_checks':coeff,'retained_price_checks':whole,
        'independent_polynomial_recurrence_checks':polynomial,'pair_covariance_exposure_checks':exposure_checks,
        'formal_derivatives_checked':len(formal)*9,'strike_integrals_checked':len(strike),'one_factor_coefficients_checked':len(coeff),
        'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'production_source_sha256':hashlib.sha256((ROOT/'taylor_pricing.py').read_bytes()).hexdigest(),
        'primary_archive_sha256':hashlib.sha256((ROOT/'results/taylor_pricing.json').read_bytes()).hexdigest(),
        'reference_archive_sha256':hashlib.sha256((ROOT/'results/continuous_quadrature.json').read_bytes()).hexdigest(),
        'protocol_sha256':hashlib.sha256((ROOT/'results/generalisation_protocol.json').read_bytes()).hexdigest(),
        'dependency_source_sha256':{'continuous_quadrature.py':hashlib.sha256((ROOT/'continuous_quadrature.py').read_bytes()).hexdigest()},
        'scope':'Arb formal-series AD, independent strike-integral derivatives and direct n2 coefficient quadrature test the implementation; the higher-dimensional guarantee additionally uses the mathematical proof.'}
    (ROOT/'results/taylor_pricing_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['formal_derivatives_checked','strike_integrals_checked','one_factor_coefficients_checked','elapsed_seconds_diagnostic']},indent=2))

if __name__=='__main__':main()
