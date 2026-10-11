"""Certified Acb quadrature for the frozen two-asset continuous price family.

The primary comparator integrates the analytic common-factor residual and
adds a separately enclosed geometric price. A direct conditional arithmetic
price, with another truncation and precision, is retained for independent
representation checking. Arb balls enclose numerical error; requested
tolerances are not treated as attained error certificates.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import flint
from flint import arb,acb,ctx

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent


def frac_arb(x):return arb(x.numerator)/x.denominator


def exact_endpoint(ball,upper):
    endpoint=ball.upper() if upper else ball.lower()
    mantissa,exponent=map(int,endpoint.man_exp())
    return Fraction(int(mantissa)*(1 << exponent)) if exponent>=0 else Fraction(int(mantissa),1 << -exponent)


def display_fraction(value,upwards=False):
    with localcontext() as c:
        c.prec=300;c.rounding=ROUND_CEILING if upwards else ROUND_FLOOR
        return str((Decimal(value.numerator)/Decimal(value.denominator)).quantize(Decimal('1e-30')))


def serialise_interval(lower,upper):
    assert lower<=upper
    return {'lower_exact_fraction':str(lower),'upper_exact_fraction':str(upper),
        'lower_decimal':display_fraction(lower),'upper_decimal':display_fraction(upper,True),
        'width_upper_decimal':display_fraction(upper-lower,True),
        'midpoint_radius_upper_decimal':display_fraction((upper-lower)/2,True)}


def model_checks(model):
    w=list(map(Fraction,model['weights']));r=[list(map(Fraction,row)) for row in model['loadings']]
    mu=list(map(Fraction,model['log_drifts']));lam=Fraction(model['lambda'])
    assert model['assets']==2 and model['outer_dimension']==1 and sum(w)==1 and min(w)>0 and lam>0
    assert all(len(row)==1 for row in r)
    rate=Fraction(model['risk_free_rate'])-Fraction(model['dividend_yield'])
    assert all(mu[i]+(r[i][0]**2+lam**2)/2==rate for i in range(2))
    rg=sum(w[i]*r[i][0] for i in range(2));geom=rg*rg
    assert geom==Fraction(model['geometric_outer_variance'])
    return w,r,mu,lam,rg


def capped_conditional(m,scale,cap):
    """Entire scalar kernel in m; branch-free normal CDF is evaluated by erf."""
    h=acb(cap)/100; strike_log=(1+h).log(); sqrt2=acb(2).sqrt()
    def phi(z):return (1+(z/sqrt2).erf())/2
    x0=m/scale;x1=(m-strike_log)/scale
    return 100*((m+scale*scale/2).exp()*(phi(x0+scale)-phi(x1+scale))
        -(phi(x0)-phi(x1))+h*phi(x1))


def evaluate(model,cap,precision_bits,truncation,absolute_tolerance,residual=True):
    ctx.prec=precision_bits
    w,r,mu,lam,rg=model_checks(model)
    wa=[frac_arb(x) for x in w];ra=[frac_arb(row[0]) for row in r]
    mua=[frac_arb(x) for x in mu];la=frac_arb(lam)
    mg0=sum(wa[i]*mua[i] for i in range(2));rg_a=frac_arb(rg)
    count={'calls':0,'analytic_verification_requests':0,'nonfinite_analytic_returns':0}
    normalisation=1/(2*arb.pi()).sqrt()
    def integrand(y,analytic):
        count['calls']+=1;count['analytic_verification_requests']+=int(analytic)
        logs=[mua[i]+ra[i]*y for i in range(2)]
        # This forwarding is required: no branch-cut assumption is inferred
        # from a real-valued integration path when Acb evaluates complex disks.
        ma=sum(wa[i]*logs[i].exp() for i in range(2)).log(analytic=analytic)
        out=capped_conditional(ma,la,cap)
        if residual:out-=capped_conditional(mg0+rg_a*y,la,cap)
        out*=(-y*y/2).exp()*normalisation
        if analytic and not out.is_finite():count['nonfinite_analytic_returns']+=1
        return out
    start=time.perf_counter()
    result=acb.integral(integrand,-truncation,truncation,
        rel_tol=arb(absolute_tolerance),abs_tol=arb(absolute_tolerance),
        eval_limit=100000,depth_limit=40)
    seconds=time.perf_counter()-start
    if not result.is_finite() or not result.imag.contains(0):
        raise ArithmeticError('A finite real integral enclosure was not obtained')
    integral=result.real
    tail=cap*(arb(truncation)/arb(2).sqrt()).erfc()
    # For real Y, 0<=residual<=H and 0<=arithmetic conditional price<=H.
    # The omitted tail interval is [0,2H Phi(-L)], not a point correction.
    correction=capped_conditional(acb(mg0),(la*la+rg_a*rg_a).sqrt(),cap).real if residual else arb(0)
    price_ball=integral+correction
    lower=exact_endpoint(price_ball,False)
    upper=exact_endpoint(price_ball,True)+exact_endpoint(tail,True)
    assert upper-lower<Fraction(1,1000)
    discount=(-frac_arb(Fraction(model['risk_free_rate'])*Fraction(model['maturity']))).exp()
    discounted_lower=exact_endpoint(frac_arb(lower)*discount,False)
    discounted_upper=exact_endpoint(frac_arb(upper)*discount,True)
    return {'case_id':f"{model['model_id']}_cap{cap}",'model_id':model['model_id'],'cap':cap,
        'representation':'common-factor residual plus exact geometric conditional price' if residual else 'direct common-factor conditional arithmetic price',
        'precision_bits':precision_bits,'normal_truncation':truncation,
        'requested_absolute_integration_tolerance':absolute_tolerance,
        'integral_ball':integral.str(45,more=True),'geometric_correction_ball':correction.str(45,more=True),
        'normal_tail_price_bound_ball':tail.str(45,more=True),
        'normal_tail_price_bound_upper_exact_fraction':str(exact_endpoint(tail,True)),
        'undiscounted_price_interval':serialise_interval(lower,upper),
        'discounted_price_interval':serialise_interval(discounted_lower,discounted_upper),
        'confidence':1,'requested_undiscounted_price_error':'.01',
        'statistical_error':0,'integration_diagnostics':count,
        'integration_seconds_diagnostic':seconds,
        'premises':'Exact frozen rational parameters; Arb/Acb outward ball arithmetic and a finite validated integral enclosure; log analytic flag forwarded; capped residual tail in [0,H].',
        'claim_limit':'Synthetic two-asset Gaussian model. No market calibration, empirical validation, quantum runtime comparison or theoretical novelty claim.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--protocol',type=Path,default=ROOT/'results/generalisation_protocol.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/continuous_quadrature.json')
    args=p.parse_args();protocol=json.loads(args.protocol.read_text())
    models=[m for m in protocol['models'] if m['assets']==2]
    start=time.perf_counter();rows=[];checks=[]
    for m in models:
        for cap in protocol['caps']:
            row=evaluate(m,cap,128,8,'1e-16',True)
            reference=evaluate(m,cap,176,10,'1e-24',False)
            a=row['undiscounted_price_interval'];b=reference['undiscounted_price_interval']
            lo=max(Fraction(a['lower_exact_fraction']),Fraction(b['lower_exact_fraction']))
            hi=min(Fraction(a['upper_exact_fraction']),Fraction(b['upper_exact_fraction']))
            assert lo<=hi,(row['case_id'],'disjoint rigorous representations')
            checks.append({'case_id':row['case_id'],'residual_and_direct_enclosures_overlap':True,
                'reference':reference,'intersection':serialise_interval(lo,hi)})
            rows.append(row)
            print('Certified continuous price',row['case_id'],a['lower_decimal'],a['upper_decimal'],flush=True)
    out={'status':'certified continuous two-asset deterministic comparator','data':'Frozen synthetic models only',
        'protocol_sha256':hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'models':models,'configurations':rows,'independent_representation_checks':checks,
        'python_flint_version':flint.__version__,'python_version':platform.python_version(),
        'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'runtime_qualification':'Already installed authorised python-flint runtime; no new package, subscription or entitlement.',
        'novelty_boundary':'Validated deterministic quadrature is established. This is a stronger matched classical comparator, not a new pricing model or proof of quantum advantage.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'prices':len(rows),'elapsed_seconds_diagnostic':out['elapsed_seconds_diagnostic'],
        'total_primary_integrand_calls':sum(r['integration_diagnostics']['calls'] for r in rows)},indent=2))


if __name__=='__main__':main()
