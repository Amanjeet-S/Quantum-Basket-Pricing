"""Decimal endpoint enclosures for the root constructor's function tables.

The normal CDF uses an alternating entire-function series, with its explicit
remainder. Decimal roundoff is bounded by a deliberately loose 1e-50 envelope:
100 significant decimal digits, <=2000 operations, and intermediates <=1e30.
The script rejects a calculation outside these predeclared limits. Pi uses
Machin's identity and alternating arctangent sums. Exp/ln use Decimal's
correctly rounded implementations. This certifies constants, not QPU faults.
"""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import json
import math
from pathlib import Path
import sys
import numpy as np

ZERO=Decimal(0); ONE=Decimal(1); HALF=Decimal('.5')
EPS=Decimal('1e-50')


def arctan_series(x):
    total=ZERO; term=x; x2=x*x; n=0
    while True:
        add=term/Decimal(2*n+1)
        total+=add
        term*=-x2; n+=1
        if abs(term/Decimal(2*n+1))<Decimal('1e-96'):
            return total


def constants():
    with localcontext() as ctx:
        ctx.prec=100
        pi=16*arctan_series(ONE/5)-4*arctan_series(ONE/239)
        return ONE/(2*pi).sqrt()


PHI0=constants()


def decimal_cdf(x):
    x=Decimal(x)
    if abs(x)>8: raise ValueError('CDF enclosure domain exceeded')
    if x==0: return HALF,ZERO
    with localcontext() as ctx:
        ctx.prec=100
        term=x; total=term; n=0; peak=abs(term)
        while True:
            nxt=-term*x*x*Decimal(2*n+1)/(2*Decimal(n+1)*Decimal(2*n+3))
            n+=1
            peak=max(peak,abs(nxt),abs(total))
            if n>32 and abs(nxt)<Decimal('1e-80'):
                assert n<250 and peak<Decimal('1e30')
                return HALF+PHI0*total, EPS+abs(PHI0*nxt)
            total+=nxt; term=nxt


def rounded_enclosure(value,error,p):
    with localcontext() as ctx:
        ctx.prec=100
        scale=Decimal(1<<p)
        lo=((value-error)*scale).to_integral_value(rounding=ROUND_HALF_EVEN)
        hi=((value+error)*scale).to_integral_value(rounding=ROUND_HALF_EVEN)
        if lo!=hi: raise AssertionError('rounding enclosure straddles a half-ulp')
        return int(lo)


def check_table(function,bits,p,radius=4):
    n=1<<bits
    lo,hi={'exp':(-radius,radius),'log':(1,2),'normal_cdf':(-8,8)}[function]
    rows=[]
    with localcontext() as ctx:
        ctx.prec=100
        for j in range(n+1):
            x=Decimal(lo)+Decimal(hi-lo)*j/n
            if function=='exp': value,error=x.exp(),EPS
            elif function=='log': value,error=x.ln(),EPS
            else: value,error=decimal_cdf(x)
            exact_round=rounded_enclosure(value,error,p)
            xf=float(x)
            reference={'exp':math.exp,'log':math.log,
                       'normal_cdf':lambda t:.5*math.erfc(-t/math.sqrt(2))}[function](xf)
            libm_round=int(round(reference*(1<<p)))
            if exact_round!=libm_round:
                raise AssertionError((function,j,exact_round,libm_round))
    return {'function':function,'segments_bits':bits,'fraction_bits':p,
            'domain':[lo,hi],
            'endpoint_count':n+1,'all_libm_rounded_constants_match_enclosed_rounding':True,
            'absolute_decimal_enclosure':str(EPS)}


def check_quantiles(bits,p):
    from statistics import NormalDist
    nd=NormalDist(); n=1<<bits; q=Decimal(1)/Decimal(1<<p)
    with localcontext() as ctx:
        ctx.prec=100
        for j in range(n):
            target=(Decimal(j)+HALF)/n
            zhat=Decimal(int(round(nd.inv_cdf(float(target))*(1<<p))))/(1<<p)
            value,error=decimal_cdf(zhat)
            bound_z=abs(zhat)+q/2
            derivative_lower=PHI0*(-bound_z*bound_z/2).exp()-EPS
            assert abs(value-target)+error <= q/2*derivative_lower
    return {'quantile_bits':bits,'coordinate_fraction_bits':p,'entry_count':n,
            'each_rounded_coordinate_within_half_ulp_of_true_quantile':True}


def geometric_correction(cap):
    """Independent positive tail-interval expression, offline coefficient."""
    with localcontext() as ctx:
        ctx.prec=100
        mu=Decimal('-.0225'); variance=Decimal('.0425'); sd=variance.sqrt()
        h=Decimal(cap)/100
        x0=mu/sd; x1=(mu-(1+h).ln())/sd
        values=[decimal_cdf(t)[0] for t in [x0,x1,x0+sd,x1+sd]]
        p0,p1,a0,a1=values
        currency=100*((mu+variance/2).exp()*(a0-a1)-(p0-p1)+h*p1)
        # CDF errors <1e-50 each; exp/log/sqrt correctly rounded at100
        # significant digits. All coefficients and arguments have magnitude
        # <101. The monetary envelope is deliberately enlarged to1e-44.
        error=Decimal('1e-44')
        value_float=float(currency)
        float_error=abs(Decimal.from_float(value_float)-currency)+error
        return {'cap':cap,'undiscounted_geometric_correction_decimal':str(currency),
                'absolute_monetary_enclosure':str(error),
                'correction_float':value_float,'float_representation_error_enclosure':str(float_error)}



def main():
    import hashlib
    output={'status':'independent offline geometric-correction enclosures',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'continuous_geometric_corrections':[geometric_correction(h) for h in [1,10,25]],
            'absolute_series_arithmetic_envelope':str(EPS)}
    path=Path(__file__).resolve().parents[1]/'results'/'geometric_constants.json'
    path.write_text(json.dumps(output,indent=2)+'\n');print(path)

if __name__=='__main__':main()
