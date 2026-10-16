"""Declared wide-domain kernels for the prespecified synthetic study.

This module provides constants and certificate data, not free quantum oracles.
All table access/interpolation gates are emitted by the existing engine.
"""
from decimal import Decimal, localcontext
import math
from encoding_resources import decimal_cdf, PI, INV_SQRT_2PI, rounded_scaled

D=Decimal
CDF_ALLOWANCE=D('6.4e-16')
CAP_ALLOWANCE=D('2e-12')
DECIMAL_ALLOWANCE=D('1e-38')


def function_span(kind):
    if kind in ('general_exp','general_log'):return 8
    if kind.startswith('general_strike:'):return 1
    if kind.startswith(('general_cdf:','general_cap:')):return 4
    raise ValueError(kind)


def clipped_cdf(x):
    # Phi(-8)<phi(8)/8<6.4e-16. Saturation is explicitly charged.
    if x<=-8:return D(0)
    if x>=8:return D(1)
    return decimal_cdf(x)


def conditional_centre(m,lam,cap):
    with localcontext() as ctx:
        ctx.prec=70
        h=D(cap)/100
        if not D('.05')<=lam<=D('.4') or not D('.01')<=h<=D('.25'):
            raise ValueError('Outside declared kernel certificate family')
        a=(1+h).ln()
        t0=m/lam;t1=(m-a)/lam
        p0,p1,a0,a1=[clipped_cdf(t) for t in (t0,t1,t0+lam,t1+lam)]
        value=((m+lam*lam/2).exp()*(a0-a1)-(p0-p1)+h*p1)/h
        # The true cap kernel is in [0,Phi(m/lambda)]. Projection to the
        # approximate upper endpoint increases error by at most its allowance.
        return min(max(value,D(0)),clipped_cdf(m/lam))


def value_centre(kind,x):
    if kind=='general_exp':return (x-4).exp()
    if kind=='general_log':return (x+D('.0625')).ln()+4
    parts=kind.split(':')
    if parts[0]=='general_strike':return (1+D(parts[1])/100*x).ln()
    lam=D(parts[1])
    if parts[0]=='general_cdf':return clipped_cdf((x-2)/lam)
    if parts[0]=='general_cap':return conditional_centre(x-2,lam,int(parts[2]))
    raise ValueError(kind)


def probability_constants(lam,h=None,cap_aware=False):
    with localcontext() as ctx:
        ctx.prec=70
        lip=INV_SQRT_2PI/lam
        curvature=INV_SQRT_2PI/D(1).exp().sqrt()/(lam*lam)
        if cap_aware:
            lip=min(lip,(1+h)/h)
            curvature=min(curvature,(1+h)*INV_SQRT_2PI/(h*lam))
        return lip+DECIMAL_ALLOWANCE,curvature+DECIMAL_ALLOWANCE


def bound_specification(kind,a,spacing,mode):
    with localcontext() as ctx:
        ctx.prec=70
        lo=value_centre(kind,a.lo);hi=value_centre(kind,a.hi)
        if kind=='general_exp':
            maximum=(min(D(8),a.hi+a.error+spacing)-4).exp()+DECIMAL_ALLOWANCE
            return lo,hi,maximum,maximum,DECIMAL_ALLOWANCE
        if kind=='general_log':
            minimum=max(D('.0625'),a.lo-a.error-spacing+D('.0625'))
            return lo,hi,1/(a.lo-a.error+D('.0625')),1/minimum**2,DECIMAL_ALLOWANCE
        parts=kind.split(':')
        if parts[0]=='general_strike':
            h=D(parts[1])/100
            return lo,hi,h+DECIMAL_ALLOWANCE,h*h+DECIMAL_ALLOWANCE,DECIMAL_ALLOWANCE
        lam=D(parts[1]);is_cap=parts[0]=='general_cap'
        h=D(parts[2])/100 if is_cap else None
        lip,curvature=probability_constants(lam,h,is_cap and mode=='cap_aware')
        return lo,hi,lip,curvature,CAP_ALLOWANCE if is_cap else CDF_ALLOWANCE


def function_table(kind,bits,q,width):
    n=1<<bits;span=D(function_span(kind))
    with localcontext() as ctx:
        ctx.prec=70
        values=[];previous=D(0)
        for j in range(n+1):
            x=span*j/n
            value=value_centre(kind,x)
            # Cumulative maximum preserves a uniform ordinate enclosure for
            # an increasing exact function. Both probability kernels remain
            # ordered after this operation and common nearest rounding.
            if kind.startswith(('general_cap:','general_cdf:')):
                value=max(previous,value);previous=value
            values.append(rounded_scaled(value,q))
    assert all(0<=v<1<<width for v in values)
    assert all(a<=b for a,b in zip(values,values[1:]))
    return [values[j]|((values[j+1]-values[j])<<width) for j in range(n)]
