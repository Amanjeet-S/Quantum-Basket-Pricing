"""Sharp continuous-support clipped-tail moment certificate.

This is a classical moment-option bound, not a new generic moment theory.
It provides a rigorous finite-oracle clipping bias bound when the oracle's
range and first/second moment bounds have been proved separately. The sharp
attaining law need not lie on the oracle's finite grid.
"""
from fractions import Fraction as F
import math


def sqrt_interval(value,bits=100):
    value=F(value)
    if value<0:raise ValueError('Negative square-root input')
    n,d=math.isqrt(value.numerator),math.isqrt(value.denominator)
    if n*n==value.numerator and d*d==value.denominator:
        exact=F(n,d);return exact,exact
    scale=1 << bits
    scaled=value*scale*scale
    root=math.isqrt(scaled.numerator//scaled.denominator)
    lo=F(root,scale);hi=F(root+1,scale)
    assert lo*lo<=value<=hi*hi
    return lo,hi


def basic_bounds(rho,u,v,r):
    """Established range, first-moment and second-moment upper bounds."""
    rho,u,v,r=map(F,(rho,u,v,r))
    if min(rho,u,v,r)<0:raise ValueError('Inputs must be nonnegative')
    if rho<=r:return {'range':F(0),'first_moment':F(0),'second_moment':F(0),'minimum':F(0)}
    result={'range':rho-r,'first_moment':u*(1-r/rho)}
    if r==0:
        result['second_moment']=sqrt_interval(v)[1]
    else:
        ratio=(rho-r)/(rho*rho) if rho<=2*r else 1/(4*r)
        result['second_moment']=v*ratio
    result['minimum']=min(result.values())
    return result


def clipped_tail_bound(rho,u,v,r):
    """Certify sup E[(X-r)+] under 0<=X<=rho, EX<=u, EX²<=v.

    All branches except an irrational constant-support root are exact
    rational formulas. That branch returns a rational enclosure of width
    at most 2**-100. 'upper' is the clipping-bias allowance to be charged.
    The dual polynomial is a pointwise majorant whenever it is present.
    """
    rho,u_raw,v_raw,r=map(F,(rho,u,v,r))
    if min(rho,u_raw,v_raw,r)<0:raise ValueError('Inputs must be nonnegative')
    u=min(u_raw,rho);v=min(v_raw,rho*u,rho*rho)
    record={'range_upper':rho,'clip_level':r,'input_mean_upper':u_raw,
            'input_second_moment_upper':v_raw,'effective_mean_upper':u,
            'effective_second_moment_upper':v,
            'scope':'sharp for the continuous support class [0,rho]; conservative for a prescribed finite grid'}
    def finish(case,lower,upper,dual,witness):
        record.update(case=case,lower=F(lower),upper=F(upper),
            dual_quadratic=dual,attaining_law=witness)
        return record
    if rho<=r or u==0 or v==0:
        if rho<=r:dual=(F(0),F(0),F(0))
        elif u==0:dual=(F(0),F(1),F(0))
        elif r>0:dual=(F(0),F(0),1/(4*r))
        else:dual=None  # EX²=0 implies X=0; no attained finite quadratic dual is claimed.
        return finish('zero',0,0,dual,((F(0),F(1)),))
    if r==0 and u*u<=v:
        return finish('mean_only',u,u,(F(0),F(1),F(0)),((u,F(1)),))
    if r==0 or (v>4*r*r and u*u>=v):
        lo,hi=sqrt_interval(v)
        assert hi>=2*r
        upper=hi-r;lower=lo-r
        dual=(hi/2-r,F(0),1/(2*hi))
        witness=((lo,F(1)),) if lo==hi else {'constant':'sqrt(effective_second_moment_upper)'}
        return finish('constant_root',lower,upper,dual,witness)
    if rho<2*r:
        result=(rho-r)*v/(rho*rho)
        p=v/(rho*rho)
        return finish('range_second_moment',result,result,
            (F(0),F(0),(rho-r)/(rho*rho)),((F(0),1-p),(rho,p)))
    if v<=4*r*r and u>=v/(2*r):
        result=v/(4*r);p=v/(4*r*r)
        return finish('second_moment_tangent',result,result,
            (F(0),F(0),1/(4*r)),((F(0),1-p),(2*r,p)))
    h=v/u;p=u*u/v
    assert 2*r<=h<=rho and 0<=p<=1
    result=u-r*u*u/v
    dual=(F(0),1-2*r/h,r/(h*h))
    return finish('mean_second_moment',result,result,dual,((F(0),1-p),(h,p)))


def polynomial(record,x):
    if record['dual_quadratic'] is None:raise ValueError('Degenerate primal inference has no finite quadratic dual record')
    a,b,c=record['dual_quadratic'];x=F(x)
    return a+b*x+c*x*x


def rational_record(value):
    if isinstance(value,F):return {'numerator':str(value.numerator),'denominator':str(value.denominator),'float':float(value)}
    if isinstance(value,dict):return {str(k):rational_record(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [rational_record(v) for v in value]
    return value
