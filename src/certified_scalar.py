"""Outward rational scalar machinery retained from the audited controller."""
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
import math
F=Fraction
INVERSE_TOLERANCE=F(1,1<<70)

def ceil_fraction(value):
    return -((-value.numerator)//value.denominator)


def floor_fraction(value):
    return value.numerator//value.denominator


@lru_cache(None)
def pi_rational_interval(terms=80):
    def atan_interval(den):
        x=F(1,den);total=F(0)
        for j in range(terms):
            total+=(-1 if j%2 else 1)*x**(2*j+1)/(2*j+1)
        remainder=x**(2*terms+1)/(2*terms+1)
        return total-remainder,total+remainder
    al,au=atan_interval(5);bl,bu=atan_interval(239)
    return 16*al-4*bu,16*au-4*bl


def decimal_fraction(value, rounding):
    with localcontext() as ctx:
        ctx.prec=80;ctx.rounding=rounding
        return Decimal(value.numerator)/Decimal(value.denominator)


@lru_cache(maxsize=20000)
def probability_interval(u):
    """Enclose f(u)=sin²(pi*u/2) for rational u in [0,1].

    Directed Decimal interval arithmetic encloses the cosine Taylor polynomial.
    Its decreasing alternating tail is enclosed by the next positive term.
    Pi itself has an exact rational Machin-series enclosure.
    """
    u=F(u)
    if not 0<=u<=1:raise ValueError('Probability input must lie in [0,1]')
    if u==0:return F(0),F(0)
    if u==1:return F(1),F(1)
    if u==F(1,2):return F(1,2),F(1,2)
    pl,pu=pi_rational_interval()
    with localcontext() as ctx:
        ctx.prec=80;ctx.rounding=ROUND_FLOOR
        xl=decimal_fraction(pl,ROUND_FLOOR)*decimal_fraction(u,ROUND_FLOOR)
        x2l=xl*xl
        ctx.rounding=ROUND_CEILING
        xu=decimal_fraction(pu,ROUND_CEILING)*decimal_fraction(u,ROUND_CEILING)
        x2u=xu*xu
        tl=tu=Decimal(1);sl=su=Decimal(1)
        for j in range(1,45):
            den=Decimal((2*j-1)*(2*j))
            ctx.rounding=ROUND_FLOOR;tl=tl*x2l/den
            ctx.rounding=ROUND_CEILING;tu=tu*x2u/den
            if j%2:
                ctx.rounding=ROUND_FLOOR;sl=sl-tu
                ctx.rounding=ROUND_CEILING;su=su-tl
            else:
                ctx.rounding=ROUND_FLOOR;sl=sl+tl
                ctx.rounding=ROUND_CEILING;su=su+tu
        # At j=45 the remaining terms decrease because |x|<pi<4.
        ctx.rounding=ROUND_CEILING
        tail=tu*x2u/Decimal(89*90)
        ctx.rounding=ROUND_FLOOR;cl=sl-tail
        ctx.rounding=ROUND_CEILING;cu=su+tail
        ctx.rounding=ROUND_FLOOR;lo=(Decimal(1)-cu)/2
        ctx.rounding=ROUND_CEILING;hi=(Decimal(1)-cl)/2
    return max(F(0),F(lo)),min(F(1),F(hi))


@lru_cache(maxsize=20000)
def inverse_probability_interval(p):
    """Enclose the inverse of increasing sin²(pi*u/2) on [0,1]."""
    p=F(p)
    if not 0<=p<=1:raise ValueError('Probability must lie in [0,1]')
    if p in (F(0),F(1),F(1,2)):return p,p
    lo,hi=F(0),F(1);eta=INVERSE_TOLERANCE
    while hi-lo>eta:
        mid=(lo+hi)/2;fl,fu=probability_interval(mid)
        if fu<p:lo=mid
        elif fl>p:hi=mid
        else:
            # |f(x)-f(y)| >= sin²(pi*|x-y|/2) >= |x-y|².
            # Therefore a probability enclosure of width <=eta² implies
            # inverse coordinate error <=eta, without dividing by a slope.
            if fu-fl>eta*eta:raise ArithmeticError('Cosine enclosure too wide')
            return max(lo,mid-eta),min(hi,mid+eta)
    return lo,hi


def amplitude_interval(lo,hi):
    return probability_interval(lo)[0],probability_interval(hi)[1]


def sqrt_upper_fraction(value,bits=100):
    """An exact rational upper square-root, without a libm assumption."""
    value=F(value)
    if value<0:raise ValueError('Negative square-root input')
    scale=1<<bits;integer=ceil_fraction(value*scale*scale)
    root=math.isqrt(integer)
    if root*root<integer:root+=1
    return F(root,scale)

