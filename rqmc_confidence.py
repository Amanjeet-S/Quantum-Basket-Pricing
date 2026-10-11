"""Rigorous finite-checkpoint confidence from independent random-shift batches.

The bound treats independent shifts as the sampling units. It assumes no
variance reduction rate in net size and uses no Student approximation. The
pair/Chernoff/Bernstein argument is self-contained; generic empirical
Bernstein RQMC inference has published antecedents and is not claimed new.
"""
from decimal import Decimal,localcontext
from fractions import Fraction as F
import math
import re

LOG_ALLOWANCE=Decimal('1e-120')
ROOT_BITS=160

def ceil_f(x):x=F(x);return -((-x.numerator)//x.denominator)
def sqrt_upper(x):
    x=F(x)
    if x<0:raise ValueError('Negative square-root argument')
    scale=1 << ROOT_BITS;n=ceil_f(x*scale*scale);root=math.isqrt(n)
    if root*root<n:root+=1
    out=F(root,scale);assert out*out>=x
    return out

def log_upper(x):
    x=F(x)
    if not 1<x<=10**40:raise ValueError('Declared log argument range is1<x<=1e40')
    with localcontext() as c:
        c.prec=160;centre=(Decimal(x.numerator)/Decimal(x.denominator)).ln()
        return F(centre+LOG_ALLOWANCE)

def exact_counter(value):
    if isinstance(value,bool):raise ValueError('Boolean counters are not integer statistics')
    if isinstance(value,int):return value
    if isinstance(value,str) and re.fullmatch(r'[0-9]+',value):return int(value)
    raise ValueError('Counters must be exact integers or digit strings')

def confidence_from_batch_statistics(stats,*,rho,alpha_variance,alpha_mean,
                                     variance_prior=None,mean_prior=None):
    rho=F(rho);av,am=F(alpha_variance),F(alpha_mean)
    if not 0<=rho<=1 or not 0<av<1 or not 0<am<1:raise ValueError('Invalid range/failure budgets')
    R=exact_counter(stats['replicates']);D=exact_counter(stats['batch_denominator']);S=exact_counter(stats['integer_batch_sum'])
    S2=exact_counter(stats['integer_batch_square_sum']);m=exact_counter(stats['independent_pairs']);Q=exact_counter(stats['integer_paired_difference_square_sum'])
    if R<2 or D<=0 or m!=R//2 or min(S,S2,Q)<0:raise ValueError('Invalid exact batch statistics')
    if not 0<=F(S,R*D)<=rho or S*S>R*S2 or F(S2,R*D*D)>rho*rho or F(Q,m*D*D)>rho*rho or F(S2)>rho*D*S or Q>S2:
        raise ValueError('Statistics violate the certified range')
    prior_var=rho*rho/4 if variance_prior is None else F(variance_prior)
    prior_mean=rho if mean_prior is None else F(mean_prior)
    if not 0<=prior_var<=rho*rho/4 or not 0<=prior_mean<=rho:raise ValueError('Invalid proved priors')
    pair=F(Q,2*m*D*D);M=rho*rho/2
    t=log_upper(1/av)/m
    # Lower-tail multiplicative Chernoff inversion for positive pair values.
    U_pair=pair+M*t+sqrt_upper(M*M*t*t+2*M*t*pair)
    U=min(prior_var,U_pair)
    ell=log_upper(2/am);b=rho*ell/(3*R)
    width=b+sqrt_upper(b*b+2*U*ell/R)
    sample=F(S,R*D);lo=max(F(0),sample-width);hi=min(prior_mean,sample+width)
    common={'replicates':R,'pairs':m,'sample_mean_fraction':str(sample),
        'paired_variance_observation_fraction':str(pair),'pair_variance_upper_fraction':str(U_pair),
        'variance_upper_fraction':str(U),'variance_prior_upper_fraction':str(prior_var),
        'unclipped_mean_halfwidth_fraction':str(width),'alpha_variance_fraction':str(av),
        'alpha_mean_fraction':str(am),'failure_budget_fraction':str(av+am),
        'log_variance_upper_fraction':str(log_upper(1/av)),
        'log_mean_upper_fraction':str(ell),'square_root_bits':ROOT_BITS,
        'confidence_premise':'independent ideal uniform digital shifts; within-batch points may be dependent'}
    if lo>hi:
        return {**common,'status':'empty_intersection','mean_lower_fraction':None,
                'mean_upper_fraction':None,'halfwidth_upper_fraction':None,
                'interpretation':'explicit confidence-failure outcome, included in the failure budget'}
    return {**common,'status':'nonempty','mean_lower_fraction':str(lo),
            'mean_upper_fraction':str(hi),'halfwidth_upper_fraction':str((hi-lo)/2),
            'estimate_fraction':str((hi+lo)/2)}
