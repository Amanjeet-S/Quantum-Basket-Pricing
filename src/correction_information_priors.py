"""Coefficient-informed priors for a fully informed first-correction opponent.

Only c1, the continuous second ratio moment, a global curvature bound and
an independently constructed financial transfer enter. No basket-price
answer, Taylor price interval, empirical option quote or simulator estimate
is an input. Every previous upper prior is retained in the intersection.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial,prod
from flint import arb,ctx
import hashlib,json
from boundary_grid_certificate import parameters,ball,serial,upper,endpoints
from boundary_derivative_certificate import payoff_norms
from taylor_pricing import compositions
ROOT=Path(__file__).resolve().parent

def up(iv):return F(iv['upper_exact_fraction'])
def low(iv):return F(iv['lower_exact_fraction'])

def continuous_second_moment(model):
    ctx.prec=512;w,d,s,mg,rg,lam=parameters(model,False);raw=[];terms=[]
    for order in range(3):
        total=arb(0)
        for alpha in compositions(order,len(w)):
            coefficient=F(factorial(order),prod(factorial(a) for a in alpha))*prod(wi**a for wi,a in zip(w,alpha))
            v=tuple(sum(alpha[i]*s[i][j] for i in range(len(w))) for j in range(len(rg)))
            constant=sum(alpha[i]*d[i] for i in range(len(w)))
            total+=ball(coefficient)*ball(constant+sum(x*x for x in v)/2).exp()
            terms.append({'order':order,'alpha':list(alpha),'coefficient':str(coefficient),'constant':str(constant),'rate':list(map(str,v))})
        raw.append(total)
    value=raw[2]-2*raw[1]+raw[0];lo,hi=endpoints(value)
    if hi<0:raise ArithmeticError('Nonnegative continuous ratio moment has a negative upper enclosure')
    return {'model_id':model['model_id'],'continuous_U2_interval':{'lower_exact_fraction':str(max(F(0),lo)),'upper_exact_fraction':str(hi)},
        'multinomial_terms':terms,'exponential_evaluations':len(terms),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Full Gaussian moment from frozen parameters. No finite grid, pricing answer or prior classical result is imported.'}

def first_probability_interval(model,cap,coefficient,financial,certificate,previous,moment):
    B=F(certificate['configuration']['monetary_signal_normaliser']);eta=F(certificate['probability_error_upper'])
    if B<=0:raise ValueError('Positive nondegenerate normaliser required')
    norm=payoff_norms(model['lambda'],cap);remainder=F(norm['M2_upper'])*up(moment['continuous_U2_interval'])/2
    c1=coefficient['c1'];beta=F(financial['financial_transfer_upper']);mass=financial['mass_interval'];ml,mh=low(mass),up(mass)
    if ml<=0 or ml>mh:raise ValueError('Certified positive financial mass required')
    continuous_lower=low(c1)-(remainder if cap is not None else F(0));continuous_upper=up(c1)+remainder
    finite_lower=max(F(0),continuous_lower-beta)/mh
    finite_upper=min(B,(continuous_upper+beta)/ml,F(previous['ideal_mean_upper']))
    if finite_lower>finite_upper:raise ArithmeticError('Independent finite mean envelopes contradict each other')
    probability_lower=max(F(0),finite_lower/B-eta)
    probability_upper=min(F(1),F(previous['probability_upper']),finite_upper/B+eta)
    if probability_lower>probability_upper:raise ArithmeticError('Implemented probability envelopes contradict each other')
    return {'probability_lower':str(probability_lower),'probability_upper':str(probability_upper),'native_positive':True,
        'finite_ideal_mean_lower':str(finite_lower),'finite_ideal_mean_upper':str(finite_upper),'ideal_mean_upper':str(finite_upper),
        'previous_probability_upper':previous['probability_upper'],'previous_ideal_mean_upper':previous['ideal_mean_upper'],
        'continuous_c1_interval':c1,'continuous_R2_mean_absolute_upper':str(remainder),'continuous_R2_nonnegative':cap is None,
        'continuous_R1_mean_lower':str(continuous_lower),'continuous_R1_mean_upper':str(continuous_upper),
        'financial_transfer_upper':str(beta),'external_mass_interval':mass,'arithmetic_probability_enlargement':str(eta),
        'information_scope':'Closed Gaussian coefficient and global curvature-moment bound plus price-independent grid transfer. No financial-price answer.',
        'continuous_moment_record_sha256':hashlib.sha256(json.dumps(moment,sort_keys=True,separators=(',',':')).encode()).hexdigest()}


def informed_complete_plan(strategy,budget,certificate,prior,profiles):
    import boundary_precision_selection as retained
    det=F(budget['total_deterministic_undiscounted_price_error_upper']);remaining=F(1,100)-det
    if remaining<=0:return retained.complete(strategy,budget,certificate,prior,profiles)
    C=F(budget['external_factor_upper'])*F(certificate['configuration']['monetary_signal_normaliser'])
    lo,hi=F(prior['probability_lower']),F(prior['probability_upper']);radius=C*(hi-lo)/2
    if radius<=remaining:
        return {'kind':'deterministic_correction_prior','T_count':0,'quantum_shots':0,'certified_correction_radius':str(radius),
            'total_price_error_upper':str(det+radius),'implemented_probability_interval':{'lower':str(lo),'upper':str(hi)},
            'implemented_probability_midpoint':str((lo+hi)/2),
            'reason':'Shared coefficient-informed first-correction interval meets monetary tolerance. Analytical decision, no quantum gain.'}
    result=retained.complete(strategy,budget,certificate,prior,profiles)
    result['coefficient_informed_probability_interval']={'lower':str(lo),'upper':str(hi)}
    return result


def shared_continuous_screen(model,cap,coefficient,moment):
    norm=payoff_norms(model['lambda'],cap);R=F(norm['M2_upper'])*up(moment['continuous_U2_interval'])/2
    base=coefficient['coefficient'];lo=low(base)-(R if cap is not None else F(0));hi=up(base)+R
    lo=max(F(0),lo)
    if cap is not None:hi=min(F(cap),hi)
    if lo>hi:raise ArithmeticError('Continuous price certificate is empty')
    radius=(hi-lo)/2
    return {'model_id':model['model_id'],'payoff':'vanilla' if cap is None else 'capped','cap':cap,
        'price_interval':{'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)},
        'price_midpoint':str((lo+hi)/2),'price_radius_upper':str(radius),'monetary_radius':'1/100',
        'stops_without_quantum':radius<=F(1,100),'scope':'Shared low-order continuous curvature screen before finite-grid or estimator cost. It uses only c0+c1, M2 and the full Gaussian second ratio moment.'}
