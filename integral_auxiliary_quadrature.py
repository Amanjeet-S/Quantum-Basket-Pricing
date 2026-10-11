"""Quadratic certificates for the unchanged Beta2 and strike registers.

All arithmetic is exact rational after reading independently certified
Gaussian norms. Finite financial moments belong to the actual outer grid;
fixed-point oracle error is a separate charge.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent


@lru_cache(None)
def norm_data():
    path=ROOT/'results/gaussian_derivative_norms.json';data=json.loads(path.read_text())
    if data['status']!='passed':raise ValueError('Certified norm inputs required')
    return data


def joint_norm(order,lam):
    row=next(r for r in norm_data()['records'] if r['derivative_order']==order and F(r['lambda'])==F(lam))
    if not row['passed']:raise ValueError('Rejected Gaussian norm')
    return F(row['norm_upper'])


def hermite_norm(order):
    row=next(r for r in norm_data()['unit_norms'] if r['unit_derivative_order']==order)
    if not row['passed']:raise ValueError('Rejected Hermite density norm')
    return F(row['norm_upper'])


def auxiliary_bias(model,cap,beta_bits,strike_bits,U_star,moment2_upper,
                   moment3_upper,moment4_upper=None,K=100):
    """Both old coupling caps and stronger quadratic caps, transparently.

    model supplies the classical common Gaussian scale. U_star must bound
    the ideal U on the same actual finite factor grid as the moment inputs.
    A uniform-box likelihood, if used, is already inside those moments.
    An external probability mass or volume factor remains outside this API.
    """
    H,lam,U,m2,m3,K=map(F,(cap,model['lambda'],U_star,moment2_upper,moment3_upper,K))
    if min(H,lam,K)<=0 or min(U,m2,m3)<0:raise ValueError('Positive financial scales and nonnegative support/moments required')
    if any(isinstance(b,bool) or not isinstance(b,int) or b<1 for b in [beta_bits,strike_bits]):raise ValueError('Positive integer midpoint widths required')
    support_m4=U*m3
    if moment4_upper is None:m4=support_m4;source='certified finite-support product U_star times moment3_upper'
    else:
        supplied=F(moment4_upper)
        if supplied<0:raise ValueError('Nonnegative fourth moment required')
        m4=min(supplied,support_m4);source='minimum of supplied fourth-moment upper bound and finite-support product'
    E=hermite_norm(1)/lam+hermite_norm(2)/lam**2
    Fsecond=hermite_norm(2)/lam+hermite_norm(3)/lam**2
    Dz3,Dz4=joint_norm(3,lam),joint_norm(4,lam)
    hZ,hS=F(1,1<<beta_bits),F(1,1<<strike_bits)
    beta_linear=H*Dz3*m3*hZ/4
    strike_linear=H*H*E*m2*hS/(4*K*lam)
    beta_quadratic=H*hZ*hZ*(Dz3*m3/12+Dz4*m4/48)
    strike_quadratic=H**3*m2*hS*hS*(Fsecond/lam**2+E/lam)/(48*K*K)
    if U==0:beta_linear=strike_linear=beta_quadratic=strike_quadratic=F(0)
    beta,strike=min(beta_linear,beta_quadratic),min(strike_linear,strike_quadratic)
    return {'schema':'qbp-quadratic-auxiliary-bias-v1','model_id':model.get('model_id'),'cap':str(H),'lambda':str(lam),
       'beta_bits':beta_bits,'strike_bits':strike_bits,'U_star':str(U),'K':str(K),
       'moment2_upper':str(m2),'moment3_upper':str(m3),'fourth_moment_upper_used':str(m4),'fourth_moment_source':source,
       'E2_upper':str(E),'F2_upper':str(Fsecond),'D3_upper':str(Dz3),'D4_upper':str(Dz4),
       'linear_beta_bias_upper':str(beta_linear),'linear_strike_bias_upper':str(strike_linear),
       'linear_total_auxiliary_bias_upper':str(beta_linear+strike_linear),
       'quadratic_beta_bias_upper':str(beta_quadratic),'quadratic_strike_bias_upper':str(strike_quadratic),
       'quadratic_total_auxiliary_bias_upper':str(beta_quadratic+strike_quadratic),
       'beta_minimum_midpoint_bias_upper':str(beta),'strike_midpoint_bias_upper':str(strike),
       'total_auxiliary_bias_upper':str(beta+strike),
       'moment_law':'Ideal financial U on the actual finite outer grid; relative likelihood weight included for uniform_box. Oracle arithmetic is charged separately.',
       'normalisation_scope':'Undiscounted monetary error before the external retained Gaussian mass or uniform-box volume factor.',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'norm_input_sha256':hashlib.sha256((ROOT/'results/gaussian_derivative_norms.json').read_bytes()).hexdigest(),
       'interpretation':'Standard weighted midpoint error applied to the emitted auxiliary laws; no change to the circuit or novelty claim for midpoint quadrature.'}


quadratic_auxiliary_bias=auxiliary_bias
