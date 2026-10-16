"""Payoff-specific scalar certificates without modifying retained integral sources.

The scalar graph uses the classical risk-neutral Gaussian conditional price.
Only the arithmetic implementing its remainder is prepared coherently.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
from flint import arb,ctx
import integral_value_certificate as old
from boundary_derivative_certificate import payoff_norms
from continuous_quadrature import frac_arb
from taylor_pricing import normal_cdf,normal_pdf

ROOT=Path(__file__).resolve().parent
PRECISION=old.PRECISION
endpoints=old.endpoints
upper=old.upper


def checked_norms(params,vanilla=False):
    p={k:F(v) for k,v in (params or {}).items()}
    lam=p['lambda'];K=p.get('K',F(100));cap=None if vanilla else p['cap']
    if lam<=0 or K!=100 or cap is not None and cap<=0:
        raise ValueError('Positive scale, strike100 and positive optional cap required')
    norm=payoff_norms(lam,cap,K)
    if not norm['passed']:raise ValueError('Unaccepted payoff norm certificate')
    return p,norm


def scalar_ball(name,x,params=None):
    if name not in ('boundary_h','boundary_density','vanilla_call','vanilla_call_prime'):
        return old.scalar_ball(name,x,params)
    p,norm=checked_norms(params,name!='boundary_h');lam=frac_arb(p['lambda']);K=frac_arb(p.get('K',F(100)))
    if name in ('boundary_h','boundary_density'):
        M2=p['M2']
        if M2<F(norm['M2_upper']):raise ValueError('Uncertified cash-Gamma normaliser')
        C=K/lam*normal_pdf(x/lam)
        if name=='boundary_density':return C/frac_arb(M2)
        H=frac_arb(p['cap']);c=1+H/K;ell=c.log()
        C-=K/lam*c*normal_pdf((x-ell)/lam)
        return (1+C/frac_arb(M2))/2
    forward=K*(x+lam*lam/2).exp()*normal_cdf((x+lam*lam)/lam)
    return forward if name=='vanilla_call_prime' else forward-K*normal_cdf(x/lam)


def scalar_bounds(name,lo,hi,params=None):
    ctx.prec=PRECISION;lo,hi=map(F,(lo,hi))
    if lo>=hi:raise ValueError('Positive scalar interval required')
    if name not in ('boundary_h','boundary_density','vanilla_call','vanilla_call_prime','cap','cap_prime'):
        return old.scalar_bounds(name,lo,hi,params)
    vanilla=name in ('boundary_density','vanilla_call','vanilla_call_prime')
    p,norm=checked_norms(params,vanilla);lam=p['lambda'];K=p.get('K',F(100))
    M2=F(norm['M2_upper']);M3=F(norm['C_r_upper']['3']);M4=F(norm['C_r_upper']['4'])
    if name in ('boundary_h','boundary_density'):
        denom=p['M2']
        if denom<M2:raise ValueError('Uncertified cash-Gamma normaliser')
        if name=='boundary_density':
            return {'range':(F(0),F(1)),
                    'first':K*old.hermite_norm(1)/(lam**2*denom),
                    'second':K*old.hermite_norm(2)/(lam**3*denom)}
        H=p['cap'];direct1=(2*K+H)*old.hermite_norm(1)/lam**2
        direct2=(2*K+H)*old.hermite_norm(2)/lam**3
        return {'range':(F(0),F(1)),
                'first':min(direct1,M3+2*M2)/(2*denom),
                'second':min(direct2,M4+5*M3+4*M2)/(2*denom)}
    if vanilla:
        # P' <= the finite-domain conditional forward, not a global bound.
        M1=upper(frac_arb(K)*(frac_arb(hi)+frac_arb(lam*lam/2)).exp())
        if name=='vanilla_call':return {'range':(F(0),M1),'first':M1,'second':M1+M2}
        return {'range':(F(0),M1),'first':M1+M2,'second':M1+3*M2+M3}
    H=p['cap'];M1=F(norm['M1_upper']);scale=1/H if bool(p.get('normalised',False)) else F(1)
    prior=old.scalar_bounds(name,lo,hi,params)
    if name=='cap':return {'range':(F(0),H*scale),'first':min(prior['first'],M1*scale),
                           'second':min(prior['second'],(M1+M2)*scale)}
    return {'range':(F(0),M1*scale),'first':min(prior['first'],(M1+M2)*scale),
            'second':min(prior['second'],(M1+3*M2+M3)*scale)}


def scalar_certificate(name,lo,span,params,q,table_bits,output_offset=0):
    """Actual N+1 unsigned ordinates and signed-magnitude interpolation.

    The interpolation formula remains base +/- floor(abs(delta)*fraction).
    Every decreasing cell pays its output floor error independently.
    """
    ctx.prec=PRECISION;lo,span,offset=map(F,(lo,span,output_offset))
    if not isinstance(q,int) or not isinstance(table_bits,int) or min(q,table_bits)<1:
        raise ValueError('Positive integer precisions required')
    fractional=q+old.power_two_exponent(span)-table_bits
    if fractional<0:raise ValueError('Table finer than actual input register')
    N=1<<table_bits;scale=1<<q;contract=scalar_bounds(name,lo,lo+span,params)
    if contract['range'][0]<offset:raise ValueError('Scalar output must remain unsigned')
    nodes=[];nodal=F(0)
    for j in range(N+1):
        x=lo+span*F(j,N);a=scalar_ball(name,frac_arb(x),params)-frac_arb(offset)
        al,au=endpoints(a);centre=(al+au)/2
        value=max(0,(centre*scale+F(1,2)).numerator//(centre*scale+F(1,2)).denominator)
        error=max(abs(F(value,scale)-al),abs(F(value,scale)-au));nodal=max(nodal,error)
        nodes.append({'input':str(x),'lower':str(al),'upper':str(au),'integer':value,'error_upper':str(error)})
    interpolation=contract['second']*(span/N)**2/8;floor=F(1,scale)
    return {'schema':'qbp-boundary-scalar-v1','name':name,'lo':str(lo),'span':str(span),
        'params':{k:str(F(v)) for k,v in (params or {}).items()},'q':q,'table_bits':table_bits,
        'fractional_remainder_bits':fractional,'output_offset':str(offset),
        'exact_function_range':list(map(str,contract['range'])),
        'first_derivative_upper':str(contract['first']),'second_derivative_upper':str(contract['second']),
        'nodes':nodes,'nodal_error_upper':str(nodal),'linear_interpolation_error_upper':str(interpolation),
        'interpolation_floor_error_upper':str(floor),'output_error_upper':str(nodal+interpolation+floor),
        'value_min_integer':min(n['integer'] for n in nodes),'value_max_integer':max(n['integer'] for n in nodes),
        'integer_width':max(1,max(n['integer'] for n in nodes).bit_length()),'precision_bits':PRECISION,
        'input_perturbation_rule':'Add declared first derivative bound times absolute input error; input clipping preserves the domain.',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'norm_source_sha256':hashlib.sha256((ROOT/'results/boundary_derivative_norms.json').read_bytes()).hexdigest()}


def gauss_constants(q):
    ctx.prec=PRECISION
    root=arb(6).sqrt();zm=(4-root)/10;zp=(4+root)/10;wm=arb(1)/2+1/(3*root)
    def rounded(a):
        lo,hi=endpoints(a);x=((lo+hi)/2)*(1<<q)+F(1,2);return F(x.numerator//x.denominator,1<<q)
    zmq,zpq,wmq=map(rounded,(zm,zp,wm));wpq=1-wmq
    if not 0<=zmq<=zpq<=1 or min(wmq,wpq)<=0:raise ValueError('Insufficient Gauss fixed precision')
    def record(a,x):
        lo,hi=endpoints(a);return {'ideal_lower':str(lo),'ideal_upper':str(hi),'actual':str(x),
                                  'error_upper':str(max(abs(x-lo),abs(x-hi)))}
    return {'q':q,'z_minus':record(zm,zmq),'z_plus':record(zp,zpq),
            'w_minus':record(wm,wmq),'w_plus':record(1-wm,wpq),
            'actual_weights_sum':'1','constant_encoding_error_charged_separately':True}


certify_graph=old.certify_graph
quantile_certificate=old.quantile_certificate
finite_input_domain=old.finite_input_domain
finite_grid_moments=old.finite_grid_moments
