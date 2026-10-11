"""Outward scalar and continuous-target certificates for integral pricing.

No main manuscript or existing result is changed. Generation requires the
existing python-flint runtime; retained rational metadata is portable.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import ceil, factorial, comb
import json, hashlib, time
from flint import arb, ctx
from continuous_quadrature import frac_arb, exact_endpoint, serialise_interval
from taylor_pricing import model_parameters, kernel_derivatives, normal_cdf, normal_pdf

ROOT=Path(__file__).resolve().parent
PRECISION=320
ENDPOINT_BITS=384

def dyadic_outward(x,upwards,bits=ENDPOINT_BITS):
    x=F(x);den=1<<bits
    scaled=x*den
    floor=scaled.numerator//scaled.denominator
    integer=floor+int(upwards and scaled.denominator!=1)
    return F(integer,den)

def endpoints(x):
    if not x.is_finite(): raise ArithmeticError("Nonfinite Arb enclosure")
    # Tiny Gaussian tails can have enormous exact dyadic denominators.
    # Coarsen outwards rather than disabling Python's integer-string guard.
    return dyadic_outward(exact_endpoint(x,False),False),dyadic_outward(exact_endpoint(x,True),True)

def upper(x):return endpoints(x)[1]
def interval(lo,hi):
    lo,hi=F(lo),F(hi)
    if lo>hi:raise ValueError("Reversed interval")
    mid=(lo+hi)/2;rad=(hi-lo)/2
    return frac_arb(mid)+arb(0,frac_arb(rad).upper())

def serial(x):return serialise_interval(*endpoints(x))

@lru_cache(None)
def norm_data():
    return json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())

def hermite_norm(j):
    rows=norm_data()['unit_norms']
    row=next(a for a in rows if a['unit_derivative_order']==j)
    if not row['passed']:raise ValueError("Uncertified Hermite norm")
    return F(row['norm_upper'])

def joint_norm(r,lam):
    row=next(a for a in norm_data()['records']
             if a['derivative_order']==r and F(a['lambda'])==F(lam))
    if not row['passed']:raise ValueError("Uncertified joint norm")
    return F(row['norm_upper'])

def joint_slope(r,lam):
    from gaussian_derivative_norm import signed_stirling
    lam=F(lam);st=signed_stirling(r)[r]
    return sum(F(abs(st[l]))/lam**l*hermite_norm(l) for l in range(1,r+1))

def joint_second(r,lam):
    from gaussian_derivative_norm import signed_stirling
    lam=F(lam);st=signed_stirling(r)[r]
    return sum(F(abs(st[l]))/lam**l*hermite_norm(l+1) for l in range(1,r+1))

def parameters(p):
    return {str(k):F(v) if k not in ('normalised',) else bool(v) for k,v in p.items()}

def scalar_ball(name,x,params=None):
    p=parameters(params or {})
    if name=='exp':return x.exp()
    if name=='log':return x.log()
    if name=='inv_square':return 1/(x*x)
    if name=='square':return x*x
    if name in ('h','h2'):
        lam=p['lambda'];D=p['D2']
        if lam<=0 or D<joint_norm(2,lam):
            raise ValueError("h requires positive scale and certified D2")
        value=1+(-1/frac_arb(lam)-x/(frac_arb(lam)**2))*normal_pdf(x)/frac_arb(D)
        return value/2 if name=='h' else value
    if name in ('cap','cap_prime'):
        cap=p['cap'];lam=p['lambda'];K=p.get('K',F(100))
        if cap<=0 or lam<=0 or K!=100:
            raise ValueError("Frozen cap kernel needs K=100,H>0,lambda>0")
        result=kernel_derivatives(x,frac_arb(lam),int(cap) if cap.denominator==1 else frac_arb(cap),
                                  1 if name=='cap_prime' else 0)[name=='cap_prime']
        return result/frac_arb(cap) if p.get('normalised',False) else result
    raise ValueError("Unsupported certified scalar "+name)

def scalar_bounds(name,lo,hi,params=None):
    """Bounds on the whole declared real interval; no sampled maxima."""
    ctx.prec=PRECISION;lo,hi=F(lo),F(hi);p=parameters(params or {})
    if lo>=hi:raise ValueError("Positive scalar interval width required")
    if name=='exp':
        a,b=frac_arb(lo).exp(),frac_arb(hi).exp()
        return {'range':(endpoints(a)[0],upper(b)), 'first':upper(b),'second':upper(b)}
    if name=='log':
        if lo<=0:raise ValueError("Logarithm domain must be strictly positive")
        return {'range':(endpoints(frac_arb(lo).log())[0],upper(frac_arb(hi).log())),
                'first':1/lo,'second':1/lo**2}
    if name=='inv_square':
        if lo<=0:raise ValueError("Inverse-square domain must be strictly positive")
        return {'range':(1/hi**2,1/lo**2),'first':2/lo**3,'second':6/lo**4}
    if name=='square':
        low=F(0) if lo<=0<=hi else min(lo*lo,hi*hi)
        return {'range':(low,max(lo*lo,hi*hi)),'first':2*max(abs(lo),abs(hi)),'second':F(2)}
    if name in ('h','h2'):
        lam=p['lambda'];D=p['D2'];factor=F(1,2) if name=='h' else F(1)
        if lam<=0 or D<joint_norm(2,lam):raise ValueError("Uncertified h normaliser")
        return {'range':(F(0),2*factor),
                'first':factor*(hermite_norm(1)/lam+hermite_norm(2)/lam**2)/D,
                'second':factor*(hermite_norm(2)/lam+hermite_norm(3)/lam**2)/D}
    if name in ('cap','cap_prime'):
        H=p['cap'];lam=p['lambda'];K=p.get('K',F(100))
        if H<=0 or lam<=0 or K!=100:raise ValueError("Invalid frozen cap")
        scaling=1/H if p.get('normalised',False) else F(1)
        d1=min(K+H,H*hermite_norm(0)/lam)
        if name=='cap':return {'range':(F(0),H*scaling),'first':d1*scaling,
                               'second':H*hermite_norm(1)/lam**2*scaling}
        return {'range':(F(0),d1*scaling),'first':H*hermite_norm(1)/lam**2*scaling,
                'second':H*hermite_norm(2)/lam**3*scaling}
    raise ValueError("Unsupported scalar "+name)

def power_two_exponent(value):
    value=F(value);n,d=value.numerator,value.denominator
    if n<=0 or n&(n-1) or d&(d-1):raise ValueError("Scalar span must be a positive power of two")
    return (n.bit_length()-1)-(d.bit_length()-1)

def scalar_certificate(name,lo,span,params,q,table_bits,output_offset=0):
    """Actual unsigned dyadic node integers and rigorous interpolation error.

    Input x=lo+k/2^q. Packed decreasing interpolation must use these
    magnitudes and a sign-selected add/subtract, rather than modular delta.
    """
    ctx.prec=PRECISION
    lo,span,offset=map(F,(lo,span,output_offset))
    if not isinstance(q,int) or not isinstance(table_bits,int) or q<1 or table_bits<1:
        raise ValueError("Positive integer precision and table width required")
    frac_bits=q+power_two_exponent(span)-table_bits
    if frac_bits<0:raise ValueError("Interpolation grid is finer than input register")
    N=1<<table_bits;scale=1<<q;contract=scalar_bounds(name,lo,lo+span,params)
    if contract['range'][0]<offset:raise ValueError("Scalar output offset does not make full domain unsigned")
    nodes=[];largest_error=F(0)
    for j in range(N+1):
        x=lo+span*F(j,N);a=scalar_ball(name,frac_arb(x),params)-frac_arb(offset)
        al,au=endpoints(a);centre=(al+au)/2
        value=(centre*scale+F(1,2)).numerator//(centre*scale+F(1,2)).denominator
        value=max(0,value)
        err=max(abs(F(value,scale)-al),abs(F(value,scale)-au))
        largest_error=max(largest_error,err)
        nodes.append({'input':str(x),'lower':str(al),'upper':str(au),
                      'integer':value,'error_upper':str(err)})
    # Convex node errors, linear interpolation M2 h²/8, one output ulp floor.
    interpolation=contract['second']*(span/N)**2/8
    arithmetic=F(1,scale)
    total=largest_error+interpolation+arithmetic
    return {'schema':'qbp-integral-scalar-v1','name':name,'lo':str(lo),'span':str(span),
      'params':{k:bool(v) if k=='normalised' else str(F(v)) for k,v in (params or {}).items()},
      'q':q,'table_bits':table_bits,'fractional_remainder_bits':frac_bits,
      'output_offset':str(offset),'exact_function_range':[str(x) for x in contract['range']],
      'first_derivative_upper':str(contract['first']),'second_derivative_upper':str(contract['second']),
      'nodes':nodes,'nodal_error_upper':str(largest_error),'linear_interpolation_error_upper':str(interpolation),
      'interpolation_floor_error_upper':str(arithmetic),'output_error_upper':str(total),
      'value_min_integer':min(n['integer'] for n in nodes),
      'value_max_integer':max(n['integer'] for n in nodes),
      'integer_width':max(1,max(n['integer'] for n in nodes).bit_length()),
      'input_perturbation_rule':'Add first_derivative_upper times absolute input error; actual argument must stay in declared domain.',
      'precision_bits':PRECISION}

@lru_cache(maxsize=24)
def quantile_certificate(L,b_y,q):
    ctx.prec=PRECISION;L=F(L)
    if L<=0 or not isinstance(b_y,int) or b_y<1 or not isinstance(q,int) or q<1:
        raise ValueError("Positive L and integer grid/precision widths required")
    N=1<<b_y;scale=1<<q
    c0=normal_cdf(-frac_arb(L));mass=normal_cdf(frac_arb(L))-c0
    nodes=[];integers=[];error=F(0)
    for j in range(N):
        prob=c0+frac_arb(F(2*j+1,2*N))*mass
        y=(2*prob-1).erfinv()*arb(2).sqrt();yl,yu=endpoints(y)
        centre=(yl+yu)/2+L
        rounded=(centre*scale+F(1,2)).numerator//(centre*scale+F(1,2)).denominator
        real=F(rounded,scale)-L
        if not -L<=real<=L:raise ArithmeticError("Rounded quantile outside retained box")
        err=max(abs(real-yl),abs(real-yu));error=max(error,err)
        nodes.append({'j':j,'probability_midpoint':str(F(2*j+1,2*N)),
                      'quantile_lower':str(yl),'quantile_upper':str(yu),
                      'offset_integer':rounded,'actual_coordinate':str(real),'node_error_upper':str(err)})
        integers.append(rounded)
    if any(a>b for a,b in zip(integers,integers[1:])):raise ArithmeticError("Nonmonotone quantile table")
    return {'schema':'qbp-integral-quantile-v1','L':str(L),'b_y':b_y,'q':q,'nodes':nodes,
            'precision_bits':PRECISION,'mass_interval':serial(mass),
            'maximum_node_error_upper':str(error),'monotone':True,
            'offset_integer_bits':max(1,max(integers).bit_length()),
            'finite_law':'Uniform labels, midpoint truncated-normal quantiles, rounded to stated dyadic coordinates; duplicate coordinates retain their labels.',
            'prepared_state':'Exact logical Hadamards plus reversible integer table; no Gaussian amplitude oracle assumed.'}

def exact_coefficients(model,cap,order):
    rows=json.loads((ROOT/'results/taylor_pricing.json').read_text())['configurations']
    row=next(x for x in rows if x['model_id']==model['model_id'] and F(x['cap'])==F(cap))
    lo=hi=F(0)
    for entry in row['orders'][:order+1]:
        c=entry['coefficient_interval'];lo+=F(c['lower_exact_fraction']);hi+=F(c['upper_exact_fraction'])
    return lo,hi

def financial_values(model,cap,y):
    w,d,s,mg,rg,lam=model_parameters(model)
    T=sum(frac_arb(wi)*(frac_arb(di)+sum(frac_arb(v)*ya for v,ya in zip(si,y))).exp()
          for wi,di,si in zip(w,d,s))
    m=frac_arb(mg)+sum(frac_arb(v)*ya for v,ya in zip(rg,y))
    base,prime=kernel_derivatives(m,frac_arb(lam),cap,1)
    arithmetic=kernel_derivatives(m+T.log(),frac_arb(lam),cap,0)[0]
    return {'T':T,'U':T-1,'mG':m,'R1':arithmetic-base,'R2':arithmetic-base-prime*(T-1)}

@lru_cache(maxsize=18)
def geometry_box(w,d,s,rg,mg,L):
    """Exact convex-vertex upper bound and derivative domains on a box."""
    from itertools import product
    ctx.prec=PRECISION;L=F(L);dimension=len(rg)
    maximum=F(1);vertices=0
    for signs in product((-1,1),repeat=dimension):
        y=[L*sg for sg in signs]
        value=sum(frac_arb(wi)*(frac_arb(di)+frac_arb(sum(v*z for v,z in zip(si,y)))).exp()
                  for wi,di,si in zip(w,d,s))
        maximum=max(maximum,upper(value));vertices+=1
    expupper=[upper(frac_arb(di+L*sum(map(abs,si))).exp()) for di,si in zip(d,s)]
    first=tuple(sum(w[i]*abs(s[i][j])*expupper[i] for i in range(len(w))) for j in range(dimension))
    second=tuple(sum(w[i]*s[i][j]**2*expupper[i] for i in range(len(w))) for j in range(dimension))
    mgspan=L*sum(map(abs,rg))
    return {'U_star':maximum-1,'T_star':maximum,'U_first':first,'U_second':second,
            'mG_lower':mg-mgspan,'mG_upper':mg+mgspan,'vertices_evaluated':vertices}

def financial_envelopes(model,cap,L,order=2):
    """Absolute derivative certificates for R_order, not generic raw P errors."""
    ctx.prec=PRECISION;H=F(cap);L=F(L);r=order
    if r not in (1,2):raise ValueError("First/second remainder only")
    w,d,s,mg,rg,lam=model_parameters(model)
    box=geometry_box(w,d,s,rg,mg,L);U=box['U_star']
    D=joint_norm(r,lam);D1=joint_norm(r+1,lam);D2=joint_norm(r+2,lam)
    E=joint_slope(r,lam);E1=joint_slope(r+1,lam);F2=joint_second(r,lam)
    remainder=H*D*U**r/factorial(r)
    prime_max=min(F(100)+H,H*hermite_norm(0)/lam)
    if r==1:remainder=min(H,remainder)
    else:remainder=min(remainder,max(H,prime_max*U))
    first=[];second=[]
    for a,up,upp in zip(rg,box['U_first'],box['U_second']):
        aa=abs(a)
        first.append(H/factorial(r)*(r*U**(r-1)*up*D+
                     U**r*(aa*E/lam+up*D1)))
        pref_second=(r*(r-1)*U**(r-2)*up*up*D if r>1 else F(0))+r*U**(r-1)*upp*D
        cross=2*r*U**(r-1)*up*(aa*E/lam+up*D1)
        chain=U**r*(aa*aa*F2/lam**2+2*aa*up*E1/lam+up*up*D2+upp*D1)
        second.append(H/factorial(r)*(pref_second+cross+chain))
    C=upper((frac_arb(2*L)/(2*arb.pi()).sqrt())**len(rg))
    weighted_second=[ss+2*L*ff+(L*L+1)*remainder for ss,ff in zip(second,first)]
    return {'order':r,'U_star':str(U),'mG_bounds':[str(box['mG_lower']),str(box['mG_upper'])],
       'U_first':[str(x) for x in box['U_first']],'U_second':[str(x) for x in box['U_second']],
       'D':str(D),'D_next':str(D1),'D_second_next':str(D2),'E':str(E),'E_next':str(E1),'F_second':str(F2),
       'remainder_absolute_upper':str(remainder),'coordinate_first':[str(x) for x in first],
       'coordinate_second':[str(x) for x in second],'gaussian_weighted_coordinate_second':[str(x) for x in weighted_second],
       'uniform_box_normaliser_upper':str(C),'uniform_box_signal_scale_upper':str(C*remainder),
       'uniform_midpoint_bias_coefficient':str(C*(2*L)**2*sum(weighted_second)/24),
       'old_quantile_bias_coefficient':str(L*sum(first)),'geometry_vertices':box['vertices_evaluated'],
       'scope':'Whole-box derivative upper envelopes for ideal coordinate midpoints. The quantile L/N bound is before dyadic quantile-node rounding; neither grid envelope includes continuous tail, auxiliary, arithmetic or statistical error. No claim that these sufficient grid-only widths are attained or optimal.'}

def midpoint_width_assessment(model,cap,L,allowance=F(1,1000)):
    """Reports sufficient widths, including unattractive ones."""
    out=[]
    for order in (1,2):
        env=financial_envelopes(model,cap,L,order)
        quant=F(env['old_quantile_bias_coefficient'])
        uniform=F(env['uniform_midpoint_bias_coefficient'])
        bq=max(1,(int(quant/allowance).bit_length()))
        while quant/F(1<<bq)>allowance:bq+=1
        bu=1
        while uniform/F(1<<(2*bu))>allowance:bu+=1
        out.append({'order':order,'allowance':str(allowance),'quantile_Lipschitz_sufficient_bits':bq,
                    'uniform_box_second_derivative_sufficient_bits':bu,'envelope':env})
    return {'model_id':model['model_id'],'cap':cap,'L':str(F(L)),'methods':out}

def finite_input_domain(model,L,b_y,q):
    """Tight tensor-box bounds covering certified and realised quantile nodes."""
    from itertools import product
    ctx.prec=PRECISION
    qu=quantile_certificate(F(L),b_y,q)
    lo=min(min(F(n['quantile_lower']),F(n['actual_coordinate'])) for n in qu['nodes'])
    hi=max(max(F(n['quantile_upper']),F(n['actual_coordinate'])) for n in qu['nodes'])
    w,d,s,mg,rg,lam=model_parameters(model)
    T_upper=F(1);vertices=0
    for y in product((lo,hi),repeat=len(rg)):
        value=sum(frac_arb(wi)*(frac_arb(di+sum(v*z for v,z in zip(si,y)))).exp()
                  for wi,di,si in zip(w,d,s))
        T_upper=max(T_upper,upper(value));vertices+=1
    mglo=mg+sum(min(v*lo,v*hi) for v in rg)
    mghi=mg+sum(max(v*lo,v*hi) for v in rg)
    return {'model_id':model['model_id'],'L_transfer':str(F(L)),
            'b_y':b_y,'q':q,'coordinate_bounds':[str(lo),str(hi)],
            'mG_bounds':[str(mglo),str(mghi)],'T_star':str(T_upper),'U_star':str(T_upper-1),
            'vertices_evaluated':vertices,'covers_quantile_nodal_enclosures_and_actual_dyadic_nodes':True,
            'scope':'Normalisation of the finite financial oracle only; continuous tail/grid transfer still uses its own box and law.'}

def assess_n2_quantile_transfer(model,cap,L,b_y,q,reference_path=None):
    """Certification by independent continuous and finite mean enclosures.

    This explicitly charges a classical-reference premise. It is validation
    of a discretisation, not free information in an advantage comparison.
    """
    ctx.prec=PRECISION;start=time.perf_counter()
    if model['outer_dimension']!=1:raise ValueError("Enumeration transfer currently supports n2 only")
    ref=Path(reference_path or ROOT/'results/continuous_quadrature.json')
    data=json.loads(ref.read_text());row=next(x for x in data['configurations']
                                            if x['model_id']==model['model_id'] and F(x['cap'])==F(cap))
    iv=row['undiscounted_price_interval'];price_lo,price_hi=F(iv['lower_exact_fraction']),F(iv['upper_exact_fraction'])
    quantiles=quantile_certificate(F(L),b_y,q)
    mass=interval(F(quantiles['mass_interval']['lower_exact_fraction']),
                  F(quantiles['mass_interval']['upper_exact_fraction']))
    sums={'R1':arb(0),'R2':arb(0)};moments={2:arb(0),3:arb(0)}
    for node in quantiles['nodes']:
        vals=financial_values(model,cap,[frac_arb(F(node['actual_coordinate']))])
        for k in sums:sums[k]+=vals[k]
        for k in moments:moments[k]+=vals['U']**k
    output={}
    for k in sums:
        c_lo,c_hi=exact_coefficients(model,cap,0 if k=='R1' else 1)
        finite=mass*sums[k]/(1<<b_y);f_lo,f_hi=endpoints(finite)
        # Compare global true correction against retained finite box correction.
        true_lo,true_hi=price_lo-c_hi,price_hi-c_lo
        bias=max(abs(f_lo-true_hi),abs(f_hi-true_lo))
        output[k]={'finite_box_mean':serial(finite),'true_continuous_remainder_interval':serialise_interval(true_lo,true_hi),
                   'financial_transfer_including_tail_and_quantile_rounding_upper':str(bias),
                   'coefficient_interval':serialise_interval(c_lo,c_hi)}
    return {'model_id':model['model_id'],'cap':cap,'L':str(F(L)),'b_y':b_y,'q':q,'kinds':output,
       'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),
       'reference_path':str(ref),'reference_premise':'Previously independently validated continuous price. Paid verification/preprocessing; not a free quantum-speedup premise.',
       'quantile_table_sha256':hashlib.sha256(json.dumps(quantiles,sort_keys=True).encode()).hexdigest(),
       'finite_conditional_U_moments':{str(k):serial(v/(1<<b_y)) for k,v in moments.items()},
       'finite_financial_evaluations':1<<b_y,'seconds_diagnostic':time.perf_counter()-start}

def finite_grid_moments(model,L,b_y,q,max_order=3,distribution='quantile'):
    """Exact product-grid moments using one-dimensional exponential sums.

    No price evaluations or continuous price are used. Uniform-box moments
    retain the relative Gaussian density inside their expectation, so they
    must be multiplied by (2L/sqrt(2pi))^p when used monetarily.
    """
    from taylor_pricing import compositions
    from math import prod
    ctx.prec=PRECISION;L=F(L);N=1<<b_y
    if max_order not in (1,2,3):raise ValueError('Only retained first three moments are supported')
    if distribution=='quantile':
        source=quantile_certificate(L,b_y,q)
        ys=tuple(F(n['actual_coordinate']) for n in source['nodes'])
    elif distribution=='uniform_box':
        if q<b_y+1:raise ValueError('Uniform coordinate midpoints must be exact')
        ys=tuple(-L+2*L*F(2*j+1,2*N) for j in range(N))
        if any((y*(1<<q)).denominator!=1 for y in ys):raise ValueError('Uniform midpoint not represented exactly')
    else:raise ValueError('Unknown finite financial law')
    w,d,s,mg,rg,lam=model_parameters(model);p=len(rg);one_dim={};requests=0
    def mgf(v):
        nonlocal requests
        if v not in one_dim:
            requests+=1
            one_dim[v]=sum(((frac_arb(v*y)-frac_arb(y*y)/2 if distribution=='uniform_box'
                             else frac_arb(v*y)).exp() for y in ys),arb(0))/N
        return one_dim[v]
    raw=[];term_counts=[]
    for j in range(max_order+1):
        total=arb(0);terms=0
        for alpha in compositions(j,len(w)):
            coefficient=F(factorial(j),prod(factorial(a) for a in alpha))
            for wi,a in zip(w,alpha):coefficient*=wi**a
            v=tuple(sum(alpha[i]*s[i][k] for i in range(len(w))) for k in range(p))
            total+=frac_arb(coefficient)*(frac_arb(sum(a*di for a,di in zip(alpha,d)))).exp()*prod(mgf(x) for x in v)
            terms+=1
        raw.append(total);term_counts.append(terms)
    centred={}
    for k in range(1,max_order+1):
        value=sum(((-1)**(k-j)*comb(k,j)*raw[j] for j in range(k+1)),arb(0))
        lo,hi=endpoints(value)
        if hi<0:raise ArithmeticError('Centred nonnegative finite moment has negative upper endpoint')
        centred[str(k)]=serialise_interval(max(F(0),lo),hi)
    return {'schema':'qbp-finite-grid-moments-v1','model_id':model['model_id'],'L':str(L),'b_y':b_y,'q':q,
      'distribution':distribution,'outer_dimension':p,'centred_U_moments':centred,
      'raw_T_moments':[serial(x) for x in raw],'multinomial_terms_by_order':term_counts,
      'distinct_one_dimensional_exponential_sums':requests,'one_dimensional_node_evaluations':requests*N,
      'finite_law':'Independent uniform labels at actual dyadic financial coordinates.',
      'scope':'Exact finite expectation enclosures by separability; includes Gaussian relative weight for uniform_box; no financial price reference is used.'}

def finite_auxiliary_bias(model,cap,beta_bits,strike_bits,moment2_upper,moment3_upper):
    """Integral R2 finite-law monetary bias before the external mass factor."""
    H,lam=F(cap),F(model['lambda']);m2,m3=F(moment2_upper),F(moment3_upper)
    if min(H,lam)>0 and min(m2,m3)>=0 and beta_bits>=1 and strike_bits>=1:
        beta=H*joint_norm(3,lam)*m3/F(1<<(beta_bits+2))
        strike=H*H*joint_slope(2,lam)*m2/(100*lam*F(1<<(strike_bits+2)))
        return {'beta_minimum_midpoint_bias_upper':str(beta),
                'strike_midpoint_bias_upper':str(strike),'total_auxiliary_bias_upper':str(beta+strike),
                'moment_law':'Finite financial input law, with relative likelihood weight included for uniform_box.',
                'scope':'Undiscounted monetary error before multiplication by retained Gaussian mass or uniform-box volume factor.'}
    raise ValueError('Invalid finite-law auxiliary bias inputs')

def rare_flag_prior(moment2_upper,U_star,probability_error_upper,external_moment_factor=F(1)):
    """Both signed-half success probabilities, using the implemented score.

    The relative Gaussian weight, when used, is already inside moment2.
    external_moment_factor is deliberately not an external monetary factor.
    """
    m,U,eta,factor=map(F,(moment2_upper,U_star,probability_error_upper,external_moment_factor))
    if min(m,eta)<0 or factor<1 or U<=0:raise ValueError('Invalid rare-flag prior inputs')
    prior=min(F(1),factor*m/(U*U)+2*eta)
    return {'positive_prior_upper':str(prior),'negative_prior_upper':str(prior),
      'exact_score_absolute_mean_upper':str(min(F(1),factor*m/(U*U))),
      'arithmetic_probability_prior_enlargement':str(2*eta),
      'scope':'Identical bound for positive and negative flags obtained from max(0,2p-1) and max(0,1-2p); their exact implemented difference is 2p-1.'}

def positive_r1_prior(model,cap,L,b_y,q,U_star,probability_error_upper,moment_record=None,distribution='quantile'):
    """Competent finite positive-remainder prior, without price evaluation.

    It retains the correlation between the geometric level and the ratio.
    The underlying-price gap bounds the payoff gap by its one-Lipschitz
    property in the underlying price, before the independent common Gaussian
    factor is integrated exactly.
    """
    ctx.prec=PRECISION;H,L,U,eta=map(F,(cap,L,U_star,probability_error_upper));N=1<<b_y
    if min(H,L,U)<=0 or eta<0:raise ValueError('Invalid positive R1 prior inputs')
    w,d,s,mg,rg,lam=model_parameters(model)
    if distribution=='quantile':ys=tuple(F(n['actual_coordinate']) for n in quantile_certificate(L,b_y,q)['nodes'])
    elif distribution=='uniform_box':
        ys=tuple(-L+2*L*F(2*j+1,2*N) for j in range(N))
        if any((y*(1<<q)).denominator!=1 for y in ys):raise ValueError('Uniform midpoint not represented exactly')
    else:raise ValueError('Unknown finite financial law')
    sums={}
    def mgf(v):
        if v not in sums:sums[v]=sum(((frac_arb(v*y)-frac_arb(y*y)/2 if distribution=='uniform_box'
                                      else frac_arb(v*y)).exp() for y in ys),arb(0))/N
        return sums[v]
    from math import prod
    weighted_T=sum((frac_arb(wi)*frac_arb(di).exp()*prod(mgf(a+b) for a,b in zip(rg,si))
                    for wi,di,si in zip(w,d,s)),arb(0))
    level=prod(mgf(a) for a in rg)
    gap=100*frac_arb(mg+lam*lam/2).exp()*(weighted_T-level)
    glo,ghi=endpoints(gap)
    if ghi<0:raise ArithmeticError('Nonnegative underlying gap has negative upper endpoint')
    moment_record=moment_record or finite_grid_moments(model,L,b_y,q,1,distribution)
    if (moment_record['model_id']!=model['model_id'] or F(moment_record['L'])!=L or
        moment_record['b_y']!=b_y or moment_record['q']!=q or moment_record['distribution']!=distribution):
        raise ValueError('Finite moment record belongs to a different prepared law')
    m1=F(moment_record['centred_U_moments']['1']['upper_exact_fraction'])
    mass=F(moment_record['raw_T_moments'][0]['upper_exact_fraction'])
    D1=joint_norm(1,lam);B=min(H,H*D1*U)
    mean=min(H*mass,H*D1*m1,ghi)
    return {'schema':'qbp-positive-r1-prior-v1','model_id':model['model_id'],'cap':cap,
      'L':str(L),'b_y':b_y,'q':q,'distribution':distribution,'normaliser':str(B),
      'correlated_underlying_gap_interval':serialise_interval(max(F(0),glo),ghi),
      'finite_U_first_moment_upper':str(m1),'finite_relative_density_mass_upper':str(mass),
      'positive_R1_mean_upper':str(mean),'probability_prior_upper':str(min(F(1),mean/B+eta)),
      'arithmetic_probability_prior_enlargement':str(eta),'distinct_one_dimensional_exponential_sums':len(sums),
      'scope':'Finite prepared law and exact common-factor expectation; no financial price reference, continuous tail or external monetary factor is inserted into the probability prior.'}

def complete_bias_budget(model,cap,kind,financial_transfer_upper,graph_certificate,
                         external_factor_interval,beta_bits=None,strike_bits=None,moments=None):
    """Total deterministic price allowance, before the statistical allowance.

    The financial transfer is supplied by its own validated certificate. It
    includes the continuous tail, ideal grid and financial-node rounding.
    The implementation uses interval midpoints for the external coefficient
    and factor. No assumption about attainability of a requested tolerance
    is made here.
    """
    if graph_certificate['status']!='passed':raise ValueError('Rejected graph cannot enter a certified budget')
    if graph_certificate['configuration'].get('output_encoding','centred')!='centred':
        raise ValueError('A complete price budget requires the centred base graph certificate; exact signed-part recombination is charged by that certificate')
    if kind not in ('integral_r2','direct_r2','positive_r1'):raise ValueError('Unknown correction kind')
    clo,chi=exact_coefficients(model,cap,0 if kind=='positive_r1' else 1)
    flo,fhi=map(F,external_factor_interval)
    if not 0<=flo<=fhi:raise ValueError('Invalid external monetary factor interval')
    B=F(graph_certificate['configuration']['monetary_signal_normaliser'])
    aux=F(0);parts={}
    if kind=='integral_r2':
        if moments is None or beta_bits is None or strike_bits is None:raise ValueError('Integral auxiliary certificate required')
        parts=finite_auxiliary_bias(model,cap,beta_bits,strike_bits,moments[2],moments[3])
        aux=F(parts['total_auxiliary_bias_upper'])
    value=F(graph_certificate['monetary_value_error_upper'])
    coefficient=(chi-clo)/2;factor=(fhi-flo)*B/2
    # The implemented probability is in [0,1], hence its monetary signed or
    # positive score is bounded by B even after arithmetic error.
    total=F(financial_transfer_upper)+coefficient+fhi*(aux+value)+factor
    return {'schema':'qbp-integral-complete-bias-v1','model_id':model['model_id'],'cap':cap,'kind':kind,
      'coefficient_midpoint':str((clo+chi)/2),'coefficient_radius_upper':str(coefficient),
      'external_factor_midpoint':str((flo+fhi)/2),'external_factor_radius_price_upper':str(factor),
      'external_factor_upper':str(fhi),'financial_transfer_including_tail_and_node_rounding_upper':str(F(financial_transfer_upper)),
      'auxiliary_before_factor_upper':str(aux),'auxiliary_components':parts,
      'arithmetic_before_factor_upper':str(value),'total_deterministic_undiscounted_price_error_upper':str(total),
      'scope':'Requires the independent continuous/finite financial transfer certificate. Statistical and estimator implementation error remains to be allocated; discount multiplies all terms.'}

def certify_graph(graph,tables,certificates=None):
    """Absolute node-error propagation with integer and signed-domain guards.

    graph.bounds are ideal semantic bounds. Implemented integer bounds are
    propagated independently and intersected with ideal +/- certified error.
    """
    import integral_value_oracle as oracle
    ctx.prec=PRECISION;scale=1<<graph.q;ulp=F(1,scale)
    certificates=certificates or graph.table_certificates
    errors={name:F(0) for name in graph.inputs}
    actual={name:(F(0),F((1<<bits)-1)) for name,bits in graph.inputs.items()}
    records=[];signed_limit=F(1<<(graph.width-1),scale);unsigned_limit=F(1<<graph.width,scale)
    failures=[]
    def absmax(pair):return max(map(abs,pair))
    def fl(x):return F((F(x)*scale).numerator//(F(x)*scale).denominator,scale)
    for name,node in graph.nodes.items():
        aa=[actual[x] for x in node.args];ee=[errors[x] for x in node.args]
        ideal=tuple(map(F,graph.bounds[name]));error=F(0);kind=node.kind
        if kind=='constant':
            value=F(oracle.rounded(node.meta*scale),scale)
            ran=(value,value);error=abs(value-node.meta)
        elif kind in ('uniform','minimum'):
            ran=ideal;error=F(0)
        elif kind=='quantile':
            cc=certificates[name];v=[n['offset_integer'] for n in cc['nodes']]
            ran=(F(min(v),scale),F(max(v),scale));error=F(0)
            if any(a>b for a,b in zip(v,v[1:])):failures.append([name,'nonmonotone quantiles'])
        elif kind=='add':
            ran=(aa[0][0]+aa[1][0],aa[0][1]+aa[1][1]);error=sum(ee)
        elif kind in ('sub','signed_sub'):
            ran=(aa[0][0]-aa[1][1],aa[0][1]-aa[1][0]);error=sum(ee)
        elif kind=='positive_sub':
            ran=(max(F(0),aa[0][0]-aa[1][1]),max(F(0),aa[0][1]-aa[1][0]));error=sum(ee)
        elif kind=='mul':
            if min(aa[0]+aa[1])<0:failures.append([name,'negative multiply operand'])
            ran=(fl(aa[0][0]*aa[1][0]),fl(aa[0][1]*aa[1][1]))
            atrue,btrue=[tuple(map(F,graph.bounds[x])) for x in node.args]
            error=min(absmax(atrue)*ee[1]+absmax(aa[1])*ee[0],
                      absmax(btrue)*ee[0]+absmax(aa[0])*ee[1])+ulp
        elif kind=='scale':
            integer,drop=oracle.scale_integer(node.meta,graph.q);implemented=F(integer,1<<drop)
            ran=(fl(aa[0][0]*implemented),fl(aa[0][1]*implemented))
            error=abs(node.meta)*ee[0]+absmax(aa[0])*abs(implemented-node.meta)+(ulp if drop else 0)
        elif kind in ('clamp','clamp_signed'):
            cap=F(oracle.rounded(node.meta*scale),scale)
            ran=(max(F(0),min(cap,aa[0][0])),max(F(0),min(cap,aa[0][1])))
            error=ee[0]+abs(cap-node.meta)
        elif kind=='function':
            cc=certificates[name];spec=node.meta
            ran=(F(cc['value_min_integer'],scale),F(cc['value_max_integer'],scale))
            error=F(cc['output_error_upper'])+F(cc['first_derivative_upper'])*ee[0]
            if aa[0][0]<0 or aa[0][1]>spec.span:
                failures.append([name,'actual scalar input outside guarded domain',list(map(str,aa[0]))])
            if max(n['integer'] for n in cc['nodes']) >=1<<graph.width:
                failures.append([name,'scalar node output exceeds register'])
            # Every decreasing interpolation subtracts at most its base word.
            vals=[n['integer'] for n in cc['nodes']]
            if min(vals)<0:failures.append([name,'negative unsigned scalar node'])
        else:raise ValueError('Unsupported graph node '+kind)
        # Intersection retains useful correlation without overlooking rounding.
        correlated=(ideal[0]-error,ideal[1]+error)
        ran=(max(ran[0],correlated[0]),min(ran[1],correlated[1]))
        if ran[0]>ran[1]:failures.append([name,'inconsistent actual and semantic bounds'])
        if kind=='signed_sub':
            if ran[0]<=-signed_limit or ran[1]>=signed_limit:
                failures.append([name,'signed modular ambiguity',list(map(str,ran))])
        elif kind not in ('minimum',):
            if ran[0]<0:failures.append([name,'unsigned operation can underflow',list(map(str,ran))])
            if ran[1]>=unsigned_limit:failures.append([name,'unsigned register overflow',list(map(str,ran))])
        if kind=='positive_sub':
            # Its discarded modulo difference remains guarded by a predicate.
            raw=(aa[0][0]-aa[1][1],aa[0][1]-aa[1][0])
            if max(map(abs,raw))>=unsigned_limit:failures.append([name,'positive-sub pre-mask magnitude overflow'])
        actual[name]=ran;errors[name]=error
        records.append({'node':name,'kind':kind,'ideal_bounds':list(map(str,ideal)),
          'actual_bounds':list(map(str,ran)),'absolute_output_error_upper':str(error)})
    out=actual[graph.output]
    if out[0]<0 or out[1]>1:failures.append([graph.output,'probability output outside[0,1]'])
    B=F(graph.configuration['monetary_signal_normaliser'])
    encoding=graph.configuration.get('output_encoding','centred')
    decode=F(1) if graph.configuration['kind']=='positive_r1' or encoding!='centred' else F(2)
    return {'schema':'qbp-integral-graph-certificate-v1','status':'passed' if not failures else 'rejected',
       'configuration':graph.configuration,'nodes':records,'failures':failures,
       'probability_error_upper':str(errors[graph.output]),
       'monetary_value_error_upper':str(decode*B*errors[graph.output]),
       'probability_actual_bounds':list(map(str,out)),
       'rounding_unit':str(ulp),'endpoint_metadata_lattice_bits':ENDPOINT_BITS,
       'decoded_observable':'Positive or negative signed part only' if encoding!='centred' else 'Full centred second correction or positive first correction',
       'scope':'Conditional on declared ideal financial ranges, independently checked model/domain premises and retained scalar/quantile tables; no tail/grid/statistical allowance is reused. Split-part monetary error concerns that part alone; a complete signed price budget requires the centred base certificate.'}
