"""Price-independent Gaussian-CDF-cell transfer certificates for n2 models.

Only frozen model parameters, Gaussian functions and actual dyadic quantile
nodes enter the construction. No continuous basket price, quadrature price
or earlier financial grid assessment is an input. Capped and vanilla first
and second residuals share the same original outer law.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from math import ceil,comb,factorial
import gzip,hashlib,json,time
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
PRECISION=512
ENDPOINT_BITS=768


def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator

def exact_endpoint(x,upper):
    y=x.upper() if upper else x.lower();m,e=map(int,y.man_exp())
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))

def outward(x,up,bits=ENDPOINT_BITS):
    x=F(x);scale=1<<bits;t=x*scale;k=t.numerator//t.denominator
    return F(k+int(up and t.denominator!=1),scale)

def endpoints(x):
    if not x.is_finite():raise ArithmeticError('Nonfinite Arb interval')
    return outward(exact_endpoint(x,False),False),outward(exact_endpoint(x,True),True)

def interval(lo,hi):
    lo,hi=F(lo),F(hi)
    if lo>hi:raise ValueError('Reversed interval')
    return ball((lo+hi)/2)+arb(0,ball((hi-lo)/2).upper())

def serial(x):
    lo,hi=endpoints(x)
    return {'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi)}

def abs_upper(x):
    lo,hi=endpoints(x);return max(abs(lo),abs(hi))

def upper(x):return endpoints(x)[1]
def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def Phi(x):return (1+(x/arb(2).sqrt()).erf())/2


def parameters(model,one_dimension=True):
    w=tuple(map(F,model['weights']));r=tuple(tuple(map(F,x)) for x in model['loadings']);mu=tuple(map(F,model['log_drifts']));lam=F(model['lambda'])
    p=model['outer_dimension'];rate=F(model['risk_free_rate'])-F(model['dividend_yield'])
    if model['spot']!=100 or model['strike']!=100 or F(model['maturity'])!=1:raise ValueError('Frozen unit-maturity K100 models required')
    if len(w)!=model['assets'] or min(w)<=0 or sum(w)!=1 or lam<=0:raise ValueError('Invalid model')
    if one_dimension and (len(w)!=2 or p!=1):raise ValueError('The complete cell certificate requires n2 and p1')
    if any(len(a)!=p for a in r):raise ValueError('Loading dimension mismatch')
    if any(mu[i]+(sum(x*x for x in r[i])+lam*lam)/2!=rate for i in range(len(w))):raise ValueError('Risk-neutral drift premise failed')
    mg=sum(a*b for a,b in zip(w,mu));rg=tuple(sum(w[i]*r[i][j] for i in range(len(w))) for j in range(p))
    d=tuple(x-mg for x in mu);s=tuple(tuple(x-y for x,y in zip(row,rg)) for row in r)
    assert sum(a*b for a,b in zip(w,d))==0
    assert all(sum(w[i]*s[i][j] for i in range(len(w)))==0 for j in range(p))
    return w,d,s,mg,rg,lam


def kernels(m,lam,cap,payoff='capped'):
    """P', P'', P''', C2=P''-P' and C2' in outward Arb arithmetic."""
    lam=ball(lam);K=arb(100);growth=K*(m+lam*lam/2).exp();x=m/lam
    if payoff=='vanilla':
        p1=growth*Phi(x+lam);C=K*phi(x)/lam;Cp=-K*x*phi(x)/(lam*lam)
    elif payoff=='capped':
        H=ball(cap);a=(1+H/K).log();xh=(m-a)/lam
        p1=growth*(Phi(x+lam)-Phi(xh+lam))
        C=(K*phi(x)-(K+H)*phi(xh))/lam
        Cp=(-K*x*phi(x)+(K+H)*xh*phi(xh))/(lam*lam)
    else:raise ValueError('Unknown payoff')
    p2=p1+C;p3=p2+Cp
    return p1,p2,p3,C,Cp


def geometry(y,model):
    w,d,s,mg,rg,lam=parameters(model);r=rg[0]
    terms=[ball(wi)*(ball(di)+ball(si[0])*y).exp() for wi,di,si in zip(w,d,s)]
    T=sum(terms,arb(0));Tp=sum((ball(si[0])*term for si,term in zip(s,terms)),arb(0));Tpp=sum((ball(si[0]**2)*term for si,term in zip(s,terms)),arb(0))
    # Jensen gives T>=1 pointwise. Intersect this proven premise with the
    # natural exponential interval so division and log keep valid domains.
    lo,hi=endpoints(T);T=interval(max(F(1),lo),hi)
    m=ball(mg)+ball(r)*y;U=T-1;mA=m+T.log()
    return T,Tp,Tpp,m,U,mA,r,lam


def residual_derivatives(y,model,cap,payoff='capped',order=2,form='stable'):
    if order not in (1,2):raise ValueError('Residual order must be one or two')
    T,Tp,Tpp,m,U,mA,r,lam=geometry(y,model);rb=ball(r)
    pG,pGG,pGGG,CG,CpG=kernels(m,lam,cap,payoff)
    pA,pAA,pAAA,CA,CpA=kernels(mA,lam,cap,payoff)
    if form=='raw':
        a=rb+Tp/T;b=Tpp/T-(Tp/T)**2
        first=pA*a-rb*pG
        second=pAA*a*a+pA*b-rb*rb*pGG
        if order==2:
            first-=Tp*pG+U*rb*pGG
            second-=Tpp*pG+2*rb*Tp*pGG+U*rb*rb*pGGG
        return first,second
    if form!='stable':raise ValueError('Unknown derivative form')
    lb=ball(lam);growth=100*(m+lb*lb/2).exp();x=m/lb;xa=mA/lb
    if payoff=='vanilla':v=growth*(Phi(xa+lb)-Phi(x+lb))
    else:
        a=(1+ball(cap)/100).log()/lb
        v=growth*((Phi(xa+lb)-Phi(xa-a+lb))-(Phi(x+lb)-Phi(x-a+lb)))
    A=Tp+rb*T;B=Tpp+2*rb*Tp+rb*rb*T
    first=A*v-rb*U*CG
    second=B*v+(A/T)**2*CA-rb*(2*Tp+rb*T)*CG-rb*rb*U*CpG
    if order==1:
        first+=Tp*pG+U*rb*pGG
        second+=Tpp*pG+2*rb*Tp*pGG+U*rb*rb*pGGG
    return first,second


@lru_cache(maxsize=4)
def gaussian_cells(L=F(6),b_y=10,q=28,precision=PRECISION):
    ctx.prec=precision;L=F(L);N=1<<b_y;lb=ball(L);cdf_low=Phi(-lb);mass=Phi(lb)-cdf_low;cell_mass=mass/N
    boundaries=[-lb]
    for j in range(1,N):boundaries.append(arb(2).sqrt()*(2*(cdf_low+mass*j/N)-1).erfinv())
    boundaries.append(lb);records=[]
    for j in range(N):
        p=cdf_low+mass*ball(F(2*j+1,2*N));node=arb(2).sqrt()*(2*p-1).erfinv();nl,nh=endpoints(node)
        sl=(nl*(1<<q)+F(1,2));sh=(nh*(1<<q)+F(1,2));il=sl.numerator//sl.denominator;ih=sh.numerator//sh.denominator
        if il!=ih:raise ArithmeticError('Quantile rounding is not uniquely certified')
        actual=F(il,1<<q);a,b=boundaries[j],boundaries[j+1];al,ah=endpoints(a);bl,bh=endpoints(b)
        if not ah<=actual<=bl:raise ArithmeticError('Rounded node not within its CDF cell')
        pA,pB=phi(a),phi(b);first=pA-pB-ball(actual)*cell_mass
        second=(1+ball(actual*actual))*cell_mass+(a-2*ball(actual))*pA-(b-2*ball(actual))*pB
        if endpoints(second)[0]<0:raise ArithmeticError('Second Gaussian cell moment not certified nonnegative')
        records.append({'cell':j,'left':serial(a),'right':serial(b),'actual_node':str(actual),'actual_node_integer':il,
            'quantile_midpoint':serial(node),'mass':serial(cell_mass),'first_moment':serial(first),'second_moment':serial(second),
            'left_cdf':serial(cdf_low+mass*j/N),'right_cdf':serial(cdf_low+mass*(j+1)/N),
            '_first':first,'_second':second,'_left':a,'_right':b,'_left_cdf':cdf_low+mass*j/N,'_right_cdf':cdf_low+mass*(j+1)/N})
    return mass,records


def tilted_outside(a,L):
    a,L=ball(a),ball(L);return (a*a/2).exp()*(Phi(a-L)+Phi(-a-L))


def outside_U_moment(model,power=2,L=6,extra_tilt=F(0)):
    w,d,s,mg,rg,lam=parameters(model);result=arb(0);terms=[]
    if power not in (1,2):raise ValueError('Only first and second moments are needed')
    for j in range(power+1):
        for i in range(j+1):
            coeff=F(comb(power,j)*(-1)**(power-j)*comb(j,i))*w[0]**i*w[1]**(j-i)
            constant=i*d[0]+(j-i)*d[1];a=i*s[0][0]+(j-i)*s[1][0]+F(extra_tilt)
            value=ball(coeff)*ball(constant).exp()*tilted_outside(a,L);result+=value
            terms.append({'power_of_T':j,'asset0_multiplicity':i,'coefficient':str(coeff),'constant':str(constant),'tilt':str(a),'term_interval':serial(value)})
    lo,hi=endpoints(result)
    if hi<0:raise ArithmeticError('Nonnegative tilted U moment has negative enclosure')
    return interval(max(F(0),lo),hi),terms


def tail_certificate(model,cap,order=2,payoff='capped',L=6):
    w,d,s,mg,rg,lam=parameters(model);gaussian_tail=Phi(-ball(L))+Phi(-ball(L))
    if order==2:
        moment,terms=outside_U_moment(model,2,L)
        fallback=ball(100 if payoff=='vanilla' else 200+F(cap))/(ball(lam)*(2*arb.pi()).sqrt())
        norm_upper=upper(fallback);norm_source='Gaussian density triangle bound'
        bound=norm_upper*upper(moment)/2
        return {'order':order,'U_moment_outside':serial(moment),'tilted_terms':terms,'boundary_M2_upper':str(norm_upper),'boundary_norm_source':norm_source,'tail_error_upper':str(bound)}
    moment,terms=outside_U_moment(model,1,L,rg[0]);unbounded=upper(100*ball(mg+lam*lam/2).exp()*moment)
    bound=min(F(cap)*upper(gaussian_tail),unbounded) if payoff=='capped' else unbounded
    return {'order':order,'weighted_U_moment_outside':serial(moment),'tilted_terms':terms,'vanilla_delta_tail_upper':str(unbounded),'tail_error_upper':str(bound)}


def coefficient_certificate(model,cap=None,payoff='vanilla',order=1):
    """Exact Gaussian geometric coefficients, independent of basket prices."""
    if order not in (0,1):raise ValueError('Only c0 and c0+c1 are exposed')
    if payoff not in ('capped','vanilla') or payoff=='capped' and (cap is None or F(cap)<=0):raise ValueError('Invalid coefficient payoff')
    ctx.prec=PRECISION;w,d,s,mg,rg,lam=parameters(model,False);variance=lam*lam+sum(v*v for v in rg);sigma=ball(variance).sqrt()
    def kernel0(m):
        growth=100*(m+sigma*sigma/2).exp();x=m/sigma
        if payoff=='vanilla':return growth*Phi(x+sigma)-100*Phi(x)
        H=ball(cap);a=(1+H/100).log();xh=(m-a)/sigma
        return growth*(Phi(x+sigma)-Phi(xh+sigma))-100*(Phi(x)-Phi(xh))+H*Phi(xh)
    def kernel1(m):
        growth=100*(m+sigma*sigma/2).exp();x=m/sigma
        if payoff=='vanilla':return growth*Phi(x+sigma)
        a=(1+ball(cap)/100).log();return growth*(Phi(x+sigma)-Phi((m-a)/sigma+sigma))
    c0=kernel0(ball(mg));c1=sum((ball(wi)*ball(di+sum(v*v for v in si)/2).exp()*kernel1(ball(mg+sum(a*b for a,b in zip(rg,si)))) for wi,di,si in zip(w,d,s)),arb(0))-kernel1(ball(mg))
    result=c0 if order==0 else c0+c1
    return {'schema':'qbp-price-independent-geometric-coefficient-v1','model_id':model['model_id'],'payoff':payoff,'cap':cap,'order':order,
        'c0':serial(c0),'c1':serial(c1),'coefficient':serial(result),'effective_variance':str(variance),'precision_bits':PRECISION,
        'source_scope':'Closed Gaussian geometric conditional-kernel coefficients; no arithmetic basket price is evaluated.'}



def transport_primitive(z):
    return ((z*z+1)*Phi(z)+z*phi(z))/2

def transport_weight_mass(u,v,endpoint,endpoint_cdf):
    ub,vb=ball(u),ball(v)
    value=transport_primitive(vb)-transport_primitive(ub)-phi(endpoint)*(vb-ub)-endpoint_cdf*(vb*vb-ub*ub)/2
    lo,hi=endpoints(value)
    if hi<0:raise ArithmeticError('Nonnegative Gaussian transport weight has negative enclosure')
    # h(z)=phi(z)+z Phi(z) is convex because h''=phi. Consequently
    # h(z)-h(e)-h'(e)(z-e) is nonnegative even outside the exact cell.
    return interval(max(F(0),lo),hi)

def financial_grid_transfer(model,cap,order=2,payoff='capped',L=6,b_y=10,q=28,max_derivative_width=F(1,20)):
    if order not in (1,2) or payoff not in ('capped','vanilla') or payoff=='capped' and (cap is None or F(cap)<=0):raise ValueError('Invalid residual or payoff')
    if not isinstance(b_y,int) or isinstance(b_y,bool) or b_y<1 or not isinstance(q,int) or isinstance(q,bool) or q<1 or F(L)<=0 or F(max_derivative_width)<=0:raise ValueError('Invalid finite grid')
    ctx.prec=PRECISION;parameters(model);mass,cells=gaussian_cells(F(L),b_y,q,PRECISION)
    firstsum=arb(0);quadratic=F(0);records=[];subdivisions=0
    for cell in cells:
        actual=F(cell['actual_node']);first,second=residual_derivatives(ball(actual),model,cap,payoff,order)
        contribution=first*cell['_first'];firstsum+=contribution
        lo=F(cell['left']['lower_exact_fraction']);hi=F(cell['right']['upper_exact_fraction']);span=hi-lo
        bound=F(0);subcells=[];weighted=F(0)
        for side,outer,e,ecdf in (('left',lo,cell['_left'],cell['_left_cdf']),('right',hi,cell['_right'],cell['_right_cdf'])):
            left,right=(outer,actual) if side=='left' else (actual,outer)
            parts=max(1,ceil((right-left)/F(max_derivative_width)))
            for j in range(parts):
                a=left+(right-left)*F(j,parts);b=left+(right-left)*F(j+1,parts)
                derivative=residual_derivatives(interval(a,b),model,cap,payoff,order)[1]
                value=abs_upper(derivative);bound=max(bound,value)
                weight=transport_weight_mass(a,b,e,ecdf);weighted_part=value*upper(weight);weighted+=weighted_part
                subcells.append({'side':side,'left':str(a),'right':str(b),'second_derivative_interval':serial(derivative),
                    'absolute_second_derivative_upper':str(value),'gaussian_transport_weight_mass':serial(weight),'remainder_upper':str(weighted_part)})
        moment_upper=F(cell['second_moment']['upper_exact_fraction']);unweighted=bound*moment_upper/2;remainder=min(unweighted,weighted)
        quadratic+=remainder;subdivisions+=len(subcells)
        record={k:v for k,v in cell.items() if not k.startswith('_')};record.update({'first_derivative_at_actual_node':serial(first),'signed_linear_contribution':serial(contribution),
            'second_derivative_cell_upper':str(bound),'whole_cell_quadratic_upper':str(unweighted),'weighted_transport_remainder_upper':str(weighted),
            'quadratic_remainder_upper':str(remainder),'derivative_subcells':subcells});records.append(record)
    tail=tail_certificate(model,cap,order,payoff,L);linear=abs_upper(firstsum);total=linear+quadratic+F(tail['tail_error_upper'])
    return {'schema':'qbp-price-independent-financial-grid-transfer-v1','status':'passed' if total<F(1,200) else 'retained allocation failure',
        'model_id':model['model_id'],'cap':cap,'payoff':payoff,'order':order,'L':str(F(L)),'b_y':b_y,'q':q,'mass_interval':serial(mass),
        'first_moment_sum_interval':serial(firstsum),'signed_first_moment_error_upper':str(linear),'quadratic_remainder_upper':str(quadratic),
        'tail_certificate':tail,'tail_upper':tail['tail_error_upper'],'financial_transfer_upper':str(total),'target_allocation':'1/200',
        'cell_count':len(cells),'derivative_subinterval_count':subdivisions,'maximum_derivative_subinterval_width':str(F(max_derivative_width)),
        'precision_bits':PRECISION,'endpoint_lattice_bits':ENDPOINT_BITS,'cells':records,
        'construction_inputs':'Frozen exact synthetic model; Gaussian CDF/density; actual rounded quantile nodes. No continuous arithmetic basket price or numerical pricing answer is an input.',
        'target':'Absolute undiscounted transfer between the continuous Gaussian residual expectation and retained mass times the uniform rounded-node expectation. Exact node rounding is included within cell moments.'}


def main():
    start=time.perf_counter();path=ROOT/'results/generalisation_protocol.json';data=json.loads(path.read_text());models=[m for m in data['models'] if m['outer_dimension']==1]
    records=[];coefficients=[]
    for model in models:
        for cap in (1,10,25):
            for order in (1,2):
                record=financial_grid_transfer(model,cap,order,'capped');records.append(record)
                print(model['model_id'],cap,order,'capped',float(F(record['financial_transfer_upper'])),flush=True)
            coefficients.append(coefficient_certificate(model,cap,'capped',1))
        for order in (1,2):
            record=financial_grid_transfer(model,None,order,'vanilla');records.append(record)
            print(model['model_id'],None,order,'vanilla',float(F(record['financial_transfer_upper'])),flush=True)
        coefficients.append(coefficient_certificate(model,None,'vanilla',1))
    result={'schema':'qbp-boundary-price-independent-grid-certificate-v1','status':'complete retained bound construction','records':records,'coefficients':coefficients,
        'unsupported_scope':[{'model_id':m['model_id'],'cap':cap,'status':'Higher-dimensional hybrid-coordinate theorem applies, but a complete monetary bound was not constructed; no ranking as inferior is implied.'} for m in data['models'] if m['outer_dimension']!=1 for cap in (1,10,25)],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':{'generalisation_protocol.json':hashlib.sha256(path.read_bytes()).hexdigest()},
        'elapsed_seconds_diagnostic':time.perf_counter()-start,'scope':'Price-independent finite-law transfer and geometric coefficient certificates. Classical setup remains a charged computational task; no global quantum advantage claim.'}
    out=ROOT/'results/boundary_financial_grid_transfer.json.gz'
    with out.open('wb') as f:
        with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(json.dumps(result,separators=(',',':')).encode())
    print(json.dumps({'status':'passed' if all(r['status']=='passed' for r in records) else 'retained failures','records':len(records),'maximum_transfer':float(max(F(r['financial_transfer_upper']) for r in records)),'seconds':time.perf_counter()-start}),flush=True)
if __name__=='__main__':main()
