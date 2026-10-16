"""Independent finite Gaussian Taylor competitor for ordinary basket calls.
Read-only production inputs; no finite-grid or continuous basket price.
"""
import sys,json,time,hashlib,os,argparse
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
from flint import arb,ctx
REPO=Path(sys.argv[sys.argv.index('--root')+1]) if '--root' in sys.argv else Path(os.environ.get('QBP_RESEARCH_ROOT',str(Path(__file__).resolve().parent)))
if not (REPO/'results/generalisation_protocol.json').exists():
    REPO=Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
sys.path.insert(0,str(REPO));sys.dont_write_bytecode=True
from taylor_pricing import model_parameters,exact_geometry_terms,dot
from boundary_derivative_certificate import payoff_norms
from boundary_grid_certificate import ball,serial,endpoints,Phi,phi

def polynomial(order,sigma,z):
    p=[arb(1)]
    for r in range(2,order):
        nxt=[arb(0)]*(len(p)+1)
        for j,c in enumerate(p):
            nxt[j]-=r*c
            nxt[j+1]-=c/sigma
            if j:nxt[j-1]+=j*c/sigma
        p=nxt
    v=arb(0)
    for c in reversed(p):v=v*z+c
    return v

def C(order,m,sigma):
    z=m/sigma
    return 100*polynomial(order,sigma,z)*phi(z)/sigma

def evaluate(model,maximum=5):
    ctx.prec=512;start=time.perf_counter();w,d,s,mg,rg,lam=model_parameters(model)
    geom=exact_geometry_terms(w,d,s,rg,maximum+1);weighted=[];raw=[]
    for row in geom:
        terms=[(ball(coef)*ball(exponent).exp(),shift) for alpha,coef,exponent,shift in row]
        weighted.append(terms);raw.append(sum((c for c,shift in terms),arb(0)))
    mu=[sum(((-1)**(k-j))*comb(k,j)*raw[j] for j in range(k+1)) for k in range(maximum+2)]
    for k,value in enumerate(mu):
        if endpoints(value)[1]<0:raise ArithmeticError('Nonnegative U moment excluded')
    sigma=ball(lam*lam+dot(rg,rg)).sqrt();m=ball(mg)
    slope=lambda x:100*(x+sigma*sigma/2).exp()*Phi(x/sigma+sigma)
    c0=slope(m)-100*Phi(m/sigma)
    c1=sum((coef*slope(m+ball(shift)) for coef,shift in weighted[1]),arb(0))-slope(m)
    coefficients=[c0,c1]
    for k in range(2,maximum+1):
        coeff=arb(0)
        for j in range(k+1):
            coeff+=((-1)**(k-j))*comb(k,j)*sum((factor*C(k,m+ball(shift),sigma) for factor,shift in weighted[j]),arb(0))
        coefficients.append(coeff/factorial(k))
    norm=payoff_norms(lam,None)['C_r_upper'];partial=c0;rows=[]
    discount=ball(-F(model['risk_free_rate'])*F(model['maturity'])).exp()
    for k in range(1,maximum+1):
        partial+=coefficients[k]
        radius=F(norm[str(k+1)])*max(F(0),endpoints(mu[k+1])[1])/factorial(k+1)
        lo,hi=endpoints(partial)
        lo=max(F(0),lo-(F(0) if k==1 else radius));hi+=radius
        dlow,dhi=endpoints(discount);prices=[x*y for x in (lo,hi) for y in (dlow,dhi)]
        plo,phi_=min(prices),max(prices);prad=(phi_-plo)/2
        rows.append({'order':k,'coefficient':serial(coefficients[k]),'partial_sum':serial(partial),'remainder_upper':str(radius),'continuous_U_moment':serial(mu[k+1]),'discounted_interval':{'lower_exact_fraction':str(plo),'upper_exact_fraction':str(phi_)},'discounted_radius_upper':str(prad),'certifies_cent':prad<=F(1,100),'geometric_kernel_terms':sum(map(len,geom[:k+1]))})
    return {'model_id':model['model_id'],'orders':rows,'first_successful_order':next((r['order'] for r in rows if r['certifies_cent']),None),'elapsed_seconds_diagnostic':time.perf_counter()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=REPO);ap.add_argument('--out',type=Path,default=REPO/'results/higherdim_classical_taylor.json');args=ap.parse_args()
    protocol=args.root/'results/generalisation_protocol.json';data=json.loads(protocol.read_text());out=[]
    for model in data['models']:
        result=evaluate(model);out.append(result)
        print(result['model_id'],result['first_successful_order'],[round(float(F(r['discounted_radius_upper'])),8) for r in result['orders']],flush=True)
    result={'schema':'qbp-independent-vanilla-taylor-screen-v1','scope':'Established finite Taylor consequence, not an original generic theorem. Continuous geometric Gaussian kernel coefficients and lognormal centred moments only. Maximum degree5, remainder order6. No financial grid or continuous basket price is an input.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256(protocol.read_bytes()).hexdigest(),'dependency_source_sha256':{f:hashlib.sha256((REPO/f).read_bytes()).hexdigest() for f in ['taylor_pricing.py','boundary_derivative_certificate.py','boundary_grid_certificate.py']},'boundary_norm_certificate_sha256':hashlib.sha256((REPO/'results/boundary_derivative_norms.json').read_bytes()).hexdigest(),'precision_bits':512,'records':out}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
