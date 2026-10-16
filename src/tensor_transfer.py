"""Price-independent tensor transfer using positive Gaussian transport MGFs.

Continuous denotes the unnormalised retained Gaussian measure throughout;
rounded denotes its actual dyadic midpoint grid with the SAME mass. Neither
measure is silently normalised inside a hybrid. No basket price is an input.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json,time
from flint import arb,ctx
from boundary_grid_certificate import (parameters,gaussian_cells,ball,interval,
    endpoints,serial,upper,abs_upper,phi,Phi)
from boundary_derivative_certificate import payoff_norms
from higherdim_preparation import get_mgf,EmpiricalMGF

ROOT=Path(__file__).resolve().parent
PRECISION=512

def positive(x,what):
    lo,hi=endpoints(x)
    if hi<0: raise ArithmeticError('Negative enclosure for nonnegative '+what)
    return interval(max(F(0),lo),hi)

class TransportMGF:
    def __init__(self,L=6,b_y=10,q=28):
        ctx.prec=PRECISION; self.L=F(L);self.b_y=b_y;self.q=q
        self.mass,self.cells=gaussian_cells(self.L,b_y,q,PRECISION)
        self.grid=get_mgf(self.L,b_y,q)
        nodes=tuple(F(c['actual_node']) for c in self.cells)
        if nodes!=self.grid.coordinates: raise ArithmeticError('Actual nodes differ')
        # This positive dominating measure avoids an ambiguous sign decision
        # for a ball enclosing a tiny centred first moment.
        first_weights=[abs_upper(c['_first']) for c in self.cells]
        self.rho1=EmpiricalMGF(nodes,first_weights)
        self.signed=EmpiricalMGF(nodes,[c['_first'] for c in self.cells])
        self.rho2mass=positive(sum((c['_second']/2 for c in self.cells),arb(0)), 'transport mass')
        self.cache_cont={};self.cache_second={};self.fallbacks=[]

    def continuous(self,rate):
        rate=F(rate)
        if rate not in self.cache_cont:
            v=ball(rate);L=ball(self.L)
            self.cache_cont[rate]=(v*v/2).exp()*(Phi(L-v)-Phi(-L-v))
        return self.cache_cont[rate]

    def rounded(self,rate): return self.mass*self.grid(F(rate))

    def second(self,rate):
        rate=F(rate)
        if rate not in self.cache_second:
            if not rate:value=self.rho2mass
            elif abs(rate)<F(1,10**20):
                # Positive support in [-L,L] gives a rigorous cancellation-safe
                # enclosure, including almost-zero exact rational loadings.
                lo,hi=endpoints(self.rho2mass);a=ball(abs(rate)*self.L)
                value=interval(endpoints(ball(lo)*(-a).exp())[0],upper(ball(hi)*a.exp()))
                self.fallbacks.append({'rate':str(rate),'reason':'tiny-rate positive support bound'})
            else:
                v=ball(rate)
                value=positive((self.continuous(rate)-self.rounded(rate)-v*self.signed(rate))/(v*v),'transport MGF')
            self.cache_second[rate]=value
        return self.cache_second[rate]

    def profile(self):
        return {'L':str(self.L),'b_y':self.b_y,'q':self.q,'cell_count':len(self.cells),
          'mass_interval':serial(self.mass),'rho1_mass_interval':serial(self.rho1(0)),
          'rho2_mass_interval':serial(self.rho2mass),'rho2_distinct_rates':len(self.cache_second),
          'continuous_distinct_rates':len(self.cache_cont),'tiny_rate_fallbacks':self.fallbacks,
          'rho1_positive_weights':'Exact rational upper bounds on absolute centred first moments',
          'node_certificate':self.grid.node_certificate,'precision_bits':PRECISION}

def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}

def scale(a,x):return {k:v*F(x) for k,v in a.items() if v*F(x)}

def mul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=tuple(x+y for x,y in zip(ka,kb));out[k]=out.get(k,F(0))+va*vb
    return {k:v for k,v in out.items() if v}

def geometry_polynomials(model,j):
    w,d,s,mg,rg,lam=parameters(model,False);n=len(w);zero=(0,)*n
    T={};Tj={};Tjj={}
    for i in range(n):
        k=tuple(int(a==i) for a in range(n));T[k]=w[i]
        Tj[k]=w[i]*s[i][j];Tjj[k]=w[i]*s[i][j]**2
    U=add(T,{zero:F(-1)});U2=mul(U,U);J2=mul(Tj,Tj)
    return {'U2':U2,'J2':J2,'U_JJ':mul(U,Tjj),'U_J2':mul(U,J2),
        'U2_JJ':mul(U2,Tjj),'U2_J2':mul(U2,J2),'U4':mul(U2,U2)}

class HybridMoments:
    def __init__(self,model,transport):
        self.model=model;self.transport=transport
        self.w,self.d,self.s,self.mg,self.rg,self.lam=parameters(model,False)
        self.p=len(self.rg);self.exp_cache={};self.counts={};self.reuse={}

    def monomial(self,alpha,j,kind):
        key=(alpha,j,kind)
        if key not in self.exp_cache:
            if alpha not in self.reuse:
                constant=sum(F(a)*d for a,d in zip(alpha,self.d))
                rates=tuple(sum(F(a)*si[k] for a,si in zip(alpha,self.s)) for k in range(self.p))
                self.reuse[alpha]=(ball(constant).exp(),rates)
            value,rates=self.reuse[alpha]
            # ν for coordinates before j, τ for coordinates after j. All
            # factors are subprobability masses, including the zero rate.
            for k,rate in enumerate(rates):
                if k==j: factor=self.transport.rho1(rate) if kind==1 else self.transport.second(rate)
                elif k<j:factor=self.transport.continuous(rate)
                else:factor=self.transport.rounded(rate)
                value=value*factor
            self.exp_cache[key]=value
        return self.exp_cache[key]

    def moment(self,poly,j,kind,name):
        value=sum((ball(c)*self.monomial(a,j,kind) for a,c in poly.items()),arb(0))
        self.counts[(j,kind,name)]=len(poly)
        return positive(value,name)

    def record(self):
        out=[]
        for j in range(self.p):
            polys=geometry_polynomials(self.model,j)
            terms={}
            for kind in (1,2):
                terms[str(kind)]={name:serial(self.moment(poly,j,kind,name)) for name,poly in polys.items()}
            out.append({'coordinate':j,'geometric_loading':str(self.rg[j]),'moments':terms})
        return {'model_id':self.model['model_id'],'coordinates':out,'distinct_exponential_monomials':len(self.reuse),
          'distinct_hybrid_monomial_evaluations':len(self.exp_cache),
          'scope':'Nonnegative geometry moments under rho1/rho2 and other-coordinate continuous/rounded hybrids. U cancellation retained in signed exponential-polynomial sums.'}

def B(rec):return interval(F(rec['lower_exact_fraction']),F(rec['upper_exact_fraction']))

def derivative_envelope(norms,row):
    M2=F(norms['M2_upper']);M3=F(norms['C_r_upper']['3']);M4=F(norms['C_r_upper']['4'])
    if norms['payoff']=='capped':
        D1=2*M2*F(norms['scalar_boundary_h_first_upper']);D2=2*M2*F(norms['scalar_boundary_h_second_upper'])
    else:
        D1=M2*F(norms['scalar_boundary_density_first_upper']);D2=M2*F(norms['scalar_boundary_density_second_upper'])
    D31=min(M4+3*M3,D2+2*D1)
    r=abs(F(row['geometric_loading']));result={}
    for kind in (1,2):
        mom={name:B(value) for name,value in row['moments'][str(kind)].items()}
        # Cauchy is valid under any finite positive measure: no mass divisor.
        uj=(mom['U2']*mom['J2']).sqrt()
        u2j_a=(mom['U2']*mom['U2_J2']).sqrt();u2j_b=(mom['U4']*mom['J2']).sqrt()
        u2j=min(upper(u2j_a),upper(u2j_b))
        if kind==1:
            components={'U_Tj':M2*upper(uj),'r_U2':r*D1*upper(mom['U2'])/2,
                'U2_Tj':M3*u2j/6}
        else:
            components={'Tj2_U_Tjj':M2*(upper(mom['J2'])+upper(mom['U_JJ'])),
                'r_U_Tj':2*r*D1*upper(uj),'U_Tj2':F(2,3)*M3*upper(mom['U_J2']),
                'r2_U2':r*r*D2*upper(mom['U2'])/2,'r_U2_Tj':r*D31*u2j/3,
                'U2_Tjj':M3*upper(mom['U2_JJ'])/6,'U2_Tj2':M4*upper(mom['U2_J2'])/12}
        result[str(kind)]={'components_upper':{k:str(v) for k,v in components.items()},
            'error_upper':str(sum(components.values(),F(0))),'cauchy_U2_Tj_upper':str(u2j),
            'cauchy_variant_a_upper':str(upper(u2j_a)),'cauchy_variant_b_upper':str(upper(u2j_b))}
    return result,{'M2':str(M2),'M3':str(M3),'M4':str(M4),'C2_prime':str(D1),
        'C2_second':str(D2),'C3_prime':str(D31)}

def tail_moment(model,L=6):
    w,d,s,mg,rg,lam=parameters(model,False);p=len(rg)
    U2=geometry_polynomials(model,0)['U2'];value=arb(0);terms=[]
    for alpha,coefficient in U2.items():
        rates=[sum(F(a)*si[j] for a,si in zip(alpha,s)) for j in range(p)]
        constant=sum(F(a)*di for a,di in zip(alpha,d));full=ball(constant+sum(v*v for v in rates)/2).exp()
        retained=full
        for rate in rates:
            retained*=Phi(ball(F(L)-rate))-Phi(ball(-F(L)-rate))
        term=ball(coefficient)*(full-retained);value+=term
        terms.append({'alpha':list(alpha),'coefficient':str(coefficient),'outside_tilt_interval':serial(full-retained)})
    return {'U2_outside_interval':serial(positive(value,'outside U2 moment')),'terms':terms}

def financial_transfer(model,cap,moments,tail):
    norms=payoff_norms(F(model['lambda']),cap);coord=[];linear=quadratic=F(0)
    for row in moments['coordinates']:
        e,n=derivative_envelope(norms,row);linear+=F(e['1']['error_upper']);quadratic+=F(e['2']['error_upper'])
        coord.append({'coordinate':row['coordinate'],'first_transport':e['1'],'second_transport':e['2']})
    tail_upper=F(norms['M2_upper'])*F(tail['U2_outside_interval']['upper_exact_fraction'])/2
    total=linear+quadratic+tail_upper
    return {'model_id':model['model_id'],'payoff':norms['payoff'],'cap':None if cap is None else str(cap),
        'order':2,'linear_transport_upper':str(linear),'quadratic_transport_upper':str(quadratic),
        'tail_upper':str(tail_upper),'financial_transfer_upper':str(total),'norms':n,
        'coordinates':coord,'target_allocation':'1/200',
        'status':'passed' if total<F(1,200) else 'retained allocation failure',
        'scope':'Price-independent full Gaussian to actual rounded tensor-grid transfer; no tensor enumeration, basket price or whole-box derivative sup.'}

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--b-y',type=int,default=10);ap.add_argument('--q',type=int,default=28)
    ap.add_argument('--out',type=Path,default=ROOT/'results/tensor_transfer.json');args=ap.parse_args()
    start=time.perf_counter();ctx.prec=PRECISION;transport=TransportMGF(6,args.b_y,args.q)
    models=json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']
    records=[];geometry=[];tails=[]
    for assets in (2,4,8):
        ms=[m for m in models if m['assets']==assets]
        g=HybridMoments(ms[0],transport).record();tail=tail_moment(ms[0],6)
        # mu differences and centred loadings are invariant across lambda in
        # the frozen protocol; prove this before reusing geometric evidence.
        p0=parameters(ms[0],False)
        for m in ms:
            pm=parameters(m,False)
            if pm[:3]!=p0[:3] or pm[4]!=p0[4]:raise ArithmeticError('Geometry reuse premise failed')
            for cap in (F(1),F(10),F(25),None):
                r=financial_transfer(m,cap,g,tail);records.append(r)
                print(m['model_id'],str(cap),float(F(r['financial_transfer_upper'])),r['status'],flush=True)
        geometry.append(g);tails.append({'assets':assets,**tail})
    result={'schema':'qbp-positive-transport-MGF-tensor-v1','status':'completed with all failures retained',
        'records':records,'geometry_moments':geometry,'tails':tails,'transport':transport.profile(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
            ('boundary_grid_certificate.py','boundary_derivative_certificate.py','higherdim_preparation.py')},
        'generation_seconds_diagnostic':time.perf_counter()-start,
        'priority':'A checked project-specific use of Taylor transport, exponential factorisation and Cauchy; no originality priority claim.'}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print('saved',args.out,'seconds',time.perf_counter()-start,flush=True)

if __name__=='__main__':main()
