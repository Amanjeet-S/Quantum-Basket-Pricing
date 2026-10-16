"""Matched first-correction transfer, sharing the positive transport theorem.

The second residual evidence is reused. Only the analytically known term
U P'(m) receives an additional transfer allowance. No price answer is read.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,time
from flint import arb,ctx
import tensor_transfer as tt
from boundary_grid_certificate import parameters,ball,endpoints,serial,upper,Phi
from boundary_derivative_certificate import payoff_norms
ROOT=Path(__file__).resolve().parent

def extra_polys(model,j):
    w,d,s,mg,rg,lam=parameters(model,False);n=len(w);zero=(0,)*n
    T={};J={};JJ={}
    for i in range(n):
        k=tuple(int(a==i) for a in range(n));T[k]=w[i]
        J[k]=w[i]*s[i][j];JJ[k]=w[i]*s[i][j]**2
    return {'mass':{zero:F(1)},'U':tt.add(T,{zero:F(-1)}),'J2':tt.mul(J,J),'JJ':JJ}

class ExtraMoments(tt.HybridMoments):
    def tilted_monomial(self,alpha,j,kind):
        key=('tilted',alpha,j,kind)
        if key not in self.exp_cache:
            constant=self.mg+self.lam*self.lam/2+sum(F(a)*d for a,d in zip(alpha,self.d))
            rates=tuple(self.rg[k]+sum(F(a)*si[k] for a,si in zip(alpha,self.s)) for k in range(self.p))
            value=ball(constant).exp()
            for k,rate in enumerate(rates):
                if k==j:factor=self.transport.rho1(rate) if kind==1 else self.transport.second(rate)
                elif k<j:factor=self.transport.continuous(rate)
                else:factor=self.transport.rounded(rate)
                value*=factor
            self.exp_cache[key]=value
        return self.exp_cache[key]

    def record_extra(self):
        out=[]
        for j in range(self.p):
            polys=extra_polys(self.model,j);kinds={}
            for kind in (1,2):
                row={}
                for name,poly in polys.items():
                    row[name]=serial(self.moment(poly,j,kind,name))
                    tilted=sum((ball(c)*self.tilted_monomial(a,j,kind) for a,c in poly.items()),arb(0))
                    row['tilted_'+name]=serial(tt.positive(tilted,'tilted '+name))
                kinds[str(kind)]=row
            out.append({'coordinate':j,'geometric_loading':str(self.rg[j]),'moments':kinds})
        return {'model_id':self.model['model_id'],'coordinates':out,'distinct_monomial_evaluations':len(self.exp_cache)}

def outside_U_tilt(model,L=6):
    w,d,s,mg,rg,lam=parameters(model,False);p=len(rg)
    poly=extra_polys(model,0)['U'];value=arb(0);terms=[]
    for alpha,c in poly.items():
        rates=[rg[j]+sum(F(a)*si[j] for a,si in zip(alpha,s)) for j in range(p)]
        constant=mg+lam*lam/2+sum(F(a)*di for a,di in zip(alpha,d))
        full=ball(constant+sum(v*v for v in rates)/2).exp();inside=full
        for v in rates:inside*=Phi(ball(F(L)-v))-Phi(ball(-F(L)-v))
        term=ball(c)*(full-inside);value+=term
        terms.append({'alpha':list(alpha),'coefficient':str(c),'tilted_outside_interval':serial(full-inside)})
    mass=(Phi(ball(L))-Phi(ball(-F(L))))**p
    return {'tilted_U_outside_interval':serial(tt.positive(value,'tilted U outside')),
        'outside_probability_interval':serial(1-mass),'terms':terms}

def transfer(model,cap,r2,extra,tail):
    norms=payoff_norms(F(model['lambda']),cap);M2=F(norms['M2_upper']);M3=F(norms['C_r_upper']['3'])
    linear=quadratic=F(0);rows=[]
    for row in extra['coordinates']:
        r=abs(F(row['geometric_loading']));out={}
        for kind in (1,2):
            m={k:tt.B(v) for k,v in row['moments'][str(kind)].items()}
            J=upper((m['mass']*m['J2']).sqrt())
            growthJ=upper((m['tilted_mass']*m['tilted_J2']).sqrt())
            if cap is not None:
                M1=F(norms['M1_upper']);P2=F(norms['scalar_cap_prime_first_upper']);P3=F(norms['scalar_cap_prime_second_upper'])
                if kind==1:components={'delta_Tj':M1*J,'delta_prime_U':r*P2*upper(m['U'])}
                else:components={'delta_Tjj':M1*upper(m['JJ']),'delta_prime_Tj':2*r*P2*J,'delta_second_U':r*r*P3*upper(m['U'])}
            elif kind==1:
                components={'forward_Tj':100*growthJ,'forward_r_U':100*r*upper(m['tilted_U']),
                    'boundary_r_U':r*M2*upper(m['U'])}
            else:
                components={'forward_Tjj':100*upper(m['tilted_JJ']),'forward_r_Tj':200*r*growthJ,
                    'boundary_r_Tj':2*r*M2*J,'forward_r2_U':100*r*r*upper(m['tilted_U']),
                    'boundary_r2_U':r*r*(3*M2+M3)*upper(m['U'])}
            amount=sum(components.values(),F(0));out[str(kind)]={'components_upper':{k:str(v) for k,v in components.items()},'additional_error_upper':str(amount)}
            if kind==1:linear+=amount
            else:quadratic+=amount
        rows.append({'coordinate':row['coordinate'],'extra_transfer':out})
    tail_upper=100*F(tail['tilted_U_outside_interval']['upper_exact_fraction'])
    if cap is not None:tail_upper=min(tail_upper,F(cap)*F(tail['outside_probability_interval']['upper_exact_fraction']))
    # The first-residual tail REPLACES the second-residual tail. Its retained
    # linear/quadratic transfer has the R2 allowances plus the U P' term.
    base_linear=F(r2['linear_transport_upper']);base_quad=F(r2['quadratic_transport_upper'])
    total=base_linear+linear+base_quad+quadratic+tail_upper
    return {'model_id':model['model_id'],'payoff':norms['payoff'],'cap':None if cap is None else str(cap),
      'order':1,'linear_transport_upper':str(base_linear+linear),'quadratic_transport_upper':str(base_quad+quadratic),
      'tail_upper':str(tail_upper),'financial_transfer_upper':str(total),'additional_U_delta_transfer':rows,
      'target_allocation':'1/200','status':'passed' if total<F(1,200) else 'retained allocation failure',
      'scope':'Matched first residual via R1=R2+U Pprime(m); its raw payoff may be unbounded but Gaussian exponential tails are charged.'}

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--b-y',type=int,default=10);ap.add_argument('--q',type=int,default=28)
    ap.add_argument('--second-input',type=Path,default=ROOT/'results/tensor_transfer.json')
    ap.add_argument('--out',type=Path,default=ROOT/'results/tensor_first_transfer.json');args=ap.parse_args()
    ctx.prec=512;start=time.perf_counter();transport=tt.TransportMGF(6,args.b_y,args.q)
    second=json.loads(args.second_input.read_text());records=[];extras=[];tails=[]
    law=second['transport']
    if F(law['L'])!=F(6) or law['b_y']!=args.b_y or law['q']!=args.q:
        raise ValueError('Second-correction input has a different retained financial law')
    models=json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']
    expected={(m['model_id'],cap) for m in models for cap in ('1','10','25',None)}
    actual=[(r['model_id'],r['cap']) for r in second['records']]
    if len(actual)!=len(expected) or set(actual)!=expected or any(r.get('order')!=2 for r in second['records']):
        raise ValueError('Second-correction input lacks unique complete frozen model/payoff/order coverage')
    for assets in (2,4,8):
        ms=[m for m in models if m['assets']==assets];extra=ExtraMoments(ms[0],transport).record_extra();tail=outside_U_tilt(ms[0])
        p0=parameters(ms[0],False)
        for model in ms:
            pp=parameters(model,False)
            if pp[:3]!=p0[:3] or pp[4]!=p0[4] or pp[3]+pp[5]**2/2!=p0[3]+p0[5]**2/2:raise ArithmeticError('Tilted geometry reuse failed')
            for cap in (F(1),F(10),F(25),None):
                c=None if cap is None else str(cap)
                r2=next(r for r in second['records'] if r['model_id']==model['model_id'] and r['cap']==c)
                record=transfer(model,cap,r2,extra,tail);records.append(record)
                print(model['model_id'],str(cap),float(F(record['financial_transfer_upper'])),record['status'],flush=True)
        extras.append(extra);tails.append({'assets':assets,**tail})
    result={'schema':'qbp-positive-transport-first-correction-v1','status':'completed with all failures retained',
       'records':records,'extra_moments':extras,'tails':tails,'transport':transport.profile(),
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'dependency_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('tensor_transfer.py','higherdim_preparation.py')},
       'second_transfer_sha256':hashlib.sha256(args.second_input.read_bytes()).hexdigest(),
       'generation_seconds_diagnostic':time.perf_counter()-start}
    args.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
