"""Executed finite-error degree selection versus complete order-eight search.

Both policies use the same validated Gaussian coefficient identities. Every
failed degree's moment work is charged. The screening policy evaluates scalar
coefficients only after its remainder bound permits a certificate, or at the
fixed order limit to retain an unresolved interval. No infinite-series or
speed-optimality claim is made.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial,prod
import hashlib,json,time
from flint import arb,ctx
from taylor_pricing import model_parameters,compositions,dot,kernel_derivatives,signed_stirling,Fourier_constants
from continuous_quadrature import frac_arb,exact_endpoint,serialise_interval
ROOT=Path(__file__).resolve().parent
ORDER_LIMIT=8;PRECISION=224;TOLERANCE=F(1,100)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Geometry:
    def __init__(self,parameters):
        self.w,self.d,self.s,self.rg=parameters;self.rows=[];self.raw=[];self.centred=[];self.exponentials=0;self.setup_seconds=0.0
    def ensure(self,degree):
        start=time.perf_counter()
        while len(self.rows)<=degree:
            j=len(self.rows);terms=[]
            for alpha in compositions(j,len(self.w)):
                weight=F(factorial(j),prod(factorial(a) for a in alpha))
                for wi,a in zip(self.w,alpha):weight*=wi**a
                exposure=tuple(sum(alpha[i]*self.s[i][v] for i in range(len(self.w))) for v in range(len(self.rg)))
                exponent=dot(alpha,self.d)+dot(exposure,exposure)/2;shift=dot(self.rg,exposure)
                factor=frac_arb(weight)*frac_arb(exponent).exp();self.exponentials+=1
                terms.append((factor,shift))
            assert len(terms)==comb(len(self.w)+j-1,j)
            self.rows.append(terms);self.raw.append(sum((x for x,_ in terms),arb(0)))
            value=sum(((-1)**(j-r))*comb(j,r)*self.raw[r] for r in range(j+1));assert exact_endpoint(value,True)>=0
            self.centred.append(value)
        self.setup_seconds+=time.perf_counter()-start

def coefficients(model,cap,k,bank):
    w,d,s,mg,rg,lam=model_parameters(model);sigma=frac_arb(lam*lam+dot(rg,rg)).sqrt();base=frac_arb(mg);sums=[];calls=0
    for j in range(k+1):
        values=[arb(0) for _ in range(k+1)]
        for factor,shift in bank.rows[j]:
            derivatives=kernel_derivatives(base+frac_arb(shift),sigma,cap,k);calls+=1
            for ell in range(k+1):values[ell]+=factor*derivatives[ell]
        sums.append(values)
    table=signed_stirling(k);out=[sums[0][0]]
    for degree in range(1,k+1):
        c=arb(0)
        for ell in range(1,degree+1):
            term=sum(((-1)**(degree-j))*comb(degree,j)*sums[j][ell] for j in range(degree+1));c+=table[degree][ell]*term
        out.append(c/factorial(degree))
    return out,calls

def radius(bank,lam,cap,k):
    bank.ensure(k+1);table=signed_stirling(k+1);bounds=Fourier_constants(k+1)
    constant=sum(abs(table[k+1][j])*bounds[j-1]/frac_arb(lam)**j for j in range(1,k+2))
    ball=cap*constant*bank.centred[k+1]/factorial(k+1)
    return max(F(0),exact_endpoint(ball,True))

def interval(partial,remainder,cap):
    lo=max(F(0),exact_endpoint(partial,False)-remainder);hi=min(F(cap),exact_endpoint(partial,True)+remainder);assert lo<=hi
    return lo,hi

def evaluate(model,cap,policy,bank):
    start=time.perf_counter();before_exp=bank.exponentials;before_setup=bank.setup_seconds
    lam=model_parameters(model)[-1];screens=[];orders=[];calls=0
    if policy=='monotone_then_screen':
        bank.ensure(1);w,d,ss,mg,rg,lam=model_parameters(model)
        phi_slope=1/(2*arb.pi()).sqrt()/frac_arb(lam)
        L_upper=min(exact_endpoint(phi_slope,True),F(100+cap,cap))
        upper=cap*frac_arb(L_upper)*bank.raw[1].log();gap=max(F(0),exact_endpoint(upper,True))
        simple=None
        # Sufficient radius test before spending even the scalar control call.
        if gap<=2*TOLERANCE:
            sigma=frac_arb(lam*lam+dot(rg,rg)).sqrt();g0=kernel_derivatives(frac_arb(mg),sigma,cap,0)[0];calls+=1
            lo=max(F(0),exact_endpoint(g0,False));hi=min(F(cap),exact_endpoint(g0,True)+gap);assert lo<=hi
            simple={'order':0,'remainder_upper':str(gap),'interval':serialise_interval(lo,hi),'certified':(hi-lo)/2<=TOLERANCE,'one_sided_monotonicity':True}
            if simple['certified']:
                return {'case_id':f"{model['model_id']}_cap{cap}",'assets':model['assets'],'cap':cap,'policy':policy,'selected':simple,'screens':[],'evaluated_orders':[simple],'kernel_calls':calls,'new_geometry_exponentials':bank.exponentials-before_exp,'geometry_setup_seconds_diagnostic':bank.setup_seconds-before_setup,'total_seconds_diagnostic':time.perf_counter()-start,'decision':'one-sided control certificate before estimation'}
        fallback=evaluate(model,cap,'remainder_screen',bank);fallback['policy']=policy;fallback['kernel_calls']+=calls;fallback['new_geometry_exponentials']=bank.exponentials-before_exp;fallback['geometry_setup_seconds_diagnostic']=bank.setup_seconds-before_setup;fallback['total_seconds_diagnostic']=time.perf_counter()-start;fallback['one_sided_preliminary_interval']=simple;fallback['one_sided_screen_upper_exact_fraction']=str(gap);return fallback
    if policy=='full_search':
        radii=[radius(bank,lam,cap,k) for k in range(ORDER_LIMIT+1)];coefs,n=coefficients(model,cap,ORDER_LIMIT,bank);calls+=n;partial=arb(0)
        for k in range(ORDER_LIMIT+1):
            partial+=coefs[k];lo,hi=interval(partial,radii[k],cap);orders.append({'order':k,'remainder_upper':str(radii[k]),'interval':serialise_interval(lo,hi),'certified':(hi-lo)/2<=TOLERANCE})
        good=[o for o in orders if o['certified']];selected=min(good or orders,key=lambda o:F(o['interval']['upper_exact_fraction'])-F(o['interval']['lower_exact_fraction']))
    elif policy=='remainder_screen':
        for k in range(ORDER_LIMIT+1):
            rem=radius(bank,lam,cap,k);qualified=rem<=TOLERANCE;screens.append({'order':k,'remainder_upper':str(rem),'remainder_qualified':qualified,'coefficient_evaluated':qualified or k==ORDER_LIMIT})
            if not qualified and k<ORDER_LIMIT:continue
            coefs,n=coefficients(model,cap,k,bank);calls+=n;partial=sum(coefs,arb(0));lo,hi=interval(partial,rem,cap);item={'order':k,'remainder_upper':str(rem),'interval':serialise_interval(lo,hi),'certified':(hi-lo)/2<=TOLERANCE};orders.append(item)
            if item['certified']:break
        good=[o for o in orders if o['certified']];selected=good[0] if good else orders[-1]
    else:raise ValueError('Unknown finite-order policy')
    return {'case_id':f"{model['model_id']}_cap{cap}",'assets':model['assets'],'cap':cap,'policy':policy,'selected':selected,'screens':screens,'evaluated_orders':orders,'kernel_calls':calls,'new_geometry_exponentials':bank.exponentials-before_exp,'geometry_setup_seconds_diagnostic':bank.setup_seconds-before_setup,'total_seconds_diagnostic':time.perf_counter()-start,'decision':'certificate before estimation' if selected['certified'] else 'order limit unresolved'}

def main():
    ctx.prec=PRECISION;protocol_path=ROOT/'results/generalisation_protocol.json';protocol=json.loads(protocol_path.read_text());result=[]
    for policy in ('full_search','remainder_screen','monotone_then_screen'):
        banks={};records=[];start=time.perf_counter()
        for model in protocol['models']:
            w,d,s,mg,rg,lam=model_parameters(model);key=(w,d,s,rg)
            bank=banks.setdefault(key,Geometry(key))
            for cap in protocol['caps']:records.append(evaluate(model,cap,policy,bank))
        result.append({'policy':policy,'execution_scope':'one shared frozen27-contract family, fresh geometry bank for this policy','configurations':records,'kernel_calls':sum(r['kernel_calls'] for r in records),'geometry_exponentials':sum(b.exponentials for b in banks.values()),'certified_cases':sum(r['selected']['certified'] for r in records),'total_seconds_diagnostic':time.perf_counter()-start,'geometry_setup_seconds_diagnostic':sum(b.setup_seconds for b in banks.values())})
    cold=[]
    for model in protocol['models']:
        if model['assets']!=8:continue
        for policy in ('full_search','remainder_screen','monotone_then_screen'):
            w,d,s,mg,rg,lam=model_parameters(model);bank=Geometry((w,d,s,rg));out=evaluate(model,25,policy,bank);out['execution_scope']='cold standalone cap25; own geometry bank; no borrowed setup';cold.append(out)
    fixed={r['case_id']:r for r in result[0]['configurations']};checks=[]
    for run in result[1:]:
     for r in run['configurations']:
         a=fixed[r['case_id']]['selected']['interval'];b=r['selected']['interval'];assert max(F(a['lower_exact_fraction']),F(b['lower_exact_fraction']))<=min(F(a['upper_exact_fraction']),F(b['upper_exact_fraction']))
         checks.append({'case_id':r['case_id'],'enclosures_overlap':True,'certified_status_preserved':r['selected']['certified']==fixed[r['case_id']]['selected']['certified']})
    out={'status':'executed finite-order selection and full-search comparison','protocol_sha256':sha(protocol_path),'source_hashes':{n:sha(ROOT/n) for n in ['taylor_order_benchmark.py','taylor_pricing.py','continuous_quadrature.py']},'maximum_order':ORDER_LIMIT,'precision_bits':PRECISION,'undiscounted_radius_target':str(TOLERANCE),'family_runs':result,'cold_standalone_presets':cold,'exact_comparison_checks':checks,'summary':{'screened_to_full_kernel_ratio':str(F(result[1]['kernel_calls'],result[0]['kernel_calls'])),'kernel_calls_saved':result[0]['kernel_calls']-result[1]['kernel_calls'],'hierarchy_to_full_kernel_ratio':str(F(result[2]['kernel_calls'],result[0]['kernel_calls'])),'monotone_screen_certificates':sum(r['selected'].get('one_sided_monotonicity',False) for r in result[2]['configurations']),'certified_status_unchanged':all(r['certified_status_preserved'] for r in checks)},'scope_limit':'standard finite error-controlled selection; no optimality, series convergence, market validation or physical quantum speedup inferred'}
    (ROOT/'results/taylor_order_benchmark.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'summary':out['summary'],'family_seconds':[(x['policy'],x['total_seconds_diagnostic']) for x in result]},indent=2))
if __name__=='__main__':main()
