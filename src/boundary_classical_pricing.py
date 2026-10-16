"""Strong certified classical sums for the identical boundary pricing target.

Only price-independent model, grid, geometric coefficient, norm and moment
certificates are inputs. No continuous arithmetic-basket price file is read.
The primary comparator uses the actual q28 two-node Gauss constants.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from datetime import datetime,timezone
import gzip,hashlib,json,time
from flint import arb,ctx
from continuous_quadrature import frac_arb,serialise_interval
from integral_value_oracle import model_parameters
from boundary_auxiliary_certificate import auxiliary_bias
from boundary_value_certificate import gauss_constants,endpoints
from boundary_derivative_certificate import payoff_norms

ROOT=Path(__file__).resolve().parent
PRECISION=512

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def ball_interval(lo,hi):
    lo,hi=map(F,(lo,hi));mid=(lo+hi)/2;rad=(hi-lo)/2
    if lo>hi:raise ValueError('Reversed interval')
    return frac_arb(mid)+arb(0,frac_arb(rad).upper())

def data_interval(row):return ball_interval(row['lower_exact_fraction'],row['upper_exact_fraction'])

def serial(a):return serialise_interval(*endpoints(a))

def midpoint_radius(row):
    lo,hi=F(row['lower_exact_fraction']),F(row['upper_exact_fraction'])
    return (lo+hi)/2,(hi-lo)/2

class Kernels:
    """Instrumented actual transcendental calls, not inferred CPU operations."""
    def __init__(self,lam,cap=None):
        self.lam=frac_arb(lam);self.cap=None if cap is None else frac_arb(cap);self.counts=Counter()
        self.counts['square_root_calls']+=2
        self.sqrt2=arb(2).sqrt();self.sqrt2pi=(2*arb.pi()).sqrt()
        self.ell=None if cap is None else self.log(1+self.cap/100)
    def exp(self,x):self.counts['exponential_calls']+=1;return x.exp()
    def log(self,x):self.counts['logarithm_calls']+=1;return x.log()
    def Phi(self,x):self.counts['Gaussian_CDF_calls']+=1;return (1+(x/self.sqrt2).erf())/2
    def phi(self,x):self.counts['Gaussian_density_calls']+=1;return self.exp(-x*x/2)/self.sqrt2pi
    def C2(self,m):
        self.counts['cash_Gamma_kernel_calls']+=1
        result=100/self.lam*self.phi(m/self.lam)
        if self.cap is not None:result-=100/self.lam*(1+self.cap/100)*self.phi((m-self.ell)/self.lam)
        return result
    def call(self,m,prime=False):
        self.counts['conditional_slope_calls' if prime else 'conditional_price_calls']+=1
        growth=100*self.exp(m+self.lam*self.lam/2);x=m/self.lam
        if self.cap is None:
            forward=growth*self.Phi(x+self.lam)
            return forward if prime else forward-100*self.Phi(x)
        xh=(m-self.ell)/self.lam
        forward=growth*(self.Phi(x+self.lam)-self.Phi(xh+self.lam))
        if prime:return forward
        lower=self.Phi(xh)
        return forward-100*(self.Phi(x)-lower)+self.cap*lower

    def price_and_slope(self,m):
        self.counts['conditional_price_calls']+=1;self.counts['conditional_slope_calls']+=1
        self.counts['shared_price_slope_kernel_calls']+=1
        growth=100*self.exp(m+self.lam*self.lam/2);x=m/self.lam
        if self.cap is None:
            prime=growth*self.Phi(x+self.lam)
            return prime-100*self.Phi(x),prime
        xh=(m-self.ell)/self.lam
        prime=growth*(self.Phi(x+self.lam)-self.Phi(xh+self.lam));lower=self.Phi(xh)
        return prime-100*(self.Phi(x)-lower)+self.cap*lower,prime


def financial_grid(model,nodes):
    """Shared across payoffs and both residual estimators for each model."""
    w,d,s,mg,rg,lam=model_parameters(model)
    if len(rg)!=1:raise ValueError('The retained complete classical study is n2')
    start=time.perf_counter();rows=[];calls=0
    for node in nodes:
        y=F(node['actual_coordinate']);yb=frac_arb(y)
        T=sum((frac_arb(wi)*(frac_arb(di+si[0]*y)).exp() for wi,di,si in zip(w,d,s)),arb(0));calls+=len(w)
        m=frac_arb(mg+rg[0]*y);U=T-1
        rows.append((node['j'],y,m,U,T))
    return rows,{'financial_nodes':len(rows),'conditional_exponential_calls':calls,
        'setup_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Actual ideal finite-grid financial T and mG, shared by all payoffs and both classical estimators for the same model.'}


def residual_sum(rows,lam,cap,kind,constants):
    kernel=Kernels(lam,cap);summands=[];total=arb(0);start=time.perf_counter()
    zm,zp=(frac_arb(F(constants[k]['actual'])) for k in ('z_minus','z_plus'))
    wm,wp=(frac_arb(F(constants[k]['actual'])) for k in ('w_minus','w_plus'))
    for j,y,m,U,T in rows:
        if kind=='stable_gauss_r2':
            partial=arb(0)
            for z,weight in ((zm,wm),(zp,wp)):
                t=1+z*U;partial+=weight*kernel.C2(m+kernel.log(t))/(t*t)
            value=U*U*partial/2
        elif kind=='raw_direct_r2':
            base,prime=kernel.price_and_slope(m)
            value=kernel.call(m+kernel.log(T))-base-U*prime
        else:raise ValueError('Unknown classical estimator')
        total+=value;summands.append({'j':j,'actual_coordinate':str(y),'residual_interval':serial(value)})
    mean=total/len(rows)
    counts=dict(kernel.counts);counts['finite_summands']=len(rows);counts['outer_mean_divisions']=1
    return mean,summands,counts,time.perf_counter()-start


def price_record(model,cap,payoff,kind,finite_mean,grid,coefficient,aux,normaliser,counts,seconds):
    ctx.prec=PRECISION
    coeff=data_interval(coefficient['coefficient']);mass=data_interval(grid['mass_interval'])
    pre=coeff+mass*finite_mean;prelo,prehi=endpoints(pre)
    transfer=F(grid['financial_transfer_upper']);masshi=F(grid['mass_interval']['upper_exact_fraction'])
    auxiliary=masshi*F(aux['bias_upper']);error=transfer+auxiliary
    enlarged=ball_interval(prelo-error,prehi+error)
    discount=(-frac_arb(F(model['risk_free_rate'])*F(model['maturity']))).exp()
    final=discount*enlarged;lo,hi=endpoints(final);point=(lo+hi)/2;radius=(hi-lo)/2
    cm,cr=midpoint_radius(coefficient['coefficient']);mm,mr=midpoint_radius(grid['mass_interval'])
    residual=serial(finite_mean);rm,rr=midpoint_radius(residual)
    pre_radius=(prehi-prelo)/2
    input_radius=cr+masshi*rr+abs(rm)*mr
    return {'model_id':model['model_id'],'payoff':payoff,'cap':cap,'kind':kind,
        'configuration':{'L':'6','b_y':10,'q':28,'finite_nodes':1024,'precision_bits':PRECISION,
            'classical_kernel_arithmetic':'Outward Arb arithmetic; q28 specifies the actual financial coordinates and Gauss constants, not a restriction to q28 kernel evaluation.',
            'shared_M2_upper':str(normaliser),'B2':str(F(normaliser)*F(grid['_U_star'])**2/2)},
        'finite_residual_mean_interval':residual,'coefficient_interval':coefficient['coefficient'],
        'retained_mass_interval':grid['mass_interval'],'auxiliary_certificate':aux,
        'financial_transfer_upper':str(transfer),'mass_times_auxiliary_upper':str(auxiliary),
        'coefficient_radius_upper':str(cr),'finite_numerical_mean_radius_upper':str(rr),
        'mass_uncertainty_contribution_upper':str(abs(rm)*mr),
        'coefficient_mass_and_finite_numerical_radius_upper':str(pre_radius),
        'input_interval_propagation_radius_upper':str(input_radius),
        'additional_outward_Arb_product_radius_upper':str(max(F(0),pre_radius-input_radius)),
        'undiscounted_point':str((prelo+prehi)/2),'undiscounted_price_interval':serial(enlarged),
        'discount_factor_interval':serial(discount),'price_interval':serialise_interval(lo,hi),
        'price_midpoint':str(point),'price_radius_upper':str(radius),'price_tolerance':'1/100',
        'status':'certified' if radius<=F(1,100) else 'retained tolerance failure',
        'operation_counts':counts,'finite_sum_seconds_diagnostic':seconds,
        'scope':'Deterministic continuous risk-neutral price interval after transfer, coefficient, retained mass, quadrature and numerical allowances. No continuous price answer is an input.'}


def run():
    start=time.perf_counter();ctx.prec=PRECISION
    dependencies=['boundary_classical_pricing.py','boundary_auxiliary_certificate.py','boundary_value_certificate.py',
        'boundary_derivative_certificate.py','integral_value_oracle.py','continuous_quadrature.py',
        'results/generalisation_protocol.json','results/original_q28_quantile_nodes.json',
        'results/boundary_financial_grid_transfer_summary.json','results/boundary_grid_independent_audit.json',
        'results/boundary_derivative_norms.json','results/boundary_finite_moments.json',
        'results/integral_finite_input_domains.json']
    provenance={p:sha(ROOT/p) for p in dependencies}
    protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text())
    original=json.loads((ROOT/'results/original_q28_quantile_nodes.json').read_text())['certificate']
    if (F(original['L']),original['b_y'],original['q'])!=(F(6),10,28):raise ValueError('Wrong finite input law')
    nodes=original['nodes']
    if len(nodes)!=1024 or any(F(n['actual_coordinate'])!=F(n['offset_integer'],1<<28)-6 for n in nodes):raise ValueError('Broken actual quantile integer encoding')
    summary=json.loads((ROOT/'results/boundary_financial_grid_transfer_summary.json').read_text())
    audit=json.loads((ROOT/'results/boundary_grid_independent_audit.json').read_text())
    if audit['status']!='passed' or summary['full_cell_certificate']['sha256']!=audit['certificate_sha256']:raise ValueError('Grid summary lacks its audited full certificate')
    moments=json.loads((ROOT/'results/boundary_finite_moments.json').read_text())
    domains=json.loads((ROOT/'results/integral_finite_input_domains.json').read_text())['domains']
    constants=gauss_constants(28);ctx.prec=PRECISION
    records=[];summand_archive={};setup=[]
    for model in protocol['models']:
        if model['outer_dimension']!=1:continue
        rows,financial_setup=financial_grid(model,nodes);financial_setup['model_id']=model['model_id'];setup.append(financial_setup)
        domain=next(d for d in domains if d['model_id']==model['model_id']);moment=next(d for d in moments['records'] if d['model_id']==model['model_id'])
        mu=moment['centred_U_moments']
        for payoff,cap in [('capped',1),('capped',10),('capped',25),('vanilla',None)]:
            grid=dict(next(r for r in summary['records'] if r['model_id']==model['model_id'] and r['payoff']==payoff and r['cap']==cap and r['order']==2));grid['_U_star']=domain['U_star']
            coefficient=next(r for r in summary['coefficients'] if r['model_id']==model['model_id'] and r['payoff']==payoff and r['cap']==cap and r['order']==1)
            norm=payoff_norms(model['lambda'],cap)
            for kind in ('stable_gauss_r2','raw_direct_r2'):
                aux=auxiliary_bias(model['lambda'],cap,domain['U_star'],mu['2']['upper_exact_fraction'],mu['3']['upper_exact_fraction'],
                    kind=kind,MU4=mu['4']['upper_exact_fraction'],MU6=mu['6']['upper_exact_fraction'])
                ctx.prec=PRECISION
                mean,summands,counts,seconds=residual_sum(rows,F(model['lambda']),cap,kind,constants)
                record=price_record(model,cap,payoff,kind,mean,grid,coefficient,aux,norm['M2_upper'],counts,seconds)
                key=f"{model['model_id']}:{payoff}:{cap}:{kind}";record['record_id']=key
                records.append(record);summand_archive[key]=summands
                print(key,record['status'],'price',float(F(record['price_midpoint'])),'radius',float(F(record['price_radius_upper'])),flush=True)
    if {p:sha(ROOT/p) for p in dependencies}!=provenance:raise RuntimeError('Input changed during classical pricing')
    grid_rows=[r for r in summary['records'] if r['order']==2]
    setup_ledger={'shared_financial_grid':setup,'new_financial_exponential_calls_total':sum(s['conditional_exponential_calls'] for s in setup),
        'retained_grid_certificate_work':{'quantile_inverse_evaluations':summary['classical_setup_ledger']['quantile_inverse_evaluations'],
            'R2_nodal_residual_derivative_evaluations':sum(r['cell_count'] for r in grid_rows),
            'R2_segment_residual_derivative_interval_evaluations':sum(r['derivative_subinterval_count'] for r in grid_rows),
            'geometric_coefficient_certificates':12,'full_grid_certificate_generation_seconds_diagnostic':summary['elapsed_seconds_diagnostic']},
        'retained_moment_work':{'conditional_exponential_evaluations':sum(r['classical_conditional_exponential_evaluations'] for r in moments['records']),
            'generation_seconds_diagnostic':sum(r['generation_seconds_diagnostic'] for r in moments['records'])},
        'scope':'Setup is retained and charged separately from the finite sum. Source-certificate construction counts include shared reusable work. Timing is diagnostic and cannot be converted to quantum T gates.'}
    result={'schema':'qbp-boundary-strong-classical-pricing-v1','status':'all24 deterministic intervals certified' if all(r['status']=='certified' for r in records) else 'contains retained tolerance failures',
        'created_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'dependency_sha256':provenance,
        'library_versions':{'python_flint':__import__('flint').__version__},'precision_bits':PRECISION,
        'gauss_constants':constants,'records':records,'summand_archive':summand_archive,'setup_ledger':setup_ledger,
        'financial_case_count':12,'estimator_count':2,'finite_sum_runs':24,
        'maximum_price_radius_upper':str(max(F(r['price_radius_upper']) for r in records)),
        'total_runtime_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Strong classical finite-sum competitors on exactly the retained quantum financial grid. Both methods price the original continuous model through price-independent transfer bounds; this is not empirical data or an externally known price sample.'}
    target=ROOT/'results/boundary_classical_pricing.json.gz'
    with gzip.open(target,'wt') as f:json.dump(result,f,separators=(',',':'));f.write('\n')
    brief={k:v for k,v in result.items() if k not in ('summand_archive','records')};brief['records']=[{k:v for k,v in r.items() if k not in ('auxiliary_certificate',)} for r in records];brief['archive_sha256']=sha(target)
    (ROOT/'results/boundary_classical_pricing_summary.json').write_text(json.dumps(brief,indent=2)+'\n')
    return result

if __name__=='__main__':run()
