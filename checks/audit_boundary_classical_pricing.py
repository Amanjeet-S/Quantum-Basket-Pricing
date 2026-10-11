"""Portable independent high-precision audit of every classical summand.

The producer's kernels and pricing functions are not imported. The reference
uses erfc CDFs, call-spread payoffs and direct quadratic-field constants.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,sys,time
from flint import arb,ctx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from continuous_quadrature import frac_arb,exact_endpoint

ROOT=Path(__file__).resolve().parents[1]
PRECISION=768

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def rational_ends(a):return F(exact_endpoint(a,False)),F(exact_endpoint(a,True))

def interval(row):
    lo,hi=F(row['lower_exact_fraction']),F(row['upper_exact_fraction']);mid=(lo+hi)/2;rad=(hi-lo)/2
    return frac_arb(mid)+arb(0,frac_arb(rad).upper())

def contained(a,row):
    lo,hi=rational_ends(a)
    assert F(row['lower_exact_fraction'])<=lo<=hi<=F(row['upper_exact_fraction'])

def Phi(x):return (-x/arb(2).sqrt()).erfc()/2

def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()

def own_geometry(model):
    w=list(map(F,model['weights']));mu=list(map(F,model['log_drifts']));rr=[list(map(F,r)) for r in model['loadings']]
    assert len(rr[0])==1 and sum(w)==1
    lam=F(model['lambda']);mg=sum(a*b for a,b in zip(w,mu));rg=sum(w[i]*rr[i][0] for i in range(len(w)))
    d=[a-mg for a in mu];s=[r[0]-rg for r in rr]
    assert sum(wi*di for wi,di in zip(w,d))==0 and sum(wi*si for wi,si in zip(w,s))==0
    for a,r in zip(mu,rr):assert a==F(model['risk_free_rate'])-F(model['dividend_yield'])-(lam*lam+r[0]*r[0])/2
    return w,d,s,mg,rg,lam

def financial(model,y):
    w,d,s,mg,rg,lam=own_geometry(model)
    T=sum((frac_arb(wi)*frac_arb(di+si*y).exp() for wi,di,si in zip(w,d,s)),arb(0))
    return frac_arb(mg+rg*y),T-1,T,lam

def cashgamma(m,lam,cap):
    L=frac_arb(lam);value=100/L*phi(m/L)
    if cap is not None:
        c=frac_arb(1+F(cap,100));value-=100/L*c*phi((m-c.log())/L)
    return value

def call_and_slope(m,lam,strike):
    L=frac_arb(lam);x=(m-frac_arb(F(strike,100)).log())/L
    slope=100*(m+L*L/2).exp()*Phi(x+L)
    return slope-frac_arb(strike)*Phi(x),slope

def payoff_pair(m,lam,cap):
    first,slope=call_and_slope(m,lam,100)
    if cap is not None:
        second,slope2=call_and_slope(m,lam,100+cap);first-=second;slope-=slope2
    return first,slope


def main():
    start=time.perf_counter();ctx.prec=PRECISION
    target=ROOT/'results/boundary_classical_pricing.json.gz';data=json.load(gzip.open(target,'rt'))
    assert data['source_sha256']==sha(ROOT/'boundary_classical_pricing.py')
    for path,digest in data['dependency_sha256'].items():assert sha(ROOT/path)==digest
    protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text());models={m['model_id']:m for m in protocol['models']}
    grid=json.loads((ROOT/'results/boundary_financial_grid_transfer_summary.json').read_text())
    nodes=json.loads((ROOT/'results/original_q28_quantile_nodes.json').read_text())['certificate']['nodes']
    moments=json.loads((ROOT/'results/boundary_finite_moments.json').read_text())['records']
    norms=json.loads((ROOT/'results/boundary_derivative_norms.json').read_text())['records']
    root=arb(6).sqrt();ideals={'z_minus':(4-root)/10,'z_plus':(4+root)/10,
        'w_minus':arb(1)/2+1/(3*root),'w_plus':arb(1)/2-1/(3*root)}
    actual={k:F(v['actual']) for k,v in data['gauss_constants'].items() if k in ideals}
    assert actual['w_minus']+actual['w_plus']==1
    for k,v in ideals.items():
        lo,hi=rational_ends(v);record=data['gauss_constants'][k]
        assert F(record['ideal_lower'])<=lo<=hi<=F(record['ideal_upper'])
        assert max(abs(actual[k]-lo),abs(actual[k]-hi))<=F(record['error_upper'])
        assert (actual[k]*(1<<28)).denominator==1
    profiles=[];summands_checked=0
    for row in data['records']:
        model=models[row['model_id']];cap=row['cap'];kind=row['kind'];payoff=row['payoff'];saved=data['summand_archive'][row['record_id']]
        total=arb(0)
        for node,record in zip(nodes,saved):
            assert node['j']==record['j'] and node['actual_coordinate']==record['actual_coordinate']
            y=F(node['actual_coordinate']);assert y==F(node['offset_integer'],1<<28)-6
            m,U,T,lam=financial(model,y);assert U.lower()>=0
            if kind=='stable_gauss_r2':
                weighted=arb(0)
                for suffix in ('minus','plus'):
                    t=1+frac_arb(actual['z_'+suffix])*U
                    weighted+=frac_arb(actual['w_'+suffix])*cashgamma(m+t.log(),lam,cap)/(t*t)
                value=U*U*weighted/2
            else:
                base,prime=payoff_pair(m,lam,cap);arithmetic=payoff_pair(m+T.log(),lam,cap)[0]
                value=arithmetic-base-U*prime
            contained(value,record['residual_interval']);total+=value;summands_checked+=1
        mean=total/len(nodes);contained(mean,row['finite_residual_mean_interval'])
        g=next(r for r in grid['records'] if r['model_id']==row['model_id'] and r['cap']==cap and r['payoff']==payoff and r['order']==2)
        coeff=next(r for r in grid['coefficients'] if r['model_id']==row['model_id'] and r['cap']==cap and r['payoff']==payoff and r['order']==1)
        assert row['financial_transfer_upper']==g['financial_transfer_upper'] and row['coefficient_interval']==coeff['coefficient']
        aux=row['auxiliary_certificate'];upperaux=F(0)
        if kind=='stable_gauss_r2':
            n=next(r for r in norms if r['payoff']==payoff and F(r['lambda'])==lam and (cap is None or F(r['cap'])==cap))
            mom=next(r for r in moments if r['model_id']==row['model_id'])['centred_U_moments']
            mu2=F(mom['2']['upper_exact_fraction']);mu3=F(mom['3']['upper_exact_fraction']);mu6=F(mom['6']['upper_exact_fraction'])
            trunc=F(n['C_r_upper']['6'])*mu6/14400
            nodecharge=sum(F(data['gauss_constants']['w_'+s]['ideal_upper'])*F(data['gauss_constants']['z_'+s]['error_upper']) for s in ('minus','plus'))
            weightcharge=sum(F(data['gauss_constants']['w_'+s]['error_upper']) for s in ('minus','plus'))
            fixed=F(n['C_r_upper']['3'])*mu3*nodecharge/2+F(n['M2_upper'])*mu2*weightcharge/2
            assert F(aux['Gauss_remainder_upper'])==trunc and F(aux['fixed_parameter_bias_upper'])==fixed
            upperaux=trunc+fixed
        assert F(aux['bias_upper'])==upperaux
        masshi=F(g['mass_interval']['upper_exact_fraction']);transfer=F(g['financial_transfer_upper'])
        # Reconstruct the full price interval from retained producer numerical
        # bounds, without importing its aggregation or interval constructor.
        mlo,mhi=map(F,(row['finite_residual_mean_interval']['lower_exact_fraction'],row['finite_residual_mean_interval']['upper_exact_fraction']))
        clo,chi=map(F,(coeff['coefficient']['lower_exact_fraction'],coeff['coefficient']['upper_exact_fraction']))
        alo,ahi=map(F,(g['mass_interval']['lower_exact_fraction'],g['mass_interval']['upper_exact_fraction']))
        products=[a*b for a in (alo,ahi) for b in (mlo,mhi)]
        prelo,prehi=clo+min(products),chi+max(products)
        penalty=transfer+masshi*upperaux
        unlo,unhi=prelo-penalty,prehi+penalty
        reported=row['undiscounted_price_interval'];assert F(reported['lower_exact_fraction'])<=unlo and unhi<=F(reported['upper_exact_fraction'])
        discount=frac_arb(-F(model['risk_free_rate'])*F(model['maturity'])).exp();contained(discount,row['discount_factor_interval'])
        dlow,dhigh=rational_ends(discount);prices=[u*d for u in (unlo,unhi) for d in (dlow,dhigh)]
        final=row['price_interval'];assert F(final['lower_exact_fraction'])<=min(prices) and max(prices)<=F(final['upper_exact_fraction'])
        lo,hi=F(final['lower_exact_fraction']),F(final['upper_exact_fraction'])
        assert F(row['price_midpoint'])==(lo+hi)/2 and F(row['price_radius_upper'])==(hi-lo)/2<=F(1,100)
        N=len(nodes);counts=row['operation_counts']
        if kind=='stable_gauss_r2':
            assert counts['cash_Gamma_kernel_calls']==2*N
            assert counts['Gaussian_density_calls']==(4 if cap is not None else 2)*N
            assert counts.get('Gaussian_CDF_calls',0)==0
            assert counts['logarithm_calls']==2*N+int(cap is not None)
            assert counts['exponential_calls']==counts['Gaussian_density_calls']
        else:
            assert counts['conditional_price_calls']==2*N and counts['conditional_slope_calls']==N
            assert counts['shared_price_slope_kernel_calls']==N and counts['exponential_calls']==2*N
            assert counts['Gaussian_CDF_calls']==(8 if cap is not None else 4)*N
            assert counts['logarithm_calls']==N+int(cap is not None)
        assert counts['finite_summands']==N and counts['outer_mean_divisions']==1
        profiles.append({'record_id':row['record_id'],'summands':N,'radius_upper':row['price_radius_upper'],'status':'passed'})
    # On the exact target the two independently bounded schemes must overlap.
    overlap=0
    for row in data['records']:
        if row['kind']!='stable_gauss_r2':continue
        other=next(r for r in data['records'] if r['model_id']==row['model_id'] and r['payoff']==row['payoff'] and r['cap']==row['cap'] and r['kind']=='raw_direct_r2')
        a,b=row['price_interval'],other['price_interval']
        assert max(F(a['lower_exact_fraction']),F(b['lower_exact_fraction']))<=min(F(a['upper_exact_fraction']),F(b['upper_exact_fraction']));overlap+=1
    result={'schema':'qbp-boundary-classical-independent-audit-v1','status':'passed','reference_precision_bits':PRECISION,
        'producer_source_sha256':data['source_sha256'],'pricing_archive_sha256':sha(target),
        'audit_source_sha256':sha(__file__),'summands_checked':summands_checked,'price_intervals_checked':len(profiles),
        'GJ2_raw_price_interval_overlaps_checked':overlap,'actual_operation_ledgers_checked':len(profiles),
        'records':profiles,'runtime_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Independent768-bit erfc/call-spread and direct Gaussian references; exact rational aggregation of every price allowance and operation-count checks. No producer pricing function or external basket price is used.'}
    (ROOT/'results/boundary_classical_pricing_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Passed',summands_checked,'summands and',len(profiles),'full price intervals')

if __name__=='__main__':main()
