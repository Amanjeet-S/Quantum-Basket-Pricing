"""Portable independent768-bit audit of continuous nine-model call screens."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
from flint import arb,ctx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from continuous_quadrature import frac_arb,exact_endpoint

ROOT=Path(__file__).resolve().parents[1]
PRECISION=768

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def ends(x):return F(exact_endpoint(x,False)),F(exact_endpoint(x,True))

def inside(x,record):
    lo,hi=ends(x)
    assert F(record['lower_exact_fraction'])<=lo<=hi<=F(record['upper_exact_fraction'])

def Phi(x):return (-x/arb(2).sqrt()).erfc()/2


def reference(model):
    w=list(map(F,model['weights']));r=[list(map(F,row)) for row in model['loadings']];mu=list(map(F,model['log_drifts']));lam=F(model['lambda'])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b));p=len(r[0]);n=len(w)
    mg=sum(wi*mi for wi,mi in zip(w,mu));rg=[sum(w[i]*r[i][j] for i in range(n)) for j in range(p)]
    residual=[tuple(a-b for a,b in zip(row,rg)) for row in r];d=[mi-mg for mi in mu]
    assert sum(w)==1 and min(w)>0 and sum(wi*di for wi,di in zip(w,d))==0
    assert all(sum(wi*si[j] for wi,si in zip(w,residual))==0 for j in range(p))
    ET=sum((frac_arb(wi)*(frac_arb(di+dot(si,si)/2)).exp() for wi,di,si in zip(w,d,residual)),arb(0))
    # Full ordered-pair expansion, independent of producer symmetry packing.
    ET2=arb(0)
    for i in range(n):
        for j in range(n):
            combined=[r[i][k]+r[j][k]-2*rg[k] for k in range(p)]
            exponent=mu[i]+mu[j]-2*mg+dot(combined,combined)/2
            ET2+=frac_arb(w[i]*w[j])*frac_arb(exponent).exp()
    U2=ET2-2*ET+1
    variance=lam*lam+dot(rg,rg);sigma=frac_arb(variance).sqrt();mg_b=frac_arb(mg)
    growth=100*(mg_b+frac_arb(variance)/2).exp()
    base_slope=growth*Phi(frac_arb((mg+variance))/sigma)
    c0=base_slope-100*Phi(mg_b/sigma)
    # Original asset covariance identity, independent of the producer's
    # residual tilt prefactor and common risk-neutral forward simplification.
    c1=-base_slope
    for wi,mi,ri in zip(w,mu,r):
        forward=100*frac_arb(mi+(lam*lam+dot(ri,ri))/2).exp()
        argument=frac_arb(mg+lam*lam+dot(rg,ri))/sigma
        c1+=frac_arb(wi)*forward*Phi(argument)
    M2=100/(frac_arb(lam)*(2*arb.pi()).sqrt())
    return ET,ET2,U2,c0,c1,M2


def main():
    start=time.perf_counter();ctx.prec=PRECISION
    path=ROOT/'results/boundary_convex_screening.json';data=json.loads(path.read_text())
    assert data['source_sha256']==sha(ROOT/'boundary_convex_screening.py')
    for p,digest in data['input_sha256'].items():assert sha(ROOT/p)==digest
    for p,digest in data['dependency_sha256'].items():assert sha(ROOT/p)==digest
    protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text());models={m['model_id']:m for m in protocol['models']}
    assert len(data['records'])==len(models)==9
    records=[];pairs=0;stops=0
    for row in data['records']:
        m=models[row['model_id']];ET,ET2,U2,c0,c1,M2=reference(m)
        for value,key in ((ET,'continuous_T_mean_interval'),(ET2,'continuous_T_square_mean_interval'),(U2,'continuous_U_square_mean_interval'),
                          (c0,'c0_interval'),(c1,'c1_interval'),(c0+c1,'coefficient_interval'),(M2,'global_M2_interval')):inside(value,row[key])
        assert U2.lower()>0 and c1.lower()>0
        upper_u=max(F(0),F(row['continuous_U_square_mean_interval']['upper_exact_fraction']))
        upper_gamma=F(row['global_M2_interval']['upper_exact_fraction'])
        assert F(row['undiscounted_convex_correction_upper'])==upper_u*upper_gamma/2
        coeff=row['coefficient_interval'];lo=F(coeff['lower_exact_fraction']);hi=F(coeff['upper_exact_fraction'])+upper_u*upper_gamma/2
        reported=row['undiscounted_price_interval'];assert F(reported['lower_exact_fraction'])==lo and F(reported['upper_exact_fraction'])==hi
        discount=frac_arb(-F(m['risk_free_rate'])*F(m['maturity'])).exp();inside(discount,row['discount_interval'])
        dl,du=map(F,(row['discount_interval']['lower_exact_fraction'],row['discount_interval']['upper_exact_fraction']))
        products=[d*x for d in (dl,du) for x in (lo,hi)];pl,pu=min(products),max(products)
        pi=row['price_interval'];assert F(pi['lower_exact_fraction'])==pl and F(pi['upper_exact_fraction'])==pu
        radius=(pu-pl)/2;assert F(row['price_radius_upper'])==radius and F(row['price_midpoint'])==(pl+pu)/2
        assert row['stops_without_sampling']==(radius<=F(1,100));stops+=row['stops_without_sampling']
        n=m['assets'];counts=row['operation_counts'];symmetric=n*(n+1)//2;pairs+=n*n
        assert counts['first_moment_terms']==n and counts['symmetric_second_moment_terms']==symmetric
        assert counts['raw_ordered_second_moment_terms_replaced_by_symmetry']==n*n and counts['Gaussian_tilt_terms']==n
        assert counts['square_root_calls']==3 and counts['Gaussian_CDF_calls']==n+2 and counts['exponential_calls']==n+symmetric+3
        assert len(row['second_moment_terms'])==symmetric and sum(t['multiplicity'] for t in row['second_moment_terms'])==n*n
        records.append({'model_id':row['model_id'],'ordered_pairs_checked':n*n,'radius_upper':row['price_radius_upper'],
            'stops_without_sampling':row['stops_without_sampling'],'status':'passed'})
    assert stops==data['tolerance_stop_count']
    result={'schema':'qbp-continuous-convex-call-screen-independent-audit-v1','status':'passed',
        'producer_source_sha256':data['source_sha256'],'screen_archive_sha256':sha(path),'audit_source_sha256':sha(__file__),
        'reference_precision_bits':PRECISION,'models_checked':len(records),'ordered_moment_pairs_checked':pairs,
        'continuous_moment_and_coefficient_intervals_checked':7*len(records),'full_price_intervals_checked':len(records),
        'tolerance_stops_verified':stops,'records':records,'runtime_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Independent full ordered-pair moments and original asset covariance coefficients with erfc CDFs; exact rational bound/discount aggregation. No producer screen function, finite grid or price answer is imported.'}
    (ROOT/'results/boundary_convex_screening_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Passed nine models and',pairs,'ordered moment pairs;',stops,'tolerance stops')

if __name__=='__main__':main()
