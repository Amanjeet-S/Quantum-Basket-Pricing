"""Independent public, synthetic classical diagnostics for frozen model family.

The retained price reference uses a positive strike integral, not the signed
closed form used to generate the quantum kernel. Digital shifts are paired
across cases. Their replicate dispersion is a numerical diagnostic, not an
error certificate or confidence guarantee. No market data are used.
"""
from fractions import Fraction as F
from pathlib import Path
from statistics import NormalDist
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
PUBLIC=ROOT/'public_data'
ND=NormalDist()


def normal_cdf(x):
    a=np.asarray(x,dtype=float)
    return np.fromiter((.5*math.erfc(-v/math.sqrt(2)) for v in a.ravel()),float,count=a.size).reshape(a.shape)


def normal_interval(lo,hi):
    """Phi(hi)-Phi(lo), avoiding cancellation in positive upper tails."""
    lo,hi=np.broadcast_arrays(lo,hi)
    assert np.all(hi>=lo-1e-13)
    positive=lo>0
    out=np.empty_like(lo)
    out[positive]=normal_cdf(-lo[positive])-normal_cdf(-hi[positive])
    out[~positive]=normal_cdf(hi[~positive])-normal_cdf(lo[~positive])
    return np.maximum(out,0)


def log_baskets(model,y):
    w=np.asarray(model['weights'],float);r=np.asarray(model['loadings'],float)
    logs=np.asarray(model['log_drifts'],float)+np.asarray(y)@r.T
    high=logs.max(axis=-1)
    ma=high+np.log(np.exp(logs-high[...,None])@w)
    mg=logs@w
    assert np.all(ma>=mg-2e-14)
    return ma,mg


def conditional_positive(ma,mg,lam,cap,order=96):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    strike_logs=np.log1p((cap/100)*(nodes+1)/2)
    ma,mg=np.asarray(ma),np.asarray(mg)
    out=np.empty_like(ma)
    # Bound working memory and make timings reproducible across array sizes.
    af,gf=ma.ravel(),mg.ravel();of=out.ravel()
    for i in range(0,len(af),2048):
        hi=(af[i:i+2048,None]-strike_logs)/lam
        lo=(gf[i:i+2048,None]-strike_logs)/lam
        of[i:i+2048]=cap*(normal_interval(lo,hi)@(weights/2))
    return out


def conditional_closed(m,lam,cap):
    h=cap/100;a=math.log1p(h);m=np.asarray(m)
    x0=m/lam;x1=(m-a)/lam
    out=100*(np.exp(m+lam*lam/2)*normal_interval(x1+lam,x0+lam)
             -normal_interval(x1,x0)+h*normal_cdf(x1))
    return np.clip(out,0,cap)


def geometric_positive(model,cap,order=192):
    w=np.asarray(model['weights'],float);mu=float(w@np.asarray(model['log_drifts'],float))
    sd=math.sqrt(float(model['lambda'])**2+float(model['geometric_outer_variance']))
    nodes,weights=np.polynomial.legendre.leggauss(order)
    return float(cap*normal_cdf((mu-np.log1p(cap/100*(nodes+1)/2))/sd)@(weights/2))


def exact_model_checks(model):
    n,d=model['assets'],model['outer_dimension']
    w=list(map(F,model['weights']));r=[list(map(F,row)) for row in model['loadings']]
    mu=list(map(F,model['log_drifts']));lam=F(model['lambda'])
    assert n==len(w)==len(r)==len(mu) and all(len(row)==d for row in r)
    assert sum(w)==1 and min(w)>0 and lam>0 and d==n-1
    cov=[[lam*lam+sum(x*y for x,y in zip(a,b)) for b in r] for a in r]
    assert cov==[list(map(F,row)) for row in model['covariance']]
    rate=F(model['risk_free_rate'])-F(model['dividend_yield'])
    assert all(mu[i]+cov[i][i]/2==rate for i in range(n))
    rg=[sum(w[i]*r[i][j] for i in range(n)) for j in range(d)]
    geo=sum(x*x for x in rg)
    dispersion=sum(w[i]*sum(x*x for x in r[i]) for i in range(n))-geo
    assert geo==F(model['geometric_outer_variance'])
    assert dispersion==F(model['weighted_residual_dispersion'])
    # Exact rational LDL without pivoting certifies positive definiteness.
    l=[[F(i==j) for j in range(n)] for i in range(n)];pivots=[]
    for j in range(n):
        pivot=cov[j][j]-sum(l[j][k]**2*pivots[k] for k in range(j))
        assert pivot>0;pivots.append(pivot)
        for i in range(j+1,n):
            l[i][j]=(cov[i][j]-sum(l[i][k]*l[j][k]*pivots[k] for k in range(j)))/pivot
    eigen=np.linalg.eigvalsh(np.array([[float(x) for x in row] for row in cov]))
    lips={str(cap):[float((100+cap)*(max(abs(row[j]) for row in r)+abs(rg[j])))
                    for j in range(d)] for cap in (1,10,25)}
    return {'model_id':model['model_id'],'assets':n,'outer_dimension':d,
            'exact_positive_weights_and_unit_sum':True,'exact_covariance_factor_identity':True,
            'exact_risk_neutral_drift_identity':True,'exact_rational_LDL_positive_pivots':True,
            'LDL_pivots_decimal_diagnostic':[float(x) for x in pivots],
            'covariance_eigenvalues_diagnostic':eigen.tolist(),
            'condition_number_diagnostic':float(eigen[-1]/eigen[0]),
            'weighted_residual_dispersion':float(dispersion),
            'dispersion_misfit_from_0.0225':float(dispersion-F('.0225')),
            'geometric_outer_variance':float(geo),'outer_price_Lipschitz_constants':lips}


def sobol_points(bits,dimension):
    """Port of the author-distributed Joe--Kuo 32-bit direction recurrence.

    See the retained BSD licence in PUBLIC/licence. This is a digital shift,
    not Owen scrambling. The first 10 three-dimensional points are checked
    against the author's published source example.
    """
    lines=(PUBLIC/'new-joe-kuo-6.21201').read_text().splitlines()
    directions=np.zeros((dimension,32),dtype=np.uint32)
    directions[0]=np.array([1<<(31-i) for i in range(32)],np.uint32)
    for dim in range(1,dimension):
        values=list(map(int,lines[dim].split()));dd,s,a=values[:3];m=values[3:]
        assert dd==dim+1
        for i in range(s):directions[dim,i]=np.uint32(m[i]<<(31-i))
        for i in range(s,32):
            v=int(directions[dim,i-s])^(int(directions[dim,i-s])>>s)
            for k in range(1,s):
                if (a>>(s-1-k))&1:v^=int(directions[dim,i-k])
            directions[dim,i]=np.uint32(v)
    indices=np.arange(1<<bits,dtype=np.uint32);gray=indices^(indices>>np.uint32(1))
    out=np.zeros((1<<bits,dimension),np.uint32)
    for j in range(bits):out^=((gray>>np.uint32(j))&1)[:,None]*directions[None,:,j]
    return out


def check_sobol_fixture():
    expected=np.array([[0,0,0],[.5,.5,.5],[.75,.25,.25],[.25,.75,.75],
                       [.375,.375,.625],[.875,.875,.125],[.625,.125,.875],
                       [.125,.625,.375],[.1875,.3125,.9375],[.6875,.8125,.4375]])
    assert np.array_equal(sobol_points(4,3)[:10].astype(float)/(1<<32),expected)
    return True


def describe_replicates(values):
    v=np.asarray(values,float)
    return {'replicate_estimates':v.tolist(),'mean':float(v.mean()),
            'replicate_sample_standard_deviation':float(v.std(ddof=1)),
            'replicate_standard_error_diagnostic':float(v.std(ddof=1)/math.sqrt(len(v)))}


def main():
    p=argparse.ArgumentParser();p.add_argument('--protocol',type=Path,default=ROOT/'results/generalisation_protocol.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/generalisation_classical_reference.json')
    p.add_argument('--bits',type=int,default=12);p.add_argument('--replicates',type=int,default=8)
    args=p.parse_args();start=time.perf_counter();protocol=json.loads(args.protocol.read_text())
    metrics=[exact_model_checks(m) for m in protocol['models']];assert check_sobol_fixture()
    directions=sobol_points(args.bits,8);rng=np.random.default_rng(20261005)
    shifted=[];plain=[]
    for rep in range(args.replicates):
        shift=rng.integers(0,1<<32,size=8,dtype=np.uint32)
        u=((directions^shift).astype(float)+.5)/(1<<32)
        shifted.append(np.fromiter((ND.inv_cdf(float(x)) for x in u.ravel()),float,count=u.size).reshape(u.shape))
        plain.append(rng.standard_normal((1<<args.bits,8)))
    rows=[]
    for model in protocol['models']:
        d=model['outer_dimension'];lam=float(model['lambda']);ys=[a[:,:d] for a in shifted]
        logstates=[log_baskets(model,y) for y in ys];mcstates=[log_baskets(model,a[:,:d]) for a in plain]
        for cap in protocol['caps']:
            correction=geometric_positive(model,cap);est48=[];est96=[];mc=[];sd=[];closed_error=0
            coordinate_corr=np.zeros(d)
            for rep,((ma,mg),(am,gm)) in enumerate(zip(logstates,mcstates)):
                r48=conditional_positive(ma,mg,lam,cap,48)
                r96=conditional_positive(ma,mg,lam,cap,96)
                mr=conditional_positive(am,gm,lam,cap,96)
                est48.append(correction+float(r48.mean()));est96.append(correction+float(r96.mean()))
                mc.append(correction+float(mr.mean()));sd.append(float(mr.std(ddof=1)))
                fc=conditional_closed(ma,lam,cap)-conditional_closed(mg,lam,cap)
                closed_error=max(closed_error,float(np.max(abs(fc-r96))))
                for j in range(d):
                    cc=np.corrcoef(r96,ys[rep][:,j])[0,1]
                    coordinate_corr[j]+=cc*cc/args.replicates
            row={'case_id':f"{model['model_id']}_cap{cap}",'model_id':model['model_id'],'cap':cap,
                 'undiscounted_geometric_correction_positive_integral':correction,
                 'shifted_Sobol_positive_96':describe_replicates(est96),
                 'plain_MC_positive_96':describe_replicates(mc),
                 'mean_plain_MC_residual_sample_standard_deviation':float(np.mean(sd)),
                 'positive_strike_quadrature_orders':[48,96],
                 'maximum_paired_replicate_price_refinement':float(np.max(abs(np.asarray(est48)-est96))),
                 'maximum_per_state_signed_formula_vs_positive_96':closed_error,
                 'first_order_linear_coordinate_correlation_squared_diagnostic':coordinate_corr.tolist()}
            if d==1:
                gh=[]
                for order in (48,96):
                    x,w=np.polynomial.hermite.hermgauss(order);ma,mg=log_baskets(model,(math.sqrt(2)*x)[:,None])
                    residual=conditional_positive(ma,mg,lam,cap,192)
                    gh.append(correction+float(residual@(w/math.sqrt(math.pi))))
                row['one_dimensional_Gauss_Hermite']={'orders':[48,96],'prices':gh,
                    'refinement_difference':abs(gh[1]-gh[0]),
                    'qualification':'convergence diagnostic, not a directed quadrature error bound'}
            rows.append(row);print('Reference completed',row['case_id'],row['shifted_Sobol_positive_96']['mean'],flush=True)
    sources=[PUBLIC/'new-joe-kuo-6.21201',PUBLIC/'sobol.cc',PUBLIC/'licence']
    out={'status':'independent positive-integral synthetic diagnostics; no confidence or advantage claim',
         'protocol_sha256':hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'model_checks':metrics,'configurations':rows,
         'design':{'seed':20261005,'points_per_replicate':1<<args.bits,'independent_digital_shifts':args.replicates,
             'shift_bits':32,'paired_cases':True,'Sobol_author_fixture_matches':True,
             'qualifications':['digital shift is not Owen scrambling','plain MC and shifted Sobol use the same sample count',
                'replicate standard errors are diagnostic, not guaranteed confidence bounds',
                'high-dimensional Gaussian state space is not exhaustively integrated',
                'risk matching does not identify pure dimension causality']},
         'source_record':{'page':'https://web.maths.unsw.edu.au/~fkuo/sobol/',
             'direction_set':'new-joe-kuo-6.21201; author page last update 16 September 2010',
             'source_algorithm':'author-distributed sobol.cc; header updated 21 October 2008',
             'licence':'author-distributed BSD licence retained alongside direction data',
             'consulted_version':'author web page, direction file, algorithm source and licence; journal papers were not read',
             'files':[{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources]},
         'elapsed_seconds':time.perf_counter()-start}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n');print(args.output)


if __name__=='__main__':main()
