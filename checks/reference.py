"""Independent finite-target and function-precision diagnostics.

No quantum-resource or novelty claim. No observed market data. Function-table
certificates assume each endpoint is rounded from its exact mathematical value;
libm calculations below are diagnostic generators, not interval enclosures.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import time
from functools import lru_cache
import numpy as np

PHI_MAX = 1 / math.sqrt(2 * math.pi)
PHI_CURVATURE = 1 / math.sqrt(2 * math.pi * math.e)
SQRT2 = math.sqrt(2)
ND = NormalDist()


def ncdf(x):
    x = np.asarray(x, dtype=float)
    return np.fromiter((0.5 * math.erfc(-t / SQRT2) for t in x.ravel()),
                       float, count=x.size).reshape(x.shape)


def interval(lo, hi):
    lo, hi = np.broadcast_arrays(lo, hi)
    return np.where(lo >= 0, ncdf(-lo) - ncdf(-hi), ncdf(hi) - ncdf(lo))


@lru_cache(maxsize=32)
def quantiles(bits):
    n = 1 << bits
    return np.array([ND.inv_cdf((j + 0.5) / n) for j in range(n)])


def payoff(a, k, h):
    return np.maximum(0.0, np.minimum(h, a-k))


@lru_cache(maxsize=32)
def model(d, rho):
    if d == 2 and rho == -1:
        # Parent's EXACT stipulated test instance; sentinel rho=-1 identifies
        # this explicitly specified covariance, not an equicorrelation input.
        return (np.array([-.015,-.03]), np.array([[-.1],[.2]]), .2,
                math.exp(-.03), np.array([[.05,.02],[.02,.08]]))
    sigma = np.linspace(.2, .3, d)
    covariance = np.outer(sigma, sigma) * ((1-rho)*np.eye(d)+rho)
    loading_matrix = np.linalg.cholesky(covariance)
    v = np.linalg.solve(loading_matrix, np.ones(d))
    v /= np.linalg.norm(v)
    c = float((loading_matrix @ v)[0])
    # Same stipulated coordinate system as existing pilot, independent kernels.
    _, _, vh = np.linalg.svd(v[None,:], full_matrices=True)
    r = loading_matrix @ vh[1:].T
    drift = .03-.02-sigma*sigma/2
    return drift, r, c, math.exp(-.03), covariance


def outer_grid(bits, dimensions):
    z = quantiles(bits)
    if dimensions == 1:
        return z[:,None]
    return np.stack(np.meshgrid(*([z]*dimensions), indexing='ij'), axis=-1).reshape(-1, dimensions)


def bases(y, drift, r):
    logs = drift + y @ r.T
    ba = 100*np.exp(logs).mean(axis=1)
    dg = logs.mean(axis=1)
    bg = 100*np.exp(dg)
    return ba, bg, dg, logs


def full_positive(ba, bg, c, k, h, order=96):
    """Positive layer-cake integral; independent of call-value subtraction."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    strikes = k + h*(nodes+1)/2
    lo = (np.log(strikes)[None,:]-np.log(ba)[:,None])/c
    hi = (np.log(strikes)[None,:]-np.log(bg)[:,None])/c
    return h/2*(interval(lo,hi) @ weights)


def full_call_subtraction(ba, bg, c, k, h):
    def call(b, strike):
        x = np.log(b/strike)/c
        return b*math.exp(c*c/2)*ncdf(x+c)-strike*ncdf(x)
    return (call(ba,k)-call(ba,k+h))-(call(bg,k)-call(bg,k+h))


def gauss_coupling(bits):
    """Exact formula evaluated diagnostically with fsum cancellation control."""
    n=1 << bits
    midpoint=quantiles(bits)
    interior=np.array([ND.inv_cdf(j/n) for j in range(1,n)])
    dens=lambda x: math.exp(-x*x/2)*PHI_MAX
    a=math.fsum(dens(x) for x in midpoint)
    b=math.fsum(dens(x) for x in interior)
    # Assumed <= 5e-15 inverse/CDF constant error is explicit, not certified libm.
    diagnostic=2*(a-b)
    return {'value': diagnostic, 'floating_aggregate_roundoff_allowance':
            8*math.ulp(max(a,b,1.0)),
            'qualification':'diagnostic, not a directed-rounding enclosure'}


def finite_case(d, rho, k, h, bits):
    drift,r,c,disc,cov=model(d,rho)
    y=outer_grid(bits,d-1)
    zz=quantiles(bits)
    n=len(zz)
    ba,bg,dg,logs=bases(y,drift,r)
    # Chunk to permit b=10 without excessive temporary memory.
    direct=[]; geom=[]; direct_var=[]; random_cont=[]; random_finite=[]
    strike=k+h*(np.arange(n)+.5)/n
    for start in range(0,len(y),256):
        bb=ba[start:start+256]; gg=bg[start:start+256]
        a=bb[:,None]*np.exp(c*zz)[None,:]
        g=gg[:,None]*np.exp(c*zz)[None,:]
        dr=payoff(a,k,h)-payoff(g,k,h)
        direct.append(dr.mean(axis=1))
        geom.append(payoff(g,k,h).mean(axis=1))
        direct_var.append((dr*dr).mean(axis=1))
        za=(np.log(strike)[None,:]-np.log(bb)[:,None])/c
        zg=(np.log(strike)[None,:]-np.log(gg)[:,None])/c
        rc=h*interval(za,zg)
        # Exact survival counts for the finite midpoint-quantile Z law.
        rf=h*(np.searchsorted(zz,zg,side='right')-
              np.searchsorted(zz,za,side='right'))/n
        random_cont.append((rc.mean(axis=1),(rc*rc).mean(axis=1)))
        random_finite.append((rf.mean(axis=1),(rf*rf).mean(axis=1)))
    direct=np.concatenate(direct); geom=np.concatenate(geom)
    direct_second=np.concatenate(direct_var)
    rc_mean=np.concatenate([a[0] for a in random_cont]); rc_second=np.concatenate([a[1] for a in random_cont])
    rf_mean=np.concatenate([a[0] for a in random_finite]); rf_second=np.concatenate([a[1] for a in random_finite])
    full=full_positive(ba,bg,c,k,h,96)
    refined=full_positive(ba,bg,c,k,h,192)
    calls=full_call_subtraction(ba,bg,c,k,h)
    geo_mu=math.log(100)+float(drift.mean())
    geo_sd=math.sqrt(float(cov.mean()))
    def scalar_call(strike):
        x=(geo_mu-math.log(strike))/geo_sd
        return math.exp(geo_mu+geo_sd*geo_sd/2)*float(ncdf(x+geo_sd))-strike*float(ncdf(x))
    continuous_geo=scalar_call(k)-scalar_call(k+h)
    residual_mean=float(direct.mean())
    row={'dimension':d,'rho':rho,'strike':k,'cap':h,'bits':bits,'grid_n':n,
         'outer_states':len(y),'loading':c,
         'direct_shared_finite_residual_mean':residual_mean,
         'full_shared_finite_residual_mean':residual_mean,
         'random_strike_shared_finite_Z_residual_mean':float(rf_mean.mean()),
         'random_strike_continuous_Z_residual_mean':float(rc_mean.mean()),
         'full_continuous_Z_residual_mean':float(full.mean()),
         'direct_finite_total_price':disc*float((direct+geom).mean()),
         'direct_residual_with_continuous_geometric_correction_price':disc*(continuous_geo+residual_mean),
         'random_residual_with_continuous_geometric_correction_price':disc*(continuous_geo+float(rc_mean.mean())),
         'continuous_geometric_correction':continuous_geo,
         'full_outer_finite_continuous_price':disc*(continuous_geo+float(full.mean())),
         'standard_deviations':{
             'direct_finite':math.sqrt(max(0,float(direct_second.mean())-residual_mean**2)),
             'random_continuous_Z':math.sqrt(max(0,float(rc_second.mean())-float(rc_mean.mean())**2)),
             'full_continuous_Z':float(full.std()),
             'full_finite_Z':float(direct.std())},
         'monetary_target_bias_bounds_undiscounted':{
             'direct_finite_Z_vs_full_continuous_Z':h/n,
             'random_continuous_Z_finite_U_vs_full_continuous_Z':h/n,
             'random_finite_Z_finite_U_vs_direct_finite_Z':h/n},
         'diagnostics':{
             'positive_quadrature_96_vs_192_max':float(np.max(abs(full-refined))),
             'stable_integral_vs_call_subtraction_max':float(np.max(abs(full-calls))),
             'inner_quantile_mean_difference':float(np.max(abs(direct-full))),
             'strike_midpoint_mean_difference':float(np.max(abs(rc_mean-full))),
             'shared_finite_strike_mean_difference':float(np.max(abs(rf_mean-direct)))}}
    lr=(k+h)*(np.linalg.norm(r,axis=1).max()+np.linalg.norm(r.mean(axis=0)))
    row['outer_lipschitz_bound']=float(lr)
    row['outer_quantile_coupling']=gauss_coupling(bits)
    row['outer_target_bias_bound_diagnostic']=float(lr*(d-1)*row['outer_quantile_coupling']['value'])
    assert row['diagnostics']['inner_quantile_mean_difference'] <= h/n+1e-9
    assert row['diagnostics']['strike_midpoint_mean_difference'] <= h/n+1e-9
    assert row['diagnostics']['shared_finite_strike_mean_difference'] <= h/n+1e-9
    return row


def continuous_reference(rho,k,h):
    """Two-asset deterministic normal integration, refinement is diagnostic."""
    drift,r,c,disc,cov=model(2,rho)
    prices=[]
    for order in [48,96,192,384]:
        nodes,weights=np.polynomial.legendre.leggauss(order)
        y=(10*nodes)[:,None]
        ba,bg,_,_=bases(y,drift,r)
        residual=full_positive(ba,bg,c,k,h,192)
        density=np.exp(-y[:,0]**2/2)*PHI_MAX
        residual_mean=10*float(weights @ (residual*density))
        # Independently integrate arithmetic conditional payoff in strike.
        u,uw=np.polynomial.legendre.leggauss(192)
        strikes=k+h*(u+1)/2
        arithmetic=h/2*(ncdf((np.log(ba)[:,None]-np.log(strikes)[None,:])/c) @ uw)
        price=disc*10*float(weights @ (arithmetic*density))
        random_h=h*interval((np.log(strikes)[None,:]-np.log(ba)[:,None])/c,
                            (np.log(strikes)[None,:]-np.log(bg)[:,None])/c)
        random_second=h*0+10*float(weights @ ((random_h*random_h) @ uw/2*density))
        full_second=10*float(weights @ (residual*residual*density))
        direct_second=10*float(weights @ (conditional_direct_second(ba,bg,c,k,h)*density))
        prices.append({'outer_order':order,'price':price,'residual_mean':residual_mean,
                       'residual_standard_deviations':{
                           'direct':math.sqrt(max(0,direct_second-residual_mean**2)),
                           'random_interval':math.sqrt(max(0,random_second-residual_mean**2)),
                           'full_conditional':math.sqrt(max(0,full_second-residual_mean**2))}})
    return {'rho':rho,'strike':k,'cap':h,'method':'positive nested Gauss-Legendre; orders and tail bound recorded',
            'values':prices,'normal_outer_tail_bound':2*h*float(ncdf(-10))*disc,
            'last_order_difference':abs(prices[-1]['price']-prices[-2]['price']),
            'qualification':'refinement diagnostics are not a rigorous quadrature-error enclosure'}


def conditional_direct_second(ba,bg,c,k,h):
    """Independent piecewise positive Z integral of the squared direct residual."""
    nodes,weights=np.polynomial.legendre.leggauss(48)
    result=[]
    for b,g in zip(ba,bg):
        cuts=sorted(set([-10.,10.]+[float(np.clip(math.log(strike/base)/c,-10,10))
                         for strike in [k,k+h] for base in [b,g]]))
        val=0.
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            z=(lo+hi)/2+(hi-lo)*nodes/2
            residual=payoff(b*np.exp(c*z),k,h)-payoff(g*np.exp(c*z),k,h)
            val+=(hi-lo)/2*float(weights @ (residual*residual*np.exp(-z*z/2)*PHI_MAX))
        result.append(val)
    return np.array(result)



def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results'/'continuous_reference.json')
    args=p.parse_args()
    start=time.perf_counter()
    finite=[finite_case(2,-1,100.,10.,b) for b in [4,6,8]]
    finite += [finite_case(2,-1,100.,h,6) for h in [1.,25.]]
    continuous=[continuous_reference(-1,100.,h) for h in [1.,10.,25.]]
    for row in finite+continuous:
        row.pop('rho',None)
        row['instance']='specified covariance [[.05,.02],[.02,.08]]; not an equicorrelation input'
    output={'status':'independent deterministic price and variance diagnostics',
            'model':'stipulated two-asset constructor', 'currency':'GBP',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'numpy_version':np.__version__,'finite':finite,'continuous_reference':continuous,
            'qualification':'continuous quadrature refinements and variance values are diagnostics, not directed error enclosures or confidence guarantees',
            'elapsed_seconds':time.perf_counter()-start}
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(args.output)

if __name__=='__main__':main()
