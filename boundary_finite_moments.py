"""Additional nonnegative finite-grid moments for the n2 boundary study.

Only the ideal dispersion U on actual dyadic input coordinates is evaluated.
No conditional call value or basket price answer enters this preprocessing.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,time
from flint import arb,ctx
from generalisation_models import models
from boundary_value_certificate import quantile_certificate,endpoints,PRECISION
from integral_value_oracle import model_parameters
from continuous_quadrature import frac_arb,serialise_interval

ROOT=Path(__file__).resolve().parent

def finite_moments(model,L=6,b_y=10,q=28,orders=(2,3,4,6)):
    if model['outer_dimension']!=1:raise ValueError('Direct nonnegative moment quadrature is presently retained for n2 only')
    ctx.prec=PRECISION;start=time.perf_counter();w,d,s,mg,rg,lam=model_parameters(model)
    nodes=quantile_certificate(L,b_y,q)['nodes'];sums={r:arb(0) for r in orders}
    for row in nodes:
        y=F(row['actual_coordinate']);U=sum((frac_arb(wi)*(frac_arb(di+si[0]*y)).exp() for wi,di,si in zip(w,d,s)),arb(0))-1
        # Weighted geometric centring proves U>=0. Interval error at U=0
        # is harmless for even moments; odd-order lower endpoints are clipped.
        for r in orders:sums[r]+=U**r
    intervals={}
    for r,total in sums.items():
        lo,hi=endpoints(total/len(nodes));intervals[str(r)]=serialise_interval(max(F(0),lo),max(F(0),hi))
    return {'schema':'qbp-boundary-n2-finite-moments-v1','status':'passed','model_id':model['model_id'],
        'L':str(F(L)),'b_y':b_y,'q':q,'distribution':'quantile','outer_dimension':1,
        'centred_U_moments':intervals,'finite_nodes':len(nodes),
        'classical_conditional_exponential_evaluations':len(nodes)*len(w),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_nodes_sha256':hashlib.sha256(json.dumps(nodes,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        'generation_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Paid deterministic expectation of powers of ideal U on the actual outer dyadic law, no financial payoff or price-answer input.'}


def main():
    records=[finite_moments(m) for m in models() if m['assets']==2]
    retained=json.loads((ROOT/'results/integral_finite_grid_moments.json').read_text())['records']
    for row in records:
        old=next(x for x in retained if x['model_id']==row['model_id'])
        for r in ('2','3'):
            a,b=row['centred_U_moments'][r],old['centred_U_moments'][r]
            assert max(F(a['lower_exact_fraction']),F(b['lower_exact_fraction']))<=min(F(a['upper_exact_fraction']),F(b['upper_exact_fraction']))
    result={'schema':'qbp-boundary-n2-moment-family-v1','status':'passed','records':records,
        'independent_factorised_M2_M3_overlap':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['generalisation_models.py','integral_value_oracle.py','integral_value_certificate.py','boundary_value_certificate.py','results/integral_finite_grid_moments.json']}}
    (ROOT/'results/boundary_finite_moments.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Retained',len(records),'direct positive moment records')

if __name__=='__main__':main()
