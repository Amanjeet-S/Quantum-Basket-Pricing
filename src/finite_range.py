"""Certified finite residual ranges from convex basket log ratios.

This bounds the actual finite oracle output. It does not assert that the
continuous Gaussian payoff has the same bounded latent support. Its existing
distribution-transfer allowance remains necessary for a continuous price.
"""
from decimal import Decimal,localcontext,ROUND_CEILING
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import math
from variance_certificate import quantile_integers
from estimation_prior import normal_density_at_zero,prior_for_row
from dyadic_plan import sharp_grid_mean_upper

ROOT=Path(__file__).resolve().parent
D=Decimal
PAD=D('1e-100')


def decimal(f):
    f=F(f);return D(f.numerator)/D(f.denominator)


@lru_cache(None)
def support(bits,q):
    values=quantile_integers(bits,q)
    assert min(values)==-max(values)
    return F(max(values),1 << q),hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()


@lru_cache(None)
def geometry(model_json,bits,q):
    m=json.loads(model_json);w=list(map(F,m['weights']));rows=[list(map(F,r)) for r in m['loadings']]
    mu=list(map(F,m['log_drifts']));p=m['outer_dimension']
    assert len(w)<=8 and p<=7 and all(x>0 for x in w) and sum(w)==1
    rg=[sum(wi*r[j] for wi,r in zip(w,rows)) for j in range(p)]
    mg=sum(wi*x for wi,x in zip(w,mu));maximum,sha=support(bits,q)
    best=D('-Infinity');winner=None;max_argument=F(0)
    with localcontext() as ctx:
        ctx.prec=220
        for signs in product((-1,1),repeat=p):
            y=[maximum*s for s in signs]
            exponents=[a-mg+sum((x-z)*yy for x,z,yy in zip(r,rg,y)) for a,r in zip(mu,rows)]
            assert sum(wi*a for wi,a in zip(w,exponents))==0
            max_argument=max(max_argument,*map(abs,exponents))
            assert max(map(abs,exponents),default=F(0))<=20
            total=sum(decimal(wi)*decimal(a).exp() for wi,a in zip(w,exponents))
            value=total.ln()
            if value>best:best,winner=value,signs
        upper=max(D(0),best+PAD)
        mg_max=mg+maximum*sum(map(abs,rg))
    return {'maximum_log_gap_upper_decimal':str(upper),
        'maximum_log_gap_centre_decimal':str(best),'maximising_vertex_signs':list(winner),
        'vertex_count':1 << p,'quantile_support_maximum_fraction':str(maximum),
        'quantile_table_sha256':sha,'maximum_exponential_argument_exact':str(max_argument),
        'geometric_log_base_maximum_exact':str(mg_max),
        'numerical_allowance_decimal':str(PAD)}


def row_certificate(model,row,moment):
    q=row['q'];cap=F(row['cap']);g=geometry(json.dumps(model,sort_keys=True),row['b_y'],q)
    with localcontext() as ctx:
        ctx.prec=220
        h=decimal(cap/F(model['strike']));lam=decimal(F(model['lambda']))
        J=D(g['maximum_log_gap_upper_decimal']);slope=normal_density_at_zero()/lam+PAD
        if row['method'] in ('full_global','full_cap'):slope=min(slope,(1+h)/h+PAD)
        if row['method']=='direct':slope=(1+h)/h+PAD
        candidates={'range':D(1),'log_gap':slope*J+PAD}
        mg_max=decimal(F(g['geometric_log_base_maximum_exact']))
        if row['method']=='direct':
            zmax,zsha=support(row['b_z'],q)
            candidates['positive_AG_gap']=D(model['spot'])*(mg_max+lam*decimal(zmax)).exp()*(J.exp()-1)/D(row['cap'])+PAD
            g={**g,'inner_quantile_support_maximum_fraction':str(zmax),'inner_quantile_table_sha256':zsha}
        elif row['method'] in ('full_global','full_cap'):
            candidates['positive_AG_gap']=D(model['spot'])*(mg_max+lam*lam/2).exp()*(J.exp()-1)/D(row['cap'])+PAD
        ideal=min(candidates.values())
        eta=D.from_float(row['arithmetic_price_error_upper'])/D(row['cap'])
        upper=min(D(1),ideal+eta+PAD)
        bound=F(upper)*(1 << q)
        maximum=bound.numerator//bound.denominator
        assert 0<=maximum<=1 << q
        rho=F(maximum,1 << q)
        legacy=prior_for_row(model,row)
        mean=min(F(legacy['finite_mean_probability_upper_decimal']),rho)
        second=min(F(moment['finite_output_second_moment_upper_decimal']),rho*mean,rho*rho)
        if maximum:mean=min(mean,sharp_grid_mean_upper(q,second,upper_integer=maximum))
        else:mean=F(0)
        second=min(second,rho*mean)
    return {**{k:row[k] for k in ('case_id','model_id','method','cap','q','exp_address_bits','log_address_bits','probability_address_bits')},
        **g,'normalised_slope_upper_decimal':str(slope),
        'ideal_grid_range_candidates_decimal':{k:str(v) for k,v in candidates.items()},
        'body_pointwise_error_upper_decimal':str(eta),'finite_output_range_upper_decimal':str(upper),
        'max_integer':maximum,'finite_lattice_range_fraction':str(rho),
        'fair_mean_upper_fraction':str(mean),'fair_second_moment_upper_fraction':str(second),
        'legacy_mean_upper_decimal':legacy['finite_mean_probability_upper_decimal'],
        'scope':'actual finite output with rounded latent values; no continuous-support or extra distribution-error inference'}


def main():
    path=ROOT/'results/architecture_resource_candidates.json';vp=ROOT/'results/classical_variance_certificates.json'
    data=json.loads(path.read_text());mom=json.loads(vp.read_text())
    models={m['model_id']:m for m in data['models']}
    key=lambda r:tuple(r[k] for k in ('case_id','method','q','exp_address_bits','log_address_bits','probability_address_bits'))
    source={key(r):r for r in data['configurations']};moment={key(r):r for r in mom['configurations']}
    assert len(source)==len(moment)==273 and set(source)==set(moment)
    rows=[row_certificate(models[r['model_id']],r,moment[k]) for k,r in source.items()]
    out={'status':'convex-vertex and arithmetic certificates for actual finite ranges',
        'configuration_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'moment_sha256':hashlib.sha256(vp.read_bytes()).hexdigest(),
        'source_hashes':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('finite_range.py','variance_certificate.py','estimation_prior.py','dyadic_plan.py')},
        'configurations':rows,'numerical_premise':'correctly rounded Decimal exp/ln at220digits; at most128vertices,8assets,7coordinates and exponent absolute bound20; aggregate error enclosed by1e-100',
        'novelty_boundary':'established convexity, arithmetic-geometric order and range normalisation; project-specific uniform finite-oracle certificate, priority unverified'}
    (ROOT/'results/finite_range_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'rows':len(rows),'strictly_below_one':sum(r['max_integer']<1 << r['q'] for r in rows),'minimum_range':min(float(F(r['finite_lattice_range_fraction'])) for r in rows)},indent=2))


if __name__=='__main__':main()
