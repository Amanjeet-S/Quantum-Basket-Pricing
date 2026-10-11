"""Freeze the financial evaluation family before resource tuning."""
from decimal import Decimal, localcontext
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

D=Decimal
ROOT=Path(__file__).resolve().parent


def models():
    out=[]
    for n,weights in [(2,['.5','.5']),(4,['.4','.3','.2','.1']),
                      (8,['.3','.2','.15','.1','.08','.07','.06','.04'])]:
        w=list(map(D,weights))
        with localcontext() as ctx:
            ctx.prec=70
            if n==2:r=[[D('-.1')],[D('.2')]]
            else:
                c=(D('.0225')/(1-sum(x*x for x in w))).sqrt()
                r=[[D(0) for _ in range(n-1)] for _ in range(n)]
                for j in range(1,n):
                    base=c/(D(j*(j+1))).sqrt()
                    for i in range(j):r[i][j-1]=+base
                    r[j][j-1]=-j*base
        for lam in map(D,['.05','.2','.4']):
            with localcontext() as ctx:
                # Frozen decimal loadings are the exact stipulated model.
                # 200 digits suffice for their finite squares and sums.
                ctx.prec=200
                rowvar=[sum(a*a for a in rr) for rr in r]
                drift=[D('.01')-(lam*lam+v)/2 for v in rowvar]
                rg=[sum(w[i]*r[i][j] for i in range(n)) for j in range(n-1)]
                dispersion=sum(w[i]*rowvar[i] for i in range(n))-sum(a*a for a in rg)
                covariance=[[lam*lam+sum(a*b for a,b in zip(ra,rb)) for rb in r] for ra in r]
                out.append({'model_id':f'n{n}_lambda{lam}','assets':n,'outer_dimension':n-1,
                    'lambda':str(lam),'weights':[str(a) for a in w],
                    'loadings':[[str(a) for a in rr] for rr in r],
                    'log_drifts':[str(a) for a in drift],
                    'covariance':[[str(a) for a in rr] for rr in covariance],
                    'weighted_residual_dispersion':str(dispersion),
                    'geometric_outer_variance':str(sum(a*a for a in rg)),
                    'spot':100,'strike':100,'maturity':1,'risk_free_rate':'.03',
                    'dividend_yield':'.02','parameter_interpretation':'frozen decimals define the exact synthetic risk-neutral model'})
    return out


def protocol():
    mm=models()
    return {'status':'prespecified before resource tuning; computational study, not empirical validation',
        'created_utc':datetime.now(timezone.utc).isoformat(),'models':mm,'caps':[1,10,25],
        'price_tolerances':[.01],'methods':['direct','interval','full_global','full_cap'],
        'bias_budget':{'arithmetic_fraction':.1,'distribution_fraction':.4,'statistical_reserved_fraction':.5,'geometric_correction_allowance':1e-12},
        'precision_search':{'q_min':18,'q_max':30,'initial_function_bits':4,'maximum_function_bits':15,
            'greedy_refinement_steps':45,'exact_count_shortlist':2,
            'full_cap_includes_global_winner_as_feasible_candidate':True},
        'compiler':{'signed_integer_bits':7,'price_register_absolute_limit':64,
            'exp_domain':[-4,4],'log_input_domain':[.0625,8.0625],
            'probability_log_domain':[-2,2],'exact_Toffoli_T_assignment':7},
        'data':'all parameters stipulated; no observed prices, provider samples or purchased access',
        'risk_matching':'weighted residual dispersion is matched to .0225 before frozen-decimal roundoff; outer variances and full covariance still differ',
        'interpretation_limits':['asset count is not a proxy for classical difficulty','no pure dimension-scaling causal inference',
            'no estimator, confidence interval, physical quantum advantage or literature novelty claim'],
        'model_generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    target=ROOT/'results'/'generalisation_protocol.json'
    if target.exists():raise SystemExit('Protocol exists; do not silently replace prespecified inputs')
    target.write_text(json.dumps(protocol(),indent=2)+'\n')
    print(target)
