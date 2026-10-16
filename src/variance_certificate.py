"""Certified finite-output second moments and classical confidence budgets.

The finite law is the actual symmetric rounded quantile constructor. Gaussian
second moments are recorded as a separate analytical comparison. No statement
about the variance of a Bernoulli quantum flag follows from these moments.
"""
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
from statistics import NormalDist
from pathlib import Path
import argparse
import hashlib
import json
import math
import time
from estimation_prior import prior_for_row, normal_density_at_zero

D = Decimal
ROOT = Path(__file__).resolve().parent
REFERENCE_Q = 48
TERMS = 64
PAD = D('1e-100')
NORMAL = NormalDist()


def rational(value):
    return Fraction(value)


def decimal(value):
    value = Fraction(value)
    return D(value.numerator)/D(value.denominator)


def quantile_integers(bits, q):
    """Exactly reproduce encoding_resources.quantile_table's integer map."""
    n = 1 << bits
    positive = [round((NORMAL.inv_cdf((j+.5)/n)+5)*(1 << q))
                -5*(1 << q) for j in range(n//2,n)]
    result = [-y for y in positive[::-1]]+positive
    assert all(a+b == 0 for a,b in zip(result,result[::-1]))
    assert max(abs(y) for y in result) < 5*(1 << q)
    return result


@lru_cache(None)
def reference_moments(bits):
    """Exact even dyadic moments, including the zeroth moment."""
    values = quantile_integers(bits,REFERENCE_Q)
    sums = [0]*(TERMS+1)
    # Symmetry halves the workload; the integer summation is exact.
    for value in values[len(values)//2:]:
        square = value*value
        power = 1
        for k in range(TERMS+1):
            sums[k] += power
            power *= square
    moments = [Fraction(2*value,(1 << bits)*(1 << (2*k*REFERENCE_Q)))
               for k,value in enumerate(sums)]
    digest = hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
    return moments, Fraction(max(values),1 << REFERENCE_Q), digest


@lru_cache(None)
def mgf_interval(bits,q,coefficient):
    """Enclose the exact MGF of each production dyadic marginal.

    For a symmetric marginal, M(t)=E cosh(tY). The exact even moments
    give a positive Taylor polynomial. Its omitted tail is geometrically
    bounded using the largest support point. The q=48 reference and a
    production table round the same finite floating seed, so their pointwise
    displacement is at most (2^-q+2^-48)/2. Multiplication by exp(+-|t|d)
    transfers the reference MGF to the exact production table.
    """
    coefficient = Fraction(coefficient)
    if coefficient == 0:
        return D(1),D(1)
    moments,maximum,_ = reference_moments(bits)
    with localcontext() as ctx:
        ctx.prec = 220
        t = decimal(abs(coefficient))
        term_coefficient,total = D(1),D(1)
        for k in range(1,TERMS+1):
            term_coefficient *= t*t/D((2*k-1)*(2*k))
            total += term_coefficient*decimal(moments[k])
        a = t*decimal(maximum)
        first = a**(2*TERMS+2)/D(math.factorial(2*TERMS+2))
        ratio = a*a/D((2*TERMS+3)*(2*TERMS+4))
        assert ratio < 1 and first/(1-ratio) < D('1e-120')
        tail = first/(1-ratio)
        displacement = D(1)/D(1 << (q+1))+D(1)/D(1 << (REFERENCE_Q+1))
        lo = (total-PAD)*(-t*displacement).exp()-PAD
        hi = (total+tail+PAD)*(t*displacement).exp()+PAD
        assert lo > 0
        return lo,hi


def gaussian_mgf(coefficient):
    coefficient = Fraction(coefficient)
    with localcontext() as ctx:
        ctx.prec = 220
        centre = decimal(coefficient*coefficient/2).exp()
        return centre-PAD,centre+PAD


def product_mgf(loadings,bits,q,continuous=False):
    lo,hi = D(1),D(1)
    for coefficient in loadings:
        a,b = gaussian_mgf(coefficient) if continuous else mgf_interval(bits,q,coefficient)
        lo,hi = lo*a-PAD,hi*b+PAD
    return lo,hi


def exponential_moment(drift,loadings,bits,q,continuous=False):
    a,b = product_mgf(loadings,bits,q,continuous)
    centre = decimal(drift).exp()
    return (centre-PAD)*a-PAD,(centre+PAD)*b+PAD


def basket_geometry(model):
    weights = list(map(rational,model['weights']))
    rows = [list(map(rational,r)) for r in model['loadings']]
    drifts = list(map(rational,model['log_drifts']))
    assert all(w > 0 for w in weights) and sum(weights) == 1
    rg = [sum(w*r[j] for w,r in zip(weights,rows))
          for j in range(model['outer_dimension'])]
    mg = sum(w*m for w,m in zip(weights,drifts))
    return weights,rows,drifts,rg,mg


def log_gap_moment(model,bits,q,continuous=False):
    """Upper bound E[(sum w_i exp(D_i)-1)^2], with signed cancellation."""
    weights,rows,drifts,rg,mg = basket_geometry(model)
    centred = [[a-b for a,b in zip(r,rg)] for r in rows]
    ds = [mu-mg for mu in drifts]
    e1lo,e1hi,e2lo,e2hi = D(0),D(0),D(0),D(0)
    for wi,di,si in zip(weights,ds,centred):
        lo,hi = exponential_moment(di,si,bits,q,continuous)
        e1lo += decimal(wi)*lo-PAD; e1hi += decimal(wi)*hi+PAD
    for i in range(len(weights)):
        for j in range(len(weights)):
            lo,hi = exponential_moment(ds[i]+ds[j],
                [a+b for a,b in zip(centred[i],centred[j])],bits,q,continuous)
            e2lo += decimal(weights[i]*weights[j])*lo-PAD
            e2hi += decimal(weights[i]*weights[j])*hi+PAD
    value = max(D(0),e2hi-2*e1lo+1+PAD)
    return {'first_moment_lower_decimal':str(e1lo),
            'first_moment_upper_decimal':str(e1hi),
            'second_moment_lower_decimal':str(e2lo),
            'second_moment_upper_decimal':str(e2hi),
            'log_gap_square_upper_decimal':str(value)}


def arithmetic_geometric_gap_square(model,bits,q,zbits=None,continuous=False):
    """Finite or Gaussian moment of A-G, inner-conditioned if zbits=None."""
    weights,rows,drifts,rg,mg = basket_geometry(model)
    aa,ag,gg = D(0),D(0),D(0)
    for i in range(len(weights)):
        for j in range(len(weights)):
            _,hi = exponential_moment(drifts[i]+drifts[j],
                [a+b for a,b in zip(rows[i],rows[j])],bits,q,continuous)
            aa += decimal(weights[i]*weights[j])*hi+PAD
        lo,_ = exponential_moment(drifts[i]+mg,
            [a+b for a,b in zip(rows[i],rg)],bits,q,continuous)
        ag += decimal(weights[i])*lo-PAD
    _,gg = exponential_moment(2*mg,[2*x for x in rg],bits,q,continuous)
    lam = rational(model['lambda'])
    if zbits is None:
        # The conditionally averaged positive gap has factor E exp(lambda Z),
        # squared here. This is exp(lambda^2), rather than E exp(2lambda Z).
        inner = decimal(lam*lam).exp()+PAD
    else:
        _,inner = gaussian_mgf(2*lam) if continuous else mgf_interval(zbits,q,2*lam)
    return D(model['spot'])**2*max(D(0),aa-2*ag+gg+PAD)*inner+PAD


def row_certificate(model,row):
    with localcontext() as ctx:
        ctx.prec = 220
        cap = D(row['cap']);method = row['method'];q = row['q'];bits = row['b_y']
        h = cap/D(model['strike']);phi0 = normal_density_at_zero()
        slope = phi0/D(model['lambda'])+PAD
        if method in ('full_global','full_cap'):
            slope = min(slope,(1+h)/h+PAD)
        gap = log_gap_moment(model,bits,q)
        gaussian_gap = log_gap_moment(model,bits,q,True)
        candidates = {'range':D(1)}
        if method != 'direct':
            candidates['log_gap'] = slope*slope*D(gap['log_gap_square_upper_decimal'])+PAD
        if method in ('full_global','full_cap','direct'):
            candidates['positive_AG_gap'] = arithmetic_geometric_gap_square(
                model,bits,q,row['b_z'] if method=='direct' else None)/cap**2+PAD
        exact_grid_upper = min(candidates.values())
        # The graph's conservative recurrence already charges quantile error.
        # With its actual rounded quantile values as ideal inputs, that first
        # error is zero, so the retained bound remains a valid body error.
        error = D.from_float(row['arithmetic_price_error_upper'])/cap
        finite_second = min(D(1),(exact_grid_upper.sqrt()+error)**2+PAD)
        legacy = prior_for_row(model,row)
        variance = min(D(legacy['finite_variance_upper_decimal']),finite_second,D('.25'))
        target = D('.005')/cap
        budgets = {}
        for label,alpha in [('single_price',D('.01')),('simultaneous_27',D(1)/2700-PAD)]:
            logarithm = (2/alpha).ln()+PAD
            bound = (2*variance+2*target/3)*logarithm/(target*target)+PAD
            n_bernstein = int(bound.to_integral_value(rounding=ROUND_CEILING))
            n_hoeffding = int((logarithm/(2*target*target)+PAD).to_integral_value(rounding=ROUND_CEILING))
            n = min(n_bernstein,n_hoeffding)
            method_name = 'Bernstein' if n_bernstein <= n_hoeffding else 'Hoeffding'
            failure = (2*(-D(n)*target*target/(2*variance+2*target/3)).exp()+PAD
                       if method_name=='Bernstein' else 2*(-2*D(n)*target*target).exp()+PAD)
            assert failure <= alpha
            budgets[label] = {'failure_allocation_decimal':str(alpha),
                'Bernstein_samples':n_bernstein,'Hoeffding_samples':n_hoeffding,
                'certified_samples':n,'selected_bound':method_name,
                'certified_failure_upper_decimal':str(failure)}
        return {**{k:row[k] for k in ('model_id','case_id','method','cap','q','b_y','b_z','b_u',
                'exp_address_bits','log_address_bits','probability_address_bits')},
            'candidate_id':row.get('candidate_id'),
            'finite_grid_log_gap_moments':gap,
            'continuous_Gaussian_log_gap_moments':gaussian_gap,
            'normalised_kernel_slope_upper_decimal':str(slope),
            'exact_grid_second_moment_candidates':{k:str(v) for k,v in candidates.items()},
            'exact_grid_second_moment_upper_decimal':str(exact_grid_upper),
            'body_pointwise_normalised_error_upper_decimal':str(error),
            'finite_output_second_moment_upper_decimal':str(finite_second),
            'legacy_mean_variance_upper_decimal':legacy['finite_variance_upper_decimal'],
            'finite_output_variance_upper_decimal':str(variance),
            'single_price_legacy_samples':legacy['classical_certified_samples'],
            'confidence_budgets':budgets,
            'scope':'exact uniformly sampled finite fractional output, not Bernoulli flag; ideal grid body error bounded by retained recurrence; no empirical variance estimate',
            'quantum_scope':'does not tighten the current flag estimator: its Bernoulli variance is p(1-p); variance-aware quantum mean estimation requires a separate construction'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config',type=Path,default=ROOT/'results/architecture_resource_candidates.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/classical_variance_certificates.json')
    a = p.parse_args();data = json.loads(a.config.read_text())
    models = {m['model_id']:m for m in data['models']}
    distinct = {}
    for row in data['configurations']:
        key = (row['case_id'],row['method'],row['q'],row['exp_address_bits'],
               row['log_address_bits'],row['probability_address_bits'])
        if key not in distinct: distinct[key] = row
    started = time.perf_counter();records = []
    for row in distinct.values():
        records.append(row_certificate(models[row['model_id']],row))
        if len(records)%25 == 0: print('Certified finite variance rows',len(records),flush=True)
    out = {'status':'analytical finite-output second moments and fixed-stage confidence budgets',
        'configuration_sha256':hashlib.sha256(a.config.read_bytes()).hexdigest(),
        'source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                         for name in ('variance_certificate.py','estimation_prior.py')},
        'reference_quantile_bits':REFERENCE_Q,'positive_cosh_terms':TERMS,
        'numerical_envelope':str(PAD),
        'reference_quantile_tables':[{'latent_bits':bits,'entries':1 << bits,
            'integer_sha256':reference_moments(bits)[2]}
            for bits in sorted({r['b_y'] for r in records}|{r['b_z'] for r in records if r['b_z']})],
        'configurations':records,'elapsed_seconds':time.perf_counter()-started,
        'novelty_boundary':'established convexity, finite-law moments, Minkowski and Bernstein; project-specific parameter certificate, originality unverified',
        'assumptions':'correctly rounded Decimal exp/ln/sqrt at220digits with1e-100 envelope; exact dyadic moment sums; the shared finite floating quantile seed is the installed constructor; no distribution-bias-to-L2 inference'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'rows':len(records),'sample_range':[min(r['confidence_budgets']['single_price']['Bernstein_samples'] for r in records),max(r['confidence_budgets']['single_price']['Bernstein_samples'] for r in records)],'elapsed_seconds':out['elapsed_seconds']},indent=2))


if __name__ == '__main__':main()
