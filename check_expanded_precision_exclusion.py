"""Close both deterministic shortcuts in the conditional high-width bound."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, isqrt
import argparse, gzip, hashlib, json
PI = F(355, 113)
POST = F(1, 10**60)
EPS = F(1, 100)
ALPHA = F(1, 3600)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p = Path(p)
    return json.load(gzip.open(p, 'rt')) if p.suffix == '.gz' else json.loads(p.read_text())
def sqrt_upper(x):
    z = x*(1 << 200)
    n = (z.numerator+z.denominator-1)//z.denominator
    k = isqrt(n)
    return F(k+(k*k < n), 1 << 100)
def failure(n, s):
    return sum((F(comb(n, k))*(1-s)**k*s**(n-k)
                for k in range((n+1)//2, n+1)), F(0))
def main(a):
    output = Path(a.output)
    assert not output.exists()
    audit = read(a.resource_review)
    assert audit['status'] == 'passed' and len(audit['conditional_high_precision_exclusion']) == 4
    factor = read(a.factor_resources)
    baselines = {r['original_model_id']: r for r in factor['records']
        if r['kind'] == 'positive_r1' and (r['original_model_id'], r['factor_basis_candidate_id']) in
        {('n2_lambda0.4_scale3', 'identity'), ('n4_lambda0.05_scale3over2', 'plane_1_2_minus'),
         ('n4_lambda0.2_scale3', 'plane_1_2_minus'), ('n8_lambda0.2_scale3', 'plane_2_3_minus')}}
    records = []
    for bound in audit['conditional_high_precision_exclusion']:
        row = baselines[bound['model_id']]
        b, p = row['complete_deterministic_budget'], row['coefficient_informed_prior']
        A = F(b['external_factor_upper'])
        d0 = F(b['coefficient_radius_upper'])+F(b['financial_transfer_upper'])+F(b['external_factor_radius'])
        eps0 = (EPS-d0)/A
        l, u = F(p['finite_ideal_mean_lower']), F(p['finite_ideal_mean_upper'])
        assert 0 <= l <= u < F(1, 2) and eps0 > 0
        radius = A*(u-l)/2
        assert d0+radius > EPS and u > 2*eps0
        M = 2
        while 2*PI*sqrt_upper(u*(1-u))/M+PI**2/M**2+POST > eps0:
            M *= 2
        assert M == bound['zero_eta_phase_length_lower']
        success = 8/PI**2
        reps = 1
        while failure(reps, success) > ALPHA:
            reps += 2
        assert reps == bound['perfect_qft_repetitions_lower'] == 25
        # The two retained prior shortcuts only worsen when eta increases:
        # the enlarged interval width and one-sided upper increase, while
        # the remaining financial tolerance decreases.
        floor = bound['minimum_fixed_nonexp_A_T']+64*(32768+15)
        complete = reps*(2*M-1)*floor
        assert complete == bound['width15_complete_T_lower'] > bound['selected_raw_T']
        records.append({'model_id': bound['model_id'], 'ideal_probability_lower': str(l),
          'ideal_probability_upper': str(u), 'zero_eta_deterministic_error': str(d0),
          'zero_eta_probability_tolerance': str(eps0), 'two_sided_radius': str(radius),
          'two_sided_shortcut_excluded': True, 'one_sided_shortcut_excluded': True,
          'minimum_M': M, 'minimum_repetitions': reps,
          'width15_complete_T_lower': complete, 'selected_raw_T': bound['selected_raw_T']})
    output.write_text(json.dumps({'schema': 'qbp-expanded-high-width-shortcut-audit-v1',
      'status': 'passed', 'source_sha256': sha(__file__),
      'resource_review_sha256': sha(a.resource_review),
      'factor_resources_sha256': sha(a.factor_resources), 'records': records,
      'scope': 'Exact closure of both retained deterministic prior shortcuts and independent zero-eta binary-M/perfect-Fourier repetition floors. The width15..22 complete cap remains conditional on the exact common-width unary-lookup and fixed graph described in the resource review.'}, indent=2)+'\n')
    print('Conditional high-width shortcut closure passed', len(records), 'cases')
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--resource-review', required=True)
    p.add_argument('--factor-resources', required=True)
    p.add_argument('--output', required=True)
    main(p.parse_args())
