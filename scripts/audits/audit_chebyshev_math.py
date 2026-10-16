"""Independent exact algebra and confidence audit of the bounded estimator.

The controller is imported only as the implementation under test. Reference
matrix arithmetic, binomial tails, coverage and resource identities below
do not call its mathematical helper functions. Retained diagnostic histories
are optional, and are checked when present. No physical randomness or
financial-price answer is inferred from these checks.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import argparse, hashlib, importlib, json, sys
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {}


def mark(name, n=1):
    COUNTS[name] = COUNTS.get(name, 0)+n


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ball(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def ends(x):
    def endpoint(y):
        m, e = map(int, y.man_exp())
        return F(m*(1 << e)) if e >= 0 else F(m, 1 << (-e))
    return endpoint(x.lower()), endpoint(x.upper())


def mprod(a, b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(2)), F(0))
                       for j in range(2)) for i in range(2))


def mpow(a, k):
    result = ((F(1), F(0)), (F(0), F(1)))
    while k:
        if k % 2:
            result = mprod(result, a)
        a = mprod(a, a)
        k //= 2
    return result


def matvec(a, v):
    return tuple(sum((a[i][j]*v[j] for j in range(2)), F(0))
                 for i in range(2))


def exact_parity_tests():
    # Rational points on the unit circle permit exact coherent calculations,
    # including both degenerate endpoint probabilities and degree 2,047.
    degrees = (1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 127, 128, 2047)
    for a, b in ((F(0), F(1)), (F(1), F(0)), (F(3, 5), F(4, 5)),
                 (F(4, 5), F(3, 5)), (F(5, 13), F(12, 13)),
                 (F(8, 17), F(15, 17))):
        assert a*a+b*b == 1
        q = ((1-2*a*a, 2*a*b), (-2*a*b, 1-2*a*a))
        cheb = [F(1), a]
        for d in range(2, max(degrees)+1):
            cheb.append(2*a*cheb[-1]-cheb[-2])
        for d in degrees:
            k = d//2
            if d % 2:
                state = matvec(mpow(q, k), (a, b))
                value = state[0]**2
            else:
                state = matvec(mpow(q, k-1), (a, b))
                # Final S_good followed by A†, observing the input zero.
                amplitude = -a*state[0]+b*state[1]
                value = amplitude**2
            assert value == cheb[d]**2, (a, d)
            mark('exact_parity_matrix_checks')


def tail(p, n, first, last):
    p = F(p)
    a, den = p.numerator, p.denominator
    numerator = sum(comb(n, j)*a**j*(den-a)**(n-j)
                    for j in range(first, last+1))
    return F(numerator, den**n)


def check_cp_record(record, k, n, failure):
    target = F(failure)/2
    lower, upper = map(F, record['lower_root_bracket'])
    assert F(record['lower']) == lower
    assert upper-lower <= F(1, 1 << record['inverse_bits'])
    if k == 0:
        assert lower == upper == 0
    else:
        assert tail(lower, n, k, n) <= target <= tail(upper, n, k, n)
    lower, upper = map(F, record['upper_root_bracket'])
    assert F(record['upper']) == upper
    assert upper-lower <= F(1, 1 << record['inverse_bits'])
    if k == n:
        assert lower == upper == 1
    else:
        assert tail(upper, n, 0, k) <= target <= tail(lower, n, 0, k)
    assert 0 <= F(record['lower']) <= F(record['upper']) <= 1
    mark('exact_binomial_root_brackets', 2)


def fixed_confidence_tests(controller):
    failure = F(1, 86400)
    intervals = []
    for k in range(101):
        record = controller.fixed_binomial_interval(k, 100, failure)
        check_cp_record(record, k, 100, failure)
        intervals.append((F(record['lower']), F(record['upper'])))
    # This finite check supports the fully general threshold-event proof.
    # It is not substituted for that proof or for physical independence.
    for numerator in range(65):
        q = F(numerator, 64)
        den = 64**100
        missed = sum(comb(100, k)*numerator**k*(64-numerator)**(100-k)
                     for k, (lo, hi) in enumerate(intervals)
                     if not lo <= q <= hi)
        assert F(missed, den) <= failure, q
        mark('exact_finite_confidence_coverage_checks')
    assert 12*failure == F(1, 7200)
    assert 12*failure+F(1, 7200) == F(1, 3600)
    mark('adaptive_and_fallback_failure_identity_checks', 2)


def scalar_branch_tests(controller):
    ctx.prec = 512
    for d in range(1, 25):
        for j in range(d):
            for v in (F(0), F(1, 8), F(1, 2), F(7, 8), F(1)):
                t = (F(j)+v)/d
                if v in (0, 1):
                    q = F(1 if (j+int(v)) % 2 == 0 else 0)
                    qlo = qhi = q
                else:
                    qlo, qhi = ends((arb.pi()*ball(F(j)+v)/2).cos()**2)
                    qlo, qhi = max(F(0), qlo), min(F(1), qhi)
                lo, hi = controller.invert_signal_interval(d, j, qlo, qhi)
                assert F(lo) <= t <= F(hi), (d, j, v)
                pl, pu = controller.arccos_to_probability_interval(t, t)
                true = (arb.pi()*ball(t)/2).cos()**2
                tl, tu = ends(true)
                # Endpoint trigonometric balls can protrude by a negligible
                # amount; their exact probabilities are known independently.
                if t == 0:
                    tl = tu = F(1)
                elif t == 1:
                    tl = tu = F(0)
                elif t == F(1, 2):
                    tl = tu = F(1, 2)
                assert pl <= tl <= tu <= pu, (d, j, v)
                mark('outward_monotone_branch_and_decode_checks')
    for p in (F(0), F(1, 10**12), F(1, 100), F(1, 2), F(99, 100),
              F(1)-F(1, 10**12), F(1)):
        lo, hi = controller.probability_to_arccos_interval(p, p)
        if p in (0, 1):
            tl = tu = 1-p
        elif p == F(1, 2):
            tl = tu = F(1, 2)
        else:
            tl, tu = ends(1-2*ball(p).sqrt().asin()/arb.pi())
        assert lo <= tl <= tu <= hi
        mark('deterministic_prior_inverse_checks')


def retained_history_tests(path):
    ctx.prec = 512
    def reference_inverse(p):
        p = F(p)
        if p in (0, 1, F(1, 2)):
            return p, p
        return ends(2*ball(p).sqrt().asin()/arb.pi())
    def reference_probability(t):
        t = F(t)
        if t in (0, 1):
            return 1-t, 1-t
        if t == F(1, 2):
            return F(1, 2), F(1, 2)
        return ends((arb.pi()*ball(t)/2).cos()**2)
    def branch(d, lo, hi):
        j = (d*lo).numerator//(d*lo).denominator
        if lo == hi == 1:
            j = d-1
        return j if j < d and d*hi <= j+1 else None
    data = json.loads(path.read_text())
    alpha = F(data['alpha_per_complete_price'])
    for case in data['cases']:
        setup = case['controller_setup']
        assert F(setup['alpha_look'])*setup['look_cap'] <= alpha/2
        assert F(setup['alpha_fallback']) == alpha/2
        fb = case['half_alpha_QPE_fallback']['T_count']
        cap = case['pilot_T_cap']
        assert cap == fb//10
        assert case['whole_history_T_cap'] == cap+fb <= F(11, 10)*fb
        mark('finite_worst_history_cap_checks')
        preparation = case['retained_financial_information']['financial_preparation_profile']
        at = preparation['primitive_counts'].get('T', 0)+preparation['primitive_counts'].get('Tdg', 0)
        rt = case['uncontrolled_zero_reflection']['T_count']
        for key, profile in case['sampler_profiles'].items():
            d = int(key)
            k = d//2
            zero = k if d % 2 else k-1
            assert profile['A_forward_calls'] == k+int(d % 2 == 1)
            assert profile['A_adjoint_calls'] == k
            assert profile['good_Z_calls'] == k
            assert profile['zero_reflection_calls'] == zero
            assert profile['per_shot_metrics']['T'] == d*at+zero*rt
            assert profile['per_shot_metrics']['measurements'] == profile['peak_logical_qubits']
            assert profile['per_shot_metrics']['reset_X'] == profile['peak_logical_qubits']
            mark('sampler_parity_cost_checks')
        for diagnostic in case['diagnostics']:
            result = diagnostic['controller_result']
            history = result['history']
            assert len(history) <= 12
            assert result['pilot_shots'] == 100*len(history)
            assert F(result['pilot_failure_upper']) == len(history)*alpha/24
            charged = 0
            for i, row in enumerate(history):
                d, j = row['degree'], row['quadrant']
                lo, hi = map(F, row['prior_angle_interval'])
                assert j <= d*lo <= d*hi <= j+1
                if i:
                    prev = history[i-1]
                    assert d == prev['degree'] or d >= 2*prev['degree']
                    assert tuple(map(F, prev['intersected_angle_interval'])) == (lo, hi)
                check_cp_record(row['signal_confidence'], row['successes'], 100, alpha/24)
                qlo, qhi = F(row['signal_confidence']['lower']), F(row['signal_confidence']['upper'])
                if j % 2:
                    invlo = reference_inverse(qlo)[0]
                    invhi = reference_inverse(qhi)[1]
                else:
                    invlo = reference_inverse(1-qhi)[0]
                    invhi = reference_inverse(1-qlo)[1]
                ia, ib = map(F, row['inverted_angle_interval'])
                assert ia <= (j+invlo)/d <= (j+invhi)/d <= ib
                newlo, newhi = map(F, row['intersected_angle_interval'])
                assert newlo == max(lo, ia) and newhi == min(hi, ib)
                mark('independent_retained_branch_inverse_checks')
                assert row['shots'] == 100
                profile = case['sampler_profiles'][str(d)]
                assert row['complete_batch_T'] == 100*profile['per_shot_metrics']['T']
                charged += row['complete_batch_T']
                assert row['cumulative_pilot_T'] == charged <= cap
                if not row.get('empty_intersection'):
                    plo, phi = map(F, row['decoded_probability_interval'])
                    assert F(row['decoded_probability_radius']) == (phi-plo)/2
                    reflo = reference_probability(newhi)[0]
                    refhi = reference_probability(newlo)[1]
                    assert plo <= max(F(setup['prior_lower']), reflo)
                    assert phi >= min(F(setup['prior_upper']), refhi)
                    if 'largest_capped_monotone_candidate' in row:
                        expected = next(x for x in range(2047, 0, -1)
                                        if branch(x, newlo, newhi) is not None)
                        assert expected == row['largest_capped_monotone_candidate']
                        assert row['candidate_adopted'] == (expected >= 2*d)
                        mark('independent_retained_degree_selection_checks')
                    mark('independent_retained_probability_decode_checks')
                mark('retained_adaptive_history_look_checks')
            assert charged == result['pilot_T_count']
            needs = result['status'] != 'certified'
            assert diagnostic['fallback_invoked_in_resource_history'] == needs
            assert diagnostic['actual_history_T_count'] == charged+int(needs)*fb
            assert diagnostic['actual_history_T_count'] <= case['whole_history_T_cap']
            if not needs:
                pl, pu = map(F, result['probability_interval'])
                assert F(result['estimate']) == (pl+pu)/2
                assert (pu-pl)/2 <= F(setup['epsilon'])
            mark('retained_terminal_probability_and_cost_checks')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result', type=Path, default=ROOT/'results/chebyshev_resources.json')
    parser.add_argument('--output', type=Path, default=ROOT/'results/chebyshev_math_audit.json')
    args = parser.parse_args()
    sys.path.insert(0, str(ROOT))
    controller = importlib.import_module('chebyshev_controller')
    own_source = sha(__file__)
    sources = {name: sha(ROOT/name) for name in
               ('chebyshev_controller.py', 'chebyshev_resources.py', 'chebyshev_estimator_proof.tex')}
    exact_parity_tests()
    fixed_confidence_tests(controller)
    scalar_branch_tests(controller)
    if args.result.exists():
        retained_history_tests(args.result)
    assert sources == {name: sha(ROOT/name) for name in sources}
    assert own_source == sha(__file__)
    out = {'schema': 'qbp-chebyshev-independent-mathematics-audit-v1',
           'status': 'passed', 'source_sha256': sources,
           'auditor_sha256': own_source, 'checks': COUNTS,
           'result_sha256': sha(args.result) if args.result.exists() else None,
           'scope': 'Independent exact two-dimensional coherent parity matrices, exact binomial root/coverage checks, outward branch inversion, adaptive look/cost bookkeeping and complete finite fallback cap. Numerical diagnostics supplement the fully stated proof. No physical execution, source priority or practical quantum advantage is inferred.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
