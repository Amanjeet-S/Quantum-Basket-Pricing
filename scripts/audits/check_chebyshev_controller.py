"""Finite arithmetic and interface checks for the bounded Chebyshev pilot.

I test all fixed-binomial count endpoints and exact coverage on a declared
probability grid. Those finite tests do not replace the full confidence
proof or certify a physical random source. Controlled deterministic counts
are interface tests, not measured quantum prices or statistical evidence.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import chebyshev_controller as c
from flint import arb, ctx


def tail(p, first, last, n=100):
    """Independent exact polynomial, without the producer comparison API."""
    p = F(p)
    return sum((F(comb(n, j))*p**j*(1-p)**(n-j)
                for j in range(first, last+1)), F(0))


def bounds(value):
    return (F(str(value.lower().fmpq())), F(str(value.upper().fmpq())))


def ball(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=ROOT/'results/chebyshev_controller_checks.json')
    args = parser.parse_args()
    ctx.prec = 256
    alpha = F(1, 3600)
    alpha_look = alpha/24
    target = alpha_look/2
    rows = [c.fixed_binomial_interval(k, 100, alpha_look) for k in range(101)]
    checks = 0
    for k, row in enumerate(rows):
        l, u = row['lower_root_bracket'], row['upper_root_bracket']
        assert 0 <= row['lower'] <= row['upper'] <= 1
        assert l[1]-l[0] <= F(1, 1 << 64)
        assert u[1]-u[0] <= F(1, 1 << 64)
        if k:
            assert tail(l[0], k, 100) <= target <= tail(l[1], k, 100)
        else:
            assert l == (F(0), F(0))
        if k < 100:
            assert tail(u[1], 0, k) <= target <= tail(u[0], 0, k)
        else:
            assert u == (F(1), F(1))
        assert row['lower'] == 1-rows[100-k]['upper']
        assert row['upper'] == 1-rows[100-k]['lower']
        if k:
            assert rows[k-1]['lower'] <= row['lower']
            assert rows[k-1]['upper'] <= row['upper']
        checks += 1

    points = sorted(set([F(i, 20) for i in range(21)]
                        + [F(1, 10**6), F(1, 1000), F(1, 100),
                           F(99, 100), F(999, 1000), F(999999, 10**6)]))
    coverage = []
    for p in points:
        excluded = [k for k, row in enumerate(rows)
                    if p < row['lower'] or p > row['upper']]
        failure = sum((F(comb(100, k))*p**k*(1-p)**(100-k)
                       for k in excluded), F(0))
        assert failure <= alpha_look
        coverage.append({'probability': str(p), 'exact_exclusion_probability': str(failure),
                         'excluded_counts': excluded})

    signal_cases = []
    for d in (1, 2, 3, 4, 7, 16, 32, 2047):
        branch_indices = range(d) if d <= 32 else (0, 1, d//2, d-2, d-1)
        for j in branch_indices:
            t = F(2*j+1, 2*d)
            signal = (ball(t)*d*arb.pi()/2).cos()**2
            ql, qu = bounds(signal)
            a, b = c.invert_signal_interval(d, j, ql, qu)
            assert a <= t <= b
            assert c.quadrant(d, F(j, d), F(j+1, d)) == j
            signal_cases.append({'degree': d, 'quadrant': j,
                                 'exact_normalised_arccos': str(t),
                                 'inverse_interval': [str(a), str(b)]})

    degree_cases = []
    for lo, hi in ((F(0), F(1)), (F(0), F(1, 8)),
                   (F(1, 5), F(21, 100)), (F(7, 10), F(71, 100)),
                   (F(9999, 10000), F(1)), (F(1, 4), F(1, 4))):
        d = c.largest_monotone_degree(lo, hi)
        cap = 2047 if lo == hi else min(2047, (1/(hi-lo)).numerator//(1/(hi-lo)).denominator)
        assert c.quadrant(d, lo, hi) is not None
        assert all(c.quadrant(v, lo, hi) is None for v in range(d+1, cap+1))
        degree_cases.append({'angle': [str(lo), str(hi)], 'degree': d, 'search_start': cap})

    costs = {d: 100*(17*d+3) for d in range(1, 2048)}
    calls = []
    zero = c.plan(F(1, 10), F(1, 5), F(1, 10**6), alpha, 0, costs)
    stopped = c.estimate(lambda d, n: calls.append((d, n)), zero)
    assert stopped['status'] == 'requires_fallback' and calls == []
    assert stopped['pilot_T_count'] == 0
    narrow = c.plan(F(1, 10), F(1001, 10000), F(1, 10000), alpha, 0, costs)
    prior_stop = c.estimate(lambda d, n: calls.append((d, n)), narrow)
    assert prior_stop['status'] == 'certified' and calls == []

    interface_runs = []
    for p in (F(0), F(1, 100), F(1, 4), F(1, 2), F(99, 100), F(1)):
        lo, hi = max(F(0), p-F(1, 20)), min(F(1), p+F(1, 20))
        setup = c.plan(lo, hi, F(1, 10000), alpha, 10**12, costs)

        def deterministic_count(d, n):
            # Finite interface exercise only, deliberately not a binomial draw.
            angle = ball(p).sqrt().acos()
            q = (d*angle).cos()**2
            ql, qu = bounds(q)
            nearest = int((ql+qu)*n/2+F(1, 2))
            return min(n, max(0, nearest))

        run = c.estimate(deterministic_count, setup)
        assert run['probability_interval'][0] <= p <= run['probability_interval'][1]
        assert run['pilot_T_count'] <= setup.pilot_T_cap
        assert len(run['history']) <= 12
        previous = None
        for row in run['history']:
            assert row['shots'] == 100
            assert row['conditional_failure_upper'] == alpha_look
            assert 1 <= row['degree'] <= 2047
            if previous is not None:
                assert row['degree'] == previous or row['degree'] >= 2*previous
            previous = row['degree']
        if run['status'] == 'certified':
            assert run['probability_radius'] <= setup.epsilon
        interface_runs.append({'exact_test_probability': str(p),
                               'test_backend': 'Deterministic nearest-mean count, not sampling evidence',
                               'run': c.rational_record(run)})

    invalid_backend_checks = 0
    generous = c.plan(F(1, 10), F(1, 5), F(1, 10**6), alpha, 10**12, costs)
    for invalid in (-1, 101, True, 2.0):
        try:
            c.estimate(lambda d, n: invalid, generous)
        except ValueError:
            invalid_backend_checks += 1
        else:
            raise AssertionError('Invalid count was accepted')

    result = {'status': 'passed', 'source_hashes': c.source_record(),
              'check_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'fixed_binomial_endpoint_cases': checks,
              'fixed_binomial_failure': str(alpha_look),
              'CP_endpoint_checks': [c.rational_record(x) for x in rows],
              'exact_fixed_binomial_coverage_grid': coverage,
              'coverage_scope': 'Finite grid check; general adaptive coverage is established in the separate proof',
              'arbitrary_parity_branch_checks': signal_cases,
              'largest_monotone_degree_checks': degree_cases,
              'zero_budget_stop': c.rational_record(stopped),
              'supplied_prior_stop': c.rational_record(prior_stop),
              'deterministic_interface_runs': interface_runs,
              'invalid_backend_checks': invalid_backend_checks,
              'Hoeffding_reference_radius': str(c.hoeffding_radius(100, alpha_look)),
              'confidence_method_selected': 'fixed_sample_Clopper_Pearson_exact_rational',
              'data_origin': 'Synthetic scalar priors and artificial interface counts; no LSEG inputs',
              'limitations': 'No physical quantum execution, empirical finance sample or random-source certification'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'fixed_binomial_endpoint_cases',
          'invalid_backend_checks', 'confidence_method_selected')}, indent=2))


if __name__ == '__main__':
    main()
