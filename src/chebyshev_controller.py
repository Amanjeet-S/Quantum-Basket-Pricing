"""Bounded, confidence-safe Chebyshev sampling at an existing pricing flag.

This is a declared finite adaptation of Rall and Fuller's published ChebAE,
not an implementation of its unchanged late-round heuristic or empirical
query coefficient. Each confidence look uses exactly 100 fresh shots. The
backend measures T_d(sqrt(p))**2, including the distinct even-degree return
projector. No financial-price answer is an input to planning or scheduling.

The controller returns a certified probability interval or explicitly asks
its frontend for the independently charged canonical-QPE fallback. It does
not execute or hide that fallback. The frontend assigns alpha/2 to each
stage and charges every sampler, observation and reset. Its supplied pilot
budget is at most one tenth of its complete fallback T count.

All endpoints are rational and outward. Fixed-binomial confidence is the
classical Clopper-Pearson inversion, with its roots enclosed by exact integer
comparisons rather than floating-point quantiles. The existing certified
scalar machinery handles sin-squared and inverse endpoints, including
stationary points. No physical quantum hardware or randomness is certified
by these arithmetic checks.
"""
from dataclasses import dataclass, asdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path
from collections.abc import Mapping

from certified_scalar import (floor_fraction, inverse_probability_interval,
                              probability_interval, sqrt_upper_fraction)

LOOK_CAP = 12
DEGREE_CAP = 2047
BATCH_SHOTS = 100
CP_BITS = 64


def digest(value):
    import json
    return sha256(json.dumps(rational_record(value), sort_keys=True,
                             separators=(',', ':')).encode()).hexdigest()


@lru_cache(None)
def _coefficients(n):
    return tuple(comb(n, j) for j in range(n+1))


def binomial_tail_compare(p, n, first, last, target):
    """Sign of the exact binomial range probability minus target.

    p=a/D gives numerator sum C(n,j)a**j(D-a)**(n-j) and
    denominator D**n. Cross multiplication proves every decision exactly.
    """
    p, target = F(p), F(target)
    if not 0 <= p <= 1 or not 0 <= target <= 1:
        raise ValueError('Invalid probability or tail target')
    if isinstance(n, bool) or not isinstance(n, int) or not 1 <= n <= 100:
        raise ValueError('Declared fixed-binomial range is1<=n<=100')
    if not 0 <= first <= last <= n:
        raise ValueError('Invalid binomial range')
    a, denominator = p.numerator, p.denominator
    b = denominator-a
    ap = [1]
    bp = [1]
    for _ in range(n):
        ap.append(ap[-1]*a)
        bp.append(bp[-1]*b)
    numerator = sum(_coefficients(n)[j]*ap[j]*bp[n-j]
                    for j in range(first, last+1))
    difference = (numerator*target.denominator
                  - denominator**n*target.numerator)
    return (difference > 0)-(difference < 0)


@lru_cache(maxsize=4096)
def fixed_binomial_interval(successes, n, failure, bits=CP_BITS):
    """Outward fixed-sample Clopper-Pearson interval and root certificates.

    A lower interior root solves Pr_q(X>=successes)=failure/2;
    an upper interior root solves Pr_q(X<=successes)=failure/2.
    The appropriate outer bracket endpoints widen the exact interval.
    The inverse brackets have width at most2**(-bits).
    """
    failure = F(failure)
    if any(isinstance(x, bool) or not isinstance(x, int)
           for x in (successes, n, bits)):
        raise ValueError('Counts and inverse precision must be integers')
    if not 1 <= n <= 100 or not 0 <= successes <= n or not 1 <= bits <= 80:
        raise ValueError('Fixed-binomial argument outside declared range')
    if not 0 < failure < 1:
        raise ValueError('Failure probability must lie strictly between0and1')
    target = failure/2

    def enclose(first, last, increasing):
        lo, hi = F(0), F(1)
        for _ in range(bits):
            midpoint = (lo+hi)/2
            comparison = binomial_tail_compare(midpoint, n, first, last, target)
            if comparison == 0:
                return midpoint, midpoint
            if (comparison < 0) == increasing:
                lo = midpoint
            else:
                hi = midpoint
        return lo, hi

    lower = ((F(0), F(0)) if successes == 0
             else enclose(successes, n, True))
    upper = ((F(1), F(1)) if successes == n
             else enclose(0, successes, False))
    return {'lower': lower[0], 'upper': upper[1],
            'lower_root_bracket': lower, 'upper_root_bracket': upper,
            'tail_target': target, 'inverse_bits': bits,
            'method': 'fixed_sample_Clopper_Pearson_exact_rational',
            'arithmetic': 'Exact binomial integer comparisons at dyadic bisection points'}


def log_upper(value, terms=80):
    """Exact upper log from positive atanh series and binary reduction."""
    value = F(value)
    if value < 1 or terms < 1:
        raise ValueError('Log argument must be at least1')
    exponent = 0
    while value >= 2:
        value /= 2
        exponent += 1

    def series(x):
        z = (x-1)/(x+1)
        return (sum((2*z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
                + 2*z**(2*terms+1)/((2*terms+1)*(1-z*z)))
    return exponent*series(F(2))+series(value)


def hoeffding_radius(n, failure):
    """Independent conservative reference, not the selected CP interval."""
    return sqrt_upper_fraction(log_upper(2/F(failure))/(2*n))


def probability_to_arccos_interval(lower, upper):
    lower, upper = F(lower), F(upper)
    if not 0 <= lower <= upper <= 1:
        raise ValueError('Invalid supplied probability prior')
    inverse_low = inverse_probability_interval(lower)
    inverse_high = inverse_probability_interval(upper)
    return 1-inverse_high[1], 1-inverse_low[0]


def arccos_to_probability_interval(lower, upper):
    lower, upper = F(lower), F(upper)
    if not 0 <= lower <= upper <= 1:
        raise ValueError('Invalid normalised arccos interval')
    return (max(F(0), 1-probability_interval(upper)[1]),
            min(F(1), 1-probability_interval(lower)[0]))


def quadrant(degree, lower, upper):
    """A certified monotone cos-squared branch; endpoints may be extrema."""
    j = floor_fraction(degree*lower)
    if lower == upper and j == degree:
        return degree-1
    return j if degree*upper <= j+1 and j < degree else None


def largest_monotone_degree(lower, upper, degree_cap=DEGREE_CAP):
    lower, upper = F(lower), F(upper)
    if not 0 <= lower <= upper <= 1 or degree_cap < 1:
        raise ValueError('Invalid angle domain or degree cap')
    width = upper-lower
    start = degree_cap if width == 0 else min(degree_cap, floor_fraction(1/width))
    for degree in range(start, 0, -1):
        if quadrant(degree, lower, upper) is not None:
            return degree
    raise ArithmeticError('The degree-one monotone domain was lost')


def invert_signal_interval(degree, j, lower, upper):
    """Outward inverse of cos²(pi*d*t/2) on one certified branch."""
    lower, upper = F(lower), F(upper)
    if not 0 <= lower <= upper <= 1 or not 0 <= j < degree:
        raise ValueError('Invalid signal interval or branch')
    if j % 2 == 0:
        left = inverse_probability_interval(1-upper)[0]
        right = inverse_probability_interval(1-lower)[1]
    else:
        left = inverse_probability_interval(lower)[0]
        right = inverse_probability_interval(upper)[1]
    return (F(j)+left)/degree, (F(j)+right)/degree


@dataclass(frozen=True)
class Plan:
    prior_lower: F
    prior_upper: F
    epsilon: F
    alpha_total: F
    alpha_pilot: F
    alpha_fallback: F
    alpha_look: F
    pilot_T_cap: int
    normalised_arccos_prior: tuple
    initial_degree: int
    batch_T_costs: tuple
    look_cap: int = LOOK_CAP
    degree_cap: int = DEGREE_CAP
    batch_shots: int = BATCH_SHOTS
    inverse_bits: int = CP_BITS
    confidence_method: str = 'fixed_sample_Clopper_Pearson_exact_rational'

    def batch_cost(self, degree):
        if not 1 <= degree <= self.degree_cap:
            raise ValueError('Sampler degree outside the declared finite range')
        return self.batch_T_costs[degree-1]


def plan(prior_lower, prior_upper, epsilon, alpha_total, pilot_T_cap,
         batch_cost_by_degree):
    """Register a complete supplied-cost pilot without observing outcomes.

    The cost mapping/callable must contain the complete100-shot T cost,
    including the parity-correct emitted sampler, observations and reset.
    It is not a query count, and the controller cannot verify those circuits.
    The frontend separately checks pilot_T_cap<=fallback_T/10.
    """
    lo, hi, epsilon, alpha = map(F, (prior_lower, prior_upper, epsilon, alpha_total))
    if not 0 <= lo <= hi <= 1 or not 0 < epsilon < 1 or not 0 < alpha < 1:
        raise ValueError('Invalid prior, tolerance or confidence')
    if isinstance(pilot_T_cap, bool) or not isinstance(pilot_T_cap, int) or pilot_T_cap < 0:
        raise ValueError('Pilot T budget must be a nonnegative integer')
    costs = []
    for d in range(1, DEGREE_CAP+1):
        cost = (batch_cost_by_degree[d] if isinstance(batch_cost_by_degree, Mapping)
                else batch_cost_by_degree(d))
        if isinstance(cost, bool) or not isinstance(cost, int) or cost < 0:
            raise ValueError('Complete batch costs must be nonnegative integers')
        costs.append(cost)
    angle = probability_to_arccos_interval(lo, hi)
    return Plan(lo, hi, epsilon, alpha, alpha/2, alpha/2,
                alpha/(2*LOOK_CAP), pilot_T_cap, angle,
                largest_monotone_degree(*angle), tuple(costs))


def estimate(count_ones, setup):
    """Run a finite pilot against fresh independent fixed-binomial shots.

    count_ones(degree,100) must supply successes for this new batch alone.
    In particular, counts are never pooled across adaptively stopped looks.
    """
    if not isinstance(setup, Plan):
        raise TypeError('A registered Plan is required')
    lo, hi = setup.normalised_arccos_prior
    degree = setup.initial_degree
    history = []
    charged = 0
    probability = (setup.prior_lower, setup.prior_upper)

    def result(status, reason):
        a, b = probability
        return {'status': status, 'reason': reason,
                'probability_interval': probability,
                'estimate': (a+b)/2 if status == 'certified' else None,
                'probability_radius': (b-a)/2,
                'normalised_arccos_interval': (lo, hi),
                'history': history, 'pilot_T_count': charged,
                'pilot_shots': len(history)*setup.batch_shots,
                'final_degree': degree,
                'pilot_failure_upper': len(history)*setup.alpha_look,
                'reserved_pilot_failure': setup.alpha_pilot,
                'reserved_fallback_failure': setup.alpha_fallback,
                'fallback_executed_here': False,
                'scope': 'Ideal independent fixed-binomial shot contract; no hardware or price answer'}

    if (probability[1]-probability[0])/2 <= setup.epsilon:
        return result('certified', 'supplied_prior_stops_without_shots')
    for index in range(setup.look_cap):
        j = quadrant(degree, lo, hi)
        if j is None:
            raise AssertionError('An updated angle interval lost its chosen monotone branch')
        cost = setup.batch_cost(degree)
        if charged+cost > setup.pilot_T_cap:
            return result('requires_fallback', 'next_complete_batch_exceeds_pilot_T_cap')
        count = count_ones(degree, setup.batch_shots)
        if isinstance(count, bool) or not isinstance(count, int) or not 0 <= count <= setup.batch_shots:
            raise ValueError('Backend must return integer successes for the fresh fixed batch')
        charged += cost
        cp = fixed_binomial_interval(count, setup.batch_shots, setup.alpha_look,
                                     setup.inverse_bits)
        a, b = invert_signal_interval(degree, j, cp['lower'], cp['upper'])
        newlo, newhi = max(lo, a), min(hi, b)
        record = {'look': index+1, 'degree': degree, 'quadrant': j,
                  'parity': 'odd_good_flag' if degree % 2 else 'even_zero_input_return',
                  'shots': setup.batch_shots, 'successes': count,
                  'conditional_failure_upper': setup.alpha_look,
                  'signal_confidence': cp, 'inverted_angle_interval': (a, b),
                  'prior_angle_interval': (lo, hi),
                  'intersected_angle_interval': (newlo, newhi),
                  'complete_batch_T': cost, 'cumulative_pilot_T': charged}
        history.append(record)
        if newlo > newhi:
            record['empty_intersection'] = True
            return result('confidence_failure', 'empty_angle_intersection')
        lo, hi = newlo, newhi
        pl, pu = arccos_to_probability_interval(lo, hi)
        probability = max(setup.prior_lower, pl), min(setup.prior_upper, pu)
        record['decoded_probability_interval'] = probability
        if probability[0] > probability[1]:
            record['empty_intersection'] = True
            return result('confidence_failure', 'empty_probability_prior_intersection')
        record['decoded_probability_radius'] = (probability[1]-probability[0])/2
        if record['decoded_probability_radius'] <= setup.epsilon:
            return result('certified', 'decoded_probability_radius_meets_tolerance')
        candidate = largest_monotone_degree(lo, hi, setup.degree_cap)
        adopted = candidate >= 2*degree
        record['largest_capped_monotone_candidate'] = candidate
        record['candidate_adopted'] = adopted
        if adopted:
            degree = candidate
        if index+1 == setup.look_cap:
            return result('requires_fallback', 'fresh_look_cap_reached')
    raise AssertionError('Finite pilot did not reach an explicit exit')


def rational_record(value):
    """Portable exact JSON records; float displays never drive decisions."""
    if isinstance(value, F):
        return str(value)
    if isinstance(value, Plan):
        out = asdict(value)
        costs = out.pop('batch_T_costs')
        out['batch_T_costs_sha256'] = digest(costs)
        out['batch_T_cost_entries'] = len(costs)
        out['confidence_scope'] = '12 fresh100-shot CP looks, no pooling, then alpha/2 QPE fallback'
        return rational_record(out)
    if isinstance(value, dict):
        return {str(k): rational_record(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [rational_record(v) for v in value]
    return value


def source_record():
    scalar = Path(__file__).with_name('certified_scalar.py')
    return {'controller_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
            'certified_scalar_sha256': sha256(scalar.read_bytes()).hexdigest(),
            'method': 'Declared bounded Chebyshev pilot with exact fixed-binomial inverses',
            'literature_status': 'Adaptation of published mechanisms; no estimator novelty claim'}
