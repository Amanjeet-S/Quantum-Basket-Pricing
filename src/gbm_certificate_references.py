"""Independent arbitrary-precision diagnostics, not outward enclosures."""
from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import math
from typing import Mapping

import mpmath as mp

def mp_number(value):
    if isinstance(value, Fraction) or (isinstance(value, str) and '/' in value):
        fraction = Fraction(value)
        return mp.mpf(fraction.numerator)/fraction.denominator
    return mp.mpf(value)


CONTRACT_FIELDS = ('spot', 'strike', 'cap', 'maturity', 'rate', 'dividend')


def numbers(contract: Mapping[str, str]):
    return tuple(mp_number(contract[name]) for name in CONTRACT_FIELDS)


def normal_cdf(x):
    return mp.erfc(-x / mp.sqrt(2)) / 2


def capped_price(contract: Mapping[str, str], sigma):
    s, k, h, t, r, q = numbers(contract)
    v = mp_number(sigma) * mp.sqrt(t)
    # Tail-valued CDFs do not form 1-Phi when the positive-tail price is small.
    def call(strike):
        d1 = (mp.log(s / strike) + (r-q) * t + v*v/2) / v
        return s * mp.exp(-q*t) * normal_cdf(d1) - strike * mp.exp(-r*t) * normal_cdf(d1-v)
    return call(k) - call(k+h)


def peak_square(contract: Mapping[str, str]):
    s, k, h, t, r, q = numbers(contract)
    return (2 * mp.log(k/s) + mp.log1p(h/k)) / t - 2 * (r-q)


def full_image(contract, lower, upper):
    lo, hi = mp_number(lower), mp_number(upper)
    cstar = peak_square(contract)
    endpoints = [capped_price(contract, lo), capped_price(contract, hi)]
    candidates = endpoints[:]
    peak = mp.sqrt(cstar) if cstar > 0 else None
    interior = peak is not None and lo < peak < hi
    if interior:
        candidates.append(capped_price(contract, peak))
    return {'lower': min(endpoints), 'upper': max(candidates),
            'endpoint_only_upper': max(endpoints), 'peak': peak,
            'interior': interior}


def digital_integral(contract, sigma):
    s, k, h, t, r, q = numbers(contract)
    sig = mp_number(sigma)
    width = mp.log1p(h/k)
    v = sig*mp.sqrt(t)
    a = (mp.log(s/k) + (r-q-sig*sig/2)*t)/v
    # Resolves sharp normal transitions in the positive normalised-strike integral.
    breaks = [mp.mpf(0), mp.mpf(1)]
    for score in range(-16, 17):
        x = (a-score)*v/width
        if 0 < x < 1:
            breaks.append(x)
    breaks = sorted(set(breaks))
    # A relative scale prevents absolute-error stopping from accepting an
    # inaccurate integral when the entire positive price is below 10**(-dps).
    # This informed refinement retains exactly the same mathematical integral.
    tail_scale = normal_cdf(a)
    f = lambda x: mp.exp(width*x)*(normal_cdf(a-width*x/v)/tail_scale)
    return mp.exp(-r*t)*k*width*tail_scale*mp.quad(f, breaks)


def chi_square_cdf(value, degrees):
    x = mp_number(value)
    return mp.gammainc(mp.mpf(degrees)/2, 0, x/2, regularized=True)


def chi_square_quantile_bisection(probability, degrees, digits=100):
    """Point approximation using a different library/solver, no interval proof."""
    with mp.workdps(digits):
        p = mp_number(probability)
        lo, hi = mp.mpf(0), mp.mpf(max(1, degrees))
        while chi_square_cdf(hi, degrees) < p:
            hi *= 2
        for _ in range(4*digits+32):
            mid = (lo+hi)/2
            if chi_square_cdf(mid, degrees) < p:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2


def exact_decimal_scatter(observations):
    """sum (x-xbar)^2 as an exact rational for finite decimal observations."""
    xs = [Fraction(Decimal(value)) for value in observations]
    total = sum(xs, Fraction())
    return sum((x*x for x in xs), Fraction()) - total*total/len(xs)


def rounded_decimal_observations(draws, places=6):
    """Round binary64 draws exactly to a decimal grid; retain no raw row files."""
    scale = 10**places
    result = []
    squared_error = Fraction()
    for x in draws:
        value = Fraction.from_float(float(x))
        # Fraction.__round__ implements ties-to-even without binary scaling loss.
        integer = round(value*scale)
        rounded = Fraction(integer, scale)
        text = format(Decimal(integer)/Decimal(scale), f'.{places}f')
        result.append(text)
        squared_error += (rounded-value)**2
    return result, squared_error


def fixture_cases():
    base = {'spot':'100', 'strike':'100', 'cap':'20', 'maturity':'1',
            'rate':'0.03', 'dividend':'0.01'}
    cases = [
      ('straddling_peak', base, '0.2', '0.6'),
      ('increasing_interval', base, '0.08', '0.2'),
      ('decreasing_interval', base, '0.6', '1.1'),
      ('singleton', base, '0.2', '0.2'),
      ('deep_in_money', {**base, 'spot':'160'}, '0.1', '1.2'),
      ('negative_rate', {**base, 'rate':'-0.03'}, '0.15', '0.75'),
      ('tiny_cap', {**base, 'cap':'1e-15', 'rate':'0', 'dividend':'0'}, '1e-9', '6e-9'),
      ('large_scale', {**base, 'spot':'1e200','strike':'1e200','cap':'2e199'}, '0.2', '0.6'),
      ('small_scale', {**base, 'spot':'1e-200','strike':'1e-200','cap':'2e-201'}, '0.2','0.6'),
      ('small_tail', {**base, 'spot':'1','rate':'0','dividend':'0'}, '0.1', '0.2'),
    ]
    with mp.workdps(180):
        peak = mp.sqrt(peak_square(base))
        quantum = mp.mpf('1e-60')
        left = mp.floor(peak/quantum)*quantum
        right = left+quantum
        cases.append(('peak_one_decimal_place', base, mp.nstr(left, 180), mp.nstr(right,180)))
        # At this spot c*=0; finite decimal perturbations give either sign.
        boundary = mp.sqrt(mp.mpf(100)*120)*mp.exp(-mp.mpf('0.02'))
        for name, sign in [('almost_zero_peak_positive',-1),('almost_zero_peak_negative',1)]:
            contract = {**base, 'spot':mp.nstr(boundary*(1+sign*mp.mpf('1e-60')), 160)}
            cases.append((name, contract, '1e-31', '1e-29'))
    return [{'name':name,'contract':dict(c),'lower':lo,'upper':hi} for name,c,lo,hi in cases]


def exact_composition_example():
    """Finite dependent law satisfying conditional beta and attaining union sum."""
    # Data A fails the parameter event, probability alpha. Numerical failure
    # occurs with conditional probability zero for A and beta for B. The two
    # events are dependent while every conditional budget remains <= beta.
    # This yields 1-alpha-beta+alpha*beta >= 1-alpha-beta.
    alpha, beta = Fraction(1,10), Fraction(1,20)
    rows = [
      {'data':'A','weight':alpha,'parameter_ok':False,'numerical_failure_given_data':Fraction(0)},
      {'data':'B','weight':1-alpha,'parameter_ok':True,'numerical_failure_given_data':beta},
    ]
    good_joint = sum((r['weight']*(1-r['numerical_failure_given_data'])
                      for r in rows if r['parameter_ok']), Fraction())
    return {'alpha':str(alpha),'beta':str(beta),'good_joint':str(good_joint),
            'union_lower_bound':str(1-alpha-beta),'bound_holds':good_joint>=1-alpha-beta}
