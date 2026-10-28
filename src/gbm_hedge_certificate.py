"""Validated population static hedging for the retained two-instrument GBM study.

Exact rational inputs define the financial target; Arb encloses its prices and
moments. Feasible rational holdings and weak-duality bounds certify the output.
The exact_fraction fields carry the certificate. Optional binary64 summaries
are rounded outwards and may be unavailable outside their conversion range.
Established Gaussian moments and convex duality are specialised here, without
a claim of theoretical priority, empirical validation or quantum advantage.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import math
from typing import Any

import flint
from flint import arb, ctx
from gbm_price_certificate import (ExactCappedGBM, rational, _arb, _bounds,
                                   _normal_cdf, _options, _serial)

F = Fraction


def _pair(value, bits):
    return _bounds(value, bits)


def _mid(value, bits):
    lower, upper = _pair(value, bits)
    return (lower + upper) / 2


def _vector(values, count, name):
    if isinstance(values, (str, bytes)):
        raise ValueError(f'{name} requires {count} rational entries.')
    try:
        result = tuple(rational(v, name) for v in values)
    except TypeError as error:
        raise ValueError(f'{name} requires {count} rational entries.') from error
    if len(result) != count:
        raise ValueError(f'{name} requires {count} rational entries.')
    return result


def _target(contract, sigma, mu):
    if not isinstance(contract, ExactCappedGBM):
        raise TypeError('contract must be ExactCappedGBM.')
    sigma, mu = rational(sigma, 'volatility'), rational(mu, 'physical_total_return')
    if sigma <= 0:
        raise ValueError('volatility must be strictly positive.')
    return sigma, mu


def _moments(c, sigma, mu):
    """Gain Gram matrix in (stock, lower call, absent upper call) basis.

    R(theta)=-theta_s*g_s+(1-theta_c)*g_K-g_J preserves the rare upper-tail
    component when a near-unit lower call almost replicates the capped claim.
    The absent upper call is an analytic representation, not a tradable hedge.
    """
    s, k, h, t, r, q, v = map(_arb, (c.spot, c.strike, c.cap, c.maturity,
                                     c.rate, c.dividend, sigma))
    j = _arb(c.strike + c.cap)
    discount, stock_factor = (-r*t).exp(), ((q-r)*t).exp()
    mean = s.log() + (_arb(mu)-q-v*v/2)*t
    sd = v*t.sqrt()

    def above(power, endpoint=None):
        coefficient = (power*mean + power*power*sd*sd/2).exp()
        if endpoint is None:
            return coefficient
        return coefficient*_normal_cdf((mean+power*sd*sd-endpoint.log())/sd)

    def call_cost(strike):
        d1 = ((s/strike).log()+(r-q+v*v/2)*t)/sd
        return s*(-q*t).exp()*_normal_cdf(d1) - strike*discount*_normal_cdf(d1-sd)

    ck, cj = call_cost(k), call_cost(j)
    premium = ck-cj
    totals = [above(power) for power in range(3)]
    lower = [above(power, k) for power in range(3)]
    upper = [above(power, j) for power in range(3)]
    emk, emj = lower[1]-k*lower[0], upper[1]-j*upper[0]
    eqk = lower[2]-2*k*lower[1]+k*k*lower[0]
    eqj = upper[2]-2*j*upper[1]+j*j*upper[0]
    exk, exj = lower[2]-k*lower[1], upper[2]-j*upper[1]
    ekj = upper[2]-(k+j)*upper[1]+k*j*upper[0]
    gs_mean = stock_factor*totals[1]-s
    gk_mean = discount*emk-ck
    gj_mean = discount*emj-cj
    ss = stock_factor*stock_factor*totals[2]-2*stock_factor*s*totals[1]+s*s
    sk = stock_factor*discount*exk-stock_factor*ck*totals[1]-s*discount*emk+s*ck
    sj = stock_factor*discount*exj-stock_factor*cj*totals[1]-s*discount*emj+s*cj
    kk = discount*discount*eqk-2*discount*ck*emk+ck*ck
    jj = discount*discount*eqj-2*discount*cj*emj+cj*cj
    kj = discount*discount*ekj-discount*cj*emk-discount*ck*emj+ck*cj
    gram = ((ss, sk, sj), (sk, kk, kj), (sj, kj, jj))
    matrix = ((ss, sk), (sk, kk))
    cross = (sk-sj, kk-kj)
    constant = kk-2*kj+jj
    return {'gram': gram, 'matrix': matrix, 'cross': cross, 'constant': constant,
            'call_cost': ck, 'upper_call_cost': cj, 'premium': premium,
            'gain_means': (gs_mean, gk_mean, gj_mean)}


def _risk(model, holdings):
    weights = (-_arb(holdings[0]), 1-_arb(holdings[1]), arb(-1))
    return sum((weights[i]*weights[j]*model['gram'][i][j]
                for i in range(3) for j in range(3)), arb(0))


def _spd(model, bits):
    m = model['matrix']
    determinant = m[0][0]*m[1][1]-m[0][1]*m[1][0]
    trace = m[0][0]+m[1][1]
    return _pair(m[0][0], bits), _pair(determinant, bits), _pair(trace, bits)


def _solve(matrix, rhs):
    determinant = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    if determinant == 0:
        raise ZeroDivisionError('Singular rational matrix.')
    return ((matrix[1][1]*rhs[0]-matrix[0][1]*rhs[1])/determinant,
            (matrix[0][0]*rhs[1]-matrix[1][0]*rhs[0])/determinant)


def _dot(first, second):
    return sum((a*b for a, b in zip(first, second)), F(0))


def _matvec(matrix, vector):
    return tuple(_dot(row, vector) for row in matrix)


def _repair(holdings, costs_upper, budget_lower):
    box = tuple(min(F(1), max(F(0), x)) for x in holdings)
    spend = _dot(costs_upper, box)
    factor = F(1) if spend <= budget_lower else budget_lower/spend
    return tuple(factor*x for x in box), factor, box


def _constraints(costs, budget):
    return ((F(-1), F(0)), (F(0), F(-1)), (F(1), F(0)),
            (F(0), F(1)), costs), (F(0), F(0), F(1), F(1), budget)


def _refine(model, costs_upper, budget_lower, bits):
    """Exact rational proposal for a conservative nearby two-dimensional QP.

    Midpoint coefficients generate proposals only. True-model risk, feasibility
    and dual bounds are checked separately, so proposal optimality is not used
    as a substitute for a certificate for the original financial problem.
    """
    m = tuple(tuple(_mid(value, bits) for value in row) for row in model['matrix'])
    h = tuple(_mid(value, bits) for value in model['cross'])
    if m[0][0] <= 0 or m[0][0]*m[1][1]-m[0][1]*m[1][0] <= 0:
        raise ArithmeticError('Midpoint proposal matrix is not positive definite.')
    a, b = _constraints(costs_upper, budget_lower)
    candidates = [_solve(m, h)]
    for normal, value in zip(a, b):
        centre = tuple(x*value/_dot(normal, normal) for x in normal)
        tangent = (-normal[1], normal[0])
        step = (_dot(tangent, h)-_dot(tangent, _matvec(m, centre))) / _dot(tangent, _matvec(m, tangent))
        candidates.append(tuple(centre[i]+step*tangent[i] for i in range(2)))
    for i, j in combinations(range(5), 2):
        try:
            candidates.append(_solve((a[i], a[j]), (b[i], b[j])))
        except ZeroDivisionError:
            continue
    feasible = [x for x in candidates if all(_dot(row, x) <= value for row, value in zip(a, b))]
    if not feasible:
        raise ArithmeticError('No rational feasible proposal was found.')
    theta = min(feasible, key=lambda x: _dot(x, _matvec(m, x))-2*_dot(h, x))
    gradient = tuple(2*(_matvec(m, theta)[i]-h[i]) for i in range(2))
    active = [i for i in range(5) if _dot(a[i], theta) == b[i]]
    proposals = []
    if gradient == (F(0), F(0)):
        proposals.append((F(0),)*5)
    for i in active:
        normal = a[i]
        scale = -_dot(normal, gradient)/_dot(normal, normal)
        if scale >= 0 and all(gradient[j]+scale*normal[j] == 0 for j in range(2)):
            lam = [F(0)]*5
            lam[i] = scale
            proposals.append(tuple(lam))
    for i, j in combinations(active, 2):
        try:
            pair = _solve(((a[i][0], a[j][0]), (a[i][1], a[j][1])), tuple(-x for x in gradient))
        except ZeroDivisionError:
            continue
        if min(pair) >= 0:
            lam = [F(0)]*5
            lam[i], lam[j] = pair
            proposals.append(tuple(lam))
    if not proposals:
        raise ArithmeticError('No exact non-negative proposal multipliers found.')
    return theta, proposals[0], len(candidates)


def _dual(model, contract, kappa, multipliers):
    m, h, c = model['matrix'], model['cross'], model['constant']
    lam = tuple(_arb(x) for x in multipliers)
    vector = (-lam[0]+lam[2]+lam[4]*_arb(contract.spot),
              -lam[1]+lam[3]+lam[4]*model['call_cost'])
    z = (h[0]-vector[0]/2, h[1]-vector[1]/2)
    determinant = m[0][0]*m[1][1]-m[0][1]*m[1][0]
    solved = ((m[1][1]*z[0]-m[0][1]*z[1])/determinant,
              (m[0][0]*z[1]-m[1][0]*z[0])/determinant)
    return c-lam[2]-lam[3]-lam[4]*_arb(kappa)*model['premium']-(z[0]*solved[0]+z[1]*solved[1])


def _float_outward(value, upwards):
    """Ancillary outward binary64 conversion; fractions remain authoritative."""
    try:
        result = float(value)
    except OverflowError:
        return None
    if not math.isfinite(result):
        return None
    exact = F.from_float(result)
    if (upwards and exact < value) or (not upwards and exact > value):
        result = math.nextafter(result, math.inf if upwards else -math.inf)
    return result if math.isfinite(result) else None


def _record_risk(interval):
    return {**_serial(interval), 'lower_float': _float_outward(interval[0], False),
            'upper_float': _float_outward(interval[1], True),
            'strictly_positive': interval[0] > 0,
            'float_summary_rule': 'Outward binary64 summaries when finite; exact fraction endpoints are authoritative.'}


def _budget_status(violation):
    if violation[0] > 0:
        return 'certified_infeasible'
    if violation[1] <= 0:
        return 'certified_feasible'
    return 'unresolved'


def certified_hedge_risk(contract: ExactCappedGBM, volatility: Any, physical_total_return: Any,
                         holdings: Any, *, absolute_tolerance: Any = '1e-50',
                         relative_tolerance: Any = '1e-12', start_bits: int = 128,
                         max_bits: int = 4096) -> dict[str, Any]:
    """Certify exact physical-law risk at arbitrary rational fixed holdings.

    Positive lower endpoints and both absolute/relative width requirements are
    attained before return. A finite budget can fail explicitly; no zero floor
    or unvalidated approximation is returned as an attained certificate.
    """
    sigma, mu = _target(contract, volatility, physical_total_return)
    theta = _vector(holdings, 2, 'holdings')
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    relative = rational(relative_tolerance, 'relative_tolerance')
    if relative <= 0:
        raise ValueError('relative_tolerance must be strictly positive.')
    bits = start_bits
    while True:
        with ctx.workprec(bits):
            model = _moments(contract, sigma, mu)
            output_bits = max(endpoint_bits, bits+16)
            bounds = _pair(_risk(model, theta), output_bits)
        if (bounds[0] > 0 and bounds[1]-bounds[0] <= tolerance
                and bounds[1]-bounds[0] <= relative*bounds[0]):
            return {'status': 'certified', 'schema': 'qbp-validated-static-hedge-risk-v1',
                    'contract_exact': contract.exact_inputs(), 'volatility_exact': str(sigma),
                    'physical_total_return_exact': str(mu),
                    'holdings_exact': [str(x) for x in theta], 'risk': _record_risk(bounds),
                    'precision_bits': bits, 'endpoint_dyadic_bits': output_bits,
                    'requested_width_exact': str(tolerance),
                    'relative_width_tolerance_exact': str(relative),
                    'python_flint_version': flint.__version__,
                    'target': 'Exact physical-law E[(Y-V-theta^T g)^2] at declared rational parameters.'}
        if bits == max_bits:
            raise ArithmeticError('Positive risk enclosure of requested width was not attained.')
        bits = min(2*bits, max_bits)


def certify_hedge(contract: ExactCappedGBM, volatility: Any, physical_total_return: Any,
                  budget_multiplier: Any, holdings: Any, multipliers: Any = None, *,
                  absolute_tolerance: Any = '1e-50', relative_tolerance: Any = '1e-12',
                  start_bits: int = 128, max_bits: int = 4096,
                  refine: bool = True) -> dict[str, Any]:
    """Audit original holdings, robust repair and an optional refined hedge.

    Candidate holdings need not satisfy the constraints. A box projection and
    a common rational budget scaling produce a separately reported feasible
    portfolio. When refine=True, proposals solve a conservative rational QP;
    true-law primal/dual inequalities independently certify their accuracy.
    All raw multipliers and negative entries remain reported; only a separately
    named non-negative multiplier vector is used for the dual lower bound.
    """
    sigma, mu = _target(contract, volatility, physical_total_return)
    kappa = rational(budget_multiplier, 'budget_multiplier')
    if kappa < 0:
        raise ValueError('budget_multiplier must be non-negative.')
    if not isinstance(refine, bool):
        raise ValueError('refine must be boolean.')
    original = _vector(holdings, 2, 'holdings')
    raw_lam = (F(0),)*5 if multipliers is None else _vector(multipliers, 5, 'multipliers')
    safe_lam = tuple(max(F(0), x) for x in raw_lam)
    tolerance, endpoint_bits = _options(absolute_tolerance, start_bits, max_bits)
    relative = rational(relative_tolerance, 'relative_tolerance')
    if relative <= 0:
        raise ValueError('relative_tolerance must be strictly positive.')
    bits = start_bits
    while True:
        try:
            with ctx.workprec(bits):
                output_bits = max(endpoint_bits, bits+16)
                model = _moments(contract, sigma, mu)
                ck = _pair(model['call_cost'], output_bits)
                cj = _pair(model['upper_call_cost'], output_bits)
                premium = _pair(model['premium'], output_bits)
                budget = (kappa*premium[0], kappa*premium[1])
                first, determinant, trace = _spd(model, output_bits)
                if min(ck[0], cj[0], premium[0], first[0], determinant[0], trace[0]) <= 0:
                    raise ArithmeticError('Required positive model quantities unresolved.')
                m_lower = determinant[0]/trace[1]
                costs = (contract.spot, ck[1])
                repaired, factor, boxed = _repair(original, costs, budget[0])
                original_risk = _pair(_risk(model, original), output_bits)
                repaired_risk = _pair(_risk(model, repaired), output_bits)
                repaired_dual = _pair(_dual(model, contract, kappa, safe_lam), output_bits)
                unrestricted = _pair(_dual(model, contract, kappa, (F(0),)*5), output_bits)
                if refine:
                    theta, lam, count = _refine(model, costs, budget[0], output_bits)
                else:
                    theta, lam, count = repaired, safe_lam, 0
                risk = _pair(_risk(model, theta), output_bits)
                dual = _pair(_dual(model, contract, kappa, lam), output_bits)
                intervals = (ck, cj, premium, original_risk, repaired_risk,
                             repaired_dual, unrestricted, risk, dual)
                if (max(pair[1]-pair[0] for pair in intervals) > tolerance
                    or min(original_risk[0], repaired_risk[0], risk[0], unrestricted[0]) <= 0
                    or any(pair[1]-pair[0] > relative*pair[0]
                           for pair in (original_risk, repaired_risk, risk, unrestricted))):
                    raise ArithmeticError('Requested positive risk/dual widths unresolved.')
                original_spend = _arb(original[0])*_arb(contract.spot)+_arb(original[1])*model['call_cost']
                original_violation = _pair(original_spend-_arb(kappa)*model['premium'], output_bits)
                budget_status = _budget_status(original_violation)
                lower_box = tuple(max(F(0), -x) for x in original)
                upper_box = tuple(max(F(0), x-1) for x in original)
                box_violation = max((*lower_box, *upper_box))
                original_status = ('certified_infeasible' if box_violation > 0
                                   else budget_status)
                spend = _arb(theta[0])*_arb(contract.spot)+_arb(theta[1])*model['call_cost']
                slack = _pair(_arb(kappa)*model['premium']-spend, output_bits)
                lower = max(F(0), dual[0], unrestricted[0])
                gap = risk[1]-lower
                if gap < 0:
                    raise ArithmeticError('Inconsistent primal/dual bounds.')
                if refine and gap > relative*risk[0]:
                    raise ArithmeticError('Refined relative optimality gap unresolved.')
                return {
                    'status': 'certified', 'schema': 'qbp-validated-static-hedge-v1',
                    'contract_exact': contract.exact_inputs(), 'volatility_exact': str(sigma),
                    'physical_total_return_exact': str(mu), 'budget_multiplier_exact': str(kappa),
                    'premium': _serial(premium), 'call_cost': _serial(ck),
                    'upper_call_cost': _serial(cj), 'budget': _serial(budget),
                    'original_holdings_exact': [str(x) for x in original],
                    'original_risk': _record_risk(original_risk),
                    'original_budget_violation': _serial(original_violation),
                    'original_budget_status': budget_status,
                    'original_lower_box_violations_exact': [str(x) for x in lower_box],
                    'original_upper_box_violations_exact': [str(x) for x in upper_box],
                    'original_box_violation_exact': str(box_violation),
                    'original_feasibility_status': original_status,
                    'boxed_original_holdings_exact': [str(x) for x in boxed],
                    'repair_common_scale_exact': str(factor),
                    'repaired_original_holdings_exact': [str(x) for x in repaired],
                    'repair_difference_exact': [str(repaired[i]-original[i]) for i in range(2)],
                    'repaired_original_risk': _record_risk(repaired_risk),
                    'repaired_original_dual': _serial(repaired_dual),
                    'raw_multipliers_exact': [str(x) for x in raw_lam],
                    'safe_multipliers_exact': [str(x) for x in safe_lam],
                    'holdings_exact': [str(x) for x in theta], 'risk': _record_risk(risk),
                    'refined_proposal': refine, 'proposal_candidates': count,
                    'multipliers_exact': [str(x) for x in lam],
                    'positive_definiteness': {'M11': _serial(first), 'determinant': _serial(determinant),
                                             'trace': _serial(trace)},
                    'minimum_eigenvalue_lower_exact': str(m_lower),
                    'unrestricted_loss': _serial(unrestricted), 'raw_dual': _serial(dual),
                    'nonnegative_risk_lower_bound_exact': '0',
                    'optimal_loss_lower_bound_exact': str(lower),
                    'optimality_gap_upper_exact': str(gap),
                    'optimality_gap_upper_float': _float_outward(gap, True),
                    'squared_distance_to_unique_optimum_upper_exact': str(gap/m_lower),
                    'robust_spend_upper_exact': str(_dot(costs, theta)),
                    'robust_budget_lower_exact': str(budget[0]), 'actual_budget_slack': _serial(slack),
                    'robust_feasibility': all(F(0) <= x <= F(1) for x in theta)
                                          and _dot(costs, theta) <= budget[0],
                    'precision_bits': bits, 'endpoint_dyadic_bits': output_bits,
                    'requested_width_exact': str(tolerance),
                    'relative_width_tolerance_exact': str(relative),
                    'python_flint_version': flint.__version__,
                    'target': 'Exact two-instrument constrained physical-law MSE with B=kappa*V, rational parameters and true transcendental prices.',
                    'scope': 'Known-law arithmetic/optimality certificate; not empirical/model validation or a quantum-advantage claim.'}
        except (ArithmeticError, ZeroDivisionError):
            if bits == max_bits:
                raise ArithmeticError('Requested hedge certificate was not attained by max_bits.') from None
            bits = min(2*bits, max_bits)
