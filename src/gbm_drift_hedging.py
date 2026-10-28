"""Continuum drift-robust funded static hedging under a known GBM volatility.

The proposal optimiser uses ordinary floating arithmetic. Every reported
primal/dual certificate is instead evaluated at exact rational parameters and
holdings, using the existing Arb gain moments. The continuum upper bound uses
a Gaussian density interpolation inequality, not an endpoint assumption.
This application of convex duality and Holder's inequality makes no claim of
theoretical priority, empirical validation or a quantum advantage.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Any
import math

import flint
from flint import arb, ctx
import numpy as np
from scipy.optimize import minimize

from gbm_price_certificate import (
    ExactCappedGBM, rational, _arb, _bounds, _serial, _integer,
    _rational_sqrt_bounds, certified_chi_square_quantile,
)
from gbm_hedge_certificate import (
    _moments, _risk, _repair, _refine, _constraints, _dual,
    _dot, _mid, _record_risk,
)

F = Fraction


def _midfloat(value, bits):
    return float(_mid(value, bits))


def _weighted(models, weights):
    def total(key, index=None):
        values = [model[key] if index is None else model[key][index]
                  for model in models]
        return sum((_arb(w)*v for w, v in zip(weights, values)), arb(0))
    matrix = tuple(tuple(sum((_arb(w)*m['matrix'][i][j]
                              for w, m in zip(weights, models)), arb(0))
                         for j in range(2)) for i in range(2))
    return {'matrix': matrix, 'cross': tuple(total('cross', i) for i in range(2)),
            'constant': total('constant'), 'call_cost': models[0]['call_cost'],
            'premium': models[0]['premium']}


def _normalise_nonnegative(values):
    exact = tuple(max(F(0), F.from_float(float(v))) for v in values)
    total = sum(exact, F(0))
    if total == 0:
        raise ArithmeticError('No positive mixture weight was proposed.')
    return tuple(x/total for x in exact)


def _proposal(models, contract, costs, budget, bits):
    """Produce feasible proposals and dual witnesses without trusting a solver."""
    count = len(models)
    indices = sorted({0, count//2, count-1})
    point_proposals = []
    for index in indices:
        theta, lam, _ = _refine(models[index], costs, budget, bits)
        point_proposals.append((theta, lam, index))

    def worst(theta):
        return max(_mid(_risk(m, theta), bits) for m in models)

    anchor, _, _ = min(point_proposals, key=lambda item: worst(item[0]))
    # An additional worst-node minimiser protects against a poor central anchor.
    worst_index = max(range(count), key=lambda i: _mid(_risk(models[i], anchor), bits))
    if worst_index not in indices:
        theta, lam, _ = _refine(models[worst_index], costs, budget, bits)
        point_proposals.append((theta, lam, worst_index))
        anchor, _, _ = min(point_proposals, key=lambda item: worst(item[0]))

    risk0 = np.array([_midfloat(_risk(m, anchor), bits) for m in models])
    scale = max(float(np.max(risk0)), np.finfo(float).tiny)
    matrices = np.array([[[_midfloat(x, bits)/scale for x in row]
                         for row in m['matrix']] for m in models])
    # Compute h-M*anchor inside Arb before conversion to avoid cancellation.
    cross = np.array([[_midfloat(m['cross'][i]-sum(
        (m['matrix'][i][j]*_arb(anchor[j]) for j in range(2)), arb(0)), bits)/scale
                       for i in range(2)] for m in models])
    constants = risk0/scale
    anchor_float = np.array([float(x) for x in anchor])
    cost_float = np.array([float(x/contract.spot) for x in costs])
    slack_float = float((budget-_dot(costs, anchor))/contract.spot)

    def losses(delta):
        return constants-2*cross@delta+np.einsum('i,nij,j->n', delta, matrices, delta)

    def constraint(x):
        return x[2]-losses(x[:2])

    def jac_constraint(x):
        derivative = -2*cross+2*np.einsum('nij,j->ni', matrices, x[:2])
        return np.column_stack((-derivative, np.ones(count)))

    financial = {'type': 'ineq', 'fun': lambda x: slack_float-cost_float@x[:2],
                 'jac': lambda x: np.array([-cost_float[0], -cost_float[1], 0.0])}
    result = minimize(lambda x: x[2], np.array([0., 0., max(1., max(constants))]),
                      jac=lambda x: np.array([0., 0., 1.]), method='SLSQP',
                      bounds=[(-anchor_float[i], 1-anchor_float[i]) for i in range(2)]
                             +[(0., None)],
                      constraints=[{'type': 'ineq', 'fun': constraint, 'jac': jac_constraint}, financial],
                      options={'ftol': 1e-12, 'maxiter': 500})
    candidate = anchor
    if np.all(np.isfinite(result.x[:2])):
        raw = tuple(anchor[i]+F.from_float(float(result.x[i])) for i in range(2))
        repaired, _, _ = _repair(raw, costs, budget)
        if worst(repaired) < worst(candidate):
            candidate = repaired

    # Sparse mixture proposals include the largest node risks. The exact dual
    # remains valid even if this numerical maximisation fails or is inaccurate.
    at_candidate = np.array([_midfloat(_risk(m, candidate), bits) for m in models])
    top = sorted(range(count), key=lambda i: at_candidate[i], reverse=True)[:8]
    selected = sorted(set(top+[0, count//2, count-1]))
    a, b = _constraints(costs, budget)
    normal = np.array([[float(x) for x in row] for row in a])
    financial_slack = np.array([float(bi-_dot(ai, anchor)) for ai, bi in zip(a, b)])
    cm = constants[selected]
    hm = cross[selected]
    mm = matrices[selected]
    size = len(selected)

    def negative_dual(x):
        weights, lam = x[:size], x[size:]
        matrix = np.einsum('n,nij->ij', weights, mm)
        z = weights@hm-normal.T@lam/2
        try:
            value = weights@cm-financial_slack@lam-z@np.linalg.solve(matrix, z)
        except np.linalg.LinAlgError:
            return 1e100
        return -float(value)

    initial = np.zeros(size+5)
    initial[selected.index(int(np.argmax(at_candidate)))] = 1.
    dual_result = minimize(negative_dual, initial, method='SLSQP',
                           bounds=[(0., 1.)]*size+[(0., None)]*5,
                           constraints=[{'type': 'eq', 'fun': lambda x: np.sum(x[:size])-1.,
                                         'jac': lambda x: np.r_[np.ones(size), np.zeros(5)]}],
                           options={'ftol': 1e-12, 'maxiter': 500})
    dual_proposals = []
    for _, lam, index in point_proposals:
        weights = [F(0)]*count
        weights[index] = F(1)
        dual_proposals.append((tuple(weights), lam, f'point_node_{index}'))
    if np.all(np.isfinite(dual_result.x)):
        try:
            sparse = _normalise_nonnegative(dual_result.x[:size])
            weights = [F(0)]*count
            for index, weight in zip(selected, sparse):
                weights[index] = weight
            lam = tuple(max(F(0), F.from_float(float(v))*F.from_float(scale))
                        for v in dual_result.x[size:])
            dual_proposals.append((tuple(weights), lam, 'normalised_floating_mixture'))
        except ArithmeticError:
            pass
    diagnostics = {'primal_solver_success': bool(result.success),
                   'primal_solver_message': str(result.message),
                   'dual_solver_success': bool(dual_result.success),
                   'dual_solver_message': str(dual_result.message),
                   'dual_selected_node_indices': selected,
                   'raw_mixture_weights_binary64_exact': [
                       str(F.from_float(float(x))) if math.isfinite(float(x)) else None
                       for x in dual_result.x[:size]],
                   'raw_financial_multipliers_solver_units_binary64_exact': [
                       str(F.from_float(float(x))) if math.isfinite(float(x)) else None
                       for x in dual_result.x[size:]],
                   'risk_scaling_binary64_exact': str(F.from_float(scale)),
                   'mixture_normalisation_rule': 'Replace negative floating weights by zero, convert each binary64 value to its exact rational, and divide by their exact positive sum. Financial multipliers are converted to exact rationals, multiplied by the exact binary64 risk scale and replaced by zero if negative.',
                   'proposal_role': 'Floating proposals only; success flags do not establish any certificate.'}
    return candidate, anchor, dual_proposals, diagnostics


def certify_drift_robust_hedge(contract: ExactCappedGBM, volatility: Any,
                              drift_lower: Any, drift_upper: Any,
                              budget_multiplier: Any, *,
                              relative_gap_tolerance: Any = '0.005',
                              cell_inflation_tolerance: Any = '0.001',
                              initial_cells: int = 8, max_nodes: int = 513,
                              precision_bits: int = 256) -> dict[str, Any]:
    """Bound the continuum minimax MSE with a feasible hedge and a mixture dual.

    The volatility is known and common to P and Q. Only the physical total
    return mu varies, so all fair prices and funding constraints remain fixed.
    An unattained gap or grid-inflation tolerance returns ``unresolved`` with
    proved bounds. An arithmetic/positive-definiteness enclosure that cannot
    be resolved at the requested precision raises ``ArithmeticError``. Neither
    form of abstention replaces the continuum by an endpoint problem.
    """
    if not isinstance(contract, ExactCappedGBM):
        raise TypeError('contract must be ExactCappedGBM.')
    sigma = rational(volatility, 'volatility')
    lower, upper = rational(drift_lower, 'drift_lower'), rational(drift_upper, 'drift_upper')
    kappa = rational(budget_multiplier, 'budget_multiplier')
    tolerance = rational(relative_gap_tolerance, 'relative_gap_tolerance')
    inflation_tolerance = rational(cell_inflation_tolerance, 'cell_inflation_tolerance')
    _integer(initial_cells, 'initial_cells', 1)
    _integer(max_nodes, 'max_nodes', 2)
    _integer(precision_bits, 'precision_bits', 128)
    if sigma <= 0 or lower > upper or kappa < 0 or tolerance <= 0 or inflation_tolerance <= 0:
        raise ValueError('Require sigma>0, ordered drift bounds, kappa>=0 and positive tolerances.')
    cells = min(initial_cells, max_nodes-1) if lower < upper else 0
    with ctx.workprec(precision_bits):
        output_bits = precision_bits+16
        while True:
            grid = ([lower+(upper-lower)*F(i, cells) for i in range(cells+1)]
                    if cells else [lower])
            width = (upper-lower)/cells if cells else F(0)
            exponent = contract.maturity*width*width/(8*sigma*sigma)
            inflation = _bounds(_arb(exponent).exp(), output_bits)
            if inflation[1]-1 <= inflation_tolerance or cells >= max_nodes-1:
                break
            cells = min(2*cells, max_nodes-1)
        models = [_moments(contract, sigma, mu) for mu in grid]
        first = models[0]
        ck = _bounds(first['call_cost'], output_bits)
        premium = _bounds(first['premium'], output_bits)
        if ck[0] <= 0 or premium[0] <= 0:
            raise ArithmeticError('Positive fair price bounds unresolved at the requested precision.')
        costs, budget = (contract.spot, ck[1]), kappa*premium[0]
        theta, anchor, duals, diagnostic = _proposal(models, contract, costs, budget, output_bits)
        theta, repair_scale, boxed = _repair(theta, costs, budget)
        risks = [_bounds(_risk(m, theta), output_bits) for m in models]
        if any(lo < 0 or hi <= 0 for lo, hi in risks):
            raise ArithmeticError('Positive node-risk bounds unresolved at the requested precision.')
        continuum_upper = max(hi for _, hi in risks)*inflation[1]
        best = None
        for weights, lam, origin in duals:
            mixed = _weighted(models, weights)
            matrix = mixed['matrix']
            determinant = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
            m11 = _bounds(matrix[0][0], output_bits)
            det = _bounds(determinant, output_bits)
            trace = _bounds(matrix[0][0]+matrix[1][1], output_bits)
            if min(m11[0], det[0], trace[0]) <= 0:
                continue
            dual = _bounds(_dual(mixed, contract, kappa, lam), output_bits)
            if best is None or dual[0] > best[0][0]:
                best = (dual, weights, lam, mixed, m11, det, trace, origin)
        if best is None:
            raise ArithmeticError('No positive-definite mixture witness resolved.')
        dual, weights, lam, mixed, m11, det, trace, origin = best
        optimum_lower = max(F(0), dual[0])
        gap = continuum_upper-optimum_lower
        if gap < 0:
            raise ArithmeticError('Inconsistent primal/dual bounds.')
        relative_gap = gap/optimum_lower if optimum_lower > 0 else None
        numerical_width = max(hi-lo for lo, hi in risks)
        certified = (relative_gap is not None and relative_gap <= tolerance
                     and inflation[1]-1 <= inflation_tolerance)
        result = {
            'schema': 'qbp-validated-drift-robust-hedge-v1',
            'status': 'certified' if certified else 'unresolved',
            'contract_exact': contract.exact_inputs(), 'volatility_exact': str(sigma),
            'drift_interval': _serial((lower, upper)), 'budget_multiplier_exact': str(kappa),
            'node_grid_exact': [str(x) for x in grid], 'grid_nodes': len(grid),
            'max_cell_width_exact': str(width), 'cell_exponent_exact': str(exponent),
            'cell_inflation': _serial(inflation),
            'max_cell_inflation_upper_exact': str(inflation[1]),
            'anchor_holdings_exact': [str(x) for x in anchor],
            'holdings_exact': [str(x) for x in theta],
            'final_repair_scale_exact': str(repair_scale),
            'final_boxed_holdings_exact': [str(x) for x in boxed],
            'call_cost': _serial(ck), 'premium': _serial(premium),
            'robust_spend_upper_exact': str(_dot(costs, theta)),
            'robust_budget_lower_exact': str(budget),
            'robust_feasibility': all(F(0) <= x <= F(1) for x in theta)
                                  and _dot(costs, theta) <= budget,
            'node_risks': [_record_risk(x) for x in risks],
            'node_maximum_lower_exact': str(max(lo for lo, _ in risks)),
            'continuum_upper_exact': str(continuum_upper),
            'mixture_weights_exact': [str(x) for x in weights],
            'financial_multipliers_exact': [str(x) for x in lam],
            'selected_dual_witness_origin': origin,
            'zero_budget_forces_zero_holdings': kappa == 0,
            'mixture_quadratic': {
                'matrix': [[_serial(_bounds(x, output_bits)) for x in row]
                           for row in mixed['matrix']],
                'cross': [_serial(_bounds(x, output_bits)) for x in mixed['cross']],
                'constant': _serial(_bounds(mixed['constant'], output_bits))},
            'mixture_positive_definiteness': {'M11': _serial(m11),
                                            'determinant': _serial(det), 'trace': _serial(trace)},
            'raw_mixture_dual': _serial(dual),
            'nonnegative_risk_lower_bound_exact': '0',
            'optimal_loss_lower_bound_exact': str(optimum_lower),
            'optimality_gap_upper_exact': str(gap),
            'relative_gap_upper_exact': None if relative_gap is None else str(relative_gap),
            'maximum_node_risk_width_exact': str(numerical_width),
            'relative_gap_tolerance_exact': str(tolerance),
            'cell_inflation_tolerance_exact': str(inflation_tolerance),
            'max_nodes': max_nodes, 'precision_bits': precision_bits,
            'endpoint_dyadic_bits': output_bits, 'python_flint_version': flint.__version__,
            'optimiser': diagnostic,
            'continuum_argument': 'Gaussian density interpolation and Holder imply L_mu <= exp(T*(u-l)^2/(8*sigma^2))*max(L_l,L_u) on every cell.',
            'target': 'inf over funded box-constrained real holdings of sup_mu E_Pmu[(Y-V-theta^T g)^2], with known common volatility and fixed Q quotes; the feasible reported holdings are exact rationals.',
            'scope': 'Known-volatility GBM conditional arithmetic/optimisation bounds; no empirical model validation, theoretical-priority or quantum-advantage claim.',
        }
        if not certified:
            result['unresolved_reason'] = ('Relative minimax gap or requested cell inflation was not attained within the declared numerical budgets.')
        return result


def certified_known_sigma_drift_interval(observations: Any, volatility: Any,
                                        observation_spacing: Any, *,
                                        dividend: Any = '0', rounding_radii: Any = None,
                                        alpha: Any = '0.05', precision_bits: int = 256) -> dict[str, Any]:
    """Gaussian known-volatility confidence interval from log-return increments.

    Independently observed latent returns satisfy R_i ~ N((mu-q-sigma^2/2)Delta,
    sigma^2 Delta). Supplied rational observations z_i obey |R_i-z_i|<=rho_i.
    A conservative chi-square(1) upper quantile yields a simultaneous outward
    enclosure of the exact normal mean interval. Delta is historical spacing;
    it need not equal the derivative's maturity.
    """
    sigma = rational(volatility, 'volatility')
    spacing = rational(observation_spacing, 'observation_spacing')
    q = rational(dividend, 'dividend')
    significance = rational(alpha, 'alpha')
    _integer(precision_bits, 'precision_bits', 128)
    if sigma <= 0 or spacing <= 0 or not 0 < significance < 1:
        raise ValueError('Require sigma>0, observation_spacing>0 and 0<alpha<1.')
    if isinstance(observations, (str, bytes)):
        raise ValueError('observations must contain log-return increments.')
    try:
        values = tuple(rational(x, 'observation') for x in observations)
    except TypeError as error:
        raise ValueError('observations must contain log-return increments.') from error
    if not values:
        raise ValueError('At least one log-return increment is required.')
    if rounding_radii is None:
        radii = (F(0),)*len(values)
    else:
        if isinstance(rounding_radii, (str, bytes)):
            raise ValueError('rounding_radii must match observations.')
        try:
            radii = tuple(rational(x, 'rounding_radius') for x in rounding_radii)
        except TypeError as error:
            raise ValueError('rounding_radii must match observations.') from error
        if len(radii) != len(values) or any(x < 0 for x in radii):
            raise ValueError('Nonnegative rounding_radii must match observations.')
    n = len(values)
    observed_mean = sum(values, F(0))/n
    mean_radius = sum(radii, F(0))/n
    centre = q+sigma*sigma/2+observed_mean/spacing
    centre_radius = mean_radius/spacing
    quantile = certified_chi_square_quantile(
        1, 1-significance, absolute_tolerance=F(1, 1 << (precision_bits//2)),
        start_bits=precision_bits, max_bits=max(precision_bits, 4096))
    quantile_upper = F(quantile['quantile']['upper_exact_fraction'])
    root = _rational_sqrt_bounds(quantile_upper/(n*spacing), precision_bits)
    half_width = sigma*root[1]
    interval = centre-centre_radius-half_width, centre+centre_radius+half_width
    return {'schema': 'qbp-validated-known-sigma-drift-interval-v1',
            'status': 'certified', 'observations': n,
            'observed_log_return_mean_exact': str(observed_mean),
            'mean_rounding_radius_exact': str(mean_radius),
            'observation_spacing_exact': str(spacing), 'volatility_exact': str(sigma),
            'dividend_exact': str(q), 'alpha_exact': str(significance),
            'nominal_coverage_exact': str(1-significance),
            'drift_centre_exact': str(centre), 'drift_centre_rounding_radius_exact': str(centre_radius),
            'gaussian_half_width_upper_exact': str(half_width),
            'drift_interval': _serial(interval), 'chi_square_quantile': quantile,
            'precision_bits': precision_bits,
            'sampling_premise': 'Independent latent Gaussian log-return increments with known common sigma and equally spaced observations; each rounding cell contains its latent increment.',
            'coverage_interpretation': 'At least 1-alpha coverage under the sampling premise. This does not validate GBM or estimate volatility.',
            'maturity_distinction': 'Historical observation spacing enters inference; derivative maturity enters the subsequent hedge problem separately.'}
