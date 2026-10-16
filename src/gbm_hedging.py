"""Population static hedging of a synthetic capped GBM claim.

Closed moments and active-set formulas implement exact-arithmetic results.
Floating-point moments, ranks, optimiser outputs and quadrature error estimates
are diagnostics, not arithmetic certificates. No eigenvalue or risk clipping.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math

import numpy as np
from scipy.integrate import quad, quad_vec
from scipy.optimize import minimize

from capped_gbm import Contract, analytic, call_analytic, normal_cdf, payoff


@dataclass(frozen=True)
class HedgeMoments:
    matrix: np.ndarray
    cross: np.ndarray
    constant: float
    gain_mean: np.ndarray
    residual_mean: float
    gain_covariance: np.ndarray
    centred_cross: np.ndarray
    residual_variance: float
    costs: np.ndarray
    premium: float


@dataclass(frozen=True)
class HedgeSolution:
    theta: np.ndarray
    loss_quadratic: float
    multipliers: np.ndarray
    active_labels: tuple[str, ...]
    candidates_considered: int
    primal_violation: float
    dual_violation: float
    stationarity_residual: float
    complementarity_residual: float
    matrix_eigenvalues: np.ndarray
    matrix_condition_number: float


def _law(c: Contract, physical_total_return: float) -> tuple[float, float]:
    if not math.isfinite(physical_total_return):
        raise ValueError("Physical total return must be finite.")
    return (math.log(c.spot)+(physical_total_return-c.dividend-.5*c.volatility**2)*c.maturity,
            c.volatility*math.sqrt(c.maturity))


def truncated_lognormal(c: Contract, physical_total_return: float,
                        power: int, lower: float) -> float:
    """E_P[X**power 1{X>lower}] for powers zero, one and two."""
    if power not in (0, 1, 2) or not math.isfinite(lower) or lower < 0:
        raise ValueError("Supported powers are 0,1,2 and the lower endpoint must be non-negative.")
    mean, sd = _law(c, physical_total_return)
    tilted_mean = mean+power*sd**2
    probability = 1.0 if lower == 0 else normal_cdf((tilted_mean-math.log(lower))/sd)
    return math.exp(power*mean+.5*power**2*sd**2)*probability


def hedge_moments(c: Contract, physical_total_return: float) -> HedgeMoments:
    """Uncentred mean-square coefficients; centred coefficients stay distinct."""
    k, upper = c.strike, c.strike+c.cap
    lower_moments = np.array([truncated_lognormal(c, physical_total_return, j, k) for j in range(3)])
    upper_moments = np.array([truncated_lognormal(c, physical_total_return, j, upper) for j in range(3)])
    interval = lower_moments-upper_moments
    ex = truncated_lognormal(c, physical_total_return, 1, 0)
    ex2 = truncated_lognormal(c, physical_total_return, 2, 0)
    call_mean = lower_moments[1]-k*lower_moments[0]
    call_second = lower_moments[2]-2*k*lower_moments[1]+k*k*lower_moments[0]
    x_call = lower_moments[2]-k*lower_moments[1]
    capped_mean = interval[1]-k*interval[0]+c.cap*upper_moments[0]
    capped_second = interval[2]-2*k*interval[1]+k*k*interval[0]+c.cap**2*upper_moments[0]
    x_capped = interval[2]-k*interval[1]+c.cap*upper_moments[1]
    call_capped = interval[2]-2*k*interval[1]+k*k*interval[0]+c.cap*(upper_moments[1]-k*upper_moments[0])
    discount = math.exp(-c.rate*c.maturity)
    stock_factor = math.exp((c.dividend-c.rate)*c.maturity)
    call_cost = call_analytic(c, c.strike)["price"]
    premium = analytic(c)["price"]
    costs = np.array([c.spot, call_cost])
    gain_mean = np.array([stock_factor*ex-c.spot, discount*call_mean-call_cost])
    claim_mean = discount*capped_mean
    matrix = np.array([
        [stock_factor**2*ex2-2*stock_factor*c.spot*ex+c.spot**2,
         stock_factor*discount*x_call-stock_factor*call_cost*ex-c.spot*discount*call_mean+c.spot*call_cost],
        [0.0, discount**2*call_second-2*discount*call_cost*call_mean+call_cost**2]])
    matrix[1, 0] = matrix[0, 1]
    cross = np.array([stock_factor*discount*x_capped-c.spot*claim_mean,
                      discount**2*call_capped-call_cost*claim_mean])-premium*gain_mean
    residual_mean = claim_mean-premium
    constant = discount**2*capped_second-2*premium*claim_mean+premium**2
    covariance = matrix-np.outer(gain_mean, gain_mean)
    centred_cross = cross-gain_mean*residual_mean
    return HedgeMoments(matrix, cross, constant, gain_mean, residual_mean,
                        covariance, centred_cross, constant-residual_mean**2, costs, premium)


def constraint_system(costs: np.ndarray, budget: float) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    p = np.asarray(costs, dtype=float)
    if p.shape != (2,) or not np.isfinite(p).all() or np.any(p <= 0):
        raise ValueError("Two strictly positive finite instrument costs are required.")
    if not math.isfinite(budget) or budget < 0:
        raise ValueError("Budget must be finite and non-negative.")
    return (np.array([[-1., 0.], [0., -1.], [1., 0.], [0., 1.], p]),
            np.array([0., 0., 1., 1., budget]),
            ("stock_lower", "call_lower", "stock_upper", "call_upper", "budget"))


def quadratic_loss(theta: np.ndarray, matrix: np.ndarray, cross: np.ndarray, constant: float) -> float:
    return float(constant-2*np.dot(cross, theta)+theta@matrix@theta)


def solve_constrained_hedge(matrix: np.ndarray, cross: np.ndarray, constant: float,
                            costs: np.ndarray, budget: float, tolerance: float = 2e-9) -> HedgeSolution:
    """Enumerate all stationary faces/vertices of a strictly convex 2D QP."""
    matrix, cross = np.asarray(matrix, dtype=float), np.asarray(cross, dtype=float)
    if (matrix.shape != (2, 2) or cross.shape != (2,) or not np.isfinite(matrix).all()
            or not np.isfinite(cross).all() or not math.isfinite(constant)
            or not np.allclose(matrix, matrix.T, atol=1e-12, rtol=1e-12)):
        raise ValueError("Finite symmetric two-dimensional quadratic coefficients are required.")
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues[0] <= 0:
        raise ValueError("The numerical matrix is not positive definite; no eigenvalue repair is performed.")
    a, b, labels = constraint_system(costs, budget)
    scales = np.maximum(1.0, np.linalg.norm(a, axis=1))
    candidates = [np.linalg.solve(matrix, cross)]
    for normal, value in zip(a, b):
        centre = normal*(value/float(normal@normal))
        tangent = np.array([-normal[1], normal[0]])
        candidates.append(centre+tangent*(float(tangent@(cross-matrix@centre))/float(tangent@matrix@tangent)))
    for first, second in combinations(range(5), 2):
        face = a[[first, second]]
        if abs(float(np.linalg.det(face))) > 1e-14*np.linalg.norm(face[0])*np.linalg.norm(face[1]):
            candidates.append(np.linalg.solve(face, b[[first, second]]))
    feasible = [theta for theta in candidates if np.max((a@theta-b)/scales) <= tolerance]
    if budget == 0.0:
        # Positive instrument costs and non-negative holdings imply the
        # singleton feasible set {0}; this is an analytic boundary case.
        feasible = [np.zeros(2)]
    if not feasible:
        raise ArithmeticError("No numerically feasible active-set candidate was found.")
    theta = min(feasible, key=lambda x: quadratic_loss(x, matrix, cross, constant))
    slack = a@theta-b
    active = tuple(i for i in range(5) if abs(slack[i]/scales[i]) <= tolerance)
    gradient = 2*(matrix@theta-cross)
    multiplier_candidates = []
    if np.max(np.abs(theta)) <= tolerance and budget == 0.0:
        budget_multiplier = max(0.0, *(2*cross/np.asarray(costs)))
        multipliers = np.array([budget_multiplier*costs[0]-2*cross[0],
                               budget_multiplier*costs[1]-2*cross[1], 0., 0., budget_multiplier])
        multiplier_candidates.append(multipliers)
    for count in (0, 1, 2):
        for support in combinations(active, count):
            multipliers = np.zeros(5)
            if support:
                solution, *_ = np.linalg.lstsq(a[list(support)].T, -gradient, rcond=None)
                multipliers[list(support)] = solution
            if multipliers.min() >= -tolerance:
                multiplier_candidates.append(multipliers)
    if not multiplier_candidates:
        raise ArithmeticError("No numerical non-negative KKT multiplier support was found.")
    multipliers = min(multiplier_candidates, key=lambda lam: np.linalg.norm(gradient+a.T@lam, ord=np.inf))
    return HedgeSolution(theta, quadratic_loss(theta, matrix, cross, constant), multipliers,
                          tuple(labels[i] for i in active), len(candidates),
                          max(0., float(np.max(slack/scales))), max(0., -float(multipliers.min())),
                          float(np.linalg.norm(gradient+a.T@multipliers, ord=np.inf)),
                          float(np.max(np.abs(multipliers*slack))), eigenvalues,
                          float(eigenvalues[-1]/eigenvalues[0]))


def slsqp_hedge(moments: HedgeMoments, budget: float) -> dict:
    """Independent general-purpose constrained solver; failures remain visible."""
    a, b, _ = constraint_system(moments.costs, budget)
    scale = max(1.0, abs(moments.constant), float(np.linalg.norm(moments.matrix)))
    result = minimize(lambda theta: quadratic_loss(theta, moments.matrix, moments.cross, moments.constant)/scale,
                      np.zeros(2), jac=lambda theta: 2*(moments.matrix@theta-moments.cross)/scale,
                      bounds=((0., 1.), (0., 1.)),
                      constraints={"type": "ineq", "fun": lambda theta: b[-1]-a[-1]@theta,
                                   "jac": lambda theta: -a[-1]},
                      method="SLSQP", options={"ftol": 1e-13, "maxiter": 1000})
    return {"theta": result.x, "success": bool(result.success), "status": int(result.status),
            "message": str(result.message), "iterations": int(result.nit),
            "loss_quadratic": quadratic_loss(result.x, moments.matrix, moments.cross, moments.constant),
            "primal_violation": max(0., float(np.max((a@result.x-b)/np.maximum(1., np.linalg.norm(a, axis=1)))))}


def _weighted_linear(mean: float, sd: float, z: float, slope: float, intercept: float) -> float:
    log_root_density = -.25*z*z-.25*math.log(2*math.pi)
    return slope*math.exp(mean+sd*z+log_root_density)+intercept*math.exp(log_root_density)


def residual_square_integral(c: Contract, physical_total_return: float,
                              theta: np.ndarray) -> tuple[float, float]:
    """Integrate a non-negative residual square independently of moment subtraction."""
    mean, sd = _law(c, physical_total_return)
    discount, stock_factor = math.exp(-c.rate*c.maturity), math.exp((c.dividend-c.rate)*c.maturity)
    premium, call_cost = analytic(c)["price"], call_analytic(c, c.strike)["price"]
    stock, call = np.asarray(theta, dtype=float)
    common = -premium+stock*c.spot+call*call_cost
    coefficients = [(-stock_factor*stock, common),
                    (discount*(1-call)-stock_factor*stock, common-discount*c.strike*(1-call)),
                    (-discount*call-stock_factor*stock, common+discount*c.cap+discount*c.strike*call)]
    boundaries = [-math.inf, (math.log(c.strike)-mean)/sd,
                  (math.log(c.strike+c.cap)-mean)/sd, math.inf]
    value, error = 0.0, 0.0
    for (slope, intercept), low, high in zip(coefficients, boundaries, boundaries[1:]):
        part, bound = quad(lambda z: _weighted_linear(mean, sd, z, slope, intercept)**2,
                           low, high, epsabs=1e-60, epsrel=3e-10, limit=300)
        value += part
        error += bound
    return value, error


def quadrature_moments(c: Contract, physical_total_return: float) -> tuple[np.ndarray, float]:
    """Independent numerical products: two gain means,M11,M12,M22,h1,h2,c."""
    mean, sd = _law(c, physical_total_return)
    discount, stock_factor = math.exp(-c.rate*c.maturity), math.exp((c.dividend-c.rate)*c.maturity)
    premium, call_cost = analytic(c)["price"], call_analytic(c, c.strike)["price"]
    boundaries = [-math.inf, (math.log(c.strike)-mean)/sd,
                  (math.log(c.strike+c.cap)-mean)/sd, math.inf]
    result = np.zeros(8)
    error = 0.0
    for region, (low, high) in enumerate(zip(boundaries, boundaries[1:])):
        def integrand(z):
            root_density = math.exp(-.25*z*z-.25*math.log(2*math.pi))
            stock = _weighted_linear(mean, sd, z, stock_factor, -c.spot)
            call = (_weighted_linear(mean, sd, z, discount, -discount*c.strike-call_cost)
                    if region else -call_cost*root_density)
            residual = (_weighted_linear(mean, sd, z, discount, -discount*c.strike-premium)
                        if region == 1 else ((discount*c.cap if region == 2 else 0)-premium)*root_density)
            return np.array([stock*root_density, call*root_density, stock*stock,
                              stock*call, call*call, stock*residual, call*residual, residual*residual])
        part, bound = quad_vec(integrand, low, high, epsabs=1e-10, epsrel=2e-10, limit=300)
        result += part
        error += float(bound)
    return result, error


def simulated_gains(c: Contract, physical_total_return: float, paths: int,
                     seed: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(paths, int) or paths < 2:
        raise ValueError("At least two independent simulated paths are required.")
    mean, sd = _law(c, physical_total_return)
    x = np.exp(mean+sd*np.random.default_rng(seed).standard_normal(paths))
    discount, stock_factor = math.exp(-c.rate*c.maturity), math.exp((c.dividend-c.rate)*c.maturity)
    gains = np.column_stack((stock_factor*x-c.spot,
                            discount*np.maximum(x-c.strike, 0)-call_analytic(c, c.strike)["price"]))
    residual = discount*payoff(x, c)-analytic(c)["price"]
    return gains, residual


def hedge_financing(c: Contract, moments: HedgeMoments, theta: np.ndarray) -> dict:
    cost = float(moments.costs@theta)
    cash = moments.premium-cost
    terminal_tail_slope = float(math.exp(c.dividend*c.maturity)*theta[0]+theta[1])
    return {"instrument_cost": cost, "initial_cash": cash,
            "initial_borrowing": max(0., -cash), "far_tail_undiscounted_wealth_slope": terminal_tail_slope,
            "far_tail_unbounded_below": terminal_tail_slope < 0,
            "stock_holding": float(theta[0]), "call_holding": float(theta[1])}
