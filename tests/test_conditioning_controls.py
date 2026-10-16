"""Meaningful geometry, independence, degeneracy and cost tests."""

from __future__ import annotations

import math
import unittest

import numpy as np

from capped_gbm import Contract, analytic, normal_pdf, synthetic_contracts
from conditioning_controls import (PilotCoefficients, conditional_population_moments,
                                    conditional_value, cost_counts, estimator_values,
                                    fit_controls, fit_pilot, gaussian_rule,
                                    raw_population_moments, residual_population_second_moment,
                                    split_observations, truncated_stock_moment)


class ConditioningTests(unittest.TestCase):
    def test_truncated_moments_against_independent_interval_quadrature(self):
        c = Contract()
        nodes, weights = np.polynomial.legendre.leggauss(400)
        mu = math.log(c.spot) + c.log_drift * c.maturity
        sd = c.volatility * math.sqrt(c.maturity)
        left = (math.log(c.strike) - mu) / sd
        right = (math.log(c.strike + c.cap) - mu) / sd
        z = 0.5 * (right - left) * nodes + 0.5 * (right + left)
        density = np.fromiter((normal_pdf(v) for v in z), dtype=float)
        for p in (0, 1, 2):
            numerical = 0.5 * (right - left) * float(weights @ (np.exp(p * (mu + sd * z)) * density))
            self.assertAlmostEqual(numerical, truncated_stock_moment(c, p, c.strike,
                                                                      c.strike + c.cap), delta=1e-8)

    def test_conditional_mean_and_projection_identity(self):
        for name, c in synthetic_contracts():
            for u in (0.25, 0.5, 0.9):
                raw = raw_population_moments(c)
                m = conditional_population_moments(c, u, 320)
                # The smooth rule resolves these registered parameters well.
                self.assertAlmostEqual(m["quadrature_mean"], analytic(c)["price"], delta=1e-7)
                self.assertAlmostEqual(m["quadrature_control_mean"], 0.0, delta=1e-9)
                residual = residual_population_second_moment(c, u, raw["oracle_beta"], 320)
                self.assertAlmostEqual(residual, m["same_beta_projection_removed_variance"], delta=2e-6)
                self.assertGreaterEqual(m["same_beta_projection_removed_variance"], -2e-6)
                self.assertLessEqual(m["oracle_control_variance"], raw["oracle_control_variance"] + 2e-6)
                excess = m["raw_oracle_coefficient_conditional_variance"] - m["oracle_control_variance"]
                self.assertAlmostEqual(excess, m["coefficient_reoptimisation_reduction"], delta=1e-10)

    def test_nested_time_variances_and_degenerate_endpoints(self):
        c = Contract()
        nodes, weights = gaussian_rule(160)
        obs0 = split_observations(c, 0.0, nodes, -nodes)
        obs1 = split_observations(c, 1.0, nodes, -nodes)
        self.assertTrue(np.all(obs0["conditional"] == analytic(c)["price"]))
        self.assertTrue(np.all(obs0["conditional_stock_control"] == 0.0))
        np.testing.assert_allclose(obs1["raw"], obs1["conditional"], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(obs1["stock_control"], obs1["conditional_stock_control"], atol=1e-13)
        fractions = (0.0, 0.25, 0.5, 0.9, 1.0)
        moments = [conditional_population_moments(c, u) for u in fractions]
        for field in ("variance", "oracle_control_variance", "raw_oracle_coefficient_conditional_variance"):
            self.assertTrue(all(a[field] <= b[field] + 1e-10 for a, b in zip(moments, moments[1:])))
        # Marginal split normals across u are not asserted to form one Brownian path.

    def test_singular_controls_moore_penrose_and_zero_control(self):
        x = np.linspace(-1.0, 1.0, 501)
        y = 2.0 * x + 0.01 * x ** 2
        single = fit_controls(y, x)
        duplicate = fit_controls(y, np.column_stack((x, x, np.zeros_like(x))))
        self.assertAlmostEqual(single[0], 2.0, places=12)
        np.testing.assert_allclose(duplicate, (1.0, 1.0, 0.0), atol=1e-12)
        np.testing.assert_allclose(fit_controls(y, np.zeros_like(x)), (0.0,), atol=0.0)
        np.testing.assert_allclose(fit_controls(1e10 * x, x, bound=1e6), (1e6,), atol=0.0)

    def test_independent_pilot_fixed_coefficient_mean(self):
        c = Contract()
        fitted = fit_pilot(c, 0.5, 8192, 840731)
        # Conditional on the pilot, any finite coefficient multiplies an exact
        # zero-mean control. A 1D deterministic rule checks that centring.
        nodes, weights = gaussian_rule(320)
        t = 0.5 * c.maturity
        stock = c.spot * np.exp(c.log_drift * t + c.volatility * math.sqrt(t) * nodes)
        m = conditional_value(c, 0.5, stock)
        control = math.exp(-c.rate * t - c.dividend * t) * stock - c.spot * math.exp(-c.dividend * c.maturity)
        for beta in (fitted.raw, fitted.conditional, -2.0, 3.0):
            self.assertAlmostEqual(float(weights @ (m - beta * control)), analytic(c)["price"], places=10)

    def test_unconditioned_control_adds_residual_stock_noise(self):
        c = Contract()
        moments = conditional_population_moments(c, 0.25)
        self.assertGreater(moments["unconditioned_control_added_variance"], 0.0)
        self.assertAlmostEqual(moments["unconditioned_control_variance"]
                               - moments["raw_oracle_coefficient_conditional_variance"],
                               moments["unconditioned_control_added_variance"], places=12)
        # This comparator remains unbiased, but retains unnecessary stock noise.

    def test_cost_counts_include_pilot_and_conditional_cdfs(self):
        raw = cost_counts("raw_pilot_stock_control", 65536, 8192)
        cond = cost_counts("conditional_own_pilot_coefficient", 65536, 8192)
        naive = cost_counts("conditional_unconditioned_stock_control", 65536, 8192)
        self.assertEqual(raw["normal_variates"], 73728)
        self.assertEqual(raw["terminal_payoff_evaluations"], 73728)
        self.assertEqual(cond["normal_cdf_evaluations"], 294912)
        self.assertEqual(cond["terminal_payoff_evaluations"], 0)
        self.assertEqual(naive["normal_variates"], 139264)

    def test_input_validation(self):
        c = Contract()
        for u in (-0.01, 1.01, math.nan):
            with self.assertRaises(ValueError):
                conditional_value(c, u, np.array([100.0]))
        with self.assertRaises(ValueError):
            fit_controls(np.array([1.0]), np.array([2.0]))
        with self.assertRaises(ValueError):
            conditional_value(c, 0.5, np.array([0.0]))

    def test_stable_extreme_stock_limits_and_range_record(self):
        c = Contract()
        diagnostics = {}
        values = conditional_value(c, 0.5, np.array([1e-300, 100.0, 1e300]), diagnostics)
        self.assertEqual(values[0], 0.0)
        self.assertAlmostEqual(values[-1], c.cap * math.exp(-c.rate * c.maturity), places=12)
        self.assertEqual(diagnostics["conditional_state_evaluations"], 3)
        self.assertEqual(diagnostics["lower_range_violations"], 0)
        self.assertEqual(diagnostics["upper_range_violations"], 0)


if __name__ == "__main__":
    unittest.main()
