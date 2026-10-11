"""Meaningful structural/regression checks for the capped-GBM companion."""

from dataclasses import replace
import math
import unittest

import numpy as np

from capped_gbm import (Contract, GREEKS, HeatGrid, analytic, analytic_bumped_greeks,
                        bumped_greeks, domain_price_bound, heat_backward_euler,
                        heat_backward_euler_thomas, monte_carlo,
                        monte_carlo_observations, mixed_gamma_absolute_bound,
                        pde_price, raw_observation_ranges,
                        study_domain, synthetic_contracts)


class CappedGBMTests(unittest.TestCase):
    def test_invalid_parameters_are_rejected(self):
        for name in ("spot", "strike", "cap", "maturity", "volatility"):
            for value in (0.0, -1.0, float("inf"), float("nan")):
                with self.assertRaises(ValueError):
                    replace(Contract(), **{name: value})
        with self.assertRaises(ValueError):
            HeatGrid(0.0, 1.0, 16, 2).validate_for(Contract())

    def test_analytic_greeks_independent_price_differences(self):
        c = Contract(spot=109.0, rate=-0.01, dividend=0.045, volatility=0.31)
        grid = HeatGrid(*study_domain(c), 32, 64)
        bumps = {"spot": 0.01, "volatility": 0.00001, "rate": 0.00001, "maturity": 0.00001}
        fd = analytic_bumped_greeks(c, grid, bumps)
        exact = analytic(c)
        for q in GREEKS:
            self.assertAlmostEqual(fd[q], exact[q], delta=2e-6 * max(1.0, abs(exact[q])))
        # The analytic values satisfy the dividend-paying Black-Scholes PDE.
        residual = (exact["theta"] + 0.5 * c.volatility ** 2 * c.spot ** 2 * exact["gamma"]
                    + (c.rate - c.dividend) * c.spot * exact["delta"] - c.rate * exact["price"])
        self.assertAlmostEqual(residual, 0.0, delta=1e-11)

    def test_callspread_financial_limits(self):
        for spot in (1e-4, 100.0, 1e8):
            c = Contract(spot=spot, rate=-0.02, dividend=0.04)
            a = analytic(c)
            self.assertGreaterEqual(a["price"], -1e-8)
            self.assertLessEqual(a["price"], c.cap * math.exp(-c.rate * c.maturity) + 1e-8)
            self.assertGreaterEqual(a["delta"], -1e-12)
        c = Contract(spot=1e8, rate=-0.02, dividend=0.04)
        self.assertAlmostEqual(analytic(c)["price"], c.cap * math.exp(-c.rate * c.maturity), delta=2e-8)
        # Cap Gamma can be negative; a test must not impose vanilla convexity.
        self.assertLess(analytic(Contract(spot=150.0))["gamma"], 0.0)

    def test_dst_matches_direct_tridiagonal_recurrence(self):
        for c in (Contract(), Contract(rate=-0.03, volatility=0.7, maturity=2.0)):
            for n in (3, 16, 63):
                for nt in (1, 7, 43):
                    grid = HeatGrid(*study_domain(c), n, nt)
                    nodes, fft = heat_backward_euler(c, grid)
                    nodes2, direct = heat_backward_euler_thomas(c, grid)
                    np.testing.assert_array_equal(nodes, nodes2)
                    np.testing.assert_allclose(fft, direct, atol=2e-11, rtol=2e-12)

    def test_unconditional_heat_maximum_principle_large_lambda(self):
        c = Contract(volatility=0.7, maturity=2.0)
        grid = HeatGrid(-0.6, 0.7, 512, 1)
        # lambda far exceeds the explicit-Euler stability restriction.
        lam = 0.5 * c.volatility ** 2 * c.maturity / grid.spacing ** 2
        self.assertGreater(lam, 10000.0)
        _, values = heat_backward_euler(c, grid)
        self.assertGreaterEqual(values.min(), -1e-11)
        self.assertLessEqual(values.max(), c.cap + 1e-11)
        self.assertEqual(values[0], 0.0)
        self.assertEqual(values[-1], c.cap)

    def test_pde_refinement_and_bump_bias_separation(self):
        c = Contract(spot=110.0)
        lo, hi = study_domain(c)
        exact = analytic(c)
        coarse = HeatGrid(lo, hi, 128, 2048)
        fine = HeatGrid(lo, hi, 2048, 524288)
        self.assertLess(abs(pde_price(c, fine) - exact["price"]),
                        abs(pde_price(c, coarse) - exact["price"]) / 20.0)
        estimate = bumped_greeks(c, fine)
        bumps = analytic_bumped_greeks(c, fine)
        tolerances = {"price": 0.002, "delta": 0.001, "gamma": 0.00015,
                      "vega": 0.03, "rho": 0.05, "theta": 0.03}
        for q in GREEKS:
            self.assertLess(abs(estimate[q] - exact[q]), tolerances[q])
            self.assertAlmostEqual(estimate[q] - exact[q],
                                   (estimate[q] - bumps[q]) + (bumps[q] - exact[q]), delta=1e-12)
        self.assertLess(domain_price_bound(c, fine), 1e-15)

    def test_exact_terminal_mc_greeks_and_bounds(self):
        c = Contract(spot=110.0)
        rows = monte_carlo(c, observations=32768, pilot_pairs=4096, seed=1451,
                           alpha=0.05 / (len(GREEKS) + 1))
        raw = [r for r in rows if r["method"].startswith("raw_exact_terminal")]
        for row in raw:
            self.assertLess(abs(row["error"]), 7.0 * row["sample_standard_error"] + 1e-10)
            self.assertLess(abs(row["error"]), row["finite_bound_half_width"])
        controls = [r for r in rows if r["method"].startswith("antithetic_stock_control")]
        for row in controls:
            self.assertLess(abs(row["error"]), 7.0 * row["sample_standard_error"] + 1e-10)
            self.assertEqual(row["production_terminal_evaluations"], 57344)
            self.assertEqual(row["pilot_terminal_evaluations"], 8192)
            self.assertEqual(row["total_terminal_evaluations"], 65536)
        # Directly check the bounded pathwise-observation ranges over a wide
        # deterministic normal grid, including both payoff kinks.
        z = np.linspace(-12.0, 12.0, 20001)
        obs = monte_carlo_observations(c, z)
        for q, (lower, upper) in raw_observation_ranges(c).items():
            idx = GREEKS.index(q)
            self.assertGreaterEqual(float(obs[:, idx].min()), lower - 1e-10)
            self.assertLessEqual(float(obs[:, idx].max()), upper + 1e-10)
        self.assertLessEqual(float(np.abs(obs[:, 6]).max()), mixed_gamma_absolute_bound(c) + 1e-10)
        mixed_rows = [r for r in raw if r["gamma_estimator"] == "mixed_pathwise_likelihood_ratio"]
        self.assertEqual(len(mixed_rows), 1)
        self.assertEqual(mixed_rows[0]["finite_bound_type"], "Hoeffding")

    def test_mixed_gamma_expectation_by_independent_gaussian_quadrature(self):
        # Integrate its smooth active-region integrand directly. This is
        # independent of both simulated observations and call differences.
        nodes, weights = np.polynomial.legendre.leggauss(160)
        for _, c in synthetic_contracts():
            v = c.volatility * math.sqrt(c.maturity)
            za = (math.log(c.strike / c.spot) - c.log_drift * c.maturity) / v
            zb = (math.log((c.strike + c.cap) / c.spot) - c.log_drift * c.maturity) / v
            z = 0.5 * ((zb - za) * nodes + zb + za)
            terminal = c.spot * np.exp(c.log_drift * c.maturity + v * z)
            integrand = (math.exp(-c.rate * c.maturity) * terminal / c.spot ** 2
                         * (z / v - 1.0) * np.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi))
            integral = 0.5 * (zb - za) * float(weights @ integrand)
            self.assertAlmostEqual(integral, analytic(c)["gamma"], delta=2e-11)

    def test_monte_carlo_seed_reproducibility(self):
        first = monte_carlo(Contract(), 128, 64, seed=321)
        second = monte_carlo(Contract(), 128, 64, seed=321)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
