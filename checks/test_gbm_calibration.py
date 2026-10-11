"""Structural and independently computed checks for Gaussian GBM calibration."""

from dataclasses import replace
import math
import unittest

import numpy as np
from scipy.stats import chi2, t

from capped_gbm import Contract, analytic
from gbm_calibration import (binomial_coverage_interval, capped_price_interval,
                             capped_volatility_peak, estimate_gaussian_returns,
                             gaussian_parameter_intervals, negative_log_likelihood,
                             trace_target_shrinkage)
from run_calibration_study import simulate_case


class GaussianCalibrationTests(unittest.TestCase):
    def setUp(self):
        self.delta = 1.0 / 252.0
        self.covariance = np.array([[0.04, 0.01, -0.006], [0.01, 0.03, 0.005], [-0.006, 0.005, 0.06]])
        self.drift = np.array([0.03, -0.02, 0.06])
        rng = np.random.default_rng(727161)
        self.returns = self.delta*self.drift + math.sqrt(self.delta)*rng.standard_normal((128, 3)) @ np.linalg.cholesky(self.covariance).T

    def test_likelihood_mle_stationarity_and_perturbations(self):
        estimate = estimate_gaussian_returns(self.returns, self.delta)
        self.assertEqual(estimate.covariance_ml_status, "full_rank_covariance_mle")
        inv = np.linalg.inv(estimate.covariance_ml_candidate)
        residual = self.returns - self.delta*estimate.log_drift
        mean_gradient = -inv @ residual.sum(axis=0)
        covariance_gradient = (estimate.observations*inv/2
                               - inv @ (residual.T @ residual) @ inv/(2*self.delta))
        np.testing.assert_allclose(mean_gradient, 0.0, atol=1e-11)
        np.testing.assert_allclose(covariance_gradient, 0.0, atol=1e-10)
        optimum = negative_log_likelihood(self.returns, self.delta, estimate.log_drift,
                                          estimate.covariance_ml_candidate)
        for sign in (-1.0, 1.0):
            self.assertGreater(negative_log_likelihood(self.returns, self.delta,
                               estimate.log_drift+sign*np.array([0.1, -0.05, 0.02]),
                               estimate.covariance_ml_candidate), optimum)
            self.assertGreater(negative_log_likelihood(self.returns, self.delta,
                               estimate.log_drift, estimate.covariance_ml_candidate*(1+sign*.05)), optimum)

    def test_ml_and_unbiased_covariance_are_distinct(self):
        estimate = estimate_gaussian_returns(self.returns, self.delta)
        np.testing.assert_allclose(estimate.covariance_unbiased,
                                   estimate.covariance_ml_candidate*128/127, rtol=1e-14)
        np.testing.assert_allclose(estimate.log_drift, self.returns.mean(axis=0)/self.delta)

    def test_singular_ml_is_flagged_and_likelihood_rejects_it(self):
        data = np.array([[.01, .02, .03], [.02, -.01, .01], [-.01, .03, .02]])
        estimate = estimate_gaussian_returns(data, self.delta)
        self.assertLess(estimate.scatter_rank, 3)
        self.assertEqual(estimate.covariance_ml_status, "singular_scatter_no_positive_definite_covariance_mle")
        # A deliberately exact singular covariance tests input rejection; the
        # candidate's near-zero eigenvalue may have either sign in float64.
        with self.assertRaises(ValueError):
            negative_log_likelihood(data, self.delta, estimate.log_drift, np.diag([.01, .02, 0.0]))
        shrunk = trace_target_shrinkage(estimate, .1)
        target = np.trace(estimate.covariance_ml_candidate)/3
        self.assertGreaterEqual(np.linalg.eigvalsh(shrunk).min(), .1*target-1e-14)
        self.assertFalse(np.allclose(shrunk, estimate.covariance_ml_candidate))
        # Marginal intervals are mathematically meaningful despite joint rank
        # deficiency; this never relabels the matrix a positive-definite MLE.
        self.assertTrue(np.all(gaussian_parameter_intervals(estimate).variance_lower > 0))
        smaller = estimate_gaussian_returns(data[:2], self.delta)
        self.assertEqual(smaller.observations, 2)
        self.assertLessEqual(smaller.scatter_rank, smaller.observations-1)
        self.assertEqual(smaller.covariance_ml_status, "singular_scatter_no_positive_definite_covariance_mle")

    def test_exact_linear_dependence_and_zero_trace_handling(self):
        first = np.arange(10.0)
        estimate = estimate_gaussian_returns(np.column_stack((first, 2*first, -first)), 1.0)
        self.assertEqual(estimate.scatter_rank, 1)
        self.assertGreater(np.linalg.eigvalsh(trace_target_shrinkage(estimate)).min(), 0)
        constant = estimate_gaussian_returns(np.ones((5, 3)), 1.0)
        with self.assertRaises(ValueError):
            trace_target_shrinkage(constant)
        with self.assertRaises(ValueError):
            gaussian_parameter_intervals(constant)

    def test_interval_quantiles_and_bonferroni_split(self):
        estimate = estimate_gaussian_returns(self.returns, self.delta)
        ordinary = gaussian_parameter_intervals(estimate)
        family = gaussian_parameter_intervals(estimate, familywise=True)
        self.assertAlmostEqual(family.marginal_alpha, .05/6)
        self.assertTrue(np.all(family.variance_lower < ordinary.variance_lower))
        self.assertTrue(np.all(family.variance_upper > ordinary.variance_upper))
        nu = 127
        variance = np.diag(estimate.covariance_unbiased)
        np.testing.assert_allclose(chi2.cdf(nu*variance/ordinary.variance_lower, nu), .975, atol=1e-13)
        np.testing.assert_allclose(chi2.cdf(nu*variance/ordinary.variance_upper, nu), .025, atol=1e-13)
        radius = (ordinary.log_drift_upper-estimate.log_drift)/np.sqrt(variance/(128*self.delta))
        np.testing.assert_allclose(t.cdf(radius, nu), .975, atol=1e-13)

    def test_normal_pivot_moments_from_independent_simulation(self):
        rng = np.random.default_rng(801172)
        repetitions, n = 12000, 63
        z = rng.standard_normal((repetitions, n))
        mean = z.mean(axis=1)
        scatter = ((z-mean[:, None])**2).sum(axis=1)
        student = mean / np.sqrt(scatter/((n-1)*n))
        self.assertLess(abs(scatter.mean()-(n-1)), .5)
        self.assertLess(abs(scatter.var(ddof=1)-2*(n-1)), 7.0)
        self.assertLess(abs(student.mean()), .04)
        self.assertLess(abs(student.var(ddof=1)-(n-1)/(n-3)), .07)
        self.assertLess(abs(float(np.corrcoef(mean, scatter)[0, 1])), .04)

    def test_volatility_peak_and_dense_image_agree(self):
        c = Contract()
        peak = capped_volatility_peak(c)
        self.assertIsNotNone(peak)
        self.assertAlmostEqual(peak**2, math.log(1.2)-.04)
        self.assertAlmostEqual(analytic(replace(c, volatility=peak))["vega"], 0.0, delta=1e-12)
        self.assertGreater(analytic(replace(c, volatility=peak*.99))["vega"], 0.0)
        self.assertLess(analytic(replace(c, volatility=peak*1.01))["vega"], 0.0)
        for low, high in ((.1, .8), (.1, .2), (.6, .8)):
            image = capped_price_interval(c, low, high)
            dense = [analytic(replace(c, volatility=float(s)))["price"] for s in np.linspace(low, high, 3001)]
            self.assertAlmostEqual(image["lower"], min(dense), delta=1e-12)
            self.assertGreaterEqual(image["upper"]+1e-12, max(dense))
            self.assertLess(image["upper"]-max(dense), 2e-6)

    def test_endpoint_price_interval_fails_at_global_peak(self):
        c = Contract()
        peak = capped_volatility_peak(c)
        true_price = analytic(replace(c, volatility=peak))["price"]
        for low, high in ((.9*peak, 1.1*peak), (.5*peak, .8*peak), (1.1*peak, 1.5*peak)):
            image = capped_price_interval(c, low, high)
            self.assertLess(image["endpoint_only_upper"], true_price)
            if low < peak < high:
                self.assertAlmostEqual(image["upper"], true_price, delta=1e-12)
            else:
                self.assertLess(image["upper"], true_price)

    def test_nonpositive_critical_square_is_monotone(self):
        c = Contract(spot=160.0)
        self.assertIsNone(capped_volatility_peak(c))
        image = capped_price_interval(c, .1, .8)
        self.assertFalse(image["interior_peak_included"])
        self.assertEqual(image["upper"], image["endpoint_only_upper"])
        self.assertGreater(analytic(replace(c, volatility=.1))["price"],
                           analytic(replace(c, volatility=.8))["price"])

    def test_physical_mean_is_not_pricing_drift(self):
        estimate = estimate_gaussian_returns(self.returns, self.delta)
        altered = estimate_gaussian_returns(self.returns+self.delta*np.array([.4, -.2, .1]), self.delta)
        np.testing.assert_allclose(altered.log_drift-estimate.log_drift, [.4, -.2, .1], atol=1e-13)
        np.testing.assert_allclose(altered.covariance_unbiased, estimate.covariance_unbiased, rtol=1e-13)
        c1 = replace(Contract(), volatility=math.sqrt(estimate.covariance_unbiased[0, 0]))
        c2 = replace(Contract(), volatility=math.sqrt(altered.covariance_unbiased[0, 0]))
        self.assertAlmostEqual(analytic(c1)["price"], analytic(c2)["price"], delta=1e-12)

    def test_binomial_boundaries_and_midpoint(self):
        self.assertEqual(binomial_coverage_interval(0, 100)[0], 0.0)
        self.assertEqual(binomial_coverage_interval(100, 100)[1], 1.0)
        lo, hi = binomial_coverage_interval(95, 100)
        self.assertLess(lo, .95)
        self.assertGreater(hi, .95)

    def test_invalid_inputs_are_rejected(self):
        for data, delta in ((np.ones((1, 3)), 1.0), (np.ones(10), 1.0),
                            (np.ones((10, 3)), 0.0), (np.array([[0.], [float("nan")]]), 1.0)):
            with self.assertRaises(ValueError):
                estimate_gaussian_returns(data, delta)
        for lo, hi in ((0.0, 1.0), (.5, .1), (.1, math.inf)):
            with self.assertRaises(ValueError):
                capped_price_interval(Contract(), lo, hi)
        estimate = estimate_gaussian_returns(self.returns, self.delta)
        for strength in (0.0, -.1, 1.1):
            with self.assertRaises(ValueError):
                trace_target_shrinkage(estimate, strength)
        for n in (1, 2, 3, 3.5):
            with self.assertRaisesRegex(ValueError, "four observations"):
                simulate_case(Contract(), n, self.delta, 10, 727161)


if __name__ == "__main__":
    unittest.main()
