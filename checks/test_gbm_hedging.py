"""Structural checks for physical-law constrained static hedging."""

import math
import unittest

import numpy as np

from capped_gbm import Contract, analytic, call_analytic, payoff, synthetic_contracts
from gbm_hedging import (constraint_system, hedge_financing, hedge_moments,
                         quadratic_loss, quadrature_moments, residual_square_integral,
                         simulated_gains, slsqp_hedge, solve_constrained_hedge)


class GBMHedgingTests(unittest.TestCase):
    def test_closed_moments_match_independent_kink_quadrature(self):
        for c in (Contract(),Contract(spot=160.0),Contract(volatility=.65,maturity=2.,cap=30.)):
            for mu in (.03,.12):
                m=hedge_moments(c,mu)
                numerical,error=quadrature_moments(c,mu)
                exact=np.array([*m.gain_mean,m.matrix[0,0],m.matrix[0,1],m.matrix[1,1],*m.cross,m.constant])
                np.testing.assert_allclose(exact,numerical,rtol=2e-9,atol=5e-9)
                self.assertLess(error,1e-4)

    def test_active_set_kkt_and_dense_feasible_comparator(self):
        matrix=np.array([[2.,.5],[.5,1.]])
        cross=np.array([1.5,.8])
        costs=np.array([100.,10.])
        for budget in (0.,2.,6.,30.,150.):
            solution=solve_constrained_hedge(matrix,cross,3.,costs,budget)
            self.assertLessEqual(solution.candidates_considered,16)
            self.assertLess(solution.primal_violation,1e-9)
            self.assertLess(solution.dual_violation,1e-9)
            self.assertLess(solution.stationarity_residual,1e-8)
            self.assertLess(solution.complementarity_residual,1e-8)
            values=[]
            for stock in np.linspace(0,1,101):
                remaining=budget-costs[0]*stock
                if remaining>=0:
                    for call in np.linspace(0,min(1.,remaining/costs[1]),101):
                        values.append(quadratic_loss(np.array([stock,call]),matrix,cross,3.))
            self.assertLessEqual(solution.loss_quadratic,min(values)+1e-9)
            if budget==0:
                np.testing.assert_array_equal(solution.theta,np.zeros(2))

    def test_budget_value_convexity_and_multiplier_subgradient(self):
        matrix=np.array([[2.,.5],[.5,1.]])
        cross=np.array([1.5,.8])
        costs=np.array([100.,10.])
        budgets=np.array([0.,1.,2.,4.,6.,8.,20.,40.,68.,100.])
        solutions=[solve_constrained_hedge(matrix,cross,3.,costs,float(b)) for b in budgets]
        values=np.array([s.loss_quadratic for s in solutions])
        self.assertTrue(np.all(np.diff(values)<=1e-10))
        secants=np.diff(values)/np.diff(budgets)
        self.assertTrue(np.all(np.diff(secants)>=-1e-10))
        for budget,solution in zip(budgets,solutions):
            supporting_values=solution.loss_quadratic-solution.multipliers[4]*(budgets-budget)
            self.assertTrue(np.all(values>=supporting_values-1e-10))
        # These budgets lie strictly within a fixed active-set regime, so the
        # envelope derivative equals minus the budget multiplier.
        for budget in (4.,20.,90.):
            centre=solve_constrained_hedge(matrix,cross,3.,costs,budget)
            step=1e-4
            low=solve_constrained_hedge(matrix,cross,3.,costs,budget-step)
            high=solve_constrained_hedge(matrix,cross,3.,costs,budget+step)
            derivative=(high.loss_quadratic-low.loss_quadratic)/(2*step)
            self.assertAlmostEqual(derivative,-centre.multipliers[4],delta=1e-8)

    def test_sherman_morrison_mse_and_variance_gap(self):
        for c in (Contract(),Contract(rate=-.015,dividend=.045,volatility=.25,spot=105.)):
            m=hedge_moments(c,.12)
            v=np.linalg.solve(m.gain_covariance,m.centred_cross)
            direct=np.linalg.solve(m.matrix,m.cross)
            inv_mean=np.linalg.solve(m.gain_covariance,m.gain_mean)
            zeta=m.gain_mean@inv_mean
            eta=m.residual_mean-m.gain_mean@v
            closed=v+inv_mean*eta/(1+zeta)
            np.testing.assert_allclose(closed,direct,rtol=1e-11,atol=1e-11)
            observed=quadratic_loss(v,m.matrix,m.cross,m.constant)-quadratic_loss(direct,m.matrix,m.cross,m.constant)
            self.assertAlmostEqual(observed,eta**2*zeta/(1+zeta),delta=1e-10)
            self.assertGreater(observed,0.)

    def test_risk_neutral_special_case_uses_each_contract_rate(self):
        for _,c in synthetic_contracts():
            m=hedge_moments(c,c.rate)
            np.testing.assert_allclose(m.gain_mean,0,atol=5e-12)
            self.assertAlmostEqual(m.residual_mean,0,delta=5e-12)
            np.testing.assert_allclose(m.matrix,m.gain_covariance,atol=5e-12)
            np.testing.assert_allclose(m.cross,m.centred_cross,atol=5e-12)

    def test_constrained_variance_is_not_general_mse_minimiser(self):
        c=Contract()
        m=hedge_moments(c,.12)
        mse=solve_constrained_hedge(m.matrix,m.cross,m.constant,m.costs,m.premium)
        variance=solve_constrained_hedge(m.gain_covariance,m.centred_cross,m.residual_variance,m.costs,m.premium)
        mse_risk,_=residual_square_integral(c,.12,mse.theta)
        variance_risk,_=residual_square_integral(c,.12,variance.theta)
        self.assertLessEqual(mse_risk,variance_risk+1e-10)
        self.assertGreater(np.max(np.abs(mse.theta-variance.theta)),1e-5)

    def test_nonnegative_residual_integral_retains_cancellation_difference(self):
        c=Contract(maturity=.05,volatility=.06)
        m=hedge_moments(c,.06)
        theta=np.array([0.,1.])
        positive,error=residual_square_integral(c,.06,theta)
        raw=quadratic_loss(theta,m.matrix,m.cross,m.constant)
        self.assertGreater(positive,0.)
        self.assertLess(positive,1e-30)
        self.assertGreaterEqual(error,0.)
        self.assertGreater(abs(raw-positive),positive)

    def test_heldout_simulation_and_dividend_financing(self):
        c=Contract(rate=-.015,dividend=.045,volatility=.25)
        m=hedge_moments(c,.06)
        gains,residual=simulated_gains(c,.06,131072,84193)
        theta=np.array([.05,.2])
        losses=(residual-gains@theta)**2
        population,_=residual_square_integral(c,.06,theta)
        se=losses.std(ddof=1)/math.sqrt(len(losses))
        self.assertLess(abs(losses.mean()-population),6*se)
        financing=hedge_financing(c,m,theta)
        self.assertAlmostEqual(financing['initial_cash']+financing['instrument_cost'],m.premium)
        self.assertFalse(financing['far_tail_unbounded_below'])
        self.assertTrue(hedge_financing(c,m,np.array([-1.,0.]))['far_tail_unbounded_below'])

    def test_external_call_spread_is_different_instrument_set(self):
        c=Contract()
        x=np.array([80.,100.,110.,120.,180.])
        np.testing.assert_allclose(np.maximum(x-c.strike,0)-np.maximum(x-c.strike-c.cap,0),payoff(x,c))
        self.assertAlmostEqual(call_analytic(c,c.strike)['price']-call_analytic(c,c.strike+c.cap)['price'],analytic(c)['price'])
        self.assertEqual(constraint_system(np.array([100.,10.]),20.)[0].shape,(5,2))

    def test_invalid_and_indefinite_inputs_are_not_repaired(self):
        for matrix in (np.array([[1.,0.],[0.,0.]]),np.array([[1.,2.],[2.,1.]])):
            with self.assertRaisesRegex(ValueError,'not positive definite'):
                solve_constrained_hedge(matrix,np.ones(2),1.,np.array([100.,10.]),20.)
        for budget in (-1.,math.inf):
            with self.assertRaises(ValueError):
                constraint_system(np.array([100.,10.]),budget)
        with self.assertRaises(ValueError):
            simulated_gains(Contract(),.06,1,1)


if __name__=='__main__':
    unittest.main()
