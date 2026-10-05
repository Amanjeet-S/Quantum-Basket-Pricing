# Research design for conditional interval encoding

## 1 Decision and estimand

I ask whether interval probability encoding of a geometric-control residual reduces the complete coherent arithmetic required for pricing a bounded arithmetic basket derivative at a fixed monetary tolerance. I do not make direction selection the contribution. I begin with the established equal-loading Gaussian factor so that the proposed evaluator has explicit interval boundaries and avoids numerical root finding.

The estimand is V = exp(-rT) E[F(A)] under a stipulated risk-neutral correlated GBM, with F(x) = min(max(x-K,0),H), positive normalised weights and K,H > 0. The initial basket value is 100. I use ex-dividend securities with stipulated dividend yields. The numerical target is model price accuracy, rather than observed exotic-option fit or hedging profitability.

## 2 Concrete construction

Let log prices be mu + RY + lambda Z, where Z is an independent standard normal coordinate and lambda is common to every component. Define B_A(Y) as the weighted arithmetic average of exp(mu+RY), and B_G(Y) as exp(weighted mean of mu+RY). Then A=B_A exp(lambda Z), G=B_G exp(lambda Z), and B_A >= B_G.

For U uniform on [0,1], set k=K+HU and z_A=log(k/B_A)/lambda, z_G=log(k/B_G)/lambda. The proposed output is H times the normal probability of the interval [z_A,z_G]. Its expectation is E[F(A)-F(G)]. Adding the analytically calculable E[F(G)] restores the original price. This output is bounded and nonnegative and avoids a difference of four conditional call values. It adds strike randomness, whose variance cost must be measured rather than assumed beneficial.

The mathematical identity is established elementary probability. I seek a resource improvement from its stable finite-precision implementation. Published comparison-register encoding and classical digital preintegration create a serious routine-combination objection. A prototype alone will not remove that objection.

## 3 Error certificate

I distinguish the continuous financial target, the intended finite Gaussian and strike representation, and the finite state actually prepared by a circuit. Continuous and discrete laws are mutually singular, so a small total variation bound cannot certify their approximation. I use a coupling or quadrature bound for that transition. Total variation may instead compare two laws on the same finite representation.

The mathematical note proves a strike-integrated residual certificate: uniform true threshold residual allowances rho_A and rho_G contribute at most rho_A+rho_G in undiscounted mean; N midpoint strikes contribute at most H/N for the controlled output. For the common-loading version no iterative root residual is needed, but logarithm, coefficient, loading and CDF errors remain. Two CDF errors each bounded by eta contribute at most 2H eta. The geometric correction needs its own allowance.

The total discounted bias must fit a prespecified allocation within monetary epsilon. The statistical estimator receives only the remaining tolerance. An analytical conservative residual standard-deviation bound is available from E(A-G) and the implemented output bias, but it applies to competent direct controlled pricing too. I will not present that certificate as a special variance advantage of the new method.

## 4 Initial experiments

I fixed twelve four-asset scenarios before observing their results: correlations 0.1, 0.6 and 0.9; strike and cap pairs (100,1), (100,10), (100,25) and (115,10). Initial asset values are 100, annualised volatilities range from 20 to 30 per cent, maturity is one year, r=0.03 and q=0.02. These are controlled synthetic scenarios, not an empirical sample.

I compare the direct payoff, the direct geometric residual, the fully conditioned payoff and controlled residual, and the random-strike conditioned payoff and controlled residual. Stable positive truncated-payoff integration supplies the full conditional comparator. This numerical reference is checked by quadrature refinement, rather than advertised as a certified reference for every input.

Eight independent batches each contain 4,096 paths. Paired inputs support diagnostic comparisons. Sample standard deviations and batch standard errors do not establish joint confidence coverage. A separate scalar case deliberately exposes variance inflation from strike randomisation when Gaussian loading is weak.

## 5 Strong comparisons required before a contribution claim

The decisive quantum opponents are an independently optimised direct geometric-control residual and stable analytical conditioning with the same control. Applicable QSP encodings and modern variance-sensitive estimators require their own compatible implementations. The candidate must be compared with their complete circuits at the same monetary tolerance and confidence, rather than with an unstable subtraction formula or an uncontrolled payoff alone.

Classical opponents will include geometric controls, conditioning, randomised quasi-Monte Carlo with appropriate factor construction, and low-dimensional quadrature where applicable. The current pilot does not implement or defeat those final classical opponents. Lower sample variance is not a wall-clock or gate-count result.

For a leading variance-sensitive screen only, if the candidate standard deviation is s_new and a direct controlled comparator has s_old, the candidate oracle overhead must remain below s_old/s_new. Actual constants, estimator overhead, deterministic bias allowances, state preparation, inverses and confidence can tighten this screen. A reported threshold is a measured budget for future circuit work, not a demonstrated quantum saving.

## 6 Construction and stop rule

The next implementation will represent the interval probability with certified finite arithmetic and a coherent bounded-output encoding. I will retain component resource counts and distinguish an executed circuit from a compositional projection. T-count is the primary initial objective; T-depth, logical qubits and classical setup time have separate units and will be reported separately.

I will retain a positive contribution claim only if the construction supplies verified resource savings against the strongest compatible, independently tuned opponent on a prespecified evaluation family. A lower envelope of existing methods, a host lookup table, a changed asset universe or a favourable subset of scenarios is insufficient. A failed construction remains recorded as a feasibility finding.

## 7 Asset choice and existing data

A general equity basket keeps the methodological claim clear. REITs can later supply an illustrative scenario using existing LSEG access after auditing prices, corporate actions, dividend treatment and available history. Historical diffusion estimates may parameterise the stipulated GBM model but do not validate market-consistent exotic prices. Direct property valuation, mortgage default modelling and commodity convenience yields are outside this initial project.

## 8 Keywords

Bounded basket derivatives; quantum mean estimation; interval probability; geometric control; finite precision; computational finance.
