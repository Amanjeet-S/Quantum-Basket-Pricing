# Research design for certified scalar conditioning

## 1 Decision and estimand

I initially asked whether interval probability encoding of a geometric-control residual reduces the complete coherent arithmetic required for pricing a bounded arithmetic basket derivative at a fixed monetary tolerance. The [stronger-comparator audit](conditioning_benchmark.md) rejects that advantage in the initial evaluation family. I now investigate how certified scalar-function specialisation changes the comparative reversible cost of direct, fully conditioned and randomised-strike evaluation. I retain the established equal-loading Gaussian factor and the same financial target. I do not infer originality from selecting the lowest existing-method cost.

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

## 8 Constructive stage update

I completed the first reversible evaluator comparison in a separate stipulated two-asset common-loading instance. The [encoding study](reversible_encoding_study.md) records its scope, precision search and exact emitted counts. I use the stronger coordinatewise bounded-variation finite-grid bounds in the [supplement](finite_grid_supplement.md), with the earlier coupling bound retained where smaller. The three constructors receive their own finite widths and precision configurations. Half the requested price tolerance is reserved for unimplemented statistical estimation. I treat the favourable assigned T counts as constructor evidence and require stronger conditioned arithmetic before claiming a complete resource contribution.

## 9 Stronger conditioning and revised decision

I implemented a cap-specific normalised conditional-payoff kernel and allowed equivalent scalar specialisation for the interval method. Full conditioning is cheaper than both interval implementations in every retained cap and tolerance case. I therefore withdraw the proposed interval improvement and retain the failed test. The new [benchmark](conditioning_benchmark.md) derives the financial derivative and interpolation certificate and retains complete lookup, arithmetic, flag, inverse and setup costs. An original positive literature contribution remains unmet. Before building a large mean estimator, I require a checked project-specific cost or precision argument across a prespecified broader financial family, tested against equally specialised opponents and compared with the published literature. I will preserve the result as a feasibility finding if that argument does not survive.

## 10 Keywords

Bounded basket derivatives; quantum mean estimation; scalar conditioning; interval probability; geometric control; finite precision; computational finance.

## 11 General certificate and estimation question

I now investigate whether a mathematically transferable financial moment certificate can reduce complete estimation cost after numerical representation and choice of encoding. Positive-asset conditional moment bounds do not require Gaussianity or a particular compiler. A joint-transform formulation gives a computable bound where the required moments exist, and a distribution-general weak-dispersion theorem identifies a difference between fractional observations and canonical flag observations. The original 27 contracts instantiate these arguments, while separate risk-neutral non-Gaussian and monitored-Asian derivations state wider mathematical scope.

The [new study](general_moment_study.md) retains exact classical prices, complete replay, published quantum confidence envelopes, a separately proved numerical controller, and all adverse comparisons. Established geometric controls, general-model transform methods, limiting approximations and bounded-variance quantum estimation remain mandatory predecessors. My candidate distinction lies in the checked finite-output moment/error certificate and its estimation consequences. A broader model list alone does not establish originality.

The next discriminating test must implement a variance-aware quantum value oracle with a fully charged estimator and compare a stronger classical continuous evaluator. Moment assumptions, numerical bias, confidence scopes, setup and precision must be matched or explicitly distinguished. I will preserve negative evidence if those tests remove the benefit or show that existing literature already supplies the argument.

## 12 Finite-value and deterministic comparison stage

The current research question is: when do financial residual support and moment certificates change the cheapest fully charged quantum encoding at a fixed monetary error and confidence, and where do deterministic classical opponents remove the apparent benefit? I constructed a coherent numerical-value oracle and fully charged dyadic bucket flags, then give every opponent the same actual finite-output range and moment certificates. Single-flag range compression and clipped flags receive a rigorously bounded excess-tail midpoint within the existing estimation budget. The comparison retains both proved scalar controllers, all 819 inherited architecture profiles, every admissible dyadic range and clipping exponent, and equal and cost-dependent dyadic allocations. Selection occurs before measurement outcomes, and the current run gives 99 per cent confidence for each individual selected price rather than a simultaneous 27-price guarantee.

I also execute exact summation for all 47 two-asset full-conditioning prefixes and validated continuous quadrature for all 9 two-asset contracts. These deterministic competitors remove the justification for inferring a low-dimensional pricing advantage from a Monte Carlo sample comparison. Their costs, exact endpoints, version records and independent checks are retained in deterministic_baseline_study.md. The new weak-dispersion range proof shows why a bounded-support single flag can recover the small residual scale; unbounded continuous laws, changing finite support and evaluator precision retain separate charges.

The study in variance_oracle_study.md records favourable, null and adverse complete resource caps. A genuine literature contribution still requires a project-specific distinction beyond established controls, moment-option bounds, dyadic mean estimation and range normalisation. The next discriminating test concerns higher-dimensional deterministic and randomised quasi-Monte Carlo pricing, sharper conditional or bucket-local moment certificates, and compatible modern quantum encodings. I will not infer novelty from increasing model dimension, changing the asset label or combining known ingredients without a substantive checked consequence.

## 13 Continuous certificates before quantum estimation

I completed the higher-dimensional comparator. A finite Taylor expansion of the smoothed capped payoff in the arithmetic–geometric ratio, with Gaussian tilts preserving covariance, supplies an explicit expected remainder. It certifies all eighteen four- and eight-asset contracts and seven two-asset controls at undiscounted radius 0.01. Validated continuous quadrature resolves the two remaining controls. A deterministic model-price answer is therefore available for all 27 frozen contracts before a quantum mean-estimation stage is required.

I also completed independently shifted finite digital-net pricing and a rigorous finite-checkpoint confidence rule. All 50,331,648 financial outputs are replayed. Twenty-four runs meet the statistical reservation and 26 meet the full undiscounted monetary tolerance. Its evaluator work exceeds equally informed fixed-IID plans in all 27 cases. These adverse results are retained rather than described as a randomised pricing improvement.

The current contribution candidate is the explicit coefficient/remainder certificate and its effect on selecting a pricing method, rather than physical quantum advantage on this family. The new source comparison checks published Taylor, Hermite, basket-bound and randomised-confidence predecessors, while recording missing final texts. The full Ju article and complete modern quantum comparisons remain open. A local 82.1 per cent T-cap reduction cannot answer those gaps or establish priority. The revised manuscript separates the mathematical result, numerical evidence and remaining claims, using British English without em dashes.
