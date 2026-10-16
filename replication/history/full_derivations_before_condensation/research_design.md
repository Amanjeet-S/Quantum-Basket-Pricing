# Research design for certified quantum basket pricing

## 1 Decision and estimand

I develop the quantum pricing method while retaining an established classical financial model. An improvement of an existing method is permitted. The undiscounted target is V=E_Q[F(A)] and the derivative price is Π=exp(-rT)V. The bounded capped payoff, positive normalised weights, dividend treatment and risk-neutral GBM are unchanged. The paper derives an equivalent martingale measure and absence of arbitrage for the stated admissible strategies. The model price, rather than observed exotic-option fit or return forecasting, supplies the accuracy target.

The revised question concerns whether a certified integral-remainder oracle reduces a complete quantum pricing budget relative to equally informed direct correction evaluation. It is not a request for quantum-mechanical asset dynamics. I distinguish the new oracle proposal and realised auxiliary law from the earlier completed financial circuits and estimators.

## 2 Selected quantum construction and completed next step

The exact order-r Taylor correction is the expectation of U^r t^(-r) A_r(ξ) φ(ξ), multiplied by H/r!, using a Beta(1,r) variable and a uniform strike. The joint Gaussian derivative norm supplies a certified finite-support normalisation. I prove signed comparison encoding and an error envelope covering tails, financial and auxiliary grids, coefficient/value arithmetic, state preparation where approximate, and statistical estimation. The first-order correction is nonnegative and has a simpler comparison flag. Auxiliary randomisation adds nonnegative variance against exact conditional correction evaluation.

I emitted and independently verified the minimum circuit preparing the finite Beta law. The complete financial value calculation is the next implementation task. The rejected signed Gaussian loader, compatible QSP specification and matched classical certificate remain opponents. The [current direction record](quantum_method_redesign.md) states the full experiment and decision rule. The remaining sections preserve earlier phases and their historical identities, assumptions and results; their error budgets are not transferred automatically to the new correction.

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

## 14 Matched estimator sensitivity and executed certificate selection

I add canonical phase-estimation amplitude estimation as a concrete stronger opponent, with exact conditional Grover signs, certified Clifford and T inverse Fourier tapes and exact rational majority confidence. Repricing the 27 old selected financial flags under identical priors, tolerances and confidence lowers every retained T upper budget. The full 819-profile encoding search is separate: range or clipping improves 18 unit-range QPE budgets, while the implemented dyadic family loses every comparison. These findings identify an estimator effect without asserting algorithm-wide superiority or physical quantum advantage.

I also execute three finite certificate policies rather than borrowing a hypothetical first-order count from the old complete search. The monotonicity-first hierarchy preserves the same 25 Taylor-resolved contracts using 1,169 Gaussian kernels and 1,552 geometric exponentials. The full search uses 120,690 and 25,080. Failed screens, scalar arithmetic width and separate cold setup are included. I prove that the particular raw absolute-derivative remainder eventually diverges whenever a centred Gaussian loading is nonzero, which provides a structural reason for finite selection. That result concerns the bound, not the actual expansion error.

The two figures answer specific comparison questions. Degree trajectories show the finite-order certificate's success and failure; matched gate caps and phase widths isolate the amplitude-estimator change. A compact table retains cap-dependent integer degree choices. Exact data and coordinates accompany the visual material. The candidate literature contribution remains the evaluable certificate, its derived finite-order limitation and its checked pricing-method consequences. Priority still requires the missing full Ju comparison, and a complete modern quantum comparison remains necessary before a stronger superiority claim. New changes remain uncommitted following completion of the earlier authorised snapshot commit.
