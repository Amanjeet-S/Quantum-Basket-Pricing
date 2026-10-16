# Higher-dimensional price certificates and the quantum comparison

## 1 Central result

I derived and implemented a finite-order continuous-price certificate for positive Gaussian baskets with a common Gaussian factor. Every finite coefficient is evaluated through correlated Gaussian exponential tilting; an explicit derivative and moment bound encloses the remainder. This gives a sufficient criterion for meeting a monetary tolerance before quantum mean estimation is required. It is a candidate mathematical contribution whose precise relation to existing Taylor pricing and basket bounds requires further proof-level comparison.

The order-zero-to-eight calculation certifies 25 of 27 contracts at undiscounted radius 0.01, including all 18 four- and eight-asset cases. The two failures are two-asset controls with common loading 0.05 and caps 10 and 25. Their valid but wide Taylor intervals are retained. Existing one-dimensional validated quadrature resolves both, supplying deterministic continuous-model price enclosures for all 27 contracts.

## 2 Derivation and financial scope

I write the arithmetic basket as the geometric basket times a ratio T>=1. Exact weighted centring proves the ordering. For a smoothed conditional capped price P_lambda, I expand g_m(t)=P_lambda(m+log t) at t=1 to a finite order. Its derivatives are falling-factorial differential operators in the log price. Bounds for Gaussian density derivatives supply a global majorant for t>=1. The expected integral remainder is therefore bounded by a computable positive centred-ratio moment.

The coefficients require the joint law of T and the geometric log price. The Gaussian tilt changes the scalar derivative evaluation point by R_G dot v_alpha; replacing the geometric log price by its mean would lose that dependence. Multinomial enumeration and the tilt yield explicit finite sums. All frozen inputs are exact rationals and all subsequent elementary operations are enclosed by Arb at 224 bits. The selected price interval includes arithmetic enclosure and the expected remainder. Reference agreement and apparent series convergence are not substituted for an error proof.

The theorem extends algebraically to unequal initial prices and monitored positive-weight averages whose joint Gaussian covariance has the required common factor. The executable remains the frozen 27-contract family. Positive weights, a positive common factor and finite Gaussian moments are substantive hypotheses. Finite-order validity does not imply convergence of an infinite Taylor series or termination at every requested tolerance.

## 3 Deterministic results and costs

| Assets | Loading | Cap | First certifying order | Narrowest retained radius |
|---:|---:|---:|---:|---:|
| 4 | 0.05 | 1 | 2 | 0.000070631608716753224735555924 |
| 4 | 0.05 | 10 | 4 | 0.000706316087167532247355559238 |
| 4 | 0.05 | 25 | 5 | 0.001765790217918830618388898095 |
| 4 | 0.2 | 1 | 1 | 1.480811710561530849752E-9 |
| 4 | 0.2 | 10 | 2 | 1.4808117105615308497513E-8 |
| 4 | 0.2 | 25 | 2 | 3.7020292764038271243781E-8 |
| 4 | 0.4 | 1 | 1 | 1.9346363194968425722E-11 |
| 4 | 0.4 | 10 | 1 | 1.93463631949684257219E-10 |
| 4 | 0.4 | 25 | 1 | 4.83659079874210643046E-10 |
| 8 | 0.05 | 1 | 2 | 0.000002281505780966665559231480 |
| 8 | 0.05 | 10 | 3 | 0.000022815057809666655592314798 |
| 8 | 0.05 | 25 | 4 | 0.000057037644524166638980786994 |
| 8 | 0.2 | 1 | 1 | 4.7832415819914946208E-11 |
| 8 | 0.2 | 10 | 1 | 4.78324158199149462072E-10 |
| 8 | 0.2 | 25 | 2 | 1.195810395497873655180E-9 |
| 8 | 0.4 | 1 | 1 | 6.24916241777908578E-13 |
| 8 | 0.4 | 10 | 1 | 6.249162417779085777E-12 |
| 8 | 0.4 | 25 | 1 | 1.5622906044447714441E-11 |

The full run evaluates every retained order, using 25,080 geometric exponentials and 120,690 scalar kernels. Its diagnostic host time is 15.365 seconds, including 6.675 seconds of geometric setup and 8.683 seconds of kernel work. These are shared-family timings and exclude process startup. The lower term count of a first certifying order is not presented as a measured standalone runtime. Quantum logical T counts, classical operation counts and CPU seconds remain distinct.

Independent verification covers 540 formal derivative comparisons, 60 strike-integral derivatives, 81 direct coefficient quadratures, 81 reference/order interval checks, all 25,080 polynomial weights and 83 covariance-exposure identities. Six targeted unit tests pass.

## 4 Randomised pricing and the stronger IID opponent

I separately execute 27 finite-oracle prices with 64 points per independently shifted digital net and eleven predeclared stopping checkpoints. The sample unit is the whole shift batch. The confidence proof uses disjoint-pair variance observations, a lower-tail Chernoff inversion and a two-sided Bernstein mean bound, with a union allocation across checkpoints. It does not assume independence within a net, a variance reduction rate in net size or Student coverage. The helper rejects non-integral counters and impossible aggregate moments. Its correctly rounded logarithms and integer square-root enclosures are charged explicitly.

The run evaluates 50,331,648 exact financial outputs under 786,432 shifts. Every output and batch sum is independently replayed, with 3,544 additional scalar graph checks. All 227 checkpoints are replayed through a separate 300-digit confidence calculation. Exact coverage checks for 66 Bernoulli laws retain 89,298 probability terms. Twenty-one shift/uncertainty unit tests pass. The ideal independent-uniform-shift assumption is stated separately from reproducible PCG execution.

Twenty-four prices meet the reserved statistical allowance. The three loading 0.05, cap 25 cases reach the fixed 65,536-shift maximum. Two still meet the full monetary tolerance because their deterministic error allowances have slack. Both complete undiscounted and discounted tolerances are attained in 26 of 27 cases; the four-asset weak-loading cap 25 case remains unresolved by this randomised procedure, at undiscounted radius 0.0111767. Its deterministic Taylor certificate succeeds. All 27 deterministic discounted price intervals are contained in their wider randomised intervals.

The fixed-IID opponent receives the same finite support and moment bounds, with 137 conditioned-prefix candidates and 99 per cent per-price confidence. The randomised evaluator work is 2.544–133.203 times its selected fixed-IID planned work across all 27 cases. These IID budgets are plans rather than new streams. All net generation, moment reduction and confidence overhead is retained separately. This adverse finding belongs to the implemented conservative randomised procedure; it is not a rejection of every published RQMC policy.

The final randomised run uses 43.685 seconds of candidate/table setup and 5.130 seconds of production execution. These diagnostic timings include the declared shared cache and are not directly converted into quantum hardware performance. The interface-hardening rerun reproduces all batch sums, stream hashes, checkpoint intervals and selections exactly.

## 5 Research claim and publication gaps

The present quantum advantage claim is unsupported on this benchmark: a deterministic answer at the required model-price accuracy is available for all 27 contracts. The earlier 17-contract T-cap reductions remain valid encoding comparisons. They establish neither physical advantage nor superiority over an algorithm that resolves the price without mean-estimation calls. Fixed-point and logical-gate accounting are numerical implementation costs, rather than measurements of a physical quantum machine.

The useful candidate is the finite coefficient-and-remainder certificate and its role in deciding whether coherent estimation is needed. The theorem and executable are a substantive financial mathematics result within their assumptions, but priority in the literature is unresolved. Published geometric conditioning, Taylor pricing, Hermite approximations, basket bounds and RQMC confidence methods are attributed in the source comparison. Ju's complete paper remains unread, and compatible complete QSP and stronger quantum mean-estimation implementations remain absent. Neither inaccessible text nor an unimplemented opponent is evidence of novelty or advantage.

I centralised these qualifications in the manuscript's research question and comparison discussion, while preserving technical restrictions at the theorems that need them. The mathematical argument, original project derivation, numerical findings and interpretation are distinguished. No new empirical data, purchase or entitlement was used. The current result is model-price certification and method comparison rather than empirical asset-price prediction.

## 6 Reproduction and preserved evidence

The Taylor, shifted-net, confidence and combined comparison sources are taylor_pricing.py, rqmc_pricing.py, rqmc_confidence.py, rqmc_work_comparison.py and higher_dimension_comparison.py. Their retained results and independent audits are under results and checks. Existing authorised NumPy and python-flint 0.9.0 runtimes are sufficient. The original 78 result files and all numerical inputs remain unchanged. The open generalisation_proofs.tex remains the sole manuscript source, with keywords below its abstract and the same author–date reference style.

## 7 Keywords

Financial mathematics; quantum pricing; certified approximation; Gaussian tilting; geometric controls; finite moments; randomised quadrature; monetary accuracy.
