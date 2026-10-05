# Constructed reversible encodings of a controlled basket call spread

## 1. Initial comparison and superseded interpretation

I initially constructed complete finite residual flag evaluators for direct geometric control, conditional interval evaluation and a fully conditioned closed-form comparator. The circuits include Gaussian quantile lookup, financial arithmetic, primitive function lookup, output comparison and arithmetic uncomputation. Against those initial comparators, the interval construction uses fewer assigned logical T gates at each of three caps and two monetary targets. The subsequent [stronger conditioning benchmark](conditioning_benchmark.md) reverses the interval/full ranking in all six cases. I retain this document as the initial construction record and withdraw its comparative interval-contribution interpretation.

For the £10 cap and a requested total pricing tolerance of £0.01, the interval flag evaluator uses 822,654 assigned T gates, compared with 1,207,388 for direct evaluation and 4,346,902 for the implemented full-conditioning formula. The corresponding reduction against direct evaluation is approximately 31.9 per cent. The interval construction requires 3,573 logical qubits in the retained schedule, compared with 3,426 for direct evaluation. Its lower T count therefore does not establish a simultaneous improvement in every resource.

I reserve half the requested tolerance for a future statistical estimator. This initial work certifies evaluator and finite-representation bias allowances; it does not produce a pricing confidence interval. Its counts remain correct for their specified circuits, but the proposed arithmetic does not survive the subsequently strengthened conditioned comparison.

## 2. Financial instance and estimands

I use the synthetic risk-neutral model

\[
S_1=100e^{-0.015-0.1Y+0.2Z},\qquad
S_2=100e^{-0.03+0.2Y+0.2Z},
\]

with independent standard normal Y and Z, one-year maturity, interest rate 0.03 and dividend yield 0.02. The covariance of log returns is

\[
\Sigma=\begin{pmatrix}0.05&0.02\\0.02&0.08\end{pmatrix}.
\]

The drifts equal \(r-q-\Sigma_{ii}/2\). I use equal weights, strike £100 and caps £1, £10 and £25. Sterling is a denomination convention for this stipulated instance, rather than a claim about an observed UK instrument. No historical observations are used.

I set \(A=(S_1+S_2)/2\), \(G=(S_1S_2)^{1/2}\) and \(F(x)=\min((x-100)^+,H)\). The price is

\[
V=e^{-0.03}\{\mathbb EF(G)+\mathbb E[F(A)-F(G)]\}.
\]

The direct evaluator retains Z. The interval evaluator retains a uniform strike coordinate U and evaluates the nonnegative normal interval described in the mathematical note. The fully conditioned evaluator integrates Z through the capped lognormal formula. These methods have different finite input laws. Their explicit finite-representation bounds connect them to the same continuous-model price; I do not equate their finite targets.

For the full comparator, I evaluate the capped conditional payoff as

\[
g(b)=b e^{\lambda^2/2}[\Phi(t_0+\lambda)-\Phi(t_1+\lambda)]
-K[\Phi(t_0)-\Phi(t_1)]+H\Phi(t_1),
\]

where \(t_0=\log(b/K)/\lambda\), \(t_1=\log(b/(K+H))/\lambda\), and \(\lambda=0.2\). The implementation uses exact integer subtraction, ordered CDF intervals, explicit product precision and final clipping. Subtraction does not lose floating-point digits in this representation. The formula's CDF error amplification is nevertheless included in its certificate. This is a certified closed-form constructor, rather than the earlier quadrature-order diagnostic.

## 3. Finite grids and monetary error allocation

I use uniform registers for Gaussian midpoint-quantile indices and the random strike. An independent review identified stronger bounded-variation certificates than the earlier global Lipschitz coupling estimate. Along any outer coordinate, a positive sum of exponentials is convex and has at most one minimum. Its composition with an increasing bounded payoff or conditional probability has total variation at most 2H. The geometric component is monotone, with variation at most H. The residual therefore has variation at most 3H.

An N-point midpoint rule has absolute integration error at most total variation divided by 2N. The outer error is consequently at most \(3H/(2N_Y)\) in this one-coordinate instance. Direct inner Gaussian discretisation and uniform-strike discretisation each cost at most \(H/N\), with their own register widths. The [finite-grid supplement](finite_grid_supplement.md) gives the proof and the tensor-coordinate extension. These are elementary applications of established variation arguments, not new quadrature theorems.

I allow at most 0.4 of the requested price tolerance for finite representation, at most 0.1 for arithmetic, and reserve 0.5 for statistical estimation. A separate £\(10^{-12}\) allowance covers the offline geometric correction. Discounting is below one, so undiscounted bounds are conservative price bounds. Every retained configuration stays within its total allocated bias. For the £10 cap at tolerance £0.01, the retained bias bounds are approximately £0.003809 for direct evaluation, £0.004460 for the interval construction and £0.004613 for full conditioning.

Each method receives its own Gaussian widths, fixed-point precision and primitive table sizes. The search refines table sizes greedily across fractional widths 18 to 30, ranks feasible configurations using a disclosed gate proxy and then emits their actual gates. This is an independent heuristic precision search for each method. I do not claim a globally optimal resource configuration, or describe the proxy as the measured T count.

The implementation uses nearest rounded table endpoints and floor-truncated products. Exp is tabulated on [-2,2], log on [0.25,4.25], and the normal CDF on [-8,8]. Inputs stay inside the declared domains; no uncharged CDF saturation is used. Piecewise linear interpolation has error at most \(\Delta^2\sup|f''|/8\), plus the endpoint and product-rounding allowances. I propagate all coefficient, quantile, cap and rescaling errors and check the signed register range. Multiplication by a known constant uses compiled set-bit arithmetic, including the factor 100 for the £1 cap, rather than placing an overflowing constant in a price register.

## 4. Reversible construction and quantum encoding

The [gate library](reversible.py) emits X, CNOT and exact Toffoli gates. Gaussian and nonlinear-function tables are one-dimensional primitive tables. I do not use a joint table mapping all financial inputs directly to a payoff. The lookup uses a depth-first selector tree, explicit data CNOTs and cleaned selector work. Table size, contents hashes and classical programme bits are retained; lookup is not a free QRAM operation.

Every financial graph computes an integer threshold q. A separate b-bit register R is uniformly prepared over \(0,\ldots,2^b-1\). The predicate \(R<q\) is XORed into a good flag, and the financial body is reversed. The threshold has b+1 bits, allowing the exact value \(q=2^b\). The resulting success probability is exactly \(q/2^b\), including probability one. Input registers are preserved and working registers return to zero.

Hadamards, reversible permutations and measurement define a genuine coherent quantum encoding. Basis-state replay establishes the finite permutation's correctness; it does not simulate a several-thousand-qubit state vector. The earlier small state-vector interference demonstration remains separately labelled. I also cost a complete Grover iterate from two flag evaluators and the appropriate reflections on the invariant subspace with clean work. The separate [phase-aware reflection circuit](reflections.py) is expanded and checked on 508 basis inputs, including its phase and clean work. A statistical estimation algorithm, controlled calls and its success-probability schedule remain unimplemented.

I assign the published exact seven-T, three-T-layer Toffoli decomposition documented by [Selinger](https://doi.org/10.1103/PhysRevA.87.042302). Gate totals are obtained by expanding each instantiated block, rather than applying an asymptotic primitive estimate. Large circuits are counted as compositions of emitted blocks and their exact inverses. Complete small tapes agree exactly with those composed counts. Reported complete depths are executable serial upper bounds, not optimised global T-depths. Uniform preparation adds its own Hadamards and one parallel layer.

## 5. Initial retained resource comparison

The table reports one clean residual flag evaluator. The monetary target includes a reserved statistical allowance which has not yet been realised. Qubit values are the upper bounds for the specified arithmetic schedule. All T counts use the same exact-Toffoli assignment.

| Cap | Requested price tolerance | Direct T count | Interval T count | Full conditioning T count | Direct qubits | Interval qubits | Full qubits |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| £1 | £0.10 | 530,334 | 455,910 | 2,318,666 | 2,981 | 3,171 | 8,106 |
| £1 | £0.01 | 794,836 | 456,414 | 4,056,528 | 3,419 | 2,846 | 10,509 |
| £10 | £0.10 | 566,818 | 460,614 | 2,342,550 | 2,988 | 2,846 | 8,109 |
| £10 | £0.01 | 1,207,388 | 822,654 | 4,346,902 | 3,426 | 3,573 | 11,028 |
| £25 | £0.10 | 635,894 | 551,670 | 2,270,114 | 2,991 | 3,200 | 8,110 |
| £25 | £0.01 | 1,663,284 | 1,347,402 | 4,552,282 | 3,428 | 4,472 | 11,030 |

The [machine-readable resource record](results/encoding_resources.json) contains precisions, error bounds, every emitted block's counts, serial T-depth bounds, qubits, complete Grover-iterate costs, table hashes and source hashes. Numerical runtime is explicitly a classical workflow diagnostic.

The interval construction reduces the assigned T count against direct evaluation by approximately 13.2 to 42.6 per cent across these six financial cases. The comparison identifies a potential precision and arithmetic benefit within this architecture. It does not identify an asymptotic improvement or establish useful physical execution times.

## 6. Initial verification and independent numerical evidence

Seventeen arithmetic test groups check expanded gate replay, boundary values, constant and variable products, input preservation, comparator flags, lookup correctness, wire reuse, clean work and streamed-versus-retained counts. Six complete financial tapes at two toy precisions agree exactly with their composed X, CNOT and Toffoli counts. Coarse finite financial circuits retain the comparison flag and reset all work; their deliberately coarse precision is not presented as a financial accuracy demonstration.

A separately implemented integer evaluator checks 3,220,200 retained-grid states across all 18 configurations and 29,120 exact crosschecks against the scalar graph. Full conditioning is exhaustive over its retained outer states. Direct and interval methods include systematic boundary and coarse pairs plus seeded random samples. These comparisons test implementation agreement and observed arithmetic errors. They do not substitute for the uniform analytical certificate.

Independent 100-digit calculations enclose all 61,786 generated function endpoints and all 58,304 retained quantile entries, using positive-half checks and exact symmetry where appropriate. The offline geometric correction is also independently enclosed. The [independent audit](results/independent_encoding_audit.json) records its source and configuration hashes. The normal CDF generator uses an alternating integral series with a first-omitted-term remainder after terms decrease and a conservative rounding enclosure; it is not treated as an exact libm call.

Independent positive-integral quadrature gives price diagnostics approximately £0.453422391, £3.745408049 and £6.740750694 for caps £1, £10 and £25. At cap £10, residual standard deviations are approximately 0.697912 for direct evaluation, 0.255281 for the interval construction and 0.254139 for full conditioning. Quadrature refinements do not certify those variance values. Full conditioning retains its statistical advantage; I do not use observed variance as a valid estimator stopping bound.

## 7. Interpretation and subsequent research gate

The favourable comparison is restricted to two assets, one stipulated covariance, terminal monitoring and the constructed arithmetic family. I have not established that this financial example is classically difficult. It is a numerical construction test, rather than evidence that adding a quantum device would improve practical pricing. REIT parameter scenarios and Asian monitoring remain later extensions.

The library's variable multiplier retains substantial arithmetic history. Its low-space alternative is implemented and separately tested, but the table above uses the history-retaining schedule. Efficient relative-phase logical ANDs, lower-space arithmetic and [published unary-iteration QROM](https://doi.org/10.1103/PhysRevX.8.041015) could reduce costs. I have not inserted their advertised counts into circuits which I did not implement.

The full comparator's certificate is conservative. Its exact conditional payoff has bounded financial sensitivity, whereas direct propagation through the displayed formula can amplify CDF and argument errors. Sharper consistent-input certificates, stable specialised functions or certified positive quadrature could materially reduce its resource count. The measured gap therefore cannot be described as inherent instability or unavoidable expense of analytical conditioning.

The next gate at this stage was to strengthen conditioning and rerun the same comparisons. The [completed audit](conditioning_benchmark.md) implements a scalar conditioned-payoff kernel and comparable interval specialisation. Full conditioning becomes cheaper in all six cases, so the general interval savings claim stops. A range-based estimator has no automatic variance gain from conditioning the same mean. A variance-sensitive estimator requires a valid bound and additional operations.

The useful project-specific result at this stage is a checked finite-grid and arithmetic construction with a positive gate comparison and retained evidence. The probability identities, variation inequality and standard gate decompositions are established methods. Whether the complete construction provides an original literature contribution remains a research question requiring comparison with the published conditioning and quantum-encoding literature.

## 8. Stronger-comparator outcome

For the £10 cap at £0.01 tolerance, the new full kernel uses 597,478 assigned T gates and 2,593 logical qubits, compared with 672,210 and 2,721 for the specialised interval kernel. Its total retained bias bound is approximately £0.004638. The new [30-configuration resource record](results/conditioning_resources.json) and [independent audit](results/conditioning_independent_audit.json) supersede the comparative interpretation of the earlier records without discarding them. The revised question concerns certified scalar conditioning and complete resource accounting; a positive original literature contribution remains unestablished.

## 9. Keywords

Quantum computational finance; capped basket options; geometric control; conditional probabilities; reversible arithmetic; bounded variation; logical resources; monetary bias certificates.
