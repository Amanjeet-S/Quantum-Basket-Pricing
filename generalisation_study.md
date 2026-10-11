# Cap-dependent certificates and the cost of coherent basket pricing

## 1 Research question and result

I investigate the following question: for geometrically controlled capped lognormal basket derivatives, when do contract-dependent slope and curvature certificates permit cheaper complete coherent arithmetic than generic Gaussian certificates, and which financial and architectural conditions prevent that saving?

The retained comparison contains 27 synthetic financial cases and 108 evaluator configurations. A cap-aware full-conditioning certificate reduces the assigned logical T-count in three paired cases, all at common loading 0.05 and cap 25. The reductions against independently selected generic full conditioning are 13.108, 3.870 and 1.033 per cent for two, four and eight assets. The other 24 cases have no saving. Curvature improves at cap 10 under the same loading, but does not cross the retained precision threshold. This is a specific contract and constructor consequence, rather than a new pricing identity or a universal quantum efficiency result.

The [proof supplement](generalisation_proofs.tex) supplies the financial model, conditional identities, dominated differentiation, integration by parts, weak-loading sharpness, shared Gaussian-grid transfer, bounded-variation alternatives, generator constraints, complete rounding rules, Born-probability argument and exact gate-count derivations. I distinguish these checked project derivations from a claim of theoretical novelty. An original positive literature contribution and an end-to-end quantum advantage remain unverified.

## 2 Financial design and prespecification

I retain the [immutable financial protocol](results/generalisation_protocol.json), generated before the full resource comparison. The family crosses two, four and eight assets with common loadings 0.05, 0.2 and 0.4 and caps 1, 10 and 25, at monetary target 0.01. Initial prices and strike are 100, maturity is one year, the risk-free rate is 0.03 and dividend yield is 0.02. Every model is a stipulated risk-neutral correlated GBM model with strictly positive weights. The residual outer Gaussian dimensions are one, three and seven. Exact rational checks establish positive-definite covariance, the factor-covariance identity and the risk-neutral drift for the frozen finite-decimal coefficients.

Before freezing the family, I replaced a proposed residual-factor construction with full-rank Helmert contrasts for four and eight assets and matched weighted residual dispersion to the two-asset anchor at 0.0225, up to a retained discrepancy below 1.8 × 10⁻⁷¹. I did not select or discard financial cases after observing resource results. Matching this dispersion does not match the entire covariance, geometric-basket variance, payoff distribution or classical difficulty. I therefore make no inference that nominal dimension itself caused a measured cost change.

After a preliminary smoke run and before the full comparison, I derived a stronger outer-coordinate sensitivity bound for the continuously inner-averaged residual. I applied it to every method by transferring the outer law before the inner grid. The [algorithm amendment](results/generalisation_protocol_amendment.json) records this chronology; the financial protocol remains unchanged. The smoke run is retained separately and is not independent validation. The signed arithmetic width is q + 7, with strict absolute safe limit 64 and realised lookup-domain checks. Historical narrower-domain counts remain in the earlier reports and are not attributed to this certificate change.

## 3 Mathematical mechanism

For relative cap h = H/100, the established cap-normalised conditional payoff is

\[
P_{h,\lambda}(m)=h^{-1}\int_1^{1+h}\Phi((m-\log k)/\lambda)\,dk.
\]

I derive the global constants

\[
L_{h,\lambda}=\min\{\phi(0)/\lambda,(1+h)/h\},\qquad
C_{h,\lambda}=\min\{\phi(1)/\lambda^2,(1+h)\phi(0)/(h\lambda)\}.
\]

The alternative curvature branch improves precisely when λ < exp(−1/2)h/(1+h). For fixed positive h, I also prove that λ‖P″‖∞ tends to (1+h)φ(0)/h as λ tends to zero, so the cap-dependent weak-loading leading constant is sharp. I have not established novelty for these elementary Gaussian-smoothing calculations.

The discrete resource consequence needs another argument. Let w be the fixed-point width, b the probability-table address width, f = q + log₂(S) − b the remainder width, and I_b the number of visited internal nodes in the actual pruned selector tree. The forward interpolation block has exactly

\[
F_b=2I_b+2f(4w+3f)+6w
\]

Toffolis. Holding the full prefix fixed, one fewer address bit saves whole-flag Toffolis exactly when

\[
I_b-I_{b-1}>4w+6f+3.
\]

The complete saving is 4[2(I_b−I_{b−1}) − (8w+12f+6)]. This includes two probability calls and their exact financial-body inverse. A smaller table increases the remainder multiplier width, so an improved approximation bound alone does not imply cheaper arithmetic. Three retained compiler diagnostics explicitly exhibit cost increases upon coarsening. The proof derives each count from the emitted majority, multiply, selector and inverse construction, with width, disjointness, support and range premises stated.

## 4 Complete prespecified results

The following counts refer to one complete clean financial flag evaluator, including every emitted forward arithmetic block, the output comparison and the exact financial-body inverse. Each exact Toffoli is assigned seven T gates. Classical table entries and programme storage, serial depth upper bounds, logical qubit bounds and explicit Grover reflection costs are retained in the [resource record](results/generalisation_resources.json). There is no estimator repetition schedule or confidence guarantee. I rank the best retained emitted configuration among two proxy-shortlisted feasible candidates at each method's declared search outcome, rather than claiming a global circuit optimum. The cap-aware search also retains the generic winner as a feasible candidate, which ensures noninferiority within this candidate set but does not establish a strict gain.

| Assets | Common loading | Cap | Direct T | Interval T | Generic full T | Cap-aware full T | Cap-aware reduction |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 0.05 | 1 | 864,290 | 692,048 | 604,968 | 604,968 | 0.000 per cent |
| 2 | 0.05 | 10 | 1,276,674 | 1,278,564 | 1,048,320 | 1,048,320 | 0.000 per cent |
| 2 | 0.05 | 25 | 1,732,402 | 1,924,384 | 1,718,948 | 1,493,632 | 13.108 per cent |
| 2 | 0.2 | 1 | 971,544 | 496,720 | 399,196 | 399,196 | 0.000 per cent |
| 2 | 0.2 | 10 | 1,384,012 | 895,944 | 784,588 | 784,588 | 0.000 per cent |
| 2 | 0.2 | 25 | 1,839,824 | 1,512,952 | 1,371,132 | 1,371,132 | 0.000 per cent |
| 2 | 0.4 | 1 | 1,360,982 | 449,568 | 396,900 | 396,900 | 0.000 per cent |
| 2 | 0.4 | 10 | 1,658,678 | 854,112 | 760,032 | 760,032 | 0.000 per cent |
| 2 | 0.4 | 25 | 2,229,094 | 1,024,072 | 913,332 | 913,332 | 0.000 per cent |
| 4 | 0.05 | 1 | 2,339,064 | 1,787,492 | 1,670,144 | 1,670,144 | 0.000 per cent |
| 4 | 0.05 | 10 | 3,740,716 | 4,152,428 | 3,943,660 | 3,943,660 | 0.000 per cent |
| 4 | 0.05 | 25 | 6,719,664 | 6,108,088 | 5,797,148 | 5,572,812 | 3.870 per cent |
| 4 | 0.2 | 1 | 2,484,286 | 1,425,284 | 1,344,924 | 1,344,924 | 0.000 per cent |
| 4 | 0.2 | 10 | 3,885,854 | 3,467,212 | 3,326,736 | 3,326,736 | 0.000 per cent |
| 4 | 0.2 | 25 | 6,864,718 | 5,079,256 | 4,893,952 | 4,893,952 | 0.000 per cent |
| 4 | 0.4 | 1 | 3,087,672 | 1,223,264 | 1,065,400 | 1,065,400 | 0.000 per cent |
| 4 | 0.4 | 10 | 4,116,420 | 2,703,708 | 2,563,736 | 2,563,736 | 0.000 per cent |
| 4 | 0.4 | 25 | 6,177,696 | 3,479,728 | 3,296,216 | 3,296,216 | 0.000 per cent |
| 8 | 0.05 | 1 | 5,887,546 | 4,615,212 | 4,490,696 | 4,490,696 | 0.000 per cent |
| 8 | 0.05 | 10 | 11,933,978 | 12,726,504 | 12,399,716 | 12,399,716 | 0.000 per cent |
| 8 | 0.05 | 25 | 31,198,538 | 33,298,104 | 32,747,736 | 32,409,468 | 1.033 per cent |
| 8 | 0.2 | 1 | 6,631,394 | 3,860,584 | 3,796,436 | 3,796,436 | 0.000 per cent |
| 8 | 0.2 | 10 | 12,677,658 | 11,337,312 | 11,137,056 | 11,137,056 | 0.000 per cent |
| 8 | 0.2 | 25 | 31,942,050 | 18,790,520 | 18,572,176 | 18,572,176 | 0.000 per cent |
| 8 | 0.4 | 1 | 7,520,996 | 3,395,224 | 3,308,200 | 3,308,200 | 0.000 per cent |
| 8 | 0.4 | 10 | 10,757,320 | 7,358,316 | 6,954,444 | 6,954,444 | 0.000 per cent |
| 8 | 0.4 | 25 | 20,387,836 | 17,547,124 | 17,692,500 | 17,692,500 | 0.000 per cent |

In the original retained shortlist, direct evaluation is cheapest in three cases, interval conditioning in one, and cap-aware full conditioning in 23, including its ties with generic full conditioning. The eight-asset, weak-loading, largest-cap cap-aware improvement still does not beat the retained direct evaluator. The table therefore does not support universal superiority of full conditioning or a general percentage saving.

## 5 Causal and matched-prefix checks

The [fixed-prefix audit](results/generalisation_precision_effects.json) holds the generic winner's q, Gaussian width, exponential width and logarithm width fixed. It finds the minimum feasible probability-table address width in the declared range for each certificate and exactly emits the resulting graph. The same three cases move from 13 to 12 probability address bits. The resulting T savings are 116,816, 116,424 and 116,032. Both the block formula and the whole-flag difference agree exactly with emitted counts. Independent retuning produces additional prefix changes, so its larger aggregate saving must not be described as entirely attributable to those fixed-prefix counts.

The [matched-interval diagnostic](results/generalisation_matched_interval.json) was added after inspecting the completed study. It preserves all original outputs and tests all 27 interval prefixes. Replacing their strike/CDF construction by the full cap kernel at identical q, Gaussian and function-address widths remains certified and is strictly cheaper in every case. At the apparent interval winner, eight assets, λ = 0.4 and H = 25, the matched full evaluator uses 17,363,612 T gates against interval's 17,547,124. This is cheaper than the original retained full winner at 17,692,500. Thus the original exception identifies heuristic search sensitivity; it does not identify an intrinsic interval advantage. I label the additional diagnostic as post-study rather than retrospectively prespecified.

The exact matched-graph difference is 2T_strike + 24w + 8(I_Ψ − I_P) Toffolis. It follows from removing one strike-log block and two argument subtractions and from the full-kernel packed nonzero support being contained in the CDF support. The generator projects full centres below corresponding CDF centres, preserves their error enclosures and takes prefix maxima before common rounding. This support argument orders selector Toffolis for this compiler. It does not order CNOT popcounts, every depth schedule or other QROM designs.

## 6 Uniform certificates and independent evidence

I reserve 0.4 of the monetary tolerance for distribution transfer, 0.1 for arithmetic and 0.5 for unimplemented statistical estimation. A 10⁻¹² geometric correction allowance is included within the deterministic half. The [independent 100-digit certificate recurrence](results/generalisation_independent_certificates.json) passes all 108 budgets, signed-safe limits and realised function-address domains, the Gaussian-grid b ≤ 18 premise, all 27 geometric constants and the CDF saturation allowance. Its internal decimal strings retain the independently recomputed bounds. Production domain-check binary floats are diagnostic summaries, not outward interval endpoints; I do not describe the implementation as a universal directed-rounding interval library.

For m in [−2,2], λ in [0.05,0.4] and h in [0.01,0.25], the normal CDF is clipped outside [−8,8]. Mills' upper bound gives CDF tail error below 6.4 × 10⁻¹⁶. The four-CDF coefficients amplify this by less than 1.169 × 10⁻¹² uniformly over the declared family. The charged 2 × 10⁻¹² full-kernel ordinate allowance includes the remaining decimal arithmetic and argument perturbations. I retain the series remainder, numerical-library assumptions and operation envelope in the proof, rather than assigning the old unclipped allowance outside its valid domain.

The [independent replay audit](results/generalisation_independent_audit.json) checks 1,619,282 distinct production-grid states, 27,648 scalar cross-checks, 348,537 function endpoints, 992,768 quantile entries, 1,045,248 primitive boundary and adjacent states, and 391 comparator, input-preservation and uncomputation checks over nine coarse complete gate tapes. The endpoint condition requires the entire independent reference enclosure to fit within the charged ordinate allowance; intersection of intervals is insufficient. Sampling checks support implementation consistency and do not establish a uniform error bound. Complete production-scale gate simulation was not performed. All 20 existing arithmetic and conditioning-contract tests pass.

The [compiler audit](results/generalisation_compiler_audit.json) checks 44 emitted interpolation blocks, three whole-graph count differences and 1,404 graph variants for non-identical intended constants aliasing after rounding. No such alias is found. The selected-premise audit finds no missing unsigned-input premise: the only lower enclosures below zero are at physically clamped direct residuals, which are nonnegative by construction.

## 7 Classical reference and financial interpretation

The [independent classical reference](results/generalisation_classical_reference.json) uses the positive strike integral, eight independent 32-bit digital shifts of 4,096 Sobol points, an equal-size plain Monte Carlo comparison and 48/96-point positive strike quadrature. The author-distributed Joe–Kuo direction data, algorithm source and BSD licence are retained under public_data with exact hashes. Digital shifting is not Owen scrambling. I consulted the author files and licence; I do not claim to have read the related journal papers. These free files require no new entitlement.

Across the 27 cases, the paired strike-quadrature refinement is at most 1.34 × 10⁻¹⁵, and the maximum per-state signed-formula versus positive-integral discrepancy is below 5.87 × 10⁻¹⁴. The largest replicate standard-error diagnostic is about 2.97 × 10⁻⁴. The two-asset 48/96-point outer Gauss–Hermite refinement is at most 3.40 × 10⁻¹⁰. These are convergence and consistency diagnostics, not rigorous quadrature remainders or pricing confidence intervals. Host workflow duration is not a quantum resource and the reference is not asserted to be the fastest classical implementation. Its prices are undiscounted; multiplication by exp(−0.03) gives the stipulated time-zero values.

All experiments remain synthetic. No LSEG observations, provider demonstrations, reconstructed private prices or publisher PDFs enter this evidence. A later REIT scenario would be a parameter application requiring a documented field and adjustment audit, rather than the source of the contribution. The derivative methodology is not restricted to real estate and does not model observed REIT jumps, changing dividends or volatility smiles.

## 8 Literature comparison and contribution status

The [literature audit](literature_audit.md) and [source ledger](SOURCES.md) distinguish final journal text, accessible author versions and abstract-only access. The closest published obstacle is Chakrabarti et al., Quantum 5, 463 (2021), whose final sections and Appendix C explicitly use derivatives, domains and arithmetic errors and anticipate component-dependent precision. Bhaskar et al. (2016) already analyse finite-precision scientific functions. Bayer et al. (2018) establish analytical basket conditioning. Stamatopoulos and Zeng (2024) establish coherent comparison and QSP payoff encodings. O'Brien and Sünderhauf (2025) optimise piecewise QSVT approximations. The published successors of Chen, Li and Neufeld and of McArdle, Gilyén and Berta were verified, with actual detailed access recorded accurately. A compound title assembled from those ingredients is insufficient to establish novelty.

The positive project-specific finding is a checked conjunction: a cap/common-loading regime can cross a monetary precision threshold, and an exact architecture condition determines whether that crossing reduces the emitted complete evaluator cost. The retained family exhibits both strict savings and null cases; the compiler supplies counterexamples to the inference that less table precision always means lower cost. I have not established that this conjunction is absent from the entire literature, or that its significance is sufficient for a top-journal paper. It is a defensible candidate contribution supported by explicit mathematics and retained evidence.

The next research barrier is an architecture comparison against stronger published QROM and known-function preparation methods, followed by a fully costed statistical estimator with a valid monetary confidence guarantee if an algorithm-level claim is pursued. A heuristic shortlist should inherit matched opponent configurations to prevent avoidable ranking artefacts. I must retain the current null and adverse results when enlarging that comparison. Neither a favourable synthetic case nor a rigorous presentation establishes end-to-end quantum advantage.

## 9 Reproduction and local status

The resource construction, fixed-prefix audits and compiler audit use the Python standard library. The classical reference and full replay audit require NumPy. With the existing requirements installed, the principal commands are:

```sh
python generalisation.py
python checks/audit_generalisation_certificates.py
python checks/classical_reference_generalisation.py
python checks/audit_generalisation.py --certify-quantiles
python checks/audit_precision_effects.py
python checks/audit_matched_interval.py
python checks/audit_generalisation_compiler.py
python -m unittest discover -s checks -p 'test_*.py' -v
```

I do not regenerate the financial protocol in place. The generator deliberately refuses to overwrite it. Its SHA-256 is ae0953d7100403118432ad16279955cd2c98ac9c75c0d5de162eb7d801e6d4af. Resource and audit files preserve their source hashes, selected parameters and scope qualifications. The complete resource workflow took about 363 seconds on the available host; this measures construction and certification work only.

The repository remains private. This step is saved locally without staging, committing or pushing changes, in accordance with the instruction to make no further commits.

## 10 Keywords

Quantum financial engineering; capped basket derivative; mathematical finance; conditional expectation; fixed-point arithmetic; Gaussian-grid transfer; constructive resource analysis.
