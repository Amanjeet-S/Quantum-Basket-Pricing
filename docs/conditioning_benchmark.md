# Strengthening analytical conditioning in a reversible basket-pricing benchmark

## 1. Finding and revised research question

I implemented a stronger fully conditioned evaluator and gave the random-strike interval evaluator comparable scalar-function specialisation. Full conditioning has the lowest assigned T count of the five constructed alternatives in all six retained cap and tolerance cases. The earlier claim that the interval construction survives its best conditioned opponent is therefore withdrawn. The initial gate counts remain valid for their particular circuits, but they do not justify the proposed comparative contribution.

For a £10 cap and a requested total pricing tolerance of £0.01, the specialised full-conditioning evaluator uses 597,478 assigned T gates. The comparably specialised interval evaluator uses 672,210, the initial interval evaluator uses 822,654 and direct geometric-residual evaluation uses 1,207,388. Full conditioning reduces the specialised interval count by 11.1 per cent and the earlier eight-CDF full formula by 86.3 per cent. These are logical evaluator costs under the selected exact-Toffoli schedule. Half the requested tolerance remains reserved for an unimplemented statistical estimator.

I now ask: under a fixed monetary bias budget, how much does certified scalar-function specialisation change the complete reversible cost of direct sampling, analytical conditioning and randomised-strike conditioning of a geometrically controlled capped basket payoff? This question supports a reproducible numerical-finance investigation. It does not presume that the interval method is preferable or that scalar tabulation is original. An original contribution remains unestablished; selecting the cheapest established constructor is insufficient.

## 2. Financial target and conditioning kernel

I retain the synthetic two-asset risk-neutral model, geometric correction, caps and monetary targets from the [initial study](reversible_encoding_study.md). No market observations or new data entitlements enter this comparison. With strike normalised to 1, relative cap \(h=H/100\), common loading \(\lambda=0.2\) and conditional log level \(m\), I define

\[
P_h(m)=\frac1h\mathbb E_Z\min\{(e^{m+\lambda Z}-1)^+,h\}
=\frac1h\int_1^{1+h}\Phi\!\left(\frac{m-\log k}{\lambda}\right)\,dk.
\tag{1}
\]

The strike integral follows by integrating the exercise indicator. It is an established layer-cake identity applied to the known lognormal conditional law. Gaussian analytical conditioning is established in [Bayer, Siebenmorgen and Tempone, Quantitative Finance 18, 491–505 (2018)](https://doi.org/10.1080/14697688.2017.1308003); I consulted the author accepted manuscript, section 3, as recorded in [SOURCES.md](SOURCES.md). I do not claim a new conditioning identity or imply that the final journal text was consulted.

The exact full residual is \(H\{P_h(m_A)-P_h(m_G)\}\), where

\[
m_A=\log\{(e^{-.015-.1Y}+e^{-.03+.2Y})/2\},
\qquad m_G=-.0225+.05Y.
\]

The arithmetic mean exceeds the geometric mean, so the exact normalised residual belongs to [0,1]. Full conditioning integrates the common Gaussian coordinate without introducing the random strike. The specialised interval residual remains

\[
H\{\Psi(m_A-L_h(U))-\Psi(m_G-L_h(U))\},
\quad \Psi(d)=\Phi(d/\lambda),\quad L_h(u)=\log(1+hu).
\tag{2}
\]

The mean of (2) over U equals the full residual. Conditional expectation also gives \(\operatorname{Var}(g_R(Y))\leq\operatorname{Var}(h_R(Y,U))\). This is an established variance inequality, not a new empirical finding or an estimator stopping certificate.

## 3. Financial error certificate

Differentiation of (1) under its finite integral gives

\[
0\leq P_h'(m)\leq L=\frac1{\lambda\sqrt{2\pi}},
\qquad |P_h''(m)|\leq C=\frac1{\lambda^2\sqrt{2\pi e}}.
\tag{3}
\]

The same bounds apply to \(\Psi\). At \(\lambda=.2\), \(L\approx1.994711402\) and \(C\approx6.049268113\). I use these intrinsic financial sensitivities in place of propagating one log-level error independently through eight CDF terms and large strike coefficients. The checked derivation and generator assumptions are retained in the [certificate review](reviews/checks/kernel_certificate_review.md).

For a scalar probability table on [-1,1] with \(2^b\) cells, spacing \(\Delta=2/2^b\), fractional step \(\delta=2^{-q}\) and endpoint-generation enclosure \(\epsilon_K\), the implemented integer interpolation obeys

\[
E_K\leq C\Delta^2/8+3\delta/2+\epsilon_K.
\tag{4}
\]

This includes nearest endpoint rounding and the floor in the interpolation product. If the computed input log errors are \(e_A,e_G\), a sufficient monetary full-evaluator allowance is

\[
\beta_{\rm full}\leq H\{L(e_A+e_G)+2E_K\}.
\tag{5}
\]

The implementation propagates the prefix's quantile, coefficient, exponential, averaging and logarithm errors before applying (3). It also conservatively propagates half-unit allowances for constants and clipping, even where the constant 1 is represented exactly. These extra allowances make the retained certificate slightly looser than (5). Exact equal-scale integer subtraction introduces no rounding. Projection to [0,1] is non-expansive relative to the exact residual, and the existing comparison-register encoding adds no second output floor.

For the specialised strike-log table, \(0\leq L_h'\leq h\) and \(|L_h''|\leq h^2\). I retain its interpolation, endpoint and product-floor errors separately. The uniform midpoint U is represented exactly at the chosen width. Both occurrences of its computed strike log are charged conservatively; I do not assume that their effects cancel.

I preserve the earlier distribution-transfer bounds: the outer midpoint Gaussian grid contributes at most \(3H/(2N_Y)\), and the direct Gaussian or random-strike coordinate contributes at most \(H/N\), with the earlier coupling alternative used when smaller. I allocate at most 0.4 of the requested tolerance to distribution representation, at most 0.1 to arithmetic and £\(10^{-12}\) to the separately enclosed geometric correction. The remaining half is reserved for statistical estimation. For the £10, £0.01 case, retained total bias bounds are approximately £0.003809 for direct evaluation, £0.004667 for the specialised interval evaluator and £0.004638 for specialised full conditioning. These are deterministic bias allowances, not confidence intervals.

## 4. Implemented scalar functions and complete resource boundary

The full evaluator computes both asset exponentials, their arithmetic average, its logarithm and the geometric log. It then evaluates \(P_h\) at the two computed scalar logs, subtracts and clips to [0,1]. The specialised interval evaluator retains the same financial prefix, computes \(1+L_h(U)\) through a cap-specific scalar table, and calls the domain-reduced \(\Psi\) table twice. It therefore receives equivalent access to known-function specialisation and range reduction.

I tabulate both probability functions on [-1,1], using an explicit bias of 1. The actual mathematical arithmetic log minimum is approximately -0.0766330123, reached at \(Y=(.015-\log2)/.3\), and its maximum for \(|Y|\leq4.75\) is approximately 0.445300. The geometric log lies in [-0.26,0.215]. I check every implemented domain including prefix errors. For a coarse interval search seed, a domain check can fail even though refined prefixes are feasible. The tuner refines that prefix before proceeding instead of discarding the method.

The cap-specific full-kernel ordinates are generated offline through the four-CDF capped-lognormal formula at 70 decimal digits. All its CDF arguments lie inside approximately [-6.116,5.2] for the declared domains and caps, within the existing [-8,8] generator certificate. Division by h amplifies the declared per-CDF enclosure. I therefore retain \(10^{-34}\) per kernel ordinate, rather than reuse \(10^{-38}\) without amplification. A separate 100-digit implementation encloses the rounded retained ordinates. The strike-log and scaled-CDF tables are independently enclosed too.

These are one-dimensional known-function tables addressed by computed scalar financial quantities. Neither evaluator uses a direct Gaussian-index-to-payoff table or a joint financial-input lookup. Every table selector, data CNOT, interpolation product, financial prefix, output flag and inverse is counted. Cap dependence is part of each table's identity and source manifest. Classical programme bits and measured incremental cached generation time are retained, without interpreting setup time as quantum runtime or assuming that it is amortised over many prices.

The lookup remains the explicitly emitted exact-Toffoli prefix-tree constructor. I assign the same [Selinger seven-T decomposition](https://doi.org/10.1103/PhysRevA.87.042302) to all alternatives; the detailed version consulted was the author manuscript recorded in SOURCES.md. I retain serial depth and qubit upper bounds. More efficient published QROM, different arithmetic schedules and synthesis could change the ranking. I do not borrow their advertised counts for unimplemented circuits.

## 5. Matched retained comparison

All five methods receive independent finite widths and heuristic precision searches under the same allocation rule. A disclosed analytical gate proxy ranks feasible configurations; actual resource totals come from emitted blocks and their inverses. I do not claim global precision optimality. The original three configurations reproduce their earlier gate totals exactly.

| Cap | Requested tolerance | Direct T | Initial interval T | Initial full formula T | Specialised interval T | Specialised full T |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| £1 | £0.10 | 530,334 | 455,910 | 2,318,666 | 391,846 | 346,542 |
| £1 | £0.01 | 794,836 | 456,414 | 4,056,528 | 385,182 | 340,690 |
| £10 | £0.10 | 566,818 | 460,614 | 2,342,550 | 385,182 | 340,522 |
| £10 | £0.01 | 1,207,388 | 822,654 | 4,346,902 | 672,210 | 597,478 |
| £25 | £0.10 | 635,894 | 551,670 | 2,270,114 | 452,634 | 399,434 |
| £25 | £0.01 | 1,663,284 | 1,347,402 | 4,552,282 | 1,303,302 | 1,193,122 |

Specialised full conditioning reduces T count against the comparably specialised interval method by approximately 8.5 to 11.8 per cent. The financial instance has not changed. The reversal arises from the implemented scalar conditioning functions and their certificates. It invalidates the inference that the initial interval/full gap reflected an unavoidable cost of full analytical conditioning.

For the £10 cap and £0.01 target, the full kernel's retained fractional width is 21, its Gaussian index width is 12, its exponential and log address widths are 10, and its probability address width is 9. The specialised interval configuration uses those same arithmetic widths and adds a 20-bit strike register and strike-log table. The selected full comparator also improves its serial resource upper bounds in this case.

| Constructor | Assigned T | Serial T-depth bound | Logical qubit bound | Classical table entries | Classical programme bits | Grover-iterate T |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Direct | 1,207,388 | 298,080 | 3,426 | 10,240 | 380,928 | 2,415,476 |
| Specialised interval | 672,210 | 127,140 | 2,721 | 7,680 | 304,128 | 1,345,148 |
| Specialised full | 597,478 | 113,550 | 2,593 | 6,656 | 248,832 | 1,195,404 |

The Grover totals include two clean flag evaluators and the explicit reflections. They exclude a mean estimator's repetitions, controlled calls, confidence rule, fault tolerance and physical runtime. The [resource record](../results/conditioning_resources.json) retains all 30 configurations, component counts, primitive domains, table hashes, setup diagnostics and source hashes.

## 6. Constructor-specific explanation

The ranking has a structural explanation within the implemented family. Since every strike in (1) is at least 1, \(P_h(m)\leq\Psi(m)\). For m in [-1,1] and h in [0.01,0.25], integrating only the upper half of the strike interval gives the positive gap

\[
\Psi(m)-P_h(m)\geq\frac52\log(1.005)\,
\phi\{5+5\log(1.25)\}>3.7583\,10^{-11}.
\tag{6}
\]

This exceeds the combined declared generator allowances. Same-width nearest rounding therefore preserves the ordering of the actual table ordinates. A packed interpolation row is zero exactly when both of its endpoint ordinates are zero. Every nonzero full-kernel row is consequently a nonzero scaled-CDF row. The exact lookup tree visits an internal selector only if its subtree contains a nonzero row, and emits two Toffolis per visited selector. Its full-kernel selector count \(Q_P\) is no greater than its scaled-CDF count \(Q_\Psi\).

At matching fractional, prefix and probability precisions, the two graphs have the same financial prefix, interpolation widths, argument subtractions, residual subtraction, final projection and output comparison. The full bias constant requires X gates. The interval method additionally evaluates the strike-log block with forward Toffoli count \(Q_{\log k}\). Including the complete arithmetic inverse, their flag-evaluator Toffoli difference is exactly

\[
Q_{\rm interval}-Q_{\rm full}
=2Q_{\log k}+4(Q_\Psi-Q_P)\geq0.
\tag{7}
\]

The full arithmetic certificate is also no larger throughout the searched widths: the interval strike-log error adds at least one fractional unit to each probability argument relative to the conservatively charged full bias constant. That allowance exceeds the difference between the two generator allowances. Full conditioning removes the strike-grid bias too. I require both perturbed-domain certificates to pass; arbitrary interval-domain feasibility does not prove full-domain feasibility in another financial model.

The [checked argument](reviews/checks/kernel_dominance_review.md) and [emitted-count audit](../results/kernel_dominance_audit.json) establish this mapping for all six retained interval configurations. This is an argument about assigned Toffoli and T counts in the selected constructor, rather than optimal quantum circuits. Integer ordering does not order binary popcounts, so it does not imply CNOT, depth, classical preprocessing or universal resource dominance. It does not determine complete estimator costs. I treat the cost identity as a checked project-specific argument and do not claim that it is a new literature theorem.

## 7. Verification

Twenty test groups pass: 17 arithmetic groups and three certificate-premise regression groups. An adversarial mutation exposed that the old certificate assumed ordered CDF operands. I now require an exact structural proof of their argument ordering before applying the nonnegative interval bound; reversed or unrelated arguments are rejected. This correction leaves the retained graphs and numerical bounds intact. Ten complete financial tapes at two coarse precisions agree exactly with the block-composed X, CNOT and Toffoli totals; their serial space and depth bounds remain conservative. Complete flag replay also checks both new variants, input preservation and zero working registers. These coarse circuits test implementation semantics, not monetary accuracy.

The [independent audit](../results/conditioning_independent_audit.json) replays 4,841,660 production-grid states and 74,496 coarse states, with 45,952 exact scalar crosschecks. It uses a separate integer interpreter and positive strike integration for its financial reference. It independently encloses 74,931 retained function endpoints and 80,000 inverse-quantile entries, checks 23,040 special primitive boundary states, and performs 171 complete new-variant flag and uncomputation checks. Full-kernel values agree with positive 64- and 192-point quadrature within approximately \(3.89\,10^{-16}\). The quadrature refinements are numerical diagnostics, not directed continuous-price or variance enclosures. The analytical certificate and endpoint enclosures supply the deterministic guarantees; sampled errors do not replace them.

## 8. Consequences for the project

The proposed interval encoding has failed the stronger-comparator gate in the initial evaluation family. I retain the code and negative result, but I withdraw the claimed interval contribution. The specialised full evaluator is an improved project artefact; its ingredients are established and its implementation alone does not demonstrate an original literature contribution.

The next research gate is a cost and approximation argument that generalises beyond this single easy model. I will first test whether a certified precision-allocation or scalar-kernel construction yields an informative, defensible improvement against equally specialised alternatives over a prespecified range of loadings, caps and basket structures. I will require literature comparison and retained evidence for that argument. If the result is only a lower envelope of existing methods, the contribution remains a reproducible feasibility finding and the positive-contribution requirement remains unmet. I will not proceed to a large estimator implementation on the premise that the discarded interval advantage survives.

The present two-asset problem is readily priced classically. It supplies a controlled numerical test and a genuine coherent quantum evaluator, rather than evidence of practical quantum advantage. REIT scenarios remain optional illustrations after the methodological question survives; changing assets cannot repair the failed comparison. No additional private data or purchased access is necessary for the current computational investigation.

## 9. Keywords

Quantum computational finance; capped basket derivatives; analytical conditioning; scalar function compilation; geometric controls; fixed-point error certificates; reversible arithmetic; comparative resources.
