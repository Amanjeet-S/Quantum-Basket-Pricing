# Deterministic classical comparators for the retained pricing target

## 1. Purpose and scope

I add two deterministic competitors before drawing a conclusion from the retained quantum and Monte Carlo confidence budgets. Exact finite summation prices the same financial interpreter outputs without statistical error. Validated continuous quadrature also prices the stipulated two-asset Gaussian models without the finite normal grid. These calculations retain the existing financial models and inherited precision pool. Exhaustive summation and validated quadrature are established methods, so their implementation is evidence for a stronger comparison, not an original contribution by itself.

The requested undiscounted price error remains 0.01 in the protocol's unspecified currency unit. I do not spend the unused statistical allowance by silently relaxing another numerical budget. Model error, market calibration and physical quantum runtime remain outside these calculations.

## 2. Exact finite summation

The inherited pool contains 273 financial prefixes. I record the exact tensor input count for each prefix. Direct pricing has 2^(d b_y+b_z) states, interval conditioning has 2^(d b_y+b_u), and full conditioning has 2^(d b_y). The prespecified enumeration limit is 131,072 states, and actual enumeration is restricted to the two-asset full_global and full_cap graphs. All 47 such prefixes are enumerated; the largest has 16,384 outer states. The total is 330,752 financial outputs.

Every output is a fixed-point integer between zero and 2^q. The sum is exact: with at most 2^17 outputs and q at most 30, its maximum is 2^47, below the unsigned 64-bit limit. Division gives an exact rational finite mean. Adding the retained geometric correction centre gives a rational price point; the previously independently checked deterministic bias certificate supplies exact lower and upper endpoints. Every price radius is below 0.005. Statistical error is zero, with confidence one conditional on the numerical certificate and the stipulated model. A pseudorandom generator and a concentration inequality are unnecessary.

The operation profile includes financial evaluator work, address assignments and increments, input and function-domain comparisons, Boolean reductions, exact output summation and rational mean construction. The selected full graphs charge 90N−8 online scalar-lane algorithm operations. These counts are neither CPU instructions nor a physical quantum cost. Table generation, table manifests, packed-word reads, unpacked-word writes, graph nodes, storage and preparation timings are reported separately. A disclosed shared cache avoids repeated setup within this run; it does not make one-off preparation free.

| Case | Enumerated states | Matched certified Monte Carlo samples | Selected online operation upper count |
|---|---:|---:|---:|
| λ=0.05, H=1 | 512 | 11,099 | 46,072 |
| λ=0.05, H=10 | 4,096 | 1,053,794 | 368,632 |
| λ=0.05, H=25 | 16,384 | 1,665,354 | 1,474,552 |
| λ=0.2, H=1 | 512 | 1,382 | 46,072 |
| λ=0.2, H=10 | 4,096 | 72,670 | 368,632 |
| λ=0.2, H=25 | 16,384 | 427,892 | 1,474,552 |
| λ=0.4, H=1 | 512 | 884 | 46,072 |
| λ=0.4, H=10 | 4,096 | 23,562 | 368,632 |
| λ=0.4, H=25 | 8,192 | 120,299 | 737,272 |

The finite-state to matched certified-Monte-Carlo sample ratio is between 0.003887 and 0.580499 across all 47 prefixes. These ratios compare financial-output evaluations; the complete enumeration count additionally charges input construction, validation and reduction. The chosen prefix minimises the retained conservative online count, with separately reported setup cells as a tie-breaker. I do not claim global precision optimisation.

The four-asset full grids contain between 2^30 and 2^45 states; the eight-asset full grids contain between 2^77 and 2^119 states. They are not enumerated under the prespecified state cap. These counts do not prove that another deterministic or sparse-grid method is infeasible.

## 3. Validated continuous quadrature

For each of the nine frozen two-asset models, I integrate the common-factor conditional residual Cλ(mA(y))−Cλ(mG(y)) against the normal density. Weighted arithmetic–geometric inequality and monotonicity of the capped payoff place this residual in [0,H]. The omitted region outside [−L,L] therefore contributes a one-sided interval [0,2H Φ(−L)]. The independently enclosed geometric capped price supplies the correction. The derivation in proof_sources/deterministic_baseline_proof.tex proves the positive strike representation, its signed analytic expression and the complete tail enclosure.

The existing authorised python-flint 0.9.0 runtime provides Arb/Acb ball arithmetic. The primary price uses 128 bits, L=8 and a requested absolute integral tolerance of 10^−16. The actual returned integral ball, geometric correction ball and outward normal-tail bound determine the certificate. A requested integration tolerance is never substituted for an achieved enclosure. All endpoints are retained as exact dyadic rationals, with outward decimal displays. A separately represented direct conditional arithmetic price uses 176 bits, L=10 and requested tolerance 10^−24. Every pair of price enclosures overlaps.

Acb integration checks complex neighbourhoods. I explicitly forward its analytic flag into the complex log of the sum of exponentials. This verifies branch-cut avoidance on each requested neighbourhood; the real integration path alone does not justify analyticity. A non-finite analytic value must lead to subdivision or failure. No clipping or heuristic error estimate enters the price certificate.

| Case | Certified undiscounted price, rounded for display | Actual primary enclosure width upper bound |
|---|---:|---:|
| λ=0.05, H=1 | 0.484579568771280 | 1.35×10^−15 |
| λ=0.05, H=10 | 2.782484533941027 | 1.26×10^−14 |
| λ=0.05, H=25 | 3.329903250937710 | 3.14×10^−14 |
| λ=0.2, H=1 | 0.467231158256496 | 1.26×10^−15 |
| λ=0.2, H=10 | 3.859472705514944 | 1.25×10^−14 |
| λ=0.2, H=25 | 6.946037115297563 | 3.12×10^−14 |
| λ=0.4, H=1 | 0.424776887288025 | 1.25×10^−15 |
| λ=0.4, H=10 | 3.833571237761709 | 1.25×10^−14 |
| λ=0.4, H=25 | 8.083183793564392 | 3.12×10^−14 |

The complete primary calculation uses 2,112 integrand calls. Primary and direct-reference calculations together take approximately 0.185 seconds on the retained host. This is a diagnostic timing, not a computational complexity claim, a portable performance promise or a quantum runtime comparison. It excludes starting the Python process and importing the library. Discounted price intervals are separately retained using an outward exponential discount factor; the comparisons above preserve the original undiscounted target.

## 4. Independent verification

The independent scalar Python-integer interpreter replays every state of every enumerated prefix and reproduces all 330,752 outputs, hashes and exact sums. The audit verifies all 273 tensor counts, rational endpoints, selection rules and online operation components. Every finite-grid price certificate contains the corresponding narrower continuous price enclosure.

The continuous audit checks 126 scalar kernel states through a separate positive strike-integral representation at 192-bit precision. Arithmetic and geometric kernels are tested at y equal to −8, −3, −1, 0, 1, 3 and 8 for every case. All signed-form and positive-integral enclosures overlap without clipping. The normal-tail bound is also checked against an outward Mills bound. These independent representations share the validated Arb library; the project does not claim independent implementation of that numerical library.

## 5. Implication for the contribution

The low-dimensional quantum interpretation must become more demanding. Improving on a certified Monte Carlo sample budget does not establish a quantum pricing advantage when exact finite summation is already cheaper in financial-output evaluations and validated continuous integration supplies narrow certificates with a few hundred scalar evaluations per case. This adverse evidence removes a weak route to claiming originality or practical advantage.

The project can still seek value in general financial residual moment bounds, certified transfer to implemented outputs, higher-dimensional oracle costs and the separation between fractional-output variance and Bernoulli-flag variance. These deterministic findings specify a domain where those results must be compared with stronger published classical methods. They do not establish a new theory, exclude higher-dimensional applications or prove a classical lower bound. Conditional smoothing and deterministic basket pricing already occur in the published literature, including Bayer, Siebenmorgen and Tempone (2018), as recorded in the project source audit.

## 6. Keywords

Deterministic pricing; validated quadrature; finite-grid enumeration; capped basket options; quantum pricing comparison; numerical error certificates.
