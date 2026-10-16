# Certified phase estimation and finite-order pricing

## 1 Research question

I test how the quantum resource comparison changes when the amplitude estimator is replaced while its financial circuit, flag, probability prior, tolerance and confidence remain fixed. I also execute a classical certificate hierarchy which stops as soon as the required financial tolerance is proved. These comparisons address two weaknesses of an encoding-only gate saving: a conservative estimator can exaggerate an apparent resource benefit, and computing every retained classical order can overstate the work needed for a particular tolerance.

The model family remains the 27 stipulated capped Gaussian baskets. No empirical exotic price, LSEG observation or purchased dataset enters this phase. The undiscounted price radius is 0.01. Each selected quantum price has 99 per cent per-price confidence; this phase does not assign a simultaneous 99 per cent guarantee to all 27 prices.

## 2 Quantum construction and guarantee

I implement the canonical phase-estimation amplitude estimator of Brassard, Høyer, Mosca and Tapp (2002), Section 4 and Theorem 12. The publication and exact consulted sections are recorded in qpe_estimator_sources.md. For a scalar flag probability p with proved prior p in [0,u], the least supported binary phase size M satisfies the published amplitude bound, using a rational upper enclosure of pi and of the maximum square root of p(1-p), together with an outward classical post-processing allowance. If the prior midpoint already meets the scalar tolerance, no quantum run is required.

The controlled Grover iterate executes financial preparation and its logical adjoint unconditionally around controlled reflections. The conditional global minus is retained as a Z on the phase control. Corrected measured ANDs are admitted through their scalar-Kraus logical-channel proof, including external references. I do not control a measurement tape or reverse its measurements. The controlled zero reflection adds two exact Toffolis to the uncontrolled reflection, hence fourteen T gates under the specified expansion. All Clifford gates, workspace and serial depth caps are retained as well.

The inverse Fourier transform has actual retained Clifford and T gate sequences. I generate candidate single-qubit sequences with unmodified pygridsynth 1.2.0, following the published Ross and Selinger synthesis convention, and independently certify their products with 384-bit complex ball arithmetic. The profile API itself uses only the standard library and the retained certificates. Tiny-angle identity replacements have nonzero norm-error charges. Aggregate scalar phases are recorded; an uncontrolled final Fourier scalar does not change its measurement law, whereas the conditional Grover sign cannot be discarded.

The compiled Fourier norm error reduces the one-run success lower bound. I round that bound outward onto a binary lattice of width 2^-80 before exact rational confidence arithmetic. The smallest supported odd repetition count has its exact binomial majority-failure upper bound at most the assigned alpha. The median then has the required scalar accuracy under independent fresh repetitions. Phase preparation, controlled calls, emitted Fourier gates, terminal phase observations, financial observations for reset and conditional resets are all charged. Peak-qubit figures and additive serial depths are conservative upper bounds.

## 3 Matched and optimised findings

For each of the 27 previously selected strong single flags, I preserve the prefix, architecture, comparison exponent, encoding, scalar prior, tolerance and confidence. The canonical implementation lowers all 27 retained iterative T upper caps. Exact QPE-to-iterative ratios range from approximately 0.103136 to 0.418573. The comparison therefore identifies an estimator effect before any further encoding optimisation. It compares proved upper budgets, rather than measured stochastic costs or a lower bound on either algorithm's complexity.

I separately search all 819 inherited profiles using canonical QPE. Relative to its selected unit-range flag, range compression improves 16 contracts and ties 11. Clipping or range selection improves 18 and ties nine, with a minimum ratio approximately 0.123282. The dyadic family loses all 27 unit-range comparisons, with ratios approximately 6.899 to 75.377. This is a result for the retained bucket partition, omitted-prefix rule and allocations. It does not exclude a different published variance-sensitive estimator or partition.

The compressed complete archive retains every profile and flag choice. The matched comparison is a separate result so that changing the estimator is not conflated with changing the financial encoding. The old archives and their source histories remain unchanged.

## 4 Executed classical stopping hierarchy

I retain three finite-order policies with the same 224-bit outward arithmetic and the same undiscounted radius target. The full search evaluates all coefficient orders zero through eight. The remainder screen computes each candidate remainder first and evaluates its coefficients only when the remainder permits certification, or at the order limit to retain a valid unresolved enclosure. Every failed screen's moment and geometric work is charged. The screen is sufficient, not a proof that its selected order is globally minimal.

The third policy first uses the established monotonicity and geometric-control argument. Since the arithmetic-to-geometric ratio T is at least one, the arithmetic price exceeds the geometric price. A uniform log-price slope bound L and Jensen's inequality give a residual mean upper bound H L log E[T]. If the resulting one-sided interval has midpoint radius at most 0.01, it answers the price without evaluating higher Taylor coefficients. Otherwise the policy falls back to the remainder screen. Arithmetic uncertainty is included in the stopping decision.

All three policies certify the same 25 contracts. Full search uses 120,690 Gaussian kernel evaluations and 25,080 geometric exponential terms. Remainder screening uses 1,183 and 1,552; the monotonicity-first hierarchy uses 1,169 and 1,552. The latter resolves the three cap-one, lambda-0.4 contracts at its one-sided preliminary stage. These counts include shared geometry within each independently executed frozen family. Separate cold eight-asset cap-25 runs retain their own setup: screened kernel counts are 495, 45 and nine at lambda 0.05, 0.2 and 0.4, with 1,287, 165 and 45 geometric terms. Each full-search cold counterpart uses 12,870 kernels and 24,310 terms.

The approximately 99.03 per cent kernel-count reduction concerns these matched implementations of a tolerance-respecting certificate. The full search often delivers a much narrower interval than requested. This is not a comparison with the best published classical pricing algorithm. Host seconds are diagnostics, not a conversion between CPU work and quantum gates. The two remaining weak-smoothing two-asset controls retain valid but wide Taylor enclosures and are already resolved by validated continuous quadrature.

## 5 Why the figures and table are included

The degree-versus-radius figure exposes a feature that a successful endpoint table conceals: at weak smoothing the bound may cease improving and fail throughout the retained degree range. It explains the finite nature of the certificate and why its failure is not a pricing discrepancy. Common axes permit comparison of the frozen covariance geometries without asserting a dimensional scaling law. The compact first-certifying-degree table has a separate role, showing how the payoff cap changes the required integer degree.

The matched quantum figure isolates the estimator change. Its complete-cap panel places all 27 matched cases relative to the equal-cap line. Its ratio-versus-phase-width panel shows the relationship with the additional phase register. Colours and shapes identify genuine contract attributes; nested hollow markers preserve near-coincident exact coordinates without jitter. Phase width is not an exact peak-qubit increase. The figures use exact retained data, keep case-labelled rows and drawing coordinates, and have self-contained native LaTeX pictures together with SVG and PNG exports.

## 6 Evidence and remaining scope

The scalar tests cover coherent Fourier and complete phase-estimation laws, including non-palindromic wire ordering and the incorrect controlled-sign counterexample. The rotation audit checks all 62 signed controlled-phase products, complete inverse DFT columns through five phase wires, literal full tapes through three wires, and all 32 supported profile counts. The controlled-oracle audit checks 2,040 reflection basis inputs and 3,276 production profile/reset combinations, plus small literal measurement-instrument diagnostics. Its algebraic channel proof supplies the unrestricted reference-system guarantee.

The full benchmark auditor independently checks all 819 profiles, 43,260 single-flag choices, 1,638 dyadic plans and 21,012 bucket records, including confidence, norm rounding, support and tail rejections, weighted budgets and selection. The classical audit checks all executed policy records and cold presets against outward coefficient and remainder calculations. Diagnostic agreements are not substituted for a proof or an interval certificate.

The positive mathematical contribution candidate remains the evaluable finite coefficient/remainder certificate and its concrete consequences for selecting a pricing method. The canonical estimator and synthesis mechanisms are established. The stopping hierarchy is a checked implementation of sufficient certificates, not a new general Taylor or Jensen theory. Source-level comparison with the full Ju article and complete compatible modern quantum encodings, including QSP, remains necessary before a priority or algorithm-wide superiority claim. The current results establish no physical quantum advantage on this synthetic family, because classical price intervals already meet its stated tolerance.

## 7 Reproduction

The retained qpe_estimator.py, qft_synthesis.py, controlled_estimation_resources.py, qpe_benchmark.py and qpe_comparison.py use standard Python and retained certificates. Run qpe_benchmark.py before qpe_comparison.py. The full old and new archives are read without alteration. Use the existing python-flint environment for taylor_order_benchmark.py, certify_qft_rotations.py, audit_qft_synthesis.py and the full benchmark audit; no NumPy is required in that environment. The separate controlled-oracle coherent audit uses the existing NumPy environment. The optional requirements-rotation-synthesis.txt pins the isolated generation packages. Regeneration is optional because the literal tapes and certificates are retained.

The figure generators reside in figures/ and resolve their default sources from the repository's results/ directory. They use the existing Matplotlib 3.10.6 environment. The main manuscript incorporates their native picture code inline so that the existing document editor does not need external graphics assets. No new financial-data access, purchase or entitlement upgrade was used.

## 8 Keywords

Financial error certification; capped Gaussian baskets; quantum amplitude estimation; controlled Grover reflections; Clifford and T synthesis; finite-order selection; geometric control; computational finance.
