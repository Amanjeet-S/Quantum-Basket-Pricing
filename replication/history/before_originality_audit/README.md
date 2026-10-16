# Certified curvature corrections for quantum basket pricing

## 1 Research question and current finding

I develop coherent derivative-pricing methods under an established classical no-arbitrage GBM model. The [main paper](generalisation_proofs.tex) contains the essential proofs and focused appendices. The [supporting index](supporting_proofs_index.md) maps full derivations, exact evidence and earlier investigations.

The payoff-boundary method removes a strike auxiliary and supplies a positive second correction for convex payoffs. The residual-specific Gaussian transport certificate keeps dispersion moments and certifies the actual rounded financial grid through scalar transforms, without a continuous basket-price answer or Cartesian enumeration. It completes the formerly missing four- and eight-asset financial transfers at unchanged numerical widths.

The matched comparison now covers all 36 original capped and ordinary-call cases. Shared low-order continuous bounds resolve 18 without quantum estimation. In the other 18 cases, second-correction methods improve complete logical T caps against the informed first correction in ten, by approximately 13.9 to 92.6 per cent. The boundary method beats raw direct second-correction evaluation by about 3.18 per cent for the eight-asset ordinary call at common scale 0.05; raw evaluation wins the other active second-correction comparisons.

Strong classical methods certify all 36 model prices. No physical quantum pricing experiment was executed. The work contains checked mathematics and original conditional computational evidence; a contribution's novelty in the literature remains unestablished. The [source comparison](certificate_source_comparison.md) identifies an exact published antecedent for the positive cell identity. The residual derivative closure and its measured implementation value remain the narrower contribution candidate; the [novelty review](novelty_scope_review.md) states the unresolved priority gaps.

The [shared threshold derivation](shared_threshold_proof.tex) develops a different positive first-correction preparation. The fixed shared threshold preparation remains cheaper after I grant the numerical rival 64 independent price-table precision allocations. The additional 48-circuit raw-payoff pool receives a separately proved smoothing-first certificate. Against the cheaper admitted numerical or smoothing-certified common-normal raw-payoff rival in each case, the fixed shared-threshold encoding reduces complete logical T caps by 10.26, 5.26, 4.99 and 2.40 per cent after the independent-table precision challenge. These are finite-policy sufficient implementation caps under synthetic contracts, not global optima or physical quantum advantage. All canonical phase resolutions and 25 confidence repetitions remain unchanged. This is a checked preparation-cost improvement, with full financial errors and synthetic evidence; it establishes neither global optimality nor practical advantage. The full [adversarial review](reviews/research_adversarial_review_2026_10_10.md) explains the strengthened baseline and the remaining contribution limits.

## 2 Financial and mathematical scope

The frozen family comprises 27 capped calls and nine ordinary calls across nine two-, four- and eight-asset models. The 144 base arithmetic graphs, all 36 complete financial cases and all original failures remain retained. The transport theorem applies to positive-weight exponential baskets with bounded conditional curvature derivatives and certified Gaussian cells; its numerical implementation is not a theorem for arbitrary market dynamics.

The first-correction competitor receives the same coefficients, derivative norms, moments and confidence, including a coefficient-informed probability interval. Both coherent lookup architectures and the unchanged integer emitter remain admitted. An exact simplification pass has zero effect on the 204 higher-dimensional production graph structures.

Inputs are reproducible synthetic model parameters. Free sources and existing authorised LSEG or institutional access remain the permitted acquisition routes. No new purchase, OptionMetrics sample or additional private dataset is assumed.

## 3 Evidence and reproduction

The full pool has 502 candidate programs and 141,676 recorded allocation decisions. Separate audits check Gaussian transforms, covariance-preserving coefficients, finite moments, whole-domain scalar contracts, actual integer outputs, clean inverse circuits, complete confidence costs and every selection tie. These scopes are distinct from quantum-device execution.

Using the existing Python-Flint environment, the current construction and checks are:

```sh
python higherdim_preparation.py
python tensor_transfer.py
python tensor_first_transfer.py
python higherdim_precision_selection.py
python complete_information_comparison.py
python checks/audit_higherdim_preparation.py
python checks/audit_tensor_transport_transforms.py
python checks/audit_tensor_transport_weight_integrals.py
python checks/audit_tensor_curvature_envelopes.py
python checks/audit_higherdim_precision.py
python checks/audit_complete_information.py
```

The expensive precision search constructs and costs a finite pool. Its runtime is research setup, not a single-price execution time. All original sources and results are preserved; the detailed audit index lists further circuit, classical and packaging checks.

## 4 Lossless archive packaging

The canonical higher-dimensional comparison is larger than GitHub's per-file limit. I retain its original compressed bytes as four ordered parts of at most 32 MiB, with a [manifest](results/higherdim_precision_comparison.json.gz.manifest.json). The portable reader verifies every part and the full compressed SHA-256 before reading it; no expanded archive is required for the auditors.

```sh
python artifact_archive.py verify results/higherdim_precision_comparison.json.gz
python artifact_archive.py count results/higherdim_precision_comparison.json.gz --key records
python checks/audit_artifact_archive.py
```

The canonical file name is ignored locally. The manifest and parts are the retained repository representation. Reassembly is optional, using the helper's `reassemble` command with an explicit output path. No proof, case or numerical precision is discarded by packaging.

## 5 Repository state

The repository remains private. No commit is created in this phase; HEAD remains 1a25d67. The preceding manuscript and revised documentation are retained verbatim in history/before_tensor_transfer. The [research design](research_design.md) fixes the information tiers, benchmark scope and unresolved contribution questions.

## 6 Matched transport-order and refinement comparison

I compare two certificates for the same second pricing residual, rather than confusing transport order with the first and second pricing corrections. The published weighted cell identity of Cerone (2002) exactly covers the positive second-order representation. The common-grid FTC adaptation shares all boundary norms, residual dispersion moments, tails and rounded nodes; neither method enumerates Cartesian prices. This source finding excludes a broad original claim for the cell identity or its transform.

At each of five block levels, second order passes 33 of 36 allocations and FTC passes 32. At the finest level, second-order allowances are about 25.7–72.3 per cent smaller. Refining the transported coordinate from one block to 1,024 improves second-order allowances by only 0.0702–0.3480 per cent and resolves no additional failures. These results favour the aggregate certificate for this cohort. They do not establish a new refinement theory, additional quantum resource savings or a general ordering of transport methods. The protocol, timing scope and additional norm-input hashes are in [results/block_certificate_protocol.json](results/block_certificate_protocol.json) and [results/block_certificate_provenance.json](results/block_certificate_provenance.json).

The [fundamental archive](results/block_certificate_fundamental.json.gz) retains every raw geometry/mass record and all 180 aggregate outcomes. Only duplicated payoff block lists are omitted; the frozen helper recovers their exact JSON values and checks 840 semantic hashes. This is semantic reconstruction, not recovery of the original compressed bytes. The [manifest](results/block_certificate_reconstruction_manifest.json) links the original archive checksum, and the unchanged summary keeps its original canonical filename and checksum. No giant expanded canonical archive is required in the repository.

```sh
python block_certificate_archive.py verify --fundamental results/block_certificate_fundamental.json.gz
python checks/audit_block_transport.py
```

Optional expansion uses the same verify command with `--expand` and an explicit local output path. The producer `block_certificate_comparison.py` reproduces the canonical experiment; the existing original result is outside the repository and its filename is ignored. Detailed proofs and the independent reconstruction audit are indexed in [supporting_proofs_index.md](supporting_proofs_index.md). The preceding manuscript and documentation are retained verbatim in `history/before_certificate_comparison`. No commit is created in this phase.

## 7 Residual-specific contribution and classical pricing audit

I have completed the next source and information test. The final Curran journal scan is now read, alongside the published polynomial-pricing analysis and an explicitly identified accepted manuscript of Caldana et al. Conditional Gaussian pricing, finite moment integration and generic positive-measure closure are established mechanisms. The [priority review](residual_priority_review.md) distinguishes them from the project's exact residual derivative calculation and tested finite arithmetic.

A new bounded classical calculation through degree three uses the same model and M2/M3/M4 constants as transport, with additional full-Gaussian moments and shifted kernel evaluations. It certifies 33 of 36 complete cent prices: eighteen stop at degree one, nine at degree two and six at degree three. All 24 higher-dimensional cases are priced before a grid or estimator is prepared. The three failures match the aggregate transfer failures and are already resolved by the retained validated classical variants. The [summary](results/residual_cost_boundary_summary.json) retains all work counts and timing exclusions.

The independently checked published geometric lower-bound/AG-upper implementation solves all 36 strike optimisations, including mixed-sign conditional covariances. Its derived separate-strike cap enclosure certifies no cent price here. This is a bounded negative finding for that particular implementation, not a defeat of every published variant. All new producer and independent audit sources are indexed; earlier quantum caps and numerical evidence remain unchanged. The current cohort establishes useful financial-oracle certification and conditional implementation results, but no additional pricing capability or confirmed original theoretical mechanism.

The complete abstract reduction, counterexamples, comparator derivation and consulted-source ledger are kept in supporting files. The main paper contains the essential reasoning and a compact table showing the classical work left by each information tier. No new commit is created.

## 8 Controlled dispersion stress and Figure 1

The new [fixed panel](results/dispersion_separation_protocol.json) extends the nine original GBM models by exact outer-loading scales 1, 1.5, 2 and 3, recomputing all risk-neutral drifts and covariances. Its 144 contracts retain the original common Gaussian scales, weights, payoffs, numerical grid and error targets. Original inputs remain unchanged. The [complete results](results/dispersion_separation_results.json) retain every interval, tail, rounded node, geometry moment, operation count and failure.

The aggregate transfer allowance passes in 81 cases, while strengthened degree-three price intervals pass in 71. Degree five resolves six of the ten provisional separating cases and leaves four. These are all capped calls with H=1. The result identifies transfer-certificate usefulness beyond the tested finite Taylor price screen; it does not establish complete quantum prices or superiority over all classical pricing methods. No original logical resource cap is carried over to a stressed financial support.

The classical opponent is strengthened by exact intersections of all available lower-order intervals and intrinsic payoff ranges. This policy was added after the first registered run began. Both earlier executed versions, including the outward-interval storage correction, are retained in phase history. The protocol records that chronology explicitly rather than describing the final opponent policy as a public or pre-outcome registration.

Figure 1 now presents the audited numerical relationship. Its fixed representative ordinary-call model has no separating sampled point, and the other panel reports every provisional and surviving separator by scale. The earlier boxed text flow repeated the neighbouring theorem and is removed. The [PNG](figures/dispersion_separation.png), [SVG](figures/dispersion_separation.svg), exact figure data and coordinate audit are retained; the main source embeds a self-contained vector figure for the native compiler. Previous paper and documentation versions are preserved in `history/before_dispersion_separation`.

## 9 Licensed-data publication and retention

I apply the university guidance supplied on 9 October 2026 through [DATA_POLICY.md](DATA_POLICY.md). Raw LSEG observations, reconstructible transformations and systematic derived datasets stay out of GitHub and third-party services, including private repositories. Code may be shared without data or credentials; limited academic summaries and figures must satisfy the stated insubstantial, non-commercial and non-systematic conditions and cite the actual data source. I do not systematically bulk download data. Underlying licensed working data are deleted when their project is complete.

The current financial inputs, figure data, simulation traces and quantum resource results are synthetic, with their producer lineage retained. Public Sobol material has its separate source licence. These outputs are not LSEG exports, and I do not add a misleading LSEG attribution. A future licensed-data analysis is reproducible through code and separately entitled local execution, without redistributing observations or a systematic derived panel.

## 10 Complete prices for the four stress candidates

The interrupted temporary pricing stage was unavailable on recovery. I regenerated the evidence with a fixed protocol before fresh execution, retaining new random traces and clearly recording that earlier outcomes were already known. The [complete classical calculation](results/four_case_classical_pricing.json) now certifies all four synthetic cap-one prices. One uses exact summation of all 1,024 nodes; three use fixed 100,000-draw bounded-score intervals with full financial transfer and tails. Their coverage assumes independent uniform random words. All sampled arithmetic is independently replayed.

The [fresh quantum comparison](results/four_survivor_quantum_resources.json) supplies complete logical plans for all sixteen method targets. Positive R1 has the smallest declared cap in every case, approximately 25.3–216.1 billion T gates. The eight-asset boundary all-12 baseline fails its monetary budget; the declared all-14 repair passes and is fully costed. This is a conditional logical-resource comparison, with no quantum hardware execution or practical advantage. It does not leave the four prices unavailable to classical methods.

Every new model, trace and table archive is generated synthetic evidence, not LSEG-derived data. [DATA_POLICY.md](DATA_POLICY.md) governs any future licensed-data work. The [supporting index](supporting_proofs_index.md) links complete proofs, independent checks, exact inputs, reconstruction records and timing scopes. The previous paper and revised project documents are preserved in `history/before_four_case_completion_and_data_policy`.

```sh
python checks/audit_four_case_classical.py
python checks/audit_four_survivor_quantum_resources.py
python checks/reproduce_four_survivor_quantum.py --output /absolute/path/to/new-empty-output --prepare-only
```

The quantum reproduction wrapper copies only the frozen source/input set into a new directory before registration and execution. Omitting `--prepare-only` repeats the synthetic construction there. Reproduction produces new timing metadata; it does not overwrite the retained experiment. Supplemental exact-arithmetic reviews can be replayed through `checks/replay_four_case_stage_reviews.py` with an explicit output directory.


## 11 Bounded Chebyshev mean-estimation development

I implement a finite reflection-only adaptation of the published Rall–Fuller Chebyshev sampling framework on the four selected positive first-correction flags. Both polynomial parities have explicit clean circuits and distinct success events. Every fresh 100-shot confidence look spends its own exact-binomial allowance; the controller has twelve looks, degree at most 2,047 and a pilot T budget of one tenth of its separately costed canonical fallback. The fallback retains the original financial prior. No sampled classical price answer enters either policy.

The [registered diagnostic results](results/chebyshev_resources.json) retain all 400 declared prior-mesh histories and their random words. The pilot certifies 111; unsuccessful paths include their full work and fallback cap. All four diagnostic medians exceed the earlier QPE cap. Equal-weight mesh means decrease in the first three cases and increase in the eight-asset case. These are descriptive synthetic histories, with no device execution, returned quantum price or expected-cost claim for the unknown financial flag. The [exact interpretation](results/chebyshev_interpretation.json) checks every charged-history average and the resource break-even identity.

The complete policy has failure allowance at most 1/3600 per price and worst-history T count at most 1.1 times its fresh half-allowance fallback, which is at most 1.188 times the earlier full-allowance QPE cap here. The [full derivation](chebyshev_estimator_proof.tex), [mathematical review](reviews/chebyshev_mathematical_review.md), source-version record and independent resource, confidence, wire-binding and break-even checks are indexed. This develops the existing quantum pricing method; it does not establish a new generic QSP theorem, defeat unchanged published ChebAE, implement an upstream financial QSP kernel or demonstrate practical advantage.

The main paper contains the essential circuit and confidence argument, a compact table distinguishing mean and median histories, and the exact break-even condition. Earlier paper and documentation versions are preserved in `history/before_bounded_chebyshev`. All new evidence is synthetic and follows DATA_POLICY.md.

```sh
python checks/check_chebyshev_controller.py --output /absolute/path/to/new-controller-check.json
python checks/audit_chebyshev_resources.py --repo /absolute/path/to/project --results /absolute/path/to/copied-results
python checks/audit_chebyshev_math.py --output /absolute/path/to/new-math-audit.json
python checks/audit_chebyshev_break_even.py
```

Fresh diagnostics can be generated in a separate new output directory with `chebyshev_resources.py --repo /absolute/path/to/project --output /absolute/path/to/new-output register`, followed by the same arguments and `run`. The existing frozen result is not overwritten. The registered policy is a bounded adaptation with fixed fresh batches, not the unchanged published early/late controller.

## 14 Exact bases, shared threshold encoding and stronger rival

I screen all 58 exact rational factor bases and rebuild complete plans for seven selected representations. The continuous no-arbitrage price is unchanged, while finite-law moments, supports, priors and arithmetic are rebuilt. The basis-only gains are modest and do not change QPE resolution. The new event encoding instead replaces two conditional-price tables with a common threshold and a clean interval predicate. Its cap normaliser remains one; the added inner normal law, strike midpoint and endpoint arithmetic are fully charged. Twelve failed numerical precision configurations remain visible.

The [colour precision map](figures/pricing_precision.png) shows every admitted and rejected rival choice. The [figure review](reviews/figure_selection_review.md) explains its purpose, the separate role of the dispersion figure and the rejection of an interpolated 3D surface or statistical cost bars. Exact figure data and self-contained native vector sources are retained. The main paper embeds them directly, so its native compilation needs no external image file.

The portable new producers require explicit repository and fresh output paths. Run registration before execution, preserving the historical protocols. The factor resource producer records the Python-source set at registration; additions to that set require a fresh reproduction registration or a directory containing its exact registered dependencies, rather than modifying a historical protocol.

```sh
python factor_basis_certificate.py --repo . --results /tmp/qbp-basis-new register
python factor_basis_certificate.py --repo . --results /tmp/qbp-basis-new run
python shared_threshold_oracle.py register --repo . --factor-resources results/factor_basis_quantum_resources.json --output /tmp/qbp-event-new
python shared_threshold_oracle.py run --repo . --factor-resources results/factor_basis_quantum_resources.json --output /tmp/qbp-event-new
python strengthened_r1_baseline.py register --repo . --factor-resources results/factor_basis_quantum_resources.json --threshold-resources results/shared_threshold_resources.json --output /tmp/qbp-rival-new
python strengthened_r1_baseline.py run --repo . --factor-resources results/factor_basis_quantum_resources.json --threshold-resources results/shared_threshold_resources.json --output /tmp/qbp-rival-new
python checks/replay_pricing_evidence.py --repo . --output /tmp/qbp-pricing-audit-new
```

The replay helper preserves repository evidence and writes all new audit output to its explicit isolated directory. Optional figure rendering uses the free dependencies in `requirements-figures.txt`. All inputs and new outputs are synthetic. No LSEG data, purchase, commit or index change is involved. These restricted finite-pool comparisons do not establish literature priority, a uniformly best encoder or a physical speed-up.

## 15 Common-normal rival and smoothing-first certificate

I cost the raw capped-payoff difference on a shared normal label, using exact monotone clamps and the rational exponential range bound 64/63. I preserve the initial broad-headroom and corrected pointwise studies. Independent review then derives a tighter mean-error certificate for the same circuits by integrating the continuous Gaussian before charging outer endpoint errors. All 48 financial budgets and 96 QPE/full-charge plans are rebuilt under a separately registered amendment. The eight-asset raw phase resolution falls from 1,024 to 512, and one four-asset winner changes its inner normal width. All winning exponential widths remain at the finite pool's upper boundary.

The [full supporting proof](raw_payoff_proof.tex) contains the clamp, CDF, variation, shifted-register, arithmetic, smoothing and monetary derivations. The [amended results](results/raw_payoff_smoothing/raw_payoff_smoothing_resources.json) and separate mathematical/resource reviews retain exact lineage. Figure 2's Event/min footers refer to the numerical pool; the main table uses the cheaper admitted rival across both pools.

The new sources accept explicit paths and preserve existing output. A fresh smoothing reproduction is:

```sh
python raw_payoff_smoothing_amendment.py register --repo . --raw-results results/raw_payoff_tight --proof proof_sources/raw_payoff_smoothing_first_certificate_fragment.tex --independent-review results/raw_payoff_independent_resource_review.json --output /tmp/qbp-raw-smoothing-new
python raw_payoff_smoothing_amendment.py run --repo . --raw-results results/raw_payoff_tight --proof proof_sources/raw_payoff_smoothing_first_certificate_fragment.tex --independent-review results/raw_payoff_independent_resource_review.json --output /tmp/qbp-raw-smoothing-new
python raw_payoff_smoothing_amendment.py audit --repo . --raw-results results/raw_payoff_tight --proof proof_sources/raw_payoff_smoothing_first_certificate_fragment.tex --independent-review results/raw_payoff_independent_resource_review.json --output /tmp/qbp-raw-smoothing-new
```

All evidence is synthetic. No LSEG observations, purchase, repository staging, commit or push is involved.

## 16 Wider precision test and paper structure

The separately registered extension adds exponential addresses eleven through fourteen at the same three normal widths, producing a combined 96-circuit raw pool. All four winners remain at exponential width ten, now inside the tested range. Independent rational-series scalar, whole-circuit and complete-charge audits preserve the 2.49 to 10.67 per cent event saving. A derived literal-selector charge floor also excludes feasible common exponential widths fifteen through twenty-two under the fixed graph, word precision, normal widths, unary lookup and canonical controller. Omitted lower widths, unequal payoff-node precisions and other arithmetic/estimators remain outside this exclusion. This is not a global quantum optimality theorem.

The main abstract has 159 words under the retained counting convention. Seven full proofs and seven standard derivation blocks are retained in its appendices; the main text keeps their structural sketches, financial assumptions, resource implications and adverse classical comparisons. Displayed equations are numbered, short asides use footnotes, JEL codes precede full-width naturally wrapped keywords, and the moment/derivative-bound notation is distinct.

## 17 Independent payoff-table precision

I separately register 300 independent PG/PA precision settings, reusing 48 exact diagonal records and emitting 252 fresh graphs. Every selected minimum is symmetric nine/nine, with normal widths ten, ten, eight and ten. The improvement is due to the previously omitted diagonal nine; asymmetric allocation supplies no additional winning gain. The strongest-rival shared-threshold savings are now 10.26, 5.26, 4.99 and 2.40 per cent. Distinct approximate tables can reverse represented ordering, so the final projection is retained and rigorously charged. The exact scalar witness and clean negative-branch replay are recorded without claiming that the witness occurs in the financial grids.

The conserved certificate uses the preceding 2E_m policy. A sharper common-shift result is prospective and uncosted. The extended literal-selector proof excludes either node at addresses fifteen through twenty-two for any feasible other node, under the retained compiler/controller. Omitted lower widths, unequal thirteen/fourteen, other word/normal precisions, shared arithmetic, other lookups and modern estimators remain unresolved. No global optimality, literature priority or practical advantage is inferred.

## 18 Common information and matched estimator transfer

I apply a common-information certificate to all 348 raw, 64 numerical and four threshold configurations. All 416 bounds tighten; twelve numerical budgets remain rejected and no canonical minimum changes. The raw common-scale term is 101E_s rather than 202E_s and the common fixed shift is charged once. A full positive-measure proof separates those operations from the finite atomic transfer.

Complete threshold preparation/reflection costs satisfy two affine inequalities against both selected rivals, yielding a sharp matched-degree comparison for every positive integer degree under the retained additive templates. Identical adaptive histories additionally require matched implemented probabilities and controller information.

The registered modern comparison includes all three selected preparations, 400 paired synthetic mesh/seed histories and 1,200 alternative histories. Numerical pricing has lower descriptive means in every case, and every median exceeds its own canonical plan. Rebuilding all 808 admitted half-alpha fallbacks across the full pool preserves phase size/Fourier charge and changes 25 repetitions to 27. Protected-cap minimisation is therefore a monotone transformation of canonical minimisation. The positive result is estimator-portable preparation saving and preserved minimum protected-cap ordering, not extra demonstrated gate reduction, adaptive-mean superiority, theoretical priority or physical advantage.

The main paper retains structural steps and interpretation; full new derivations are in raw_payoff_proof.tex and the main appendix. Recorded results are under results/common_shift and results/matched_chebyshev. A fresh portable matched run registers matched_chebyshev.py with explicit --repo, --output and --inputs results/common_shift/amended_inputs.json, then runs with the same repo/output. Original registered protocols retain their original path/version records. No licensed observations are used.
