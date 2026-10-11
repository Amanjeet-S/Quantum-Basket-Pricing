# Supporting proofs and retained evidence

## 1. Main-paper scope and source preservation

The condensed paper presents the current finite-correction quantum methods, its matched direct second-correction and positive first-correction opponents, and their complete financial accuracy contract. I retain the detailed derivations and earlier negative studies as supporting evidence. Condensation changes the presentation; it does not remove a failed configuration, weaken a proof premise or replace the original financial law.

The complete pre-condensation source is preserved verbatim in [history/full_derivations_before_condensation/generalisation_proofs.tex.txt](history/full_derivations_before_condensation/generalisation_proofs.tex.txt), with SHA-256 `d822fa3b91f6dcffdb9075a509d31d786f7fd2a84f74e8bfe0920d7b7db84774`. The current main source remains [generalisation_proofs.tex](generalisation_proofs.tex). Its focused appendices retain the necessary proof details; the files below preserve the extended derivations and earlier experiments.

## 2. Mathematical derivations

| Main-paper object | Detailed source | Retained mathematical scope |
| --- | --- | --- |
| Classical GBM, equivalent pricing probability and admissible-strategy no-arbitrage argument | `integral_remainder_quantum_proof.tex`, subsection “The classical no-arbitrage target”; pre-condensation source | Full volatility row rank, constant physical drift and dividend yields, Brownian filtration, equivalent density, discounted gain martingales and the stipulated price functional. |
| Exact Gaussian tilt coefficients and covariance shift | Pre-condensation source, sections “A finite-order deterministic certificate for Gaussian baskets” and “Finite Gaussian evaluation of every coefficient and remainder moment”; `taylor_pricing.py` | Multinomial expansion, tilted Gaussian convolution, retained `R_G` shift, finite moments and two-sided finite-order price intervals. |
| Combined Gaussian derivative norm and sharpness over positive threshold measures | `joint_norm_proof.tex`; `gaussian_derivative_norm.py`; `norm_literature_comparison.md` | Signed derivative polynomial, positive-measure norm bound, point-mass attainment, complete rational root isolation and outward envelopes. Operator sharpness is distinguished from capped-call remainder attainment. |
| Integral Taylor correction and exact Beta law | `integral_remainder_quantum_proof.tex`; `beta_minimum.py` | Integral remainder, positive ratio, bounded signed density polynomial, clean reversible integer minimum and continuous zero extension. |
| Full realised financial arithmetic and absolute error propagation | `integral_value_proof.tex`; `integral_value_certificate.py`; `integral_value_oracle.py` | Actual dyadic financial coordinates, finite vertex normaliser, decreasing/nonmonotone interpolation, integer guards, signed final clamp, probability endpoints and the monetary budget. |
| Quadratic auxiliary error | `integral_auxiliary_quadrature_proof.tex`; `integral_auxiliary_quadrature.py` | Exact Beta cell masses, first and second centred cell moments, quadratic Beta and strike bounds, same finite financial law and optional `M4 <= U_star M3`. |
| Exact product-grid moment priors and signed-part identity | `integral_value_proof.tex`; `integral_value_certificate.py`; `integral_value_oracle.py` | Weighted finite exponential sums, actual label multiplicities, zero-order weighted moment, correlated positive first-remainder bound and exact integer signed-part recombination. |
| Canonical estimator, confidence and Fourier synthesis | Pre-condensation source, sections “Certified Clifford and T synthesis of the inverse Fourier transform” and “A complete canonical phase-estimation comparator”; `integral_estimation_extension_proof.tex`; `integral_pricing_budget.py` | Controlled signs, phase law, operator-norm allowance, exact binomial confidence, omitted distant phases, complete preparations and resets. |
| Financial tail and grid transfer | `integral_value_proof.tex`; `continuous_quadrature.py`; `integral_value_certificate.py` | Gaussian tilt tails, physical-coordinate quantile coupling, density-weighted midpoint error and paid reference-based one-dimensional validation. Grid-only bounds are not complete price certificates. |

The main paper contains the essential proofs. These files provide expanded derivations and implementation-specific details, not additional claims of theoretical priority. Source versions and checked assumptions are retained in `SOURCES.md`, `norm_literature_comparison.md`, `quantum_method_sources.md` and `qsp_source_ledger.md`.

## 3. Finite certificates and verification records

| Evidence | File or record | Interpretation |
| --- | --- | --- |
| Fixed-precision complete financial graphs | `results/integral_graph_certificate_assessment.json` | All 81 method/model/cap configurations retain node ranges and arithmetic errors. Passing guards alone does not establish a financial tolerance. |
| Full outward scalar and quantile nodes | `results/integral_scalar_table_archive.json.gz` | Deterministic compressed archive keyed by canonical certificate hashes. The prototype archive contains 60 certificates and 242,747 nodes. |
| Tight actual-input domains | `results/integral_finite_input_domains.json` | Tensor-box vertices cover ideal quantile enclosures and actual dyadic coordinates; their normalisation scope is finite. |
| Product-grid moments | `results/integral_finite_grid_moments.json` | Exact finite prepared-law moments from one-dimensional sums; no continuous financial price is required. |
| Competent positive first-remainder priors | `results/integral_positive_r1_priors.json` | Level-ratio correlation is retained in the underlying-price gap. |
| Paid continuous-to-finite financial transfer | `results/integral_financial_grid_assessment.json` | Nine two-asset reference comparisons include tails, financial grid and actual node rounding. The continuous reference already resolves the price. |
| Higher-dimensional sufficient-width diagnostics | `results/integral_grid_width_envelopes.json` | Grid-only envelopes retain impractical widths and failures; they do not imply complete accepted quantum prices. |
| Matched method protocol | `results/integral_oracle_protocol.json` | Prespecified monetary target, confidence, architecture, support, method controls and zero-quantum screening. |
| Full production evaluator profiles | `results/integral_production_profiles.json` | Logical emitted-body counts and clean controlled composition. The final resource study must link its selected precisions to their financial budgets and independent checks. |
| Estimator extension | `results/integral_extended_estimator.json` | Actual compiled Fourier profiles through 64 phase bits and their error allowances; the 32-bit profile is preserved. |

Independent checks inspect alternate scalar formulas, inverse-normal CDF bracketing, table magnitudes, clean reversible workspaces, phase conventions and probability endpoints. The installed study should link the exact final audit filenames and source hashes in its manifest. Prototype circuit replay and large-width count composition have distinct scopes; neither is reported as a device experiment.

## 4. Earlier investigations retained outside the condensed argument

| Investigation | Retained source | Reason for retaining it |
| --- | --- | --- |
| Signed Gaussian loading and its positive-mixture opponent | `signed_gaussian_loading_proof.tex`; `coherent_tilt_derivation.tex`; `positive_envelope_comparator.md` | Full continuous and finite identities, success probabilities, precision amplification, coherent state-specific loading and the adverse complete loader comparison. |
| Quantum signal-processing applicability | `qsp_applicability_proof.tex`; `qsp_source_ledger.md` | Exact antecedent assumptions and the limits of a source comparison; an unimplemented compatible QSP opponent is not declared defeated. |
| Earlier numerical-value and dyadic controller studies | Pre-condensation source; `variance_oracle_comparison.py`; `complete_estimation_study.py` and associated results | Prior finite-target resource evidence and rejected comparison hypotheses remain reproducible. They are historical supporting studies, not evidence for the present integral oracle. |
| Randomised finite-grid and higher-dimensional studies | Pre-condensation source; `rqmc_work_comparison.py`; `higher_dimension_comparison.py` and associated results | Retained stochastic confidence calculations, realised financial-output replay and adverse evaluator-work comparisons. |

These records remain available for checking the development and its boundaries. They are not inserted into the condensed main paper merely to increase its length or apparent scope.

## 5. Final independent precision and resource evidence

The complete finite-pool comparison is results/integral_precision_comparison.json.gz, with 97 programs and 11,872 recorded decisions. Its selected local bit maps and complete monetary/confidence budgets are independently checked by checks/audit_precision_selector.py and checks/audit_integral_precision_budget.py. The scalar ordinate subset archive is results/integral_precision_scalar_archive.jsonl.gz. Full-coherent fixture, production profile, scalar, Fourier and auxiliary audits have separate records; none is an execution of quantum pricing hardware.


## 6 Current curvature-correction argument

The current main paper develops the strike-integrated boundary method, stable Gauss-Jacobi opponent, ordinary-call positive correction and price-independent financial grid. The preceding concise main paper and its documentation are preserved verbatim in history/before_boundary_correction. Sections1–5 of this index retain historical evidence and its original limitations; they are not the current resource ranking.

| Mathematical object | Detailed supporting source | Evidence |
| --- | --- | --- |
| Sparse piecewise-affine payoff curvature, global capped M1/M2 and higher derivative norms | boundary_derivative_proof.tex; boundary_derivative_certificate.py | results/boundary_derivative_norms.json; results/boundary_derivative_audit.json |
| Positive/signed curvature preparation, stable two-node opponent and dyadic constants | boundary_value_proof.tex; boundary_value_oracle.py; boundary_value_certificate.py | results/boundary_production_profiles.json.gz; results/boundary_constructor_assessment.json.gz |
| Actual-law moments through degrees4 and6; paid auxiliary allowances | boundary_finite_moments.py; boundary_auxiliary_certificate.py | results/boundary_finite_moments.json; results/boundary_auxiliary_independent_audit.json |
| Signed Gaussian cell transport, stable residual derivatives, tails and coefficient enclosures | boundary_grid_transfer_proof.tex; boundary_grid_certificate.py | results/boundary_financial_grid_transfer.json.gz; results/boundary_financial_grid_transfer_summary.json; results/boundary_grid_independent_audit.json |
| Full frozen precision pool and complete monetary/confidence caps | boundary_precision_selection.py | results/boundary_precision_comparison.json.gz; results/boundary_precision_scalar_archive.jsonl.gz; results/boundary_precision_provenance.json |
| Continuous convex screening before grids or sampling | boundary_convex_screening.py; main-paper convex-screen corollary | results/boundary_convex_screening.json; results/boundary_convex_screening_audit.json |
| Strong price-independent classical GJ/raw sums | boundary_classical_pricing.py | results/boundary_classical_pricing.json.gz; results/boundary_classical_pricing_audit.json |

The supporting proof fragments are read with the main definitions and its author–date bibliography. Full cell traces and scalar ordinate archives preserve every enclosure without expanding the main paper. The financial family comprises36 cases,144 accepted arithmetic graphs and12 complete two-asset pricing comparisons. Unsupported higher-dimensional cases remain explicit.

Independent repository-relative checks are checks/audit_boundary_derivative_certificate.py, checks/audit_boundary_value_oracle.py, checks/audit_boundary_production_profiles.py, checks/audit_boundary_grid.py, checks/audit_boundary_auxiliary.py and checks/audit_boundary_classical_pricing.py. The final precision/monetary audit passes all 173 candidates, 730 retained profiles, 23,772 decisions and 570 full QPE flag plans. Its repository-relative script and final record are linked below. The continuous convex-screen auditor is checks/audit_boundary_convex_screening.py. Earlier audit snapshots remain retained. Source versions and attribution are recorded in boundary_sources.md.

The final reproducible allocation and monetary check is [checks/audit_boundary_precision.py](checks/audit_boundary_precision.py), with [results/boundary_precision_independent_audit.json](results/boundary_precision_independent_audit.json). It checks retained profile aggregates and complete charges; literal fresh production-body expansion is separately recorded in results/boundary_production_profiles_audit.json. The original audit source and record are preserved as snapshots.

Frozen audit-version filenames and their exact hash resolution are recorded in [results/boundary_audit_linkage.json](results/boundary_audit_linkage.json). The original grid audit identified by the financial summary is retained as boundary_grid_original_audit.json; the portable replay is separately retained as boundary_grid_independent_audit.json. Numerical construction inputs and the executed comparison remain unchanged. The final source-preservation and compilation record is [results/boundary_project_validation.json](results/boundary_project_validation.json).

## 7. Completed tensor financial transfer and informed comparison

The current main argument completes all 36 frozen financial cases. Earlier sections of this index preserve the historical studies and their original limitations; the resource rankings are superseded by the fully informed comparison below. The preceding paper and updated documentation are retained verbatim in history/before_tensor_transfer.

| Mathematical or numerical object | Detailed supporting source | Evidence |
| --- | --- | --- |
| Actual rounded Gaussian nodes, factorial empirical-MGF series and moments through degree six | higherdim_preparation.py; higherdim_preparation_proof.tex | results/higherdim_preparation.json; results/higherdim_preparation_audit.json; results/higherdim_preparation_costs.json |
| Positive Gaussian transport, exact exponential transforms, tensor hybrid masses, residual derivative envelopes and tails | tensor_transfer.py; tensor_transfer_proof.tex; tensor_transport_independent_review.tex | results/tensor_transfer.json; results/tensor_transport_transform_audit.json; results/tensor_transport_weight_audit.json; results/tensor_curvature_envelope_audit.json |
| Matched first-correction transfer and covariance-preserving forward tilt | tensor_first_transfer.py; tensor_transfer_proof.tex | results/tensor_first_transfer.json; results/tensor_parameter_provenance_audit.json |
| Complete higher-dimensional greedy pool and charged estimator plans | higherdim_precision_selection.py | results/higherdim_precision_comparison.json.gz.manifest.json and ordered parts; results/higherdim_precision_independent_audit.json; results/higherdim_precision_scalar_archive.jsonl.gz |
| Exact integer simplification and its zero production effect | boundary_integer_optimisation.py; boundary_integer_optimisation_proof.tex | results/boundary_integer_optimisation_audit.json; results/boundary_integer_constant_table_audit.json; results/boundary_integer_optimisation_production_structure.json |
| Strong continuous ordinary-call Taylor opponent | higherdim_classical_taylor.py; higherdim_classical_taylor_proof.tex | results/higherdim_classical_taylor.json; results/higherdim_classical_taylor_audit.json |
| Coefficient-informed R1 interval, shared low-order screen and final 36-case selection | correction_information_priors.py; correction_information_proof.tex; complete_information_comparison.py | results/complete_information_comparison.json.gz; results/complete_information_independent_audit.json |
| Lossless archive access and byte-exact reconstruction | artifact_archive.py | results/higherdim_precision_packaging_audit.json; results/artifact_archive_audit.json; results/higherdim_precision_auditor_archive_access.json |

The final information audit checks 113 first-correction intervals, 2,562 complete QPE flag plans, 988 analytical midpoint plans and all 36 Gaussian coefficient sets. The higher-dimensional pool audit checks 329 candidates, 117,904 decisions and 324,687 fresh scalar ordinates. It replays retained profile accounting; literal small-width clean-circuit checks and the established primitive proofs have a separate scope. Neither is reported as a quantum-device experiment.

Repository-relative auditors include checks/audit_higherdim_preparation.py, checks/audit_tensor_transport_transforms.py, checks/audit_tensor_transport_weight_integrals.py, checks/audit_tensor_curvature_envelopes.py, checks/audit_tensor_parameter_provenance.py, checks/audit_higherdim_precision.py, checks/audit_complete_information.py, checks/audit_boundary_integer_optimisation.py and checks/audit_artifact_archive.py. The compressed canonical filename is a logical archive reference: the portable loader resolves verified repository shards when that oversized file is absent. Numerical data and its original compressed hash are unchanged.

The strongest candidate contribution, closest antecedents and unresolved originality comparison are documented in [novelty_scope_review.md](novelty_scope_review.md). Correctness, empirical logical-resource findings and literature priority are distinct claims.

The final installed-source, archive-portability, preservation and native compilation checks are recorded in [results/tensor_project_validation.json](results/tensor_project_validation.json). The prior phase sources and results remain unchanged except the eight snapshotted documentation and ignore files; no commit is created.

## 9 Matched transport certificates and source-level comparison

The current source mapping is [certificate_source_comparison.md](certificate_source_comparison.md), with exact consulted versions in [certificate_source_access.json](certificate_source_access.json). The Gaussian positive cell identity in the earlier transport derivation is a published specialisation of Cerone (2002). Its original source is retained for reproducibility; this attribution governs the current contribution claim. The remaining candidate concerns the residual-specific moment closure and its computational use.

| Object | Supporting source | Retained evidence and scope |
| --- | --- | --- |
| Weighted FTC adaptation, positive measures, exact Gaussian transforms, raw-mass hybrid factorisation and Cauchy refinement | `block_transport_review.tex` | Complete supporting derivations, established mechanisms explicitly attributed, cost and numerical-monotonicity qualifications. |
| Same-grid R2 certificates at five nested partitions | `block_certificate_comparison.py`; `results/block_certificate_protocol.json` | `results/block_certificate_summary.json`: 180 outcomes, all failures, same nodes/norms/tails and precisely scoped joint timing records. |
| Transform and geometry audit | `checks/audit_block_transport.py` | `results/block_transport_math_audit.json`: production transforms, positive interval quadrature, one-block agreement and independent small-grid polynomial moments. |
| Portable fundamental representation | `block_certificate_archive.py` | `results/block_certificate_fundamental.json.gz`; `results/block_certificate_reconstruction_manifest.json`: all raw evidence retained, omitted duplicated lists restored by exact semantic hashes; no original gzip-byte identity claim. |
| Independent reconstruction | `checks/audit_block_fundamental_archive.py` | `results/block_fundamental_independent_audit.json`: 157,476 block records, 2,047,188 components and 840 exact semantic matches. |
| Source and nested-input provenance | `results/block_certificate_provenance.json` | Original canonical checksum, frozen helper/producer/auditor hashes and additional boundary/Gaussian norm-input hashes. |

The joint timing covers geometry plus both certificate orders for twelve payoffs. It excludes common prefix preparation and serialisation, and higher levels reuse geometry formed at the first level. No order-specific cost reduction, cold-run speed ratio, new accepted pricing case or change to the previous quantum resource caps is inferred. The previous manuscript and revised documentation are preserved in `history/before_certificate_comparison`.

The preceding main-paper audit and setup paragraphs remain verbatim in the phase snapshot. Their records retain the 483 scalar contracts, 324,687 fresh ordinates, 46,441 program certificates, 2,828 profiles, all 117,904 precision decisions and 1,744 QPE plans; the final information audit additionally checks all 36 coefficient sets, nine continuous ratio moments, 113 priors, 2,562 QPE plans and 988 analytical midpoint plans. The MGF construction record retains 512 inverse-error-function evaluations, 30,720 empirical-power recurrences, 38,247 coefficient operations, 627 remainder exponentials, 3,227 multinomial/gap exponentials and 1,056 vertex exponentials. These full ledger counts are kept outside the concise main argument and do not change the new block-comparison timing scope.

## 10 Residual-specific priority, general reduction and classical pricing

| Object | Detailed source | Evidence and scope |
| --- | --- | --- |
| General positive-measure finite-transform reduction and exact degree-four derivation | `abstract_transfer_reduction.tex`; `independent_reduction_review.md` | Minimal product-measure and outside-tail assumptions; a reduction to established tools, not a new generic theory. |
| No uniform ordering of first/second transport | `checks/audit_transport_counterexample.py` | `results/transport_counterexample_audit.json`: explicit smooth Gaussian-cell counterexample and independent padded quadrature. |
| Bounded continuous Taylor price test | `residual_same_information_taylor.py`; `classical_price_bounds_proof.tex` | `results/residual_same_information_taylor.json`: degree cap three, moment cap four, every failed order and full work counts. |
| Independent Taylor reconstruction | `checks/audit_continuous_degree3.py` | `results/continuous_degree3_independent_audit.json`: direct polynomial expansion and 1,024-bit alternate Gaussian kernels; all 36 cases. |
| Published geometric bounds and certified mixed-sign optimisation | `caldana_conditioning_benchmark.py`; `classical_price_bounds_proof.tex` | `results/caldana_conditioning_benchmark.json`: all 36 strike and contract records, brackets, bisections, covariance signs and failed radii. |
| Independent geometric-bound check | `checks/audit_conditioning_bounds.py` | `results/conditioning_bounds_independent_audit.json`: original-covariance parameters, minima/roots, candidate values and final price propagation. |
| Same-model information and scope | `residual_priority_review.md`; `results/residual_cost_boundary_summary.json` | Continuous moments require additional queries; matching the grid-admission set is observed, not a universal equivalence. |
| Source versions and successful Safari access | `classical_antecedent_comparison.md`; `classical_antecedent_sources.json`; `residual_priority_access.json` | Final Curran scan distinguished from the accepted Caldana body and author-preprint supplements. |

Financial audit replay uses `--stage . --out` followed by an explicit local output filename. Supporting fragments use the main paper's definitions and bibliography and do not require a replacement main document. Existing producers/results remain unchanged; the preceding paper and revised documentation are preserved in `history/before_residual_priority`. The phase provenance links source hashes, original inputs and installed evidence.

## 11 Dispersion stress, retained history and scientific Figure 1

| Object | Detailed source | Retained scope |
| --- | --- | --- |
| Exact scaled model and pricing conditions | `dispersion_scaling_proof.tex`; `checks/audit_dispersion_models.py` | Recomputed rational drifts/covariances, constant asset forwards, positive-weight ratio and unchanged raw Gaussian law. |
| Complete fixed panel | `dispersion_separation_study.py`; `results/dispersion_separation_protocol.json` | `results/dispersion_separation_results.json`; `results/dispersion_separation_summary.json`: all 144 contracts, exact intersections, degree-five selection and every failed allowance. |
| Independent full comparison | `checks/audit_dispersion_certificates.py`; `results/dispersion_certificate_audit.json` | Alternate moment/transform constructions, all price coefficients and tails, every intersection and all 36 scale-one comparisons. |
| Figure selection and exact data | `figures/dispersion_separation_plot_data.json`; `figures/dispersion_separation_coordinates.json` | Fixed representative with no favourable separator; complete-panel counts. Log coordinates are display values, not decision arithmetic. |
| Scientific plot and native source | `figures/plot_dispersion_separation.py`; `figures/dispersion_separation_figure.tex` | Matching PNG/SVG and self-contained vector source embedded in the existing main document. The redundant text-flow figure is removed. |
| Registration and strengthening chronology | `history/registered_run1_before_intersection`; subsequent phase history | Earlier protocol/producer/results preserved, including initial metadata correction and exact-endpoint storage improvement. Final policy is not falsely presented as pre-outcome. |

The main paper states the stress identities and final counts; detailed proofs, exact case evidence, visual audit and prior source versions remain separate. Figure 1 has a specific interpretive role: it shows how the two certificate targets change with dispersion and how the stronger classical check eliminates six provisional separators. The earlier concise main source and revised documentation remain verbatim in `history/before_dispersion_separation`.

## 12 Complete four-case prices and licence-safe evidence

| Object | Detailed source | Retained evidence and scope |
| --- | --- | --- |
| Direct capped-price derivative, cell, tail and confidence proof | `four_case_classical_proof.tex`; main direct-price appendix | Global full-loading bounds, exact integer labels, fixed-sample Hoeffding confidence and full continuous-price decoding. |
| Fresh complete classical prices | `four_case_classical_pricing.py` | `results/four_case_classical_protocol.json`, `results/four_case_classical_pricing.json`, all fresh random words and q40 lower/upper arrays. The two-asset exact-sum arrays are also retained. |
| Complete score replay | `checks/audit_four_case_classical.py` | `results/four_case_classical_independent_audit.json`: all 400,000 realised scores plus 1,024 exact-sum scores using alternate 256-bit call formulas and fresh 1,024-bit cell reconstruction. |
| New direct first-residual allowance | `four_survivor_direct_r1_proof.tex` | Full loading and geometric derivative terms, raw mass and bounded tail. Numerical records separately compare its valid minimum with the older allowance. |
| Complete quantum plans and failed baseline | `four_survivor_quantum_resources.py` | `results/four_survivor_quantum_protocol.json`, full resources and summary: all sixteen targets, eighty profiles and 104 flag-confidence plans. No hardware execution. |
| Literal ordinates and width reconstruction | `results/four_survivor_quantum_scalar_ordinates.jsonl.gz`; `results/four_survivor_quantum_table_width_reconstruction.json.gz` | Exact underlying ordinates; three width-specific packed tables restored without changing executed source or results. |
| Independent resource audit | `checks/audit_four_survivor_quantum_resources.py` | `results/four_survivor_quantum_independent_audit.json`: profiles, actual Fourier tapes, controlled charges, budgets and confidence. |
| Additional independent mathematical reviews | `reviews/four_quantum_mathematical_review.md`; two review auditors | Global derivatives, coefficients, priors, finite moments and auxiliary allowances rebuilt separately. |
| Exact trace/resource algebra and data provenance | `reviews/four_case_stage_audits`; `results/four_case_data_publication_review.json` | Independent binary arithmetic, confidence, resource algebra and verified synthetic origin. A portable replay wrapper supplies the original stage layout. |
| Licence and retention rules | `DATA_POLICY.md` | No LSEG observations or systematic derived datasets to GitHub or other third parties; code and permitted limited academic summaries only. Current numerical evidence is synthetic. |

The source and timing records distinguish fresh reconstruction from the earlier lost temporary files. The old reported stochastic observations are not presented as retained evidence. New files are kept in persistent local project storage until installation. Earlier installed inputs remain byte-exact, and revised documents have exact phase snapshots. The selected classical results include no LSEG data. Publisher PDFs remain outside repository distribution.


## 13 Bounded modern-estimator comparison

| Argument or evidence | Retained artifact | Scope |
| --- | --- | --- |
| Exact odd/even circuits, conditional confidence, termination and break-even proof | `chebyshev_estimator_proof.tex`; `history/before_shared_threshold_and_figures/chebyshev_main_detail.tex.txt` | Minimal parity-correct Chebyshev signals, exact-binomial roots, every-look allocation, full fallback and expected-charge condition. |
| Outward controller and exhaustive count checks | `chebyshev_controller.py`; `checks/check_chebyshev_controller.py`; `results/chebyshev_controller_checks.json` | Two-sided prior, 101 count intervals, decoded stopping, all finite limits and explicit fallback statuses. |
| Registered complete profiles and histories | `chebyshev_resources.py`; `results/chebyshev_protocol.json`; `results/chebyshev_resources.json`; `results/chebyshev_resource_summary.json` | All four retained R1 preparations, both parities, fresh half-allowance QPE, 400 declared mesh histories, random words and unsuccessful charges. |
| Independent full resource replay | `checks/audit_chebyshev_resources.py`; `results/chebyshev_resource_audit.json` | Literal phase opcodes, reflection tapes, financial preparation linkage, all 673 batches, 67,300 pseudorandom words and observation/reset costs. |
| Independent confidence and branch replay | `checks/audit_chebyshev_math.py`; `results/chebyshev_math_audit.json`; `reviews/chebyshev_mathematical_review.md` | Exact coherent matrices, binomial roots/coverage, outward inverses and all retained scheduling decisions. |
| Additive break-even proof check | `checks/audit_chebyshev_break_even.py`; `results/chebyshev_break_even_audit.json` | Reconstructs the exact earlier proof hash and checks the added expectation identity without changing its earlier audit record. |
| Exact descriptive interpretation | `chebyshev_interpretation.py`; `results/chebyshev_interpretation.json` | All mean/median charges, strict threshold comparisons and rejection of an unknown-flag expected-cost inference. |
| Concrete sampler and scratch binding | `emit_chebyshev_sampler_templates.py`; `checks/audit_chebyshev_sampler_templates.py`; `results/chebyshev_sampler_templates.json`; `results/chebyshev_sampler_template_audit.json` | Actual success wires, excluded endpoint-high-zero wires, mapped core/reflection work, opcode adjoints and unchanged emitted costs. |
| Source and data provenance | `chebyshev_source_access.json`; `results/chebyshev_phase_provenance.json` | Final published material consulted, declared adaptation and synthetic-only input lineage. |

The sampler wire-binding supplement is separately linked in the phase registry. The main diagnostic table has a specific purpose: successful cheap branches coexist with adverse medians and a larger guaranteed cap. A figure would duplicate these four comparisons without supplying a further scientific relationship. Earlier documents have exact snapshots in `history/before_bounded_chebyshev`.

The Chebyshev supporting proof now compiles as a standalone document in the native editor. Its explicitly marked proof body is byte-identical to the previously checked fragment; the presentation wrapper supplies the document setup and its existing Rall and Fuller reference. The layout repair also places the main paper's geometric variance and covariance definitions in a multiline display. The earlier mathematical and numerical audit records remain unchanged. The presentation lineage is retained in `results/chebyshev_standalone_layout_validation.json`; `results/chebyshev_standalone_break_even_audit.json` replays the existing exact break-even checks with the marked body.

## 14 Exact bases and a stronger positive pricing encoding

| Object | Supporting derivation or source | Checked evidence |
| --- | --- | --- |
| Exact covariance and no-arbitrage invariance, changed finite laws and complete direct transfer | `factor_basis_proof.tex`; `factor_basis_certificate.py` | `results/factor_basis_certificates.json`; independent certificate and financial audits |
| Complete rebuilt quantum costs for seven bases, all methods and failures | `factor_basis_quantum_resources.py` | `results/factor_basis_quantum_resources.json`; `results/factor_basis_quantum_independent_audit.json` |
| General positive threshold coupling, CDF/strike/arithmetic certificate and clean quantum preparation | `shared_threshold_proof.tex`; `shared_threshold_oracle.py` | `results/shared_threshold_resources.json`; scalar archive; independent predicate, financial and full-depth audits |
| Exact classical conditional-variance comparison and general-payoff exclusions | `positive_threshold_coupling_proof.tex` | `reviews/shared_threshold_derivation.md`; final application attack |
| Strengthened numerical first-correction opponent, all 64 allocations and twelve rejections | `strengthened_r1_baseline.py` | `results/strengthened_r1_resources.json`; `results/strengthened_r1_baseline_audit.json` |
| Registered negative normaliser/auxiliary screen | Historical unchanged `positive_r1_feasibility.py.txt`; portable replay | `results/positive_r1_feasibility.json`; portable numerical replay record |
| Colour discrete precision map and source-preserving stress figure | `figures/plot_pricing_precision.py`; `figures/plot_dispersion_separation_v2.py` | Exact plot data, native vector fragments and figure audit |

The numerical-only comparison is retained historically. Against the cheaper admitted numerical or smoothing-certified raw-payoff rival, the fixed event encoding now saves 2.49--10.67 per cent. This supersedes the weaker all-twelve-address comparison for the four stressed cases, while leaving the original 36-contract study and its scope intact. No query-resolution reduction is claimed. Complete proofs are standalone sources; the main paper retains the identities and essential error composition. The preceding documents and visual assets are preserved in `history/before_shared_threshold_and_figures`. The installation and native-compiler record is `results/shared_threshold_project_validation.json`. Reproduction writes only to a new isolated output directory through `checks/replay_pricing_evidence.py`.

The detailed Chebyshev table and preceding main-paper appendix are retained verbatim in the phase snapshot and condensed in the main argument. The standalone Chebyshev proof remains complete and unchanged.

## 17 Common-normal raw-payoff and smoothing amendment

| Claim or evidence | Source | Checked scope |
| --- | --- | --- |
| Complete payoff, CDF, arithmetic and monetary proof | `raw_payoff_proof.tex` | Exact clamps, normal Kolmogorov transfer, shifted real register identities, intrinsic scalar errors, conditional smoothing and final price guarantee. |
| Frozen 48-circuit raw pool | `raw_payoff_oracle_tight.py`; `results/raw_payoff_tight/` | Both coherent architectures, all scalar ordinates, domain guards, priors and original pointwise plans. Earlier broad-headroom evidence remains in `results/raw_payoff_*`. |
| Proof-bound certificate amendment | `raw_payoff_smoothing_amendment.py`; `results/raw_payoff_smoothing/` | All 48 unchanged circuits and 96 rebuilt full QPE charges, declared known predecessors and every selected precision. |
| Independent mathematics and resources | `results/raw_payoff_independent_mathematics_review.json`; `results/raw_payoff_independent_resource_review.json` | Fresh coefficient/scalar/CDF/eta recomputation and own-formula confidence/charge checks. |
| Changed winning clean preparation | `check_raw_payoff_smoothing_winner_preparations.py`; `results/raw_payoff_smoothing/changed_winner_preparation_audit.json` | Both architectures of the new normal-eight winner, actual good flag, reflected core and excluded high-zero wire. |
| Exact variation refinement | `proof_sources/raw_payoff_exact_variation_fragment.tex` | Complete conditional variation formula; not used for amended costs without a further recorded rebuild. |

Against the cheaper admitted numerical or smoothing-certified common-normal raw-payoff rival in each case, the fixed shared-threshold encoding reduces complete logical T caps by 10.67, 5.26, 5.18 and 2.49 per cent. These are finite-policy sufficient implementation caps under synthetic contracts, not global optima or physical quantum advantage.

## 18 Wider precision exclusion and main-paper appendices

The separately registered extension adds exponential addresses eleven through fourteen at the same three normal widths, producing a combined 96-circuit raw pool. All four winners remain at exponential width ten, now inside the tested range. Independent rational-series scalar, whole-circuit and complete-charge audits preserve the 2.49 to 10.67 per cent event saving. A derived literal-selector charge floor also excludes feasible common exponential widths fifteen through twenty-two under the fixed graph, word precision, normal widths, unary lookup and canonical controller. Omitted lower widths, unequal payoff-node precisions and other arithmetic/estimators remain outside this exclusion. This is not a global quantum optimality theorem.

The main abstract has 159 words under the retained counting convention. Seven full proofs and seven standard derivation blocks are retained in its appendices; the main text keeps their structural sketches, financial assumptions, resource implications and adverse classical comparisons. Displayed equations are numbered, short asides use footnotes, JEL codes precede full-width naturally wrapped keywords, and the moment/derivative-bound notation is distinct.

The new producer is `expanded_raw_payoff.py`; resources and scalar ordinates are in `results/expanded_raw_payoff/`. Exact failed registration/source/log records are in `history/expanded_raw_payoff_failed/`. Independent mathematical/resource/shortcut checkers and records are retained. The full exclusion derivation is in the appendix of `raw_payoff_proof.tex`, with exact fragment in `proof_sources/raw_payoff_wider_precision_exclusion_fragment.tex`. Earlier evidence and its bound hashes remain unchanged.

## 19 Projected independent-width preparation

I separately register 300 independent PG/PA precision settings, reusing 48 exact diagonal records and emitting 252 fresh graphs. Every selected minimum is symmetric nine/nine, with normal widths ten, ten, eight and ten. The improvement is due to the previously omitted diagonal nine; asymmetric allocation supplies no additional winning gain. The strongest-rival shared-threshold savings are now 10.26, 5.26, 4.99 and 2.40 per cent. Distinct approximate tables can reverse represented ordering, so the final projection is retained and rigorously charged. The exact scalar witness and clean negative-branch replay are recorded without claiming that the witness occurs in the financial grids.

The conserved certificate uses the preceding 2E_m policy. A sharper common-shift result is prospective and uncosted. The extended literal-selector proof excludes either node at addresses fifteen through twenty-two for any feasible other node, under the retained compiler/controller. Omitted lower widths, unequal thirteen/fourteen, other word/normal precisions, shared arithmetic, other lookups and modern estimators remain unresolved. No global optimality, literature priority or practical advantage is inferred.

`asymmetric_raw_payoff.py` uses explicit repository, proof, output and frozen-protocol paths. New resources and scalar tables are in `results/asymmetric_raw_payoff/`. Independent mathematics, literal resources, branch/inverse and either-side exclusion checkers and records are retained. Both proof fragments are embedded in the standalone supporting appendix and stored exactly in `proof_sources/`; no frozen earlier evidence is overwritten.
