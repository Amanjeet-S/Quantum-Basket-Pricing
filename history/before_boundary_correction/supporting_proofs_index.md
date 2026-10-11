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
