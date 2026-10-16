# Certified curvature corrections for quantum basket pricing

## 1 Research question and current finding

I develop coherent derivative-pricing methods under an established classical no-arbitrage GBM model. The [main paper](generalisation_proofs.tex) contains the essential proofs and focused appendices. The [supporting index](supporting_proofs_index.md) maps full derivations, exact evidence and earlier investigations.

The payoff-boundary method removes a strike auxiliary and supplies a positive second correction for convex payoffs. The residual-specific Gaussian transport certificate keeps dispersion moments and certifies the actual rounded financial grid through scalar transforms, without a continuous basket-price answer or Cartesian enumeration. It completes the formerly missing four- and eight-asset financial transfers at unchanged numerical widths.

The matched comparison now covers all 36 original capped and ordinary-call cases. Shared low-order continuous bounds resolve 18 without quantum estimation. In the other 18 cases, second-correction methods improve complete logical T caps against the informed first correction in ten, by approximately 13.9 to 92.6 per cent. The boundary method beats raw direct second-correction evaluation by about 3.18 per cent for the eight-asset ordinary call at common scale 0.05; raw evaluation wins the other active second-correction comparisons.

Strong classical methods certify all 36 model prices. No physical quantum pricing experiment was executed. The work contains checked mathematics and original conditional computational evidence; a contribution's novelty in the literature remains unestablished. The [source comparison](certificate_source_comparison.md) identifies an exact published antecedent for the positive cell identity. The residual derivative closure and its measured implementation value remain the narrower contribution candidate; the [novelty review](novelty_scope_review.md) states the unresolved priority gaps.

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
