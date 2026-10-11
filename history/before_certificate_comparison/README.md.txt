# Certified curvature corrections for quantum basket pricing

## 1 Research question and current finding

I develop coherent derivative-pricing methods under an established classical no-arbitrage GBM model. The [main paper](generalisation_proofs.tex) contains the essential proofs and focused appendices. The [supporting index](supporting_proofs_index.md) maps full derivations, exact evidence and earlier investigations.

The payoff-boundary method removes a strike auxiliary and supplies a positive second correction for convex payoffs. The new Gaussian transport certificate keeps dispersion moments and certifies the actual rounded financial grid through scalar transforms, without a continuous basket-price answer or Cartesian enumeration. It completes the formerly missing four- and eight-asset financial transfers at unchanged numerical widths.

The matched comparison now covers all 36 original capped and ordinary-call cases. Shared low-order continuous bounds resolve 18 without quantum estimation. In the other 18 cases, second-correction methods improve complete logical T caps against the informed first correction in ten, by approximately 13.9 to 92.6 per cent. The boundary method beats raw direct second-correction evaluation by about 3.18 per cent for the eight-asset ordinary call at common scale 0.05; raw evaluation wins the other active second-correction comparisons.

Strong classical methods certify all 36 model prices. No physical quantum pricing experiment was executed. The work contains checked mathematics and original conditional computational evidence; a contribution's novelty in the literature remains unestablished. The [novelty review](novelty_scope_review.md) identifies close antecedents and the exact comparison still needed.

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
