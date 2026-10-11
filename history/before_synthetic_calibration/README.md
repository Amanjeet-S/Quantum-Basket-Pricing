# Quantum Basket Pricing

## 1 Research purpose and status

My primary portfolio purpose is to demonstrate mathematical and computational preparation for quantitative finance graduate study. I study coherent pricing of basket derivatives under a stipulated classical, no-arbitrage lognormal model. Gaussian conditioning, geometric controls and payoff-boundary corrections connect the financial target to finite quantum preparations. I charge tails, rounding, interpolation, state preparation, confidence, reflections and terminal resets.

The [preparation audit](PREPARATION.md) maps retained proofs and code to current entry requirements and published curricula, with specific remaining gaps and a [complete public module/descriptor register](results/preparation_review/module_descriptors.json), including unavailable and prior-cohort source limits. A novel quantum theory or practical quantum advantage is not required for this preparation purpose; mathematical understanding and reproducible financial computation are the priority.

The [main paper](generalisation_proofs.tex) is a **working paper**. Its main text contains intuition and structural steps; appendices and [supporting proofs](supporting_proofs_index.md) retain full derivations. The present work does not establish an original quantum theory.

## 2 Current evidence

The current comparison gives shared-threshold, raw-payoff and numerical conditional-price preparations the same financial model, information, monetary radius 0.01 and per-price failure allowance 1/3600. Every applicable method receives the same cancellation-aware optimisation.

| Result | Checked scope |
| --- | --- |
| Shared-threshold logical T-cap reduction | 2.34–10.00 per cent against the cheaper equally optimised implemented rival in each of four stressed contracts |
| Terminal cleanup deferral | 0.097–0.195 per cent against cancellation-optimised clean QPE |
| Resource comparison | 416 settings; 808 admitted architecture profiles; 1,616 phase plans; 3,232 complete variant charges; twelve retained rejected budgets |
| Matched Chebyshev diagnostics | 1,200 alternative histories; the numerical rival has lower descriptive mesh means in all four cases |
| Classical pricing | Strong classical methods certify all 36 original contracts and four stressed contracts |

These are synthetic, finite-policy logical resource comparisons. They establish neither practical quantum advantage nor adaptive expected-cost superiority. Original generic smoothing, reflection and cancellation mechanisms are excluded from the novelty claim. The [originality assessment](originality_assessment.md) identifies the guarded financial certificate and verified comparisons as the contribution candidate; literature priority remains unestablished.

Sharp Gaussian component moduli are proved and interval-checked, but have not yet been selected in the registered resource costs. The [capped-GBM companion](capped_gbm_companion.tex) now connects replication, PDE stability/convergence and analytic/PDE/Monte Carlo prices and five Greeks. Its seven synthetic cases include a stronger bounded mixed-gamma comparator and adverse control rankings. The companion additionally derives L² projection/singular-control geometry and compares conditional prices and independent pilots with explicit evaluator work. The [current plan](quantum_method_redesign.md) records remaining preparation work; a fair sharp-bound cost study remains a separate research extension.

## 3 Navigation

| Purpose | Entry point |
| --- | --- |
| Financial model, results and main argument | [Working paper](generalisation_proofs.tex) |
| Mathematical preparation, curricula and gaps | [Preparation audit](PREPARATION.md) |
| Classical one-asset pricing/PDE/Greek benchmark | [Companion paper](capped_gbm_companion.pdf) and [synthetic study](results/capped_gbm/study_manifest.json) |
| Conditional expectation and jointly projected controls | [Source/runner](run_conditioning_study.py) and [conditioning evidence](results/capped_gbm/conditioning/study_manifest.json) |
| Complete proofs and retained evidence | [Supporting index](supporting_proofs_index.md) |
| Published antecedents and versions actually consulted | [Source register](SOURCES.md) |
| Exact current costs | [Cancellation-aware summary](results/deferred_cleanup/deferred_cleanup_summary.json) |
| Current quantum-controller diagnostics | [Matched summary](results/matched_deferred_chebyshev/matched_chebyshev_summary.json) |
| Latest research validation | [Validation record](results/originality_project_validation.json) |
| Earlier investigations | [Historical documentation](history/before_documentation_consolidation/supporting_proofs_index.md) |

The main paper embeds its current colour figures as self-contained vector source. The earlier pricing_precision.png/svg exports are historical. Earlier protocols, failures, source versions and derivations remain preserved; they are not current resource rankings.

## 4 Reproduce the current independent checks

Use Python 3.13 and python-flint 0.9.0, specified in [requirements-quadrature.txt](requirements-quadrature.txt). NumPy in [requirements.txt](requirements.txt) supports earlier numerical stages. From the repository root, create a fresh output directory and run:

    qbp_audit_dir=$(mktemp -d /tmp/qbp-audit.XXXXXX)
    python check_sharp_gaussian_moduli.py --output "$qbp_audit_dir/sharp.json"
    python check_deferred_resources_independent.py --repo . --results results/deferred_cleanup --proof proof_sources/cancellation_aware_preparation_fragment.tex --output "$qbp_audit_dir/resources.json"
    python check_matched_deferred_chebyshev.py --repo . --results results/matched_deferred_chebyshev --deferred-resources results/deferred_cleanup/deferred_cleanup_resources.json.gz --deferred-review "$qbp_audit_dir/resources.json" --input-evidence results/deferred_cleanup/matched_B_inputs.json --output "$qbp_audit_dir/matched.json"

The relocation argument verifies the supplied-input file against the original frozen checksum. It does not rewrite the protocol or substitute new evidence. These checks reconstruct certified mathematics and logical charges; they execute no quantum hardware.

## 5 Reproduce the capped-GBM companion

With Python 3.13 and NumPy 2.3.5, run from the repository root:

    python3.13 -m unittest discover -s checks -p test_capped_gbm.py
    qbp_gbm_audit_dir=$(mktemp -d /tmp/qbp-gbm.XXXXXX)
    python3.13 run_capped_gbm_study.py --output "$qbp_gbm_audit_dir"

The retained study contains 336 PDE comparisons, 147 Monte Carlo comparisons and explicit analytic bump bias. The analytic call spread is the strongest baseline for this one-asset model. Exact-arithmetic proofs and ideal-sampling confidence bounds are distinguished from floating-point numerical evidence. [The source ledger](results/capped_gbm/source_ledger.json) identifies the published versions actually consulted. [The colour figure](results/capped_gbm/figures/capped_gbm_convergence.svg) displays non-monotone errors and gamma bias rather than hiding them.

## 6 Reproduce conditional valuation and controls

Using the same Python/NumPy environment:

    python3.13 -m unittest discover -s checks -p test_conditioning_controls.py
    qbp_conditioning_audit_dir=$(mktemp -d /tmp/qbp-conditioning.XXXXXX)
    python3.13 run_conditioning_study.py --output "$qbp_conditioning_audit_dir"

The study retains 21 case/splits, 147 comparisons and 672 independent pilot-risk diagnostics. Conditional moments are floating-point quadrature diagnostics, not interval certificates. The population variance order applies to jointly projected controls; fitted coefficients and evaluator costs remain separate. [The source ledger](results/capped_gbm/conditioning/source_ledger.json) records the published material actually consulted. [Current companion validation](results/capped_gbm/conditioning/validation.json) covers this extension. Earlier capped-GBM comparison files remain unchanged; their original 14-page paper/source hashes resolve to [the stage snapshot](history/before_projection_controls/capped_gbm_companion.pdf).

## 7 Data and preservation

Current numerical results use synthetic model inputs and permitted public reference material. No LSEG observations are included. [DATA_POLICY.md](DATA_POLICY.md) implements the university's restrictions: raw, reconstructible and systematic licensed data must not enter GitHub, including private repositories.

Registered evidence and hashes remain immutable. Concision changes navigation and presentation, with earlier documents retained byte-for-byte in the [documentation snapshot](history/before_documentation_consolidation/README.md).

