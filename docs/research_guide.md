# Quantum Basket Pricing: detailed research guide

Historical study commands in Sections 5–8 run inside the isolated workspace created by `python3.13 replication/reproduce.py --output /tmp/qbp-replay` from the organised repository root. Section 4 prepares its own fresh audit workspace. The reconstructed `README.md` preserves original instructions and paths; current extensions identify organised-tree commands separately.

## 1 Research purpose and status

My primary portfolio purpose is to demonstrate mathematical and computational preparation for quantitative finance graduate study. I study coherent pricing of basket derivatives under a stipulated classical, no-arbitrage lognormal model. Gaussian conditioning, geometric controls and payoff-boundary corrections connect the financial target to finite quantum preparations. I charge tails, rounding, interpolation, state preparation, confidence, reflections and terminal resets.

The [preparation audit](PREPARATION.md) maps retained proofs and code to current entry requirements and published curricula, with specific remaining gaps and a [complete public module/descriptor register](../results/preparation_review/module_descriptors.json), including unavailable and prior-cohort source limits. A novel quantum theory or practical quantum advantage is not required for this preparation purpose; mathematical understanding and reproducible financial computation are the priority.

The [main paper PDF](../paper/generalisation_proofs.pdf) is the final project paper for the defined MSc-portfolio scope ([LaTeX source](../paper/generalisation_proofs.tex)). Its main text contains intuition and structural steps; appendices and [supporting proofs](supporting_proofs_index.md) retain full derivations. The [finalisation record](project_status.md) states the checked scope. The manuscript remains unpublished and does not establish an original quantum theory.

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

Sharp Gaussian component moduli are proved and interval-checked, but have not yet been selected in the registered resource costs. The [capped-GBM companion](../paper/capped_gbm_companion.tex) now connects replication, PDE stability/convergence and analytic/PDE/Monte Carlo prices and five Greeks. Its seven synthetic cases include a stronger bounded mixed-gamma comparator and adverse control rankings. The companion additionally derives L² projection/singular-control geometry and compares conditional prices and independent pilots with explicit evaluator work. Synthetic Gaussian-return inference adds covariance-rank and exact parameter-interval proofs, including an interior-volatility price maximum where endpoint-only uncertainty propagation fails. Physical estimation remains separate from risk-neutral pricing. The [current plan](quantum_method_redesign.md) records completed project work and broader preparation gaps; a fair sharp-bound cost study remains a separate research extension.

The separate [hedging companion](../paper/constrained_hedging_companion.tex) connects physical mean-square loss, funded trading, convexity/KKT proofs and budget shadow prices. It compares constrained holdings with cash-only, unrestricted and same-constraint variance baselines, retaining the original arithmetic-floor diagnostics and adverse heuristic rankings. The [validated extension](validated_hedging.md) resolves positive risks with outward feasibility and primal-dual bounds.

## 3 Navigation

| Purpose | Entry point |
| --- | --- |
| Financial model, results and main argument | [Final project paper PDF](../paper/generalisation_proofs.pdf) and [LaTeX source](../paper/generalisation_proofs.tex) |
| Mathematical preparation, curricula and gaps | [Preparation audit](PREPARATION.md) |
| Classical one-asset pricing/PDE/Greek benchmark | [Companion paper](../paper/capped_gbm_companion.pdf) and [synthetic study](../results/capped_gbm/study_manifest.json) |
| Gaussian-return estimation and price uncertainty | [Calibration runner](../scripts/run_calibration_study.py) and [synthetic evidence](../results/capped_gbm/calibration/study_manifest.json) |
| Conditional expectation and jointly projected controls | [Source/runner](../scripts/run_conditioning_study.py) and [conditioning evidence](../results/capped_gbm/conditioning/study_manifest.json) |
| Constrained financial hedging | [Hedging companion](../paper/constrained_hedging_companion.pdf) and [registered study](../results/capped_gbm/hedging/study_manifest.json) |
| Complete proofs and retained evidence | [Supporting index](supporting_proofs_index.md) |
| Published antecedents and versions actually consulted | [Source register](SOURCES.md) |
| Exact current costs | [Cancellation-aware summary](../results/deferred_cleanup/deferred_cleanup_summary.json) |
| Current quantum-controller diagnostics | [Matched summary](../results/matched_deferred_chebyshev/matched_chebyshev_summary.json) |
| Latest research validation | [Validation record](../results/originality_project_validation.json) |
| Earlier investigations | [Historical documentation](../replication/history/before_documentation_consolidation/supporting_proofs_index.md) |

The main paper embeds its current colour figures as self-contained vector source. The earlier pricing_precision.png/svg exports are historical. Earlier protocols, failures, source versions and derivations remain preserved; they are not current resource rankings.

## 4 Reproduce the current independent checks

Use Python 3.13 and python-flint 0.9.0. Run the following from the organised repository root; it creates a new reconstructed workspace and independent output directory:

```sh
qbp_project_root=$(pwd)
qbp_replay_root=$(mktemp -d /tmp/qbp-replay.XXXXXX)
qbp_audit_dir=$(mktemp -d /tmp/qbp-audit.XXXXXX)
python replication/reproduce.py --output "$qbp_replay_root"
cd "$qbp_replay_root"
python check_sharp_gaussian_moduli.py --output "$qbp_audit_dir/sharp.json"
python check_deferred_resources_independent.py --repo . --results results/deferred_cleanup --proof proof_sources/cancellation_aware_preparation_fragment.tex --output "$qbp_audit_dir/resources.json"
python "$qbp_project_root/replication/prepare_matched_replay.py" --repository "$qbp_project_root" --workspace "$qbp_replay_root" --deferred-review "$qbp_audit_dir/resources.json" --output "$qbp_audit_dir/matched_bindings"
python check_matched_deferred_chebyshev.py --repo . --results "$qbp_audit_dir/matched_bindings" --deferred-resources results/deferred_cleanup/deferred_cleanup_resources.json.gz --deferred-review "$qbp_audit_dir/resources.json" --input-evidence "$qbp_audit_dir/matched_bindings/matched_B_inputs.json" --output "$qbp_audit_dir/matched.json"
```

The preparation helper restores exact original protocol/input bindings from the [checked recovery archive](../replication/recovery/matched_replay_bindings_manifest.json). It verifies the declared path/hash transformations and all nine scientific-source premises; released redacted exports remain unchanged. The supplied input's original checksum remains authoritative for this audit. The unchanged auditor then checks financial allowances, literal logical charges and paired histories. Preparation is not itself a passed audit, and these checks execute no quantum hardware.

## 5 Reproduce the capped-GBM companion

With Python 3.13 and NumPy 2.3.5, run from the repository root:

    python3.13 -m unittest discover -s checks -p test_capped_gbm.py
    qbp_gbm_audit_dir=$(mktemp -d /tmp/qbp-gbm.XXXXXX)
    python3.13 run_capped_gbm_study.py --output "$qbp_gbm_audit_dir"

The retained study contains 336 PDE comparisons, 147 Monte Carlo comparisons and explicit analytic bump bias. The analytic call spread is the strongest baseline for this one-asset model. Exact-arithmetic proofs and ideal-sampling confidence bounds are distinguished from floating-point numerical evidence. [The source ledger](../results/capped_gbm/source_ledger.json) identifies the published versions actually consulted. [The colour figure](../results/capped_gbm/figures/capped_gbm_convergence.svg) displays non-monotone errors and gamma bias rather than hiding them.

## 6 Reproduce conditional valuation and controls

Using the same Python/NumPy environment:

    python3.13 -m unittest discover -s checks -p test_conditioning_controls.py
    qbp_conditioning_audit_dir=$(mktemp -d /tmp/qbp-conditioning.XXXXXX)
    python3.13 run_conditioning_study.py --output "$qbp_conditioning_audit_dir"

The study retains 21 case/splits, 147 comparisons and 672 independent pilot-risk diagnostics. Conditional moments are floating-point quadrature diagnostics, not interval certificates. The population variance order applies to jointly projected controls; fitted coefficients and evaluator costs remain separate. [The source ledger](../results/capped_gbm/conditioning/source_ledger.json) records the published material actually consulted. [Conditioning-stage validation](../results/capped_gbm/conditioning/validation.json) covers this extension. Earlier capped-GBM comparison files remain unchanged; their original 14-page paper/source hashes resolve to [the stage snapshot](../replication/history/before_projection_controls/capped_gbm_companion.pdf).

## 7 Reproduce Gaussian estimation and price uncertainty

Use Python 3.13 with [NumPy 2.3.5 and SciPy 1.16.3](../configs/environments/requirements-calibration.txt):

    python3.13 -m unittest discover -s checks -p test_gbm_calibration.py
    qbp_calibration_audit_dir=$(mktemp -d /tmp/qbp-calibration.XXXXXX)
    python3.13 run_calibration_study.py --output "$qbp_calibration_audit_dir"

Eight calibration cases and five frequency/horizon regimes retain 13,000 independent synthetic datasets, with aggregate results only. At the price's volatility maximum, the full volatility-interval image covers 948/1,000 datasets; endpoint-only propagation covers 0/1,000. The exact-model proof is distinct from the floating-point quantiles and prices. The [validated scalar bridge](validated_calibration.md) now implements rounding-aware scatter, rigorous chi-squared brackets and outward price-image extrema. Successful-return risk and abstentions are distinguished; no market-model validation is claimed. The original calibration experiment remains unchanged, and its companion source/export hashes resolve to [the preserved pre-certificate companion](../replication/history/before_validated_calibration/). The [source ledger](../results/capped_gbm/calibration/source_ledger.json) records the actual consulted versions. [Current companion validation](../results/capped_gbm/calibration/validation.json) covers all three companion stages. The earlier 22-page paper/source hashes resolve to [the preserved conditioning-stage paper](../replication/history/before_synthetic_calibration/capped_gbm_companion.pdf).

## 8 Reproduce constrained static hedging

Use the Python/NumPy/SciPy environment above, from the repository root:

    python3.13 -m unittest discover -s checks -p test_gbm_hedging.py
    qbp_hedging_audit_dir=$(mktemp -d /tmp/qbp-hedging.XXXXXX)
    python3.13 run_hedging_study.py --output "$qbp_hedging_audit_dir"

The [protocol](../configs/hedging_protocol.json) registers seven contracts, three physical drifts and six budgets before outcomes. Initial costs use risk-neutral prices; loss is the physical expectation of squared discounted terminal liability error, including bias. The study retains 126 constrained mean-square hedges, strong population comparisons and 420 matched Monte Carlo diagnostics. Static holdings use reinvested-dividend stock and one lower-strike call. The unrestricted baseline is a formal relaxation; the exact call spread changes the instrument set. [Validation](../results/capped_gbm/hedging/validation.json) and [source versions](../results/capped_gbm/hedging/source_ledger.json) retain the mathematical/numerical distinction.

The colour funding/risk figure has a separate supplemental grid. To reproduce its exports with Matplotlib 3.10.6:

    python3.13 plot_hedging.py --output "$qbp_hedging_audit_dir/figures"

## 9 Data and preservation

Current numerical results use synthetic model inputs and permitted public reference material. No LSEG observations are included. [DATA_POLICY.md](DATA_POLICY.md) implements the university's restrictions: raw, reconstructible and systematic licensed data must not enter GitHub, including private repositories.

Registered evidence and hashes remain immutable. Concision changes navigation and presentation, with earlier documents retained byte-for-byte in the [documentation snapshot](../replication/history/before_documentation_consolidation/README.md).

## 10 Physical-drift uncertainty

The [drift guide](drift_hedging.md) supplies organised-tree commands for the new extension. It connects known-volatility Gaussian inference to a compact-drift minimax hedge, a continuous-risk certificate and same-budget baselines. Worst-case protection and performance at the true drift are reported separately; original registered results remain unchanged.
