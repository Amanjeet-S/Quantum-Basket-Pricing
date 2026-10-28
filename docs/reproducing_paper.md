# Reproducing the studies

## 1. Scope

The organised repository separates [papers](../paper/README.md), [numerical modules](../src/), [experiment entry points](../scripts/README.md), [tests](../tests/), [specifications](../configs/README.md) and [retained results](../results/). Frozen scientific files preserve their original bytes and execution paths through the [source layout](../replication/source_layout.json).

Current financial experiments use synthetic inputs. These commands require no paid financial dataset or provider credentials. The [data policy](DATA_POLICY.md) applies separately to any authorised observations held locally.

## 2. Environment and repository checks

Use Python 3.13 in a separate environment. From the organised repository root:

```sh
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python scripts/check_repository.py
python -m unittest discover -s tests -p 'test_reproduce.py'
```

The checker verifies retained hashes, the paper inventory and current navigation links. The tests exercise the reconstruction tool. These operations do not regenerate the scientific results. Optional study environments remain in [configs/environments](../configs/environments/).

## 3. Restore the original execution layout

From the same repository root:

```sh
python replication/reproduce.py --verify-only
python replication/reproduce.py --output /tmp/qbp-replay
cd /tmp/qbp-replay
```

Use a new directory outside the repository for each reconstruction. The tool verifies every mapped file before copying and checks each copy. Its output preserves the original README and relative paths. Execute the frozen study commands there, rather than invoking moved scripts directly from the organised tree.

## 4. Companion regression checks

In the reconstructed workspace, with the numerical environment active:

```sh
python -m unittest discover -s checks -p 'test_capped_gbm.py'
python -m unittest discover -s checks -p 'test_conditioning_controls.py'
python -m unittest discover -s checks -p 'test_gbm_calibration.py'
python -m unittest discover -s checks -p 'test_gbm_hedging.py'
```

These four suites contain 40 implementation checks for pricing, conditional controls, calibration and hedging. They test selected cases and boundary behaviour; they do not substitute for the proofs or a full rerun of the retained experiments.

## 5. Full study reproduction

Follow the reconstructed README for each study's commands, dependencies and evidence requirements. Preserve the retained outputs and select new destinations for fresh runs. Timing measurements depend on the executing machine; source agreement alone does not guarantee identical floating-point results.

The [paper guide](../paper/README.md) identifies the current PDFs and sources. The [research guide](research_guide.md) provides the detailed study map. Rebuilding a PDF, auditing a released summary and rerunning an experiment are different operations and should be reported separately.

## 6. Native pricing and exact trace audit

These new studies run directly in the organised repository, independently of the historical reconstruction. First follow the [C++ build instructions](../cpp/README.md) to create the two executables in `/tmp/qbp-native-build`. From the repository root, with Python 3.13 active:

```sh
python -m pip install -r configs/environments/requirements-native.txt
python -m unittest discover -s tests -p 'test_cpp_residual_certificate.py'
python scripts/run_cpp_pricing_study.py --executable /tmp/qbp-native-build/qbp_gbm --output /tmp/qbp-native-comparison
python scripts/certify_cpp_heat_trace.py --executable /tmp/qbp-native-build/qbp_gbm --output /tmp/qbp-native-traces
```

Use new external output directories. The runners retain protocols, exact inputs, failures, environment information and source/executable hashes. [Native evidence](../results/capped_gbm/cpp/comparison/study_manifest.json) separates implementation checks from the four retained stronger-Python derivative diagnostics. The trace runner exits with status 1 for four retained strict-range overshoots; its exact residual bounds cover those excesses. The [trace evidence](../results/capped_gbm/cpp/traces/study_manifest.json) certifies the published-coefficient discrete system; it does not certify the original continuous financial price. The [mathematical note](cpp_pricing.md) proves this scope and records the informed design amendments. The earlier scientific records remain unchanged; the current companion exports and their preserved predecessors are identified in the [paper manifest](../paper/manifest.json).

## 7. Validated scalar parameter-to-price bridge

The [certificate guide](validated_calibration.md) gives the exact input semantics, Gaussian/rounding premises and reproduction commands. Install [its environment](../configs/environments/requirements-price-certificate.txt), run `python -m unittest discover -s tests -p 'test_gbm_price_certificate.py'`, then run `python scripts/run_gbm_certificate_study.py --output /tmp/qbp-validated-calibration` from this organised tree. The output must be a new or empty external directory. Numerical/reference failures cause a nonzero exit status; true-parameter containment misses are retained statistical diagnostics. [The evidence](../results/capped_gbm/validated_calibration/validation.json) records all attempted datasets, abstentions, source hashes and independent checks. The updated companion's prior bytes remain in [the snapshot](../replication/history/before_validated_calibration/); the historical layout still reproduces those frozen bytes.

## 8. Validated static hedging

Use the same price-certificate environment from the organised repository root:

```sh
python -m unittest discover -s tests -p 'test_gbm_hedge_certificate.py'
python scripts/run_hedge_certificate_study.py --output /tmp/qbp-validated-hedging
```

The output must be new or empty and outside the repository. The [guide](validated_hedging.md) distinguishes stored, repaired and refined holdings and the exact fair-price target from the archived rounded-price calculation. The [retained evidence](../results/capped_gbm/validated_hedging/validation.json) contains all cases and the actual verification status.

## 9. Known-volatility drift uncertainty

Use the same certificate environment from the organised tree:

```sh
python -m unittest discover -s tests -p 'test_gbm_drift_hedging.py'
python scripts/run_drift_hedge_study.py --output /tmp/qbp-drift-hedging
```

The [guide](drift_hedging.md) gives the financial target, protocol and result boundaries. The [figure generator](../scripts/plot_drift_hedging.py) accepts `--results /tmp/qbp-drift-hedging --output /tmp/qbp-drift-figures`; these must be existing study evidence and a new external figure destination respectively. Its native LaTeX fragment and CSV preserve the plotted information.

## 10. Final verification and matched-protocol restoration

The [finalisation record](project_status.md) distinguishes completed checks from optional research. Full historical test discovery also imports the retained tilt experiment, which requires the separately pinned `mpmath==1.4.1`; the price-certificate environment deliberately uses `mpmath==1.3.0`. Use separate environments for these two suites. The first full-suite attempt with the latter version is retained as an environment failure; the correctly versioned historical suite passes all 193 tests.

The released matched protocol/input contain provenance rebindings, whereas the authoritative result commits to their original bytes. Follow the explicit recovery commands in [research guide Section 4](research_guide.md#4-reproduce-the-current-independent-checks). The [preparation helper](../replication/prepare_matched_replay.py) checks exact original/exported hashes, the declared provenance-only transformations, all nine scientific sources and a fresh independent resource review before writing a new external folder. It preserves both released exports and the unchanged auditor. Preparation alone is not a successful audit.
