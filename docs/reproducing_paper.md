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
