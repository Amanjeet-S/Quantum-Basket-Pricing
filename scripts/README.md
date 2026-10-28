# Study and verification entry points

## 1. Frozen research scripts

This directory contains study runners, plotting commands and other executable research entry points. [Independent audits](audits/) are separate from the automated [tests](../tests/). Reusable numerical modules are retained in [src/](../src/).

These files retain their original bytes, imports and relative paths. Their organised locations are for navigation; execute the original commands in a reconstructed workspace using the [reproduction guide](../docs/reproducing_paper.md). Moving a file does not make its dependencies importable from its new directory.

## 2. Current repository check

The new repository checker runs directly from the organised repository root:

```sh
python scripts/check_repository.py
```

It checks retained hashes, the paper inventory and current navigation links. It does not rerun the financial experiments or establish the mathematical claims.

## 3. Native pricing studies

The new [C++ comparison](run_cpp_pricing_study.py) and [exact heat-trace auditor](certify_cpp_heat_trace.py) run directly from this organised tree. They load the existing Python engine without modifying its frozen bytes. Build the [native engine](../cpp/README.md), install the [small validation environment](../configs/environments/requirements-native.txt), then follow the [native reproduction commands](../docs/reproducing_paper.md#6-native-pricing-and-exact-trace-audit). Both runners require a fresh output destination outside the repository and preserve failures. They do not overwrite retained evidence.

## 4. Validated parameter-to-price study

The [scalar certificate runner](run_gbm_certificate_study.py) loads the exact-arithmetic [engine](../src/gbm_price_certificate.py) and [independent references](../src/gbm_certificate_references.py) directly from this tree. The [guide](../docs/validated_calibration.md) supplies the environment, commands, outcome interpretation and retained evidence. Study outputs require a new or empty external destination. The [hedging certificate runner](run_hedge_certificate_study.py) uses the same environment; its [guide](../docs/validated_hedging.md) explains exact feasibility, risk/optimality bounds and independent diagnostics. The [drift runner](run_drift_hedge_study.py) and [guide](../docs/drift_hedging.md) extend that calculation to a complete physical-drift interval and statistical inference.
