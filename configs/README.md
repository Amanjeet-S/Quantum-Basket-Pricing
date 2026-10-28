# Study specifications and environments

## 1. Registered specifications

The [hedging specification](hedging_protocol.json), [figure protocol](figure1_protocol.json) and [initial figure selection](figure1_protocol_initial_selection.json) retain the existing registered designs. Their organised locations do not change the paths recorded in earlier studies; the reproduction tool restores those original paths.

## 2. Dependencies

The root [dependency lock](../requirements.lock) specifies the main numerical environment. [Study-specific requirement files](environments/) retain the earlier dependency declarations, including optional rotation-synthesis environments. These environments serve different calculations and need not be combined.

The [reproduction guide](../docs/reproducing_paper.md) distinguishes direct repository checks from frozen study execution. Neither a dependency declaration nor a successful installation establishes numerical reproducibility across different hardware and arithmetic libraries.

The [native validation environment](environments/requirements-native.txt) pins the high-precision and binary64 comparators. The [native amendment record](cpp_pricing_amendment_history.json) retains the initial protocol and informed stress/baseline changes; this development study is not described as an untouched prospective protocol.

The [price-certificate environment](environments/requirements-price-certificate.txt) pins the validated-arithmetic and independent diagnostic libraries. The [certificate amendment record](gbm_certificate_amendments.json) preserves the informed numerical corrections without changing any retained earlier experiment.

## 3. Physical-drift uncertainty

The [protocol](drift_hedge_protocol.json) fixes the compact-drift cases, budgets, confidence diagnostics and first-three-dataset inference bridge. It records informed development and the bridge amendment. The [guide](../docs/drift_hedging.md) gives the exact statistical/model scope.
