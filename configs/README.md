# Study specifications and environments

## 1. Registered specifications

The [hedging specification](hedging_protocol.json), [figure protocol](figure1_protocol.json) and [initial figure selection](figure1_protocol_initial_selection.json) retain the existing registered designs. Their organised locations do not change the paths recorded in earlier studies; the reproduction tool restores those original paths.

## 2. Dependencies

The root [dependency lock](../requirements.lock) specifies the main numerical environment. [Study-specific requirement files](environments/) retain the earlier dependency declarations, including optional rotation-synthesis environments. These environments serve different calculations and need not be combined.

The [reproduction guide](../docs/reproducing_paper.md) distinguishes direct repository checks from frozen study execution. Neither a dependency declaration nor a successful installation establishes numerical reproducibility across different hardware and arithmetic libraries.
