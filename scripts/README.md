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
