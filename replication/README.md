# Reproduction

## 1. Source preservation

The [source layout](source_layout.json) maps each retained file to its original execution path and SHA-256 digest. [Earlier snapshots](history/) preserve the historical records. The reconstruction tool checks the declared inventory before copying it; it does not alter recorded source hashes or run an experiment.

## 2. Verify or reconstruct

Run these commands from the organised repository root:

```sh
python replication/reproduce.py --verify-only
python replication/reproduce.py --output /tmp/qbp-replay
```

Choose a new output directory outside the repository. The reconstructed workspace restores the original flat source layout and original README. Follow that README for the frozen study commands. The [reproduction guide](../docs/reproducing_paper.md) provides installation and representative regression checks.

## 3. Interpretation

Hash verification checks agreement with the supplied inventory. Reconstruction restores paths; neither operation independently validates proofs, regenerates numerical findings or establishes a quantum speed-up. Fresh experiments require their own outputs and should preserve the retained evidence.
