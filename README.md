# Quantum Basket Pricing

## 1. Research

I study how Gaussian conditioning, approximation error and implementation costs affect the valuation of capped basket derivatives under an explicitly specified classical, no-arbitrage lognormal model. I use the project to demonstrate mathematical and computational preparation for quantitative finance graduate study, connecting derivative pricing to probability, numerical analysis, statistical inference, constrained optimisation and quantum computation.

The [quantum pricing working paper](paper/generalisation_proofs.pdf), with its [LaTeX source](paper/generalisation_proofs.tex), compares certified pricing constructions under declared information, error and resource assumptions. Its retained synthetic contracts show logical T-cap reductions of 2.34–10.00 per cent against the cheaper equally optimised implemented rival; the comparisons with classical exact summation and validated integration establish neither practical quantum advantage nor a general superiority result.

The [capped-GBM companion](paper/capped_gbm_companion.pdf) develops analytic/PDE/Monte Carlo prices and Greeks, conditional controls and Gaussian parameter inference. The [constrained-hedging companion](paper/constrained_hedging_companion.pdf) derives a funded mean-square objective, convexity and KKT conditions, with constrained, unrestricted and cash-only comparisons. The [paper guide](paper/README.md) and [manifest](paper/manifest.json) connect the documents to their retained evidence. All current financial experiments are synthetic; the main paper remains a working paper.

## 2. Reproduction

Use Python 3.13 and install the declared environment, check the organised repository, then reconstruct the original execution layout in a fresh directory:

```sh
python -m pip install -r requirements.lock
python scripts/check_repository.py
python replication/reproduce.py --output /tmp/qbp-replay
cd /tmp/qbp-replay
python -m unittest discover -s checks -p 'test_capped_gbm.py'
```

The repository check verifies layout and retained-file invariants; it does not rerun every pricing experiment or certify every research claim. The [reproduction guide](docs/reproducing_paper.md) and [replication guide](replication/README.md) explain the reconstruction, dependencies and individual study commands. Run each study in a fresh workspace because historical evidence uses its original paths.

The [implementation](src/), [script guide](scripts/README.md) and [tests](tests/) distinguish reusable research modules, executable studies and checks. The [configuration guide](configs/README.md) records protocols and specialised environments; the [data guide](data/README.md) identifies permitted public reference material. These research modules are not an installed software package.

## 3. Research record

The [preparation audit](docs/PREPARATION.md) maps demonstrated mathematical work to the reviewed programmes and identifies remaining gaps. The [source register](docs/SOURCES.md), [supporting index](docs/supporting_proofs_index.md) and [originality assessment](docs/originality_assessment.md) distinguish established methods, project-specific derivations, numerical findings and unresolved literature comparisons.

The [results](results/) retain favourable and adverse findings, checksums and source-access records. The [history](replication/history/) preserves earlier research artefacts; the [layout record](replication/source_layout.json) preserves their execution-path mapping. No LSEG observations or publisher article texts are distributed. Public reference material retains its source notices, and a private repository does not extend third-party redistribution rights.
