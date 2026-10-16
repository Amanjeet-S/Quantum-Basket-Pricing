# Quantum Basket Pricing

## 1 Purpose and status

I use this project to demonstrate mathematical and computational preparation for quantitative finance graduate study. It connects derivative valuation under a stipulated classical, no-arbitrage lognormal model to Gaussian conditioning, numerical analysis, statistical inference, constrained hedging and quantum pricing resource estimates.

The main quantum paper is a working paper. All current financial experiments use synthetic inputs. The resource comparisons do not establish practical quantum advantage or an original quantum theory.

## 2 Papers and research navigation

| Purpose | Entry point |
| --- | --- |
| Quantum basket pricing: financial model, certificates and comparisons | [Working paper PDF](papers/generalisation_proofs.pdf) · [LaTeX source](papers/generalisation_proofs.tex) |
| Classical pricing, PDEs, Greeks, conditional controls and Gaussian inference | [Capped-GBM companion PDF](papers/capped_gbm_companion.pdf) · [LaTeX source](papers/capped_gbm_companion.tex) |
| Constrained static hedging, convexity and KKT conditions | [Hedging companion PDF](papers/constrained_hedging_companion.pdf) · [LaTeX source](papers/constrained_hedging_companion.tex) |
| Mathematical preparation, programme curricula and remaining gaps | [Preparation audit](docs/PREPARATION.md) |
| Source versions and material actually consulted | [Source register](docs/SOURCES.md) |
| Full derivations and evidence map | [Supporting index](docs/supporting_proofs_index.md) · [Supporting papers](papers/supporting/) · [Proof fragments](proof_sources/) |
| Current claims and remaining research work | [Originality assessment](docs/originality_assessment.md) · [Research plan](docs/quantum_method_redesign.md) |
| Detailed results, limitations and study commands | [Research guide](docs/research_guide.md) |
| Data restrictions and permitted publication | [Data policy](docs/DATA_POLICY.md) |

## 3 Repository layout

| Directory | Contents |
| --- | --- |
| [src/](src/) | Pricing engines, study runners, numerical checks and figure generators |
| [checks/](checks/) | Tests and independent mathematical/numerical checks |
| [results/](results/) | Retained synthetic evidence, validation records and source ledgers |
| [protocols/](protocols/) | Registered study and figure protocols |
| [requirements/](requirements/) | Dependencies for the individual study environments |
| [docs/](docs/) | Research documentation and provenance records |
| [history/](history/) | Earlier sources, papers and checksum-bound study snapshots |
| [reviews/](reviews/) | Retained mathematical and editorial assessments |
| [public_data/](public_data/) | Permitted public reference material and its source licences |
| [tools/](tools/) | Layout mapping and reconstruction of the original execution paths |

The files in `src/` are research scripts, rather than a conventional installed Python package. Earlier studies retain their original relative paths and checksums; I use the reconstruction tool below to reproduce them from this organised layout.

## 4 Reproduce the studies

Use Python 3.13 and the dependencies specified for each study in [requirements/](requirements/). Reconstruct a fresh execution workspace, then run a representative pricing test:

```sh
python3.13 tools/reproduce.py --output /tmp/qbp-replay
cd /tmp/qbp-replay
python3.13 -m unittest discover -s checks -p 'test_capped_gbm.py'
```

The tool reconstructs the original file layout from [the recorded mapping](tools/repository_layout.json), preserving the paths used by frozen evidence. Run the detailed study commands in the restored `README.md`; the same research instructions are available in [the research guide](docs/research_guide.md). Choose a new output directory for each reconstruction and each study rerun.

## 5 Evidence and preservation

The reorganisation changes navigation and file locations, with the recorded mapping preserving the original research artefacts. Existing proofs, numerical findings, adverse comparisons, source-access limits and historical records remain available. No LSEG observations are distributed in this repository.
