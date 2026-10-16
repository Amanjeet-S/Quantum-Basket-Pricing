# Certified curvature corrections for quantum basket pricing

## 1 Research question and result

I develop coherent quantum derivative-pricing methods under an established classical no-arbitrage GBM model. The [main paper](generalisation_proofs.tex) contains the essential proofs and focused appendices. The [supporting index](supporting_proofs_index.md) maps full derivations, exact archives, audits and earlier investigations.

The current method evaluates payoff-boundary curvature analytically, removes the strike auxiliary and supplies a positive second-correction encoding for convex piecewise-affine payoffs, including ordinary basket calls. I certify shared sharp capped-payoff norms and a financial-grid transfer without a continuous basket price as an input. Four complete quantum evaluators receive the same controls, moments, financial error and confidence.

Raw direct second-correction pricing is the cheapest quantum second-correction implementation in all eight active cases. It beats the positive first correction in five cases, with reductions of approximately 13.9 to 91.2 per cent. The new boundary method also wins against that first-correction opponent in those five cases but loses to raw direct evaluation. All twelve prices are resolved by strong classical evaluators. A continuous convex screen also resolves four ordinary-call models, including four- and eight-asset cases, without a grid or sampling. These are conditional logical-resource results, without practical quantum advantage or established theoretical priority.

## 2 Financial scope

The frozen family has 27 original capped basket contracts and 9 ordinary calls across nine models, giving 36 cases. All 144 base arithmetic graphs pass. Complete financial transfers and pricing plans cover 12 two-asset cases; the other 24 cases retain arithmetic evidence and explicit financial limitations. Convex finite piecewise-affine payoffs are covered by the analytical boundary theorem, without claiming that every payoff has an implemented circuit or empirical validation.

All data are exact synthetic inputs. Free data and existing authorised LSEG or institutional access remain the only permitted acquisition routes. No new purchase, OptionMetrics sample or paid entitlement is assumed.

## 3 Reproduction and retained evidence

The frozen comparison retains 173 candidate programs, 23,772 allocation decisions, 338 scalar certificates and 538 complete strategies. The independent audit records distinguish mathematical certificates, small-width quantum circuit fixtures, large-width emitted counts and classical price evaluation. No physical quantum pricing experiment was executed.

From this repository, using the existing Python-Flint environment for outward arithmetic:

```sh
python boundary_precision_selection.py
python boundary_classical_pricing.py
python checks/audit_boundary_derivative_certificate.py
python checks/audit_boundary_value_oracle.py
python checks/audit_boundary_production_profiles.py
python checks/audit_boundary_grid.py
python checks/audit_boundary_auxiliary.py
python checks/audit_boundary_classical_pricing.py
python checks/audit_boundary_convex_screening.py
python checks/audit_boundary_precision.py
```

The final precision and confidence audits are listed in the supporting index. Existing fundamental compiler, estimator, Fourier and prior evidence remains unchanged. Setup for coefficients, quantile cells, derivative certificates, moments and scalar tables is charged separately. The [source ledger](boundary_sources.md) records published status and the versions actually consulted; the [research design](research_design.md) fixes the comparisons and admissibility policy.

## 4 Repository state

The repository remains private. No new commit is created in this phase. HEAD remains 1a25d67 after the previously authorised removal of the unpublished latest commit. The concise paper and its preceding documentation are snapshotted verbatim in history/before_boundary_correction; all earlier sources and results are retained.
