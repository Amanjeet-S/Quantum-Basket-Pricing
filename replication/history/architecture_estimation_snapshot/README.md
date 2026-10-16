# Certified quantum encoding of bounded basket derivatives

## 1 Research question and status

I investigate when contract-dependent slope and curvature certificates permit cheaper complete coherent arithmetic for geometrically controlled capped basket derivatives, and which financial and architectural conditions prevent that saving. The proposed interval advantage has failed the stronger-comparator test. The probability identities, geometric control, Gaussian conditioning, derivative-based precision selection and interpolation are established ingredients. A subsequent complete-estimation comparison now implements a quantum confidence controller and stronger coherent lookup modes. An original positive literature contribution and quantum advantage remain unverified.

I replaced the earlier direction and precision selector with one concrete evaluator. For positive normalised weights, arithmetic basket A and geometric basket G satisfy A >= G. I write the capped payoff as F(x) = min(max(x - K, 0), H), retain the analytical correction E[F(G)], and estimate the nonnegative residual E[F(A) - F(G)]. A uniformly sampled strike turns that residual into a conditional Gaussian interval probability. Under the established common-loading factor, both interval boundaries are logarithms and no iterative root solver is required.

The initial financial object is a terminal equity basket call spread under stipulated correlated GBM. REITs remain a possible later LSEG parameter scenario. Asian monitoring is an extension only after the core evaluator survives its benchmarks. Neither asset substitution nor increased nominal dimension is a contribution.

## 2 Completed feasibility work

The [complete architecture and estimation study](architecture_estimation_study.md) now retains 819 circuit profiles from 273 independently certified precision candidates across three lookup modes, plus 60 rejected transfers. I derive a financial residual-amplitude prior, apply it to both classical and quantum estimation, and prove a modified iterative quantum controller with explicit monetary accuracy, confidence and worst-case resource caps. The [expanded proof supplement](generalisation_proofs.tex) includes the financial derivations, corrected lookup Kraus branches, complete gate costs, adaptive coverage and finite classical postprocessing. Its references use the author–date style of the paper in “Map quant finance prerequisites”.

Under the measured four-T lookup mode, cap-aware evaluator savings remain in two of 27 cases: 5.840 and 0.647 per cent. The former four-asset cap-aware saving disappears. Complete quantum-budget selection also retains two paired reductions, 5.838 and 0.647 per cent, against the generic full evaluator. Direct evaluation still wins three cases; interval conditioning wins none in the inherited pool. The prior improves the same controller's complete worst-case T budget by factors 1.82 to 36.15 across all 819 profiles. These are ideal logical bounds, rather than physical runtimes or quantum advantage. The equally supplied classical prior yields certified finite-output sampling budgets. The original 108 evaluator rows remain intact, with all same-prefix count regressions passing.


The [broader certificate study](generalisation_study.md) retains 27 frozen synthetic financial cases and 108 complete evaluator configurations. A cap-aware full-conditioning certificate reduces the assigned T-count in three paired cases by 13.108, 3.870 and 1.033 per cent; all other 24 cases retain no saving. The [proof supplement](generalisation_proofs.tex) derives the financial and distribution bounds, weak-loading sharpness and exact condition under which one fewer lookup address bit actually reduces emitted arithmetic cost. A coarser table can increase multiplier cost, so a tighter error bound alone does not establish a saving. The [resource record](results/generalisation_resources.json), [independent uniform certificates](results/generalisation_independent_certificates.json), [production replay audit](results/generalisation_independent_audit.json) and [classical reference](results/generalisation_classical_reference.json) preserve the evidence.

The original broader shortlist has one interval winner. A separately labelled [post-study matched-prefix check](results/generalisation_matched_interval.json) constructs a certified cheaper full kernel at that same prefix and at every other retained interval prefix. The original 108 rows remain intact. This identifies heuristic search sensitivity, rather than an intrinsic interval advantage. Financial cases were frozen before the run; a shared distribution-bound amendment and post-study diagnostics are transparently recorded.

The [stronger conditioning benchmark](conditioning_benchmark.md) reports five complete constructed evaluators for three caps at two monetary targets. Specialised full conditioning has the lowest assigned T count in all six cases. At cap £10 and target £0.01, it uses 597,478 T gates against 672,210 for the comparably specialised interval evaluator and 1,207,388 for direct evaluation. These are constructor-specific logical costs with half the tolerance reserved for unimplemented statistical estimation. The [current resource record](results/conditioning_resources.json) and [independent audit](results/conditioning_independent_audit.json) preserve the evidence. The [initial study](reversible_encoding_study.md) and its results remain as the superseded comparison.

The [research design](research_design.md) sets the comparison and decision rules. The [mathematical note](mathematical_note.md) checks the identities and monetary bias certificates. The [source record](SOURCES.md) distinguishes published versions from the versions actually consulted. The [pilot report](results/pilot_report.md) describes twelve fixed synthetic scenarios and an adverse variance example. The [machine-readable results](results/pilot.json) retain inputs, seeds, diagnostics and a script hash.

The small state-vector check verifies finite-output Born probabilities and Grover interference. Its host-evaluated probability truth table is explicitly recorded. It does not implement scalable reversible arithmetic, demonstrate a financial accuracy target or establish a speedup. The pilot standard deviations and batch standard errors are screening diagnostics, not certified variance bounds or pricing confidence intervals.

## 3 Reproduction

Use Python 3.11 or later with NumPy. In a virtual environment, install [requirements.txt](requirements.txt), then run:

```sh
python pilot.py
python audit_certificates.py
python -m unittest discover -s checks -p "test_*.py" -v
python resource_study.py --certify-quantiles --small-checks --output results/conditioning_resources.json
python checks/reference.py
python reflections.py
python checks/audit_composition.py
python checks/shadow_audit.py --config results/conditioning_resources.json --output results/conditioning_independent_audit.json --samples 262144 --certify-quantiles --new-variant-gates
python checks/audit_kernel_dominance.py
```

The new resource study uses only the Python standard library. Its independent reference and replay audit require NumPy. It emits complete primitive-block gates, checks function domains and price-bias budgets, and independently encloses every retained quantile and function constant. The two-asset resource instance is distinct from the earlier four-asset pilot.

For the broader frozen family, run:

```sh
python generalisation.py
python checks/audit_generalisation_certificates.py
python checks/classical_reference_generalisation.py
python checks/audit_generalisation.py --certify-quantiles
python checks/audit_precision_effects.py
python checks/audit_matched_interval.py
python checks/audit_generalisation_compiler.py
```

For the subsequent complete-estimation layer, run:

```sh
python precision_candidates.py
python architecture_study.py
python checks/audit_architecture.py
python checks/audit_primitive_formula.py
python estimation_prior.py
python checks/audit_estimation_prior.py
python complete_estimation_study.py
python checks/audit_complete_estimation.py
python build_estimator_study.py --simulation
```

The intermediate architecture selector minimises Grover-iterate T count. The complete selector then examines every feasible prefix at the specified confidence and retains all method ties. Simultaneous confidence applies to 27 deterministically selected prices under one chosen implementation, rather than every candidate output jointly. The historical assigned three-layer Toffoli depth remains explicitly separate from the audited expansion's conservative five-layer bound. Classical and quantum cost units remain separate.

The financial protocol is immutable in ordinary reproduction. The model generator refuses to overwrite it. Free author-distributed Sobol directions, algorithm source and the BSD licence are retained in public_data. The reference uses digital shifts, with replicate dispersion labelled as a diagnostic. It does not require SciPy or new data access. Retained historical resource hashes refer to their original source state; the broader outputs record the extended sources.

The default pilot run uses eight independent batches of 4,096 paths for each of twelve scenarios. It writes results/pilot.json. Runtime is recorded separately from the numerical outputs and is not a quantum cost estimate. The exact script and NumPy version are retained in the output.

## 4 Next decision

I have completed a finite comparison of financial certificates, three coherent lookup constructions and a fully specified confidence estimator. The original precision gain is architecture dependent: a cheaper lookup removes one of the three earlier strict cap-aware savings. The original interval winner disappears when every method inherits the same certified opponent prefixes. These findings support a precise project-specific cost argument and expose search sensitivity. They do not establish a general optimal algorithm or theoretical priority.

I will next implement a credible complete classical conditional-pricing comparison and examine stronger published quantum estimators or known-function constructions under matched monetary error. The source ledger identifies the mathematical antecedents and the exact versions consulted. A simulation coverage rate, a new implementation or a useful bound alone does not demonstrate a positive literature contribution. No market validation or physical quantum speedup has been established.

## 5 Data and access

All current experiments use explicitly stipulated synthetic parameters. No LSEG observations, credentials, purchased data or publisher PDFs are included. Existing LSEG access may later support market-informed scenarios after a field and adjustment audit. Raw and reconstructible provider observations remain local and excluded from Git, including when the remote repository is private.

## 6 Keywords

Quantum computational finance; basket derivatives; geometric control variates; conditional expectation; reversible arithmetic; monetary error certification.
