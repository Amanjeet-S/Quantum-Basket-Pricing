# Financial certificates, coherent lookup and complete estimation cost

## 1 Research question and keywords

I investigate the following question: at absolute price error 0.01 and 99 per cent confidence, when do contract-dependent conditioning certificates and a proved positive-residual prior reduce the complete logical quantum cost of capped basket pricing after each financial opponent receives the same certified precision candidates under stronger lookup architectures?

Keywords: mathematical finance; quantum amplitude estimation; capped basket derivatives; conditional expectation; coherent lookup; confidence coverage; finite precision.

I retain the 27 frozen synthetic contracts, the original 108 configurations and their historical outputs. I change the computational comparison by inheriting opponent prefixes, implementing stronger lookup modes and completing the estimator. This is a subsequent study, rather than a relabelling of the original search as exhaustive. I use no new market data, purchased access or provider demonstration sample.

## 2 A financial prior with an explicit finite-output transfer

For positive normalised weights, arithmetic basket A, geometric basket G and payoff F_H(x)=min((x−K)_+,H), the residual F_H(A)−F_H(G) is nonnegative and bounded by H. With weighted residual loading R_G and dispersion D=Σ_i w_i‖R_i−R_G‖², matched risk-neutral drifts give E[A]=K exp(r−d) and E[G]=K exp(r−d−D/2). Monotonicity and the payoff's unit Lipschitz constant therefore bound the residual mean by B_AG=K exp(r−d)(1−exp(−D/2)).

The common-factor conditioned payoff gives an additional bound. The logarithmic basket gap J=m_A−m_G is nonnegative. Jensen's inequality and Gaussian exponential moments give E[J]≤J_0=D/2+log(Σ_i w_i exp(−(R_i−R_G)·R_G)). If L_cap is the proved global slope bound of the normalised conditioned kernel, the currency residual mean is at most B_*=min(H,B_AG,H L_cap J_0). The [proof supplement](../papers/generalisation_proofs.tex) derives both forms, verifies integrability and transfers the bound to each implemented finite flag probability. I use u=min(1,(B_*+B_arithmetic+B_distribution)/H), with explicit numerical enclosures, rather than a fitted mean or pilot variance. The separately added geometric correction enters final price error, rather than the residual flag probability.

For a finite output X in [0,1], E[X²]≤E[X]=p, so Var(X)≤p(1−p). The prior gives variance at most u(1−u) for u≤1/2 and at most 1/4 otherwise. I apply this same prior to the classical comparator. A fully derived Bernstein bound and an independently checked Hoeffding alternative supply certified classical sample budgets. Those budgets count finite-function evaluations; they do not represent CPU time or quantities commensurate with quantum T gates.

These inequalities are elementary consequences of established financial moments, Jensen's inequality and concentration theory. Geometric-average variance reduction has published antecedents, including Kemna and Vorst (1990), while variance-aware quantum financial estimation is established in Rebentrost et al. (2018). I do not infer novelty of the precise residual bound from a keyword search or from having written its derivation.

## 3 A coherent lookup construction with a branchwise proof

I emit three lookup modes: the original exact seven-T compute/erase pair, a four-T clean logical AND followed by its unitary adjoint, and a four-T clean AND followed by corrected measurement erasure. Their respective costs per visited selector node are κ=14, 8 and 4 T gates. The last construction follows the published unary-iteration and temporary-AND foundations of Babbush et al. (2018) and Gidney (2018). I do not assign four-T costs by changing a label on the original Toffoli tape.

The proof gives every gate of a phase-free clean AND isometry. Measuring its target in the X basis contributes phase (−1)^(zab), and conditional CZ cancels that phase before conditional reset. Each corrected Kraus branch is exactly 1/√2 times deletion of the computed AND, including external entanglement. Induction through the pruned prefix traversal proves that each complete record branch is 2^(−I/2) times the same logical dirty-output XOR, with input-independent record probabilities. A fresh instrument implements the logical adjoint. This establishes the coherent channel required for amplification, rather than only a correct basis-state output distribution.

The architecture changes the precision question. At fixed q, register width w, function span and remainder width f, a one-bit table coarsening saves lookup nodes ΔI but increases multiplication work. Its complete two-kernel flag saving is 4[κΔI−7(8w+12f+6)] T gates. It is positive exactly when κΔI>7(8w+12f+6). At κ=4 this requires ΔI>14w+21f+10.5. Cheap lookup therefore makes a coarser function table harder to justify. This project-specific conjunction is an architecture-dependent argument, rather than a new generic lookup theorem.

I separately charge measurements, worst-case conditional CZ, reset, uniform preparation, complete arithmetic inverses, comparisons, reflections, table-programme storage and logical qubits. Arithmetic retains the exact seven-T Toffoli convention. The old assigned three-layer depth is distinguished from a separately emitted phase-audited sequence and its conservative five-layer depth bound. Neither T depth nor measurement depth is a physical elapsed-time prediction. Constant-prefix specialisation, general arithmetic improvements and alternative published QROM trade-offs could reduce costs further.

## 4 A fully specified quantum estimator

I implement a conservative iterative amplitude estimator for the actual finite flag probability. The established amplification identity of Brassard et al. (2002) gives p_k=sin²((2k+1)θ) for p=sin²θ. The estimator uses fresh independent shots, exact rational quadrant decisions, outward scalar inverse calculations, and a separately proved three-or-five power rule. Grinko et al. (2021) is the closest published iterative mechanism; its author Appendix B already discusses factors three, five and seven. I therefore do not present those factors as an invention.

The local proof handles stationary endpoints and inversion padding, derives conditional Hoeffding coverage under adaptive powers, and gives deterministic round and Grover-query caps. An empty confidence intersection is explicitly recorded as failure. It is not repaired into an apparently certified price. The financial prior also permits a largest safe initial odd multiplier without a pilot and bounds the amplitude derivative over every retained interval. I prove the resulting tightened caps rather than transferring the published IQAE coefficient to a modified schedule.

The statistical currency tolerance is exactly 0.005. Together with the certified deterministic currency allowance of at most 0.005, the midpoint estimate has discounted absolute price error at most 0.01 with probability at least 0.99 under the stipulated positive rate and ideal logical-shot assumptions. The separate simultaneous scope allocates failure probability 0.01/27 to each of 27 preselected prices. It does not certify all retained candidate outputs jointly at 99 per cent. The selection uses deterministic model and cost records, rather than observed shot outcomes.

If S is the worst-case shot cap and Q is the Grover-call cap, the complete T bound is S C_A+Q(2C_A+C_reflection). Lookup measurements and corrections are multiplied by S+2Q, and flag measurements are added separately. Sequential shots reuse qubits. No QFT, arbitrary quantum rotation or synthesis bias is hidden in this construction. The declared ideal gate model excludes physical faults; a physical backend needs its own accuracy and correlation argument.

## 5 Retained numerical comparison

<!-- BEGIN NUMERICAL RESULTS -->
I retain 273 feasible financial prefixes, 60 rejected transfers and all 819 architecture profiles. Every original seven-T winner appears in the inherited pool, and all 108 same-prefix count regressions pass. The intermediate selector minimises emitted Grover-iterate T count. The complete selector independently minimises the proved full confidence budget across every retained prefix, preserving exact ties.

The cap-aware and generic full-conditioning comparison is architecture dependent. At λ=0.05 and H=25, the table shows evaluator savings followed by complete single-price confidence-budget savings, each in per cent. The three columns charge κ=14, 8 and 4 T gates per visited lookup selector node. All unlisted cases have zero paired saving.

| Assets | Exact7T lookup pair: evaluator / complete | Unitary8T lookup pair: evaluator / complete | Measured4T lookup: evaluator / complete |
| ---: | ---: | ---: | ---: |
| 2 | 13.108 / 13.106 | 9.952 / 9.950 | 5.840 / 5.838 |
| 4 | 3.870 / 3.869 | 0.000 / 0.000 | 0.000 / 0.000 |
| 8 | 1.033 / 1.033 | 0.892 / 0.892 | 0.647 / 0.647 |

The measured lookup eliminates the former four-asset strict cap-aware saving. Its two surviving paired reductions do not imply two unique overall wins: the eight-asset case is still cheaper by direct evaluation. Across 27 complete selections, exact7T has 22 generic/cap-aware ties, two unique cap-aware wins and three direct wins. Unitary8T and measured4T each have 23 generic/cap-aware ties, one unique cap-aware win and three direct wins. Interval conditioning attains no minimum in any of these complete selections. The 81 matched full/interval cost identities also pass, with a minimum same-prefix gap of 42,408 T gates. That dominance concerns interval versus full conditioning, rather than full versus direct evaluation.

Across all 819 profiles, the financial prior reduces the same controller's complete worst-case T bound by factors 1.82 to 36.15 relative to prior one. The classical comparator receives the same prior and requires between 6,744 and 11,495,538 certified finite-function samples for the single-price target. These quantities have different units. For the best measured-lookup implementation in each case, complete single-price T bounds range from 6,494,832,000 to 79,664,194,275,306. Large logical budgets remain even after the proved reductions; no physical speedup follows. The 540 ideal-binomial diagnostic runs all covered their stipulated probabilities, but this finite observation neither proves the confidence theorem nor represents quantum hardware execution.

<!-- END NUMERICAL RESULTS -->

## 6 Verification and reproduction

The [resource candidates](../results/architecture_resource_candidates.json), [complete estimation budgets](../results/complete_estimation_resources.json) and [independent complete audit](../results/complete_estimation_audit.json) retain source hashes, monetary and confidence allocations, all candidate rows and separate resource units. The [lookup audit](../results/architecture_audit.json) checks complex Kraus matrices for small complete lookup maps, selected branches for a larger map, dirty outputs, coherent enables and coarse complete financial compositions. Its general branchwise induction supports the channel claim beyond those finite checks. The [estimator checks](../results/estimator_checks.json) replay small complete preparation, Grover and controlled-Grover state vectors and verify the exact Toffoli phases. Numerical simulation uses ideal binomial draws and is explicitly labelled; its observed coverage cannot replace the proof.

I reproduce the new layer with Python and the existing NumPy dependency:

```sh
python precision_candidates.py
python architecture_study.py
python checks/audit_architecture.py
python checks/audit_primitive_formula.py
python estimation_prior.py
python checks/audit_estimation_prior.py
python complete_estimation_study.py
python checks/audit_complete_estimation.py
python -m unittest discover -s checks -p 'test_*.py' -v
```

The candidate generator first checks the immutable historical source state and re-certifies target-method arithmetic and distribution budgets independently. The reference and quantum state-vector checks require NumPy. The primary resource constructors and confidence controller otherwise use the Python standard library. The native LaTeX editor compiles the standalone proof source independently of a terminal TeX installation.

## 7 Contribution boundary and next research decision

The concrete argument is the interaction between financial residual bounds, confidence-certified amplification and architecture-dependent precision costs. The retained implementation can test that argument without a paid dataset. Its project-specific derivation and checked construction establish a precise finite comparison. I have not established the absence of its components or their combination from prior literature, nor superiority over the strongest published estimator.

I must still compare a complete classical conditional-pricing routine and credible stronger quantum estimation or known-function implementations. O'Brien and Sünderhauf (2025), McArdle et al. (2026) and Stamatopoulos and Zeng (2024) motivate relevant unimplemented alternatives. Their different state-preparation or approximation assumptions prevent assigning their reported function-only costs directly to this computed payoff graph. A complete classical cost comparison must include setup, numerical certification and sampling or deterministic integration on matched financial accuracy. Logical T gates and classical sample counts alone cannot demonstrate quantum advantage.

The financial models remain stipulated GBM scenarios. REIT parameters from existing authorised LSEG access could later provide a transparently labelled scenario, after a field and adjustment audit, but could not turn synthetic pricing into independent validation of real-estate dynamics. Asset substitution and more nominal dimensions do not themselves create a positive literature contribution. I will retain null and adverse findings if stronger comparisons remove the apparent improvement.

## 8 References

Babbush, R., Gidney, C., Berry, D. W., Wiebe, N., McClean, J., Paler, A., Fowler, A. and Neven, H. (2018). Encoding electronic spectra in quantum circuits with linear T complexity. *Physical Review X*, 8(4), 041015. [DOI](https://doi.org/10.1103/PhysRevX.8.041015).

Brassard, G., Høyer, P., Mosca, M. and Tapp, A. (2002). Quantum amplitude amplification and estimation. *Contemporary Mathematics*, 305, 53–74. [DOI](https://doi.org/10.1090/conm/305/05215).

Gidney, C. (2018). Halving the cost of quantum addition. *Quantum*, 2, 74. [DOI](https://doi.org/10.22331/q-2018-06-18-74).

Grinko, D., Gacon, J., Zoufal, C. and Woerner, S. (2021). Iterative quantum amplitude estimation. *npj Quantum Information*, 7, 52. [DOI](https://doi.org/10.1038/s41534-021-00379-1).

Kemna, A. G. Z. and Vorst, A. C. F. (1990). A pricing method for options based on average asset values. *Journal of Banking & Finance*, 14(1), 113–129. [DOI](https://doi.org/10.1016/0378-4266(90)90039-5).

McArdle, S., Gilyén, A. and Berta, M. (2026). Quantum state preparation without coherent arithmetic. *Physical Review Letters*, 136, 240603. [DOI](https://doi.org/10.1103/ntvs-c48s).

O'Brien, L. and Sünderhauf, C. (2025). Quantum state preparation via piecewise QSVT. *Quantum*, 9, 1786. [DOI](https://doi.org/10.22331/q-2025-07-03-1786).

Rebentrost, P., Gupt, B. and Bromley, T. R. (2018). Quantum computational finance: Monte Carlo pricing of financial derivatives. *Physical Review A*, 98, 022321. [DOI](https://doi.org/10.1103/PhysRevA.98.022321).

Stamatopoulos, N. and Zeng, W. J. (2024). Derivative pricing using quantum signal processing. *Quantum*, 8, 1322. [DOI](https://doi.org/10.22331/q-2024-04-30-1322).
