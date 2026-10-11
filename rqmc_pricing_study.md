# 1 Scope and frozen targets

I price the 18 frozen four-asset and eight-asset contracts and retain the nine two-asset contracts as controls. I use the same finite financial functions, arithmetic certificates, rounded quantile tables and continuous-to-finite distribution allowances as the existing quantum comparison. My candidate pool comprises all 137 fully conditioned prefixes, covering both the generic and cap-aware financial kernels. I select one prefix per contract before generating production shifts. The selection minimises a prospective Bernstein replicate count times 64 times the conservative financial scalar-lane count. Its global range and moment priors are deterministic certificates. Neither an observed replicate mean nor an observed replicate variance enters this choice.

I retain the full candidate inputs, chosen graph nodes, table manifests, source hashes and deterministic tie-breaking rule in `results/rqmc_pricing_protocol.json`. The production seed is 2026100707. The distinct-seed preliminary implementation check uses 2026100701 and does not enter the production inference or graph selection.

# 2 Exact shifted batches and confidence

Each replicate contains 64 deterministic digital-net points. I use the previously retained author-distributed Joe–Kuo directions, truncate the integer words to each finite input address width, and apply one uniform bitwise exclusive-or shift per coordinate. A uniform shift makes every individual point exactly uniform on the finite input group, because translation by exclusive-or is a bijection. The points within a replicate remain dependent. Independent shifts make the replicate batch means independent and identically distributed, with the same expectation as the finite quantum target. I prove these statements and the resulting range and moment priors in `rqmc_engine_proof.tex`.

Every financial output is the exact integer (a), normalised as (a2^{-q}). A 64-point sum is at most (2^{36}). I retain exact batch sums, sum of batch squares and sum of squared differences between disjoint pairs. The paired quantity estimates the variance of a whole shifted batch. The confidence helper first obtains an upper variance bound from the lower tail of these bounded nonnegative paired observations, then applies a Bernstein mean bound. It does not assume a variance reduction proportional to the number of points. It uses no Student approximation and derives no coverage claim from the realised sample variance.

The checkpoints are 64, 128, 256, 512, 1,024, 2,048, 4,096, 8,192, 16,384, 32,768 and 65,536 shifts. At every checkpoint I allocate (1/2200) to the variance event and (1/2200) to the mean event. A union bound gives a total failure allowance of one per cent for one selected price over all eleven checkpoints. The coverage statement therefore permits stopping within this prescribed family. It does not give 99 per cent joint confidence for all 27 prices. The ideal exactly uniform independent shifts are a mathematical assumption. PCG64DXSM provides a reproducible execution, rather than a proof of that assumption.

# 3 Results and independent reconstruction

The production execution uses 786,432 random shifts and 50,331,648 financial evaluations. All intervals are nonempty. The reserved statistical allowance is 0.005 in undiscounted price units, and the retained deterministic allowance is at most 0.005. The realised experiment meets the reserved statistical allowance for 24 of 27 prices. The three contracts with common loading 0.05 and cap 25 reach the maximum without meeting this reservation. I retain those outcomes without extending the checkpoint family.

The final monetary interval uses the actual retained deterministic allowance, its geometric correction and the discount factor. It has radius at most 0.01 for 26 prices. Two of the three failed statistical reservations still meet this complete monetary tolerance because their deterministic allowance is smaller than the reservation. The four-asset loading-0.05 cap-25 contract remains wider, with radius 0.0108463574802773. I distinguish a passed statistical reservation from a passed complete monetary interval throughout the retained results.

| Assets | Prices | Statistical reservation passed | Complete monetary tolerance passed | Financial evaluations |
|---:|---:|---:|---:|---:|
| 2 | 9 | 8 | 9 | 15,564,800 |
| 4 | 9 | 8 | 8 | 16,596,992 |
| 8 | 9 | 8 | 9 | 18,169,856 |
| Total | 27 | 24 | 26 | 50,331,648 |

The final execution time is 5.130 seconds on the retained local backend. All-candidate preparation, including unique table generation and evaluator preparation, takes 43.685 seconds and is recorded separately. I repeat the same protocol after strengthening the confidence helper's input checks; all financial batch sums, stream hashes, confidence checkpoints and monetary intervals reproduce exactly. Timings are implementation diagnostics. They are not comparable units to a logical T count or a physical quantum runtime.

The independent audit uses the earlier separate integer interpreter and a sequential Sobol recurrence, rather than importing the new execution engine or its vector Gray-code implementation. It replays all 50,331,648 financial outputs, all 786,432 exact batch sums, every retained checkpoint moment, both stream hashes and the prefix selection. It checks 3,544 financial points separately using the Python arbitrary-precision scalar graph and reconstructs monetary interval endpoints at 300 digits. All checks pass. Seeded reconstruction verifies these computations, not statistical independence or coverage.

# 4 Matched fixed-budget independent sampling

I also calculate fixed-budget independent-sampling plans over the same 137 prefixes from the same certified range and global moments. These plans use the minimum of the range-aware Hoeffding and Bernstein budgets with a one per cent per-price failure allowance and the same statistical reservation. Their selection uses only deterministic retained priors. I compute the plans after the digital-shift experiment, but no realised digital-shift output enters the calculation. They are guaranteed planning counts under the ideal independent-address premise, and are not newly executed independent-sampling prices.

The realised financial evaluator work of the shifted-net protocol is larger than the selected independent-sampling plan in every case, by factors 2.5436 to 133.2033. This is an adverse result for this particular fixed-size shifted-net confidence construction. It does not prove that randomised quasi-Monte Carlo is inferior in general. Its rigorous confidence allowance and checkpoint penalty can overwhelm the observed reduction in replicate variance. I do not present this experiment as a classical speedup or use its apparent numerical agreement to certify missing cases.

`results/rqmc_work_comparison.json` retains all fixed-budget candidate plans and matched minima. I charge independent-address generation, exact chunked mean reduction, digital-net construction, shift generation, point exclusive-or operations and checkpoint moment reductions separately. Financial scalar-lane operations, random word requests, memory bytes, special-function calls and arbitrary-precision calls have different units. I do not add those units into a fictitious physical runtime. The existing 108 independent-sampling executions remain unchanged.

# 5 Interpretation and source scope

This phase provides a complete independently replayed comparator for the higher-dimensional finite quantum targets and retains its failures. I do not identify the finite-group unbiasedness proof, generic bounded-replicate confidence argument or use of digital shifts as a new literature contribution. Any project-specific contribution must arise from checked financial arguments and comparison with established methods, rather than from this implementation or the volume of its numerical evidence.

I use the retained direction file `new-joe-kuo-6.21201`, author-distributed `sobol.cc` dated 21 October 2008, and BSD licence, copyright Frances Y. Kuo and Stephen Joe 2008. Exact hashes appear in the frozen protocol. The previously consulted author page is <https://web.maths.unsw.edu.au/~fkuo/sobol/>. I do not imply that the associated journal papers were read for this phase. No new data purchase, subscription or entitlement is used, and every financial contract in this experiment is synthetic.
