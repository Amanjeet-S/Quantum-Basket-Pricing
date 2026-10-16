# Independent confidence review

## 1. Mathematical findings

The pair-based variance inversion and the mean Bernstein radius are correct. The factor one half in the squared pair difference is necessary and is retained. The pair variables are independent because the pairs are disjoint and the random shifts are independent. The observed variance estimate and the sample mean can share observations. Coverage follows on the intersection of a variance confidence event and a mean deviation event constructed with the true variance, so independence between these two events is unnecessary.

The deterministic prior is a bound on the variance of an entire shifted batch. Jensen's inequality transfers a global second-moment bound to this batch. It supplies no variance reduction rate in the 64-point batch size. The implementation correctly treats R, rather than 64R, as the number of independent observations. All eleven prescribed checkpoints can share earlier observations because the final simultaneous claim uses a union bound. The exact allocation 22/2200 is 0.01 for one selected price. It does not establish a joint 99 per cent statement for all 27 prices.

The correctly rounded 160-digit Decimal quotient and logarithm premise is sufficient for the stated 10^-120 enclosing increment on the declared argument range. The integer square-root enclosure is exact and does not overflow because its calculations use arbitrary-precision integers. Empty intersections are explicitly retained as confidence failures and provide no estimate.

These mechanisms have established statistical antecedents. The self-contained derivation supplies a checked certificate for the experiment, rather than priority for a generic confidence theorem.

## 2. Independent retained evidence

The auditor imports no project confidence helper. It reconstructs every sum, square sum and paired difference-square sum from the retained integer batch sums. It checks all 27 selected prices and 227 stored checkpoints, including stopping decisions, failure allocations, exact rational endpoints and the absence of an earlier successful checkpoint. A separate 300-digit evaluation of each logarithm and square root confirms that the retained enclosures are outward. All checks passed.

An additional exact audit enumerates every aggregated pair type (00, mixed, 11) for Bernoulli batches at R=64 and R=128, for p=a/32 with a=0,...,32. The multiplicity 2^mixed restores both orientations of each mixed pair. Its probability arithmetic uses integer numerators and the exact denominator 32^R. All 66 laws normalise exactly. All pair failure probabilities satisfy 1/2200, and all union failure probabilities satisfy 1/1100. There are 89,298 law and pair-type probability terms. Sixty-two laws have a strictly positive final interval-failure probability, and 30 have a strictly positive pair-bound failure probability, so the examples check nontrivial tail events. The largest observed pair failure probability is approximately 3.302e-8 and the largest pair-or-interval failure probability is approximately 2.008e-6. These finite examples support the implementation, rather than prove coverage for arbitrary bounded batch distributions.

Fourteen portable unit and property tests pass. They check exact square roots across tiny, near-square and very large inputs, the logarithm enclosure margin, binary range scaling, correct use of the variance prior, the zero-range case, explicit empty intersections, odd replicate counts, large integer statistics, invalid ranges and priors, the finite checkpoint allocation, rejection of non-integral and boolean counters, necessary aggregate moment inequalities, and preservation of exact digit-string statistics.

## 3. Verified input-contract strengthening

The retained production statistics are valid. The revised confidence function now rejects non-integer counters before conversion, including floating-point, fractional and boolean counters. It requires the exact necessary inequalities S2 <= rho D S and Q <= S2. These checks strengthen error reporting for malformed external inputs without changing valid retained outputs. The range inequality follows from T_j^2 <= rho D T_j. The pair inequality follows from (T_a-T_b)^2 <= T_a^2+T_b^2 for nonnegative integer batch sums. Four additional tests verify these guards and ensure valid digit-string statistics yield the same output as integer statistics.

The final retained checkpoint audit checks every source hash against the actual files and against both the protocol and results ledgers. No repository files or Git state were changed in this independent review.
