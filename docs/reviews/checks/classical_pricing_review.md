# Independent review of complete classical finite-output pricing

## 1. Scope

I review the frozen deterministic pricing protocol and its 108 completed classical outputs. The protocol compares 273 inherited finite configurations, then chooses one prefix for each of 27 cases and four methods. Its scope is 99 per cent confidence for each price under ideal independent uniform finite input addresses. It does not assert that all 108 realised outputs have simultaneous 99 per cent coverage. The executed seeded NumPy streams are reproducible pseudorandom evidence; replay does not prove the independent random-bit premise of the confidence theorem.

The independent script is `checks/audit_classical_pricing.py`. It regenerates the tables, recomputes every candidate budget and operation profile, verifies deterministic prefix selection, and replays all selected sample streams through the earlier independent interpreter `audit_generalisation.vector_integer`. It accumulates first and second moments using NumPy object arrays whose products and reductions are Python integer operations. This moment reduction is independent of the production limb accumulator. The completed audit record is `results/classical_pricing_audit.json`.

## 2. Deterministic choice and fixed confidence budgets

For each candidate, I recompute the Bernstein and Hoeffding counts at 300 decimal digits from the certified finite-output variance upper bound, the exact normalised tolerance $t=0.005/H$ and the exact failure allocation $\alpha=0.01$. I verify the retained integer count and bound choice and require the entire independently computed failure bound to lie below the emitted failure upper bound and the allocation. The pointwise arithmetic and distribution certificates remain separate from that sampling guarantee.

The selection objective is the certified sample count multiplied by the retained conservative scalar-lane operation count. I verify its integer arithmetic directly. I then reconstruct every case-and-method candidate pool, identify the exact minimum objective, verify every recorded tied key, and check the stated tie order using table bytes and the full prefix key. These choices use deterministic model and compiler calculations before samples are generated. No realised sample variance or pricing observation enters the prefix or sample-count choice. The objective describes this implementation in the finite inherited pool; it is neither a global classical optimum nor a physical speed comparison with quantum gates.

## 3. Exact replay and sample moments

The replay reconstructs the seed sequence from the frozen common seed, the lexicographic case index and the fixed method index. It uses the recorded batch size and the graph's declared input insertion order. For every batch, it reproduces the name-tagged little-endian input byte stream and the little-endian output byte stream. Full SHA-256 comparisons supplement equality of exact accumulated sample moments. Four selected states from the first batch of every output are also evaluated by the scalar Python-integer graph.

For $N$ samples with output integers $a_i$, the independent reductions retain $S=\sum_i a_i$ and $Q=\sum_i a_i^2$ without a fixed-width squared-sum accumulator. I verify the exact reported mean $S/(N2^q)$, the exact sample-variance numerator $NQ-S^2$ and its denominator $N(N-1)2^{2q}$. The realised variance is a diagnostic of the recorded run. It does not replace the deterministic upper bound in the confidence budget.

## 4. Monetary interval derivation

Let $p$ be the exact finite fractional-output mean and $\overline X=S/(N2^q)$ the exact sample mean. Let $C$ be the continuous geometric correction, with retained centre $c$ and proved allowance $e_C$. Let $a$ and $b$ be the retained arithmetic and distribution currency bounds. On the sampling event, $H|\overline X-p|\leq0.005$. The deterministic certificate then gives

$$\left|V_{\rm undiscounted}-(c+H\overline X)\right|\leq a+b+e_C+0.005.$$

The code sets $e_C$ to at least the exact decimal $10^{-12}$, rather than using a potentially lower binary-float rendering of that allowance. The budget $a+b+e_C\leq0.005$ gives the intended total absolute currency tolerance before discounting. The nonnegative capped payoff lies in $[0,H]$, so intersecting the undiscounted interval with that range preserves coverage. Multiplication by the positive discount factor $e^{-r}$ preserves endpoint order. The retained models have nonnegative interest rates and unit maturity, so discounting cannot enlarge the intended absolute tolerance.

I compare the emitted lower and upper normalised sample-mean decimals with the exact rational mean. I reproduce the geometric correction centre from the model and cap in the producer's upward-rounding context. I separately recompute it in the default nearest-rounding context and retain their difference, which must be smaller than $10^{-60}$. The local context inherited by the correction generator can change its last retained digits; comparing its contexts as exact identical strings would be too strong. Its continuous accuracy rests on the separately proved correction allowance; an uncertified quadrature price is not used as a coverage test. I independently evaluate the endpoint formula using the actual retained centre at 300 digits and require the emitted endpoints to enclose that complete reference formula with its numerical allowance.

The displayed decimal midpoint need not equal the exact rational midpoint of the displayed endpoints. I therefore compare all displayed endpoint and midpoint decimals as exact fractions, and verify that the larger midpoint-to-endpoint distance is at most the reported radius. This explicitly checks the rounding of the displayed midpoint. Every reported radius must then be at most the exact fraction $1/100$. A realised interval's width and reproducible construction do not prove that it contains the unknown continuous price; that coverage statement follows from the stated sampling and deterministic certificates.

## 5. Cost interpretation and claim boundary

The audit charges a cold shared table regeneration separately from candidate preparation, random-stream regeneration, independent financial evaluation and Python object moment reduction. These timings are diagnostics for the stated runtime and hardware execution. The selected output's table-generation total is reconstructed from the same first-generation table records with each required table counted once. That accounting is different from repeatedly charging a shared bank to every price, and it is different from amortising setup away entirely.

The producer's selection-and-setup timer starts after loading its variance certificates. It consequently excludes generating those certificates and the earlier inherited precision certificates. The variance record retains a separate generation time. A complete workflow discussion must report that work separately; the pricing setup field alone is not the complete research or certification cost. The sample-count-times-operation objective explicitly excludes heterogeneous setup costs and does not purport to minimise full wall-clock work.

The finite laws, candidate scopes, all selected sample counts, full stream hashes and exact moment sums are auditable. The fixed confidence guarantee is a mathematical statement under the ideal independent input-law premise. The work demonstrates an implemented and checked classical comparator for the same finite target. It does not establish that this is the strongest continuous classical numerical pricing method, a theoretical novelty or a physical quantum advantage.

## 6. Completed audit findings

All 273 candidate profiles and budgets passed. The stated tie rules reproduced all 108 selected prefixes. I regenerated and checked 283 unique tables, replayed every one of the 49,788,913 sample outputs, matched all input and output stream hashes, and reproduced every exact integer sample aggregate. All 432 scalar cross-checks passed. All 108 monetary endpoint and midpoint-radius checks passed.

The independent cold table regeneration took approximately 50.0 seconds, the earlier integer interpreter replay approximately 13.3 seconds, and the independent Python object moment reduction approximately 2.83 seconds. Total audit time was approximately 69.6 seconds. These measurements describe this verification run, not an alternative pricing throughput promise or a physical quantum comparison. No unresolved mathematical or execution blocker remains in the reviewed pricing phase.
