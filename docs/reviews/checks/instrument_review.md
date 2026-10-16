# Independent logical-instrument and financial-prior review

## 1. Scope and versions

I reviewed /private/tmp/qbp-next-step/architecture_precision_proof.tex, financial_estimation_proof.tex, and /private/tmp/qbp-next-architecture/architecture.py on 5 October 2026. I also read the architecture audit implementation and retained results. Keywords: clean logical AND, quantum instrument, coherent QROM, Kraus branches, financial residual, precision cost.

The reviewed implementation hash is c9fbf6169b72d335b588a626300ce2caea14055e2587eb0870e8074412ec751a. The architecture proof hash is 6ef1a78825421050d83ddf69e99eeea5ae431c30a9026be963129b86edd551db and the financial proof hash is 6d5f99e97d0c9d0d6c26f488e50df51913792418f8c4bba950e1d9b7847adc55. The retained architecture_audit.json hash is 494478973a0d9b65ae327d290b7c32609db108c7393e948df7200408a553f661.

I found no blocking mathematical error in the logical-instrument, financial-prior or architecture-coarsening arguments, subject to their declared clean-work and ideal logical-operation assumptions. One retention issue remains: architecture_audit.json records implementation hash 5981e5c190a0a56916c91a97b179fce0cc16b2c3630102f89021ad84105e0998, preceding the currently reviewed implementation. The visible change shares pure-arithmetic cache keys across modes while retaining architecture-specific keys for lookup-containing blocks. It appears benign because the pure blocks contain the same exact gates and no measurements, but the audit must be regenerated against the frozen source before claiming a matching retained audit.

## 2. Four-T AND and corrected erasure

I independently checked the chronological sequence H, T, CX from b, T†, CX from a, T, CX from b, T†, H, S†. With a clean target, the three control pairs other than 11 finish at target zero with phase one. For 11, the state immediately before the final Hadamard is proportional to (i,−i), so the Hadamard gives i times target one and S† cancels that phase. The proof's intermediate vectors are consistent with the gate order. The clean-target map is therefore exactly |a,b,0⟩ to |a,b,ab⟩, including coherence with a reference. This does not establish a four-T Toffoli on a dirty target. The reversed adjoint sequence is an exact four-T inverse on the computed image.

For erasure, H followed by outcome z supplies the phase (−1)^(zab) and scalar 1/√2. Conditional CZ supplies the same phase again, and conditional X resets the measured target. Thus each corrected Kraus branch is precisely 1/√2 times removal of the encoded AND. Orthogonality of the preserved controls gives input-independent probability one half. The proof handles arbitrary data and external references, rather than only classical basis probabilities.

## 3. Recursive lookup, coherent enable and fresh adjoint

The recursive plan computes child=prefix AND address_bit. Its positive subtree preserves that child. The intervening CNOT changes it to prefix AND NOT address_bit, and the second CNOT restores the original conjunction after the negative subtree. Lower selectors start and end clean. Addresses and incoming prefixes are unchanged. Leaf writes are CNOTs into the arbitrary dirty output. The induction therefore satisfies the erasure premises at every visited internal node.

Pruning depends on the fixed classical table, so the number and placement of erasures do not depend on the quantum address or enable. A coherent enable is simply the preserved incoming root prefix. Even its zero component passes through the same corrected measurements. Every record branch is consequently 2^(−I/2) times the same enabled XOR map, with clean private selectors and no relative phase. This equality of complete Kraus maps proves the channel statement on externally entangled inputs.

The logical XOR map is self-inverse on the clean-work subspace. Reinvoking the complete instrument with fresh outcomes implements that logical adjoint. The outer inverse tape emits a fresh lookup macro rather than reversing measurements. The unitary eight-T mode uses a compute and its literal inverse within each selector pair and also realises the same logical XOR map. These are sufficient ideal logical premises for replacing the lookup inside preparation and Grover amplification. Physical gate errors, measurement errors and feedback latency remain outside this proof.

## 4. Resource identities and limitations

Each visited selector node emits two exact seven-T Toffolis, a four-T compute and four-T inverse, or a four-T compute with corrected Clifford erasure. Hence κ=14, 8 or 4 is correct for the emitted plans. The primitive counters separately retain measurements, conditional CZ and conditional X. They are not obtained by recosting the old Boolean Toffolis.

At fixed endpoint width w, input precision and span, one-bit table coarsening increases remainder width f by one. The exact multiply-count increment is 8w+12f+6 Toffolis. Combining it with the lookup-node difference gives the stated forward block count and the complete-flag saving 4[κ(I_b−I_(b−1))−7(8w+12f+6)]. The factor four accounts for two scalar kernels and their logical inverses. The κ=4 strict threshold simplifies correctly to I_b−I_(b−1)>14w+21f+10.5. The matched interval/full gap also accounts correctly for two extra forward subtractions and their inverses, giving 168w T gates, together with twice the strike block and the lookup difference.

These equalities require the specified emitter, fixed prefix and precision, no relevant rounded-constant aliases, valid function domains and an accurate visited-node count. The implementation continues to compute a selector even where a constant root could admit a cheaper specialised construction. Thus κ describes this explicit constructor; it is not a T-optimality or strongest-possible-QROM claim. The retained audit includes complete complex matrices for small lookups and selected branches for the larger lookup, with its scope recorded. The general induction, rather than the selected numerical branches, supports the larger coherent-channel claim.

## 5. Financial prior and interpretation

I checked D_R=Σw_i||R_i−R̄||²≥0 and the matched-drift moments E[A]=K exp(r−d), E[G]=K exp(r−d−D_R/2). Capped-payoff monotonicity and its one-Lipschitz property then give the arithmetic–geometric mean bound. The conditional logarithmic gap is non-negative and integrable. Jensen's inequality and Gaussian exponential moments give the stated J₀; substituting the matched drifts produces its alternative expression. The global conditioned-kernel Lipschitz bound therefore gives the stronger cap-dependent mean certificate.

Adding the residual discretisation and arithmetic bias before dividing by H correctly transfers that continuous mean bound to the implemented finite probability. The separately added geometric-price correction does not affect this probability. Since the finite output lies in [0,1], its variance is at most p(1−p), and the piecewise prior-based variance bound is valid. The classical Bernstein and Hoeffding derivations are sound and correctly permit the same prior to benefit classical estimation.

The fresh lookup instrument implements exactly the same finite function under ideal operations, so it does not add a new deterministic arithmetic bias to that financial certificate. It does add measurements and feedback that must be multiplied by all preparation and adjoint invocations. This review supports the precise constructor-specific identities and channel assumptions. It establishes neither a new generic lookup method, estimator superiority, physical runtime advantage nor a positive literature contribution by itself.
