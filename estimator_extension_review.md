# Extended canonical estimator review

## 1. Fourier extension

For a controlled phase of angle theta, replacing diag(1,1,1,exp(i theta)) by the identity has operator norm error |exp(i theta)-1| <= |theta|. At distance d, the inverse Fourier transform contains m-d such gates with |theta| = pi/2^d. The total omitted-gate allowance is therefore the sum of (m-d) pi_upper/2^d for d=32,...,m-1. Unitary telescoping combines this allowance with every retained scalar-rotation certificate. Omitted gates contribute no retained scalar correction. The emitted widths 1 through 32 remain unchanged. Wider registers reuse certified distances through 31 and explicitly omit later distances; the mathematical proof does not require new scalar rotation synthesis.

The retained aggregate scalar multiplies the final, uncontrolled Fourier transform. It therefore has no effect on the outcome probabilities. The independent audit nevertheless reconstructs that scalar from the raw rotation records and checks small direct Fourier matrices after applying it, including the declared bit order.

## 2. Probability and confidence transfer

If two unitary transforms differ in operator norm by at most eta, their resulting normalised pure states differ in Euclidean norm by at most eta. The trace distance between these pure states is at most that Euclidean distance. In detail, for overlap z, 1-|z|^2 <= 2(1-|z|) <= 2(1-Re z). Consequently, any measurement event has probability discrepancy at most eta. Subtracting an upward norm enclosure from the established 8/pi^2 ideal amplitude-estimation success lower bound is justified. The outward rounding to a 2^-80 grid can only decrease the success bound. Exact binomial majority amplification then controls the prescribed failure probability for independent fresh repetitions.

The amplitude error uses the Bernoulli variance envelope p(1-p), bounded by u(1-u) when u<=1/2 and by 1/4 otherwise. A moment bound on a numerical financial output cannot replace this flag variance without another encoding argument. The width is the least supported power of two meeting the conservative amplitude error criterion, including the retained postprocessing allowance. Required widths beyond 64 are explicitly rejected.

## 3. Complete controlled Grover accounting

The controlled Grover implementation applies a controlled good-state phase, an unconditional A adjoint, the controlled zero reflection, unconditional A, and the conditional global-minus Z. With control zero, the unconditional A and its adjoint cancel. With control one, the operation is exactly the required Grover iterate. The conditional minus cannot be dropped in phase estimation. The primitive count is the count of A plus its literal adjoint plus the emitted reflections. The adjoint swaps T with Tdg and S with Sdg. Twice the total T count is valid, but twice each individual gate-type count is generally incorrect.

The zero reflection on c core wires and one control uses c-1 clean conjunction wires. It contains 2(c-1) exact Toffolis, two controlled-Z gates after the separate good predicate is included, one conditional-minus Z, and 2c X gates. Its literal seven-T expansions therefore use 14(c-1) T or Tdg gates. The conservative serial depths charge four or five layers per exact Toffoli as explicitly declared.

Every allocated wire is terminally observed and conditionally reset for a fresh repetition. Phase observations are already included in the inherited measurement metric. The remaining observations add wires-phase_bits, while conditional resets add wires. The inherited expanded-operation metric includes the phase Hadamards but omits the phase measurement operations themselves. Its final complete-operation increment must therefore be 2*wires, not 2*wires-phase_bits. The independent audit checks the complete per-shot count and its product with the prescribed repetitions.

This is a resource and finite-confidence certificate. It does not assert that a physical quantum experiment was performed, that its noise obeys the logical-unitary model, or that logical gate counts determine practical speed.
