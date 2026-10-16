# 1. What is established

I derived the signed Gaussian-tilt preparation for the finite Taylor remainder, its exact continuous and finite component normalisations, the postselection probability, the signed financial correction identity, a stable integral representation of the observable, and the resulting finite-precision interface. The accompanying proof retains all assumptions. Established linear-combination and amplification ingredients require attribution and the complete method still requires a literature comparison before any novelty claim.

I emitted a genuine quantum signal circuit for a finite classical terminal model with two independent Gaussian-density coordinate restrictions. It uses eight qubits and 60 exact Toffolis. Every arbitrary rotation is synthesised and retained: the signal preparation costs 4,416 T gates and 11,379 primitive gates. Eleven distinct rotation tapes pass independent Acb384 certificates against angles rebuilt from the model, with all 42 uses charged. The preparation norm-error upper bound is approximately 2.215×10⁻⁸. Conditional component phases are retained, and the total physical phase cancels through adjoint pairs.

At dispersion δ=1/4, the signal succeeds with probability 0.00153885. Nineteen ordinary amplification iterations cost 173,288 T gates and ten qubits, including the full zero reflection. Their certified aggregate norm-error upper bound is approximately 8.637×10⁻⁷. Literal primitive-gate replay returns success probability 0.99836037 and clean work. The finite signal identities agree at about floating roundoff; that replay is diagnostic, while the independent matrix certificates supply rotation rigour.

# 2. The principal competitive problem

The displayed signed preparation suffers cancellation. For a nondegenerate weak-dispersion family, its success probability is Θ(δ^(4s)). Amplification reduces the preparation-query penalty to Θ(δ^(−2s)), but the resulting ideal estimation order is δ^(2s)/ε. An equally informed direct remainder range on the certified finite support has order δ^(4s)/ε and uses the base preparation directly. The importance preparation is therefore asymptotically inferior in that regime when observable costs are comparable. This is a conclusion about this representation and regime, not a lower bound for every quantum importance method.

The toy confirms why the direct opponent matters. Its range weight exceeds the importance moment weight by a factor about 12.769, while amplified preparation requires 39 signal calls. Even before evaluating the financial observable, this is about 3.054 times as many preparation calls after allowing for the narrower importance weight. The same directly prepared base state costs 736 T gates; the corresponding signal-only T ratio is about 18.44. These are matched signal diagnostics and must not be presented as complete pricing-cost ratios.

# 3. Observable and amplification requirements

The quotient R/U^(2s) can be evaluated stably through the integral Taylor remainder. It becomes a Beta(1,2s) expectation of the derivative divided by its norm certificate. This bounded construction removes a singular numerical division, but requires an emitted derivative oracle, Beta preparation or certified quadrature, and the associated representation error.

Ordinary integer amplification leaves a failure component. A coherent pricing calculation can define the observable to vanish outside the good selector and divide its expectation by the certified ideal good probability κ. Its bias must include preparation error and uncertainty in κ. Treating the amplified state as exactly the importance distribution is incorrect. Rejection sampling requires measurement and reset, while exact final amplification phases require extra synthesised rotations. The finite-moment normalisation must be certified tightly enough that small cancellation probabilities are not rounded to an unverified value.

# 4. Classical model and no-arbitrage

The finite model is defined using a positive classical measure and its asset means are normalised to the stipulated forward. The exact common discounted-expectation functional is linear and positive, prices cash correctly, and assigns zero price to the forward cash flows. The quantum weighting supplies a computational identity for that functional; it is not a substitute risk-neutral measure. The toy verifies the same financial identity and a vanilla call curve under that measure.

The finite Gaussian restriction is a stated toy target, not the paper's exact continuous GBM model and not empirical evidence. Its forward normalisation avoids a drift mismatch within the toy. A production application must retain the continuous-model target or explicitly charge the finite approximation. Independently rounded estimated derivative prices can violate joint no-arbitrage inequalities despite valid individual tolerances; joint feasibility or the exact underlying expectation functional is needed for an across-contract claim.

# 5. Work still required

This package implements and certifies the signal preparation and amplification. The financial remainder oracle, signed payoff flag, full deterministic price-error allocation, estimation controller, observations and reset remain unimplemented. There is no complete quantum pricing algorithm or advantage result yet. The present candidate is worth retaining as a rigorously checked construction and a rejection test for an attractive but costly importance proposal. A positive complete-cost result must identify an input regime or emitted observable saving which defeats the direct classical-control and direct signed-range opponents.
