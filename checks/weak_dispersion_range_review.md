# Independent review of the bounded-support range argument

## 1. Reviewed version and finding

I reviewed proof_sources/weak_dispersion_range_proof.tex with SHA-256 48e84837e0656499e9afbd53d90310cb0b5690fd4c634cbf1fd232cb22cd7974. The final version incorporates explicit H>0, finite nonnegative δ, finite nonnegative Lipschitz constant, the same-input pointwise arithmetic error, the implemented range [0,1], the prior-midpoint equality case, the inverse-angle bootstrap argument, the retained rational upper enclosure of π and the scalar implementation's logarithm guard. I found no remaining mathematical defect in the stated finite-range and outer-query conclusion.

The derivation establishes a sufficient upper bound under its assumptions. It makes no priority claim, quantum lower bound, physical runtime claim or claim of global encoding optimality.

## 2. Uniform log-gap and numerical transfer

The exponential integral remainder is exact for every real argument. Bounding e^(t h) by e^δ, integrating 1−t on [0,1] and summing with the positive weights gives the stated second-order exponential remainder. Exact weighted centring cancels the linear term. Jensen's inequality supplies a lower bound of one on the weighted exponential sum. Monotonicity of the logarithm and log(1+z)≤z then give both claimed bounds on J. A global Lipschitz constant controls P(m+J)−P(m) uniformly even if the base m depends on the supported input or on ε.

For bounded v and fixed B,M,C,L,H, δ=BM|ε|+Cε² is O(|ε|), so the residual range is O(ε²). This conclusion requires those constants to remain uniform and H to remain bounded away from zero. The paper's stated fixed-constant example meets that condition. Weighted centring must hold exactly; a neglected first-order centring discrepancy could destroy the conclusion.

The same-input error inequality adds η to the ideal pointwise upper bound. The independent [0,1] constraint justifies taking the minimum with one. If η=O(ε²), the numerical range retains that order. Given an outward certified upper bound ρ and grid spacing 2^−q, flooring 2^qρ is legitimate because the implemented output numerator is an integer. Flooring an uncertified floating approximation would not be legitimate; the proof explicitly requires outward certification first.

## 3. Dyadic endpoints and scalar controller bounds

For a positive certified integer upper bound a, the exact integer choice c=(a−1).bit_length() is the least c satisfying a≤2^c. Consequently a≤2^c<2a, including powers of two. This gives the stated strict factor-two range bound. At a=1, c=0; at a=2^q, c=q. At a=0, every nonnegative grid output is zero and no logarithm or range division is needed. A certified a need not be attained, so the range is least among dyadic choices justified by that certificate, without claiming it is the smallest possible range for the actual financial function.

For an admissible range R, the range flag's Born probability is E[hat X]/R. Tolerance e_R=t/R therefore gives finite-mean error t. When e_R is at least half the prior upper bound, the prior midpoint has worst-case error at most e_R, including equality. Otherwise the prior exceeds 2e_R and hence exceeds a fixed positive constant under the stated scaling assumptions. The inverse normalised angle is then bounded away from zero, so its reciprocal bootstrap multiplier is bounded.

I checked the two retained planning implementations. Their derivative enclosure is at most one because the variance surrogate is at most 1/4 and its binary upper square-root cannot exceed 1/2. The three/five target is bounded by ceil(1/(5e_R)). The full-batch multiplier uses the actual rational upper enclosure Π_+ of π, so its cap is bounded by Π_+/[4(e_R−10^−60)]. On the supported domain e_R≥10^−6, the endpoint guard cannot approach the denominator's zero. Bounded initial and maximum multipliers give bounded rounds, finite shot tables, shot caps and Grover-call caps at fixed supported confidence. The α=0.01 setting lies well within the retained logarithm guard. These claims inherit the scalar controllers' ideal Born-probability and fresh independent-shot assumptions.

## 4. Retained illustrative checks

The accompanying check source verifies 6,113 exact integer dyadic endpoint cases for q from one to 64, including powers of two and their neighbours. It separately checks the zero-prior zero-shot branch for both controllers.

I also retain twenty bounded-support illustrations with ε=2^−n, n=1,…,20, uniform input addresses producing v=−1,0,0,1, weights 1/2, and centred perturbations h=(εv,−εv). The mathematical test function is P(m)=m at m=0 with H=L=1; this is not a financial pricing model or a benchmark. The log gap is log cosh(εv). For |ε|≤1/2, e^|ε|<2, so the proposition supplies J≤ε². With q=2n+3 and nearest-grid rounding, η=ε²/16. The conservative upper bound ρ=17ε²/16 gives certified integer maximum eight, c=3 and R=ε² exactly. The attained positive numerator is four in the high-precision illustrative calculation; the conservative maximum need not be attained.

Taking t=ε²/32 fixes the range probability tolerance at 1/32. The checked conservative controller caps with prior one and α=0.01 are 6,000 shots and 28,000 Grover calls for three/five, and 2,468 shots and 9,402 Grover calls for full-batch, independently of n. The finite output has mean ε²/4 and variance ε⁴/16, while its canonical flag variance is (ε²/4)(1−ε²/4). These exact rational identities illustrate the observation distinction and range normalisation. They do not establish a query separation, since this deliberately simple function can also be analysed directly. The transcendental rounding values are checked at 220-digit Decimal precision as diagnostics; the conservative range certificate follows from the analytic inequality and nearest-rounding bound.

## 5. Limits of application

The result concerns the outer estimator factor for a bounded numerical input family. Registers, financial evaluation, arithmetic certification and support construction may grow as ε decreases. The illustrative q grows with n; the twenty executed cases do not claim that the project's finite compiler and experiment domain implements an infinite asymptotic sequence.

An unbounded Gaussian law lacks the fixed support premise. A finite quantile grid restores boundedness for each instance, but increasing its address width changes its support. Continuous-law distribution transfer and the deterministic pricing budget must be charged separately. A finite-mean tolerance alone does not certify the corresponding continuous financial price. The proof preserves these distinctions and therefore supports the narrow conclusion that canonical flag variance does not establish a universal quantum lower bound or an unavoidable need for dyadic estimation.

## 6. Keywords

Weak dispersion; bounded support; range normalisation; quantum mean estimation; numerical error; comparator fairness.
