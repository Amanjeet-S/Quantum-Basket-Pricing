# Native pricing: mathematics and numerical evidence

## 1. Purpose and financial target

I implement a small C++17 capped-European-call engine to connect financial derivations, numerical analysis and executable programming. It develops the established [capped-GBM companion](../paper/capped_gbm_companion.pdf); it is neither a new financial model nor a quantum-theory claim. The [native interface and build instructions](../cpp/README.md), [comparison protocol](../results/capped_gbm/cpp/comparison/protocol.json) and [amendment history](../configs/cpp_pricing_amendment_history.json) identify the implementation and evidence. All contracts are synthetic. No new academic sources were consulted for this implementation; the companion retains the methodological attributions and source record.

For positive spot, strike, cap, maturity and volatility, risk-neutral GBM gives the same no-arbitrage value as the companion. Writing $D=\exp(-rT)$, $v=\sigma\sqrt T$ and $b=r-q-\sigma^2/2$, the pointwise payoff identity implies

$$
\min\{(X-K)^+,H\}=\int_K^{K+H}\mathbf1_{\{X>k\}}\,dk,
\qquad V=D\int_K^{K+H}\Phi\!\left(\frac{\log(S/k)+bT}{v}\right)dk. \tag{1}
$$

Tonelli applies because the integrand is non-negative and its strike interval has finite length. Thus $0\leq V\leq DH$. Under the stated positive parameters, the strike interval is compact and bounded away from zero. Each price derivative used here is continuous and bounded on that interval in a sufficiently small parameter neighbourhood, so differentiation under the integral is justified by dominated convergence. The call-spread formula and its Greeks provide another representation of exactly this value, rather than an independent market model.

Direct subtraction of two calls can lose a small capped value. I retain a normal-interval formula using `erfc`, `log1p`, `expm1` and fused multiplication, with separately evaluated endpoint coordinates for broad strike intervals. Severe price cancellation invokes a positive digital-strike quadrature. Its adaptive error indicator is diagnostic, not an enclosure. Extremely small transformed transition scales are explicitly rejected when the fallback cannot resolve them in the supported floating-point range. This rejection includes $S=K=10^{-100},H=10^{300},\sigma=10^{-8}$ and $S=K=1,H=10^{300},\sigma=10^{-100}$ with $T=1,r=q=0$, which previously exposed silent zero prices. Additional range checks reject underflowed normal terms that may recover representable scale after multiplication by large financial prefactors, and positive discount/cash intermediates lost to underflow. These checks exclude the reviewed $S=10^{300},K=10^{308},H=10^{307},\sigma=0.46,T=1,r=q=0$ counterexample; they do not supply a general arithmetic certificate. Ordinary tail underflow remains a numerical limitation. I make no all-finite-input accuracy claim.

## 2. Discrete stability and convergence

The companion transforms valuation to heat flow on a fixed domain containing both payoff kinks. With $M$ intervals, $N$ steps, $h=(U-L)/M$, $k=T/N$ and $\lambda=\sigma^2 k/(2h^2)$, the interior recurrence is

$$
A u^{n+1}=u^n+f,\qquad
A=\operatorname{tridiag}(-\lambda,1+2\lambda,-\lambda), \tag{2}
$$

where $f$ includes both fixed Dirichlet boundary contributions. The left boundary is zero and the right boundary is $H$; when there is one interior node, both boundaries belong to its equation. The native solver factors this matrix once and uses Thomas elimination at each step. The existing Python DST solver evaluates the same recurrence by a different algorithm.

For $\lambda\geq0$, $A^{-1}$ is non-negative. To prove this, suppose $Az\geq0$ but $z$ has a negative minimum at node $j$. Extend $z$ by zero at the two boundary nodes. Its neighbours are at least $z_j$, so $(Az)_j=z_j+\lambda(2z_j-z_{j-1}-z_{j+1})\leq z_j<0$, a contradiction. Applying this also to $Az=0$ and $A(-z)=0$ proves a trivial null space; a finite square matrix with trivial null space is invertible. Since $A\mathbf1\geq\mathbf1$, multiplication by its non-negative inverse gives $A^{-1}\mathbf1\leq\mathbf1$. Consequently

$$
\|A^{-1}\|_\infty\leq1. \tag{3}
$$

Differences with identical boundaries therefore contract at each step. Comparison with the constant full-grid subsolution zero and supersolution $H$ proves preservation of $[0,H]$ for initial values in that interval. No CFL restriction is required in exact arithmetic. The elimination pivots obey $p_1=1+2\lambda$ and $p_j=1+2\lambda-\lambda^2/p_{j-1}\geq1+\lambda>0$, by induction, which excludes exact zero pivots. These results alone do not prove floating-point accuracy.

The companion's theorem ‘Unconditional stability and a nonsmooth-payoff error bound’ supplies convergence for the actual kinked payoff by mollification, rather than assuming four derivatives of that payoff. In its notation, the nodal bound is $2L_fm_1\eta+T(ah^2/12+a^2k/2)L_fc_3\eta^{-3}$, where $a=\sigma^2/2$ and $\eta$ stays inside the two boundary margins. Choosing $\eta$ proportional to $(h^2+k)^{1/4}$ gives the stated conservative convergence bound. Finite-domain, interpolation, parameter-formation and arithmetic errors require separate treatment. Numerical Greek refinement can be non-monotone: a price error bound $E$ is amplified by at most $E/d$ for a central first difference and $4E/d^2$ for a central second difference, before Taylor bias. The study therefore separates PDE error at identical bumps from analytic bump bias and does not fit a convergence order.

## 3. Exact residual certificate for exported traces

I additionally check the native solver's exported trace with [exact rational arithmetic](../src/cpp_residual_certificate.py). Every decimal JSON value is first parsed as binary64 and then converted exactly to a rational number. The target is equation (2) defined by that published $\lambda$, cap and initial vector. It is not the exact scheme formed directly from the original real-valued financial parameters.

Let $\widehat u^n$ be the published interior vectors and let $u^n$ solve this exact discrete target, with $u^0=\widehat u^0$. Define

$$
\rho^{n+1}=A\widehat u^{n+1}-\widehat u^n-f,
\qquad R_{n+1}=\|\rho^{n+1}\|_\infty. \tag{4}
$$

Subtracting (2) from (4), $e^{n+1}=A^{-1}(e^n+\rho^{n+1})$. Equation (3) and the triangle inequality give $\|e^{n+1}\|_\infty\leq\|e^n\|_\infty+R_{n+1}$. Induction from $e^0=0$ yields

$$
\|\widehat u^N-u^N\|_\infty\leq B_N:=\sum_{n=1}^N R_n. \tag{5}
$$

The implementation evaluates every residual, maximum and sum exactly as an integer fraction. It retains numerator and denominator, and rounds the displayed binary64 bound upwards by checking its exact fraction and using the next representable number if necessary. If no finite binary64 upper bound exists, it reports that fact and retains the exact fraction. It never clips an overshoot into the payoff range. If the published initial vector lies in $[0,H]$, the exact target does too. Applying (5) at every intermediate step gives a range excess bounded by $B_n\leq B_N$; the auditor checks this implication separately from the strict range property.

For separately supplied binary64 interpolation weight $0\leq\theta\leq1$, discount $D\geq0$ and reported scalar $\widehat p$, set $\widehat I=D[(1-\theta)\widehat u_j^N+\theta\widehat u_{j+1}^N]$. Constant boundary errors are zero. Convex interpolation and (5) imply

$$
\left|\widehat p-D[(1-\theta)u_j^N+\theta u_{j+1}^N]\right|
\leq|\widehat p-\widehat I|+DB_N. \tag{6}
$$

Both terms on the right are computed exactly. This covers scalar evaluation mismatch for the explicitly recorded interpolation data. It does not certify original coordinate formation, the initial payoff evaluation, coefficient formation, discount evaluation, continuum discretisation, domain truncation, model adequacy or any Greek. Those exclusions are essential before composing a complete financial price bound. The argument is a project-specific application of the established discrete maximum principle, without a theoretical-priority claim.

## 4. Retained comparisons and interpretation

The [comparison manifest](../results/capped_gbm/cpp/comparison/study_manifest.json) retains 120 price/Greek outputs over 20 contracts, 126 PDE outputs over three coupled refinements, 21 complete-grid comparisons with Python Thomas/DST and 100 financial identities. Native comparisons pass the declared mixed absolute/relative thresholds against 80-digit digital-strike integration cross-checked with 100-digit call spreads. Those references are high-precision numerical checks, not outward-rounded certificates. Seven ordinary contracts precede disclosed stress additions; the [history](../configs/cpp_pricing_amendment_history.json) records informed changes and earlier failures. This is not an untouched prospective protocol.

The stronger ordinary-precision Python digital-integral baseline is retained alongside the original subtractive baseline. It often matches or improves on the native price error. Four derivative diagnostics fail their fixed tolerance: small-volatility gamma, and delta/gamma/rho for the extremely broad strike interval. Signed-integrand cancellation and coordinate subtraction explain those limitations. The native broad-interval gamma has absolute error about $5.82\times10^{-11}$ against magnitude $3.99\times10^5$, which illustrates why absolute and relative errors must be read together. The study reports that not all compared methods pass; it does not identify a language-wide accuracy ranking. Non-monotone PDE gamma refinement is retained.

The [trace study](../results/capped_gbm/cpp/traces/study_manifest.json) checks 21 native traces, including one-interior-node meshes, with exact residual bounds and seven independent exact one-node target comparisons. Four traces fail the deliberately strict floating-point range condition; those failures remain recorded and unclipped. Their exact range excesses lie within the accumulated residual bounds. The largest final-grid bound is approximately $1.65\times10^{-13}$, and the largest scalar-target bound is $1.56\times10^{-13}$; neither is a continuous-price certificate. The [nine Python certificate tests](../tests/test_cpp_residual_certificate.py) additionally check perturbed multi-step targets by rational Gaussian elimination, upward rounding, malformed traces and retained overshoots. Native assertions cover analytic identities, singular limits, supported-range rejection and solver behaviour; assertion counts are not counts of independent financial experiments.

The tested Apple arm64 compiler gives `long double` only 53 mantissa bits. Direct compiler builds are verified; CMake/CTest were unavailable. The Thomas implementation requires $O(MN)$ arithmetic, whereas the existing DST implementation costs $O(M\log M)$ for this constant-coefficient problem. No timing comparison, general C++ speedup, empirical model validation or practical quantum advantage is claimed. The positive outcome is a reproducible native financial calculation with proof/code connections and explicit numerical boundaries. Broader analysis/probability, empirical diagnostics and personal C++ proficiency still require separate evidence.
