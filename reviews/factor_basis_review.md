# Independent review of exact Gaussian factor-basis selection

## 1. Mathematical scope

This is a read-only review of the retained financial mathematics and a staging implementation of an independent rational auditor. It does not certify resource improvements, numerical financial prices or theoretical novelty. No external source was accessed by this reviewer. The literature leads in Section 7 require primary-source verification by the coordinating researcher.

Let the retained financial specification be

\[
\frac{dS_i(t)}{S_i(t)}=(r-q_i)\,dt+\lambda\,dW_0(t)+R_i\,dW(t),
\]

with standard independent Brownian factors, deterministic loadings, a common scalar factor and unit maturity. The stipulated synthetic drift is
\(\mu_i=r-q_i-(\lambda^2+\|R_i\|_2^2)/2\). The current frozen models have common dividend yield, unit maturity and normalised log returns from spot 100. If \(O^\top O=I_p\), the replacement \(R\mapsto RO\) leaves
\(RO(RO)^\top=RR^\top\), every row variance and every risk-neutral drift exactly unchanged. It defines the same multivariate lognormal law at all times, because \(OW\) is again a standard Brownian motion. Equivalently, coupling the transformed model to \(W'=O^\top W\) gives \(RO W'=RW\) pathwise. Discounted dividend-reinvested asset values remain martingales under the same risk-neutral measure. There is no new financial asset model or change in derivative target.

With \(r_G=\sum_iw_iR_i\), centred rows \(s_i=R_i-r_G\), centred drifts \(d_i=\mu_i-\sum_kw_k\mu_k\) and geometric drift \(m_G=\sum_iw_i\mu_i\), the transformed quantities are exactly \(r_GO\), \(s_iO\), \(d_i\) and \(m_G\). The geometric variance \(\lambda^2+\|r_G\|^2\), \(s_i\cdot s_k\), \(r_G\cdot s_i\), weighted residual dispersion and all continuous Gaussian exponential-moment arguments remain identical. Thus the continuous coefficients \(c_0,c_1\) and continuous ratio moments can be reused only after these exact equalities are proved. The invariant continuous Taylor screen cannot be made successful solely by an orthogonal factor rotation.

The non-linear conditional-price functions satisfy \(F_{RO}(y)=F_R(Oy)\), which implies equality of full Gaussian expectations. This pointwise identity is not \(F_{RO}(y)=F_R(y)\).

## 2. What does not remain invariant

A product truncated Gaussian box is not rotationally invariant. Its mass \(a^p\) is unchanged by using the same truncation parameters in each basis, but the retained region in the original coordinates becomes a rotated cube. The associated unnormalised retained integrals differ in general. The product law on rounded equal-mass Gaussian quantile nodes is likewise not rotationally invariant. Orthogonality preserves covariance, not higher cumulants of a finite product law.

The auditor includes a simple exact counterexample. Let the two independent finite factors take values -1 and 1 with equal probabilities and let the original loading be \((1,0)\). Under a rational rotation with cosine \(3/5\) and sine \(4/5\), its finite fourth moment changes from 1 to \(1201/625\). The full Gaussian fourth moment remains 3 in both bases. This counterexample does not use pricing approximations and rules out any inference of a common finite target solely from continuous covariance equality.

Consequently every candidate requires a fresh finite grid moment record, correlated finite underlying-gap envelope, finite domain and normaliser, R2 hybrid geometry, R2 tail, additional first-residual hybrid geometry and its tail, exact arithmetic certificate, implemented flag prior and complete statistical resource charge. Identical numerical values may occur, but must be demonstrated rather than assumed. Loading-dependent cache keys must include the transformed geometry. The direct capped-R1 excluded-mass bound \(H(1-a^p)\) is basis-independent; this does not make the older tilted tails invariant.

The rounded-grid first moments and second transport-measure mass are unchanged when the same one-dimensional node law is used in each coordinate. Their geometric loading multipliers change. All tensor hybrid measures remain unnormalised, each other-coordinate factor having mass \(a\), so direct R1 transfer carries the essential multiplier \(a^{p-1}\).

## 3. Direct derivative certificate

For full transformed rows \(R'_i=R_iO\), let \(A_j=\max_i|R'_{ij}|\) and \(r'_j=(r_GO)_j\). If \(\ell(y)=\log\sum_iw_i\exp(\mu_i+R'_i\cdot y)\), its softmax weights \(\pi_i\) satisfy

\[
\ell_j=\sum_i\pi_iR'_{ij},\qquad
\ell_{jj}=\sum_i\pi_i(R'_{ij})^2-\ell_j^2.
\]

For \(F=P(\ell)-P(m_G+r'_G\cdot y)\), the first bound is
\(|F_j|\le M_1(A_j+|r'_j|)\). Retaining the cancellation gives

\[
(P\circ\ell)_{jj}=(P''-P')(\ell)\ell_j^2
 +P'(\ell)\sum_i\pi_i(R'_{ij})^2,
\]

and hence
\(|F_{jj}|\le(M_1+M_2)(A_j^2+(r'_j)^2)\). Here \(M_1\) bounds \(|P'|\) and \(M_2\) bounds \(|P''-P'|\). A centred-row maximum cannot replace \(A_j\) in this proof. The auditor independently checks the complete gradient and Hessian conjugation as well as this cancellation and the inequalities for exact rational softmax laws and derivative values. These finite checks accompany the universal proof; they are not substitutes for it.

Writing \(b_1=\sum_c|\mu_{1,c}|\) and \(b_2=\sum_c\mu_{2,c}/2\), the direct allowance is

\[
a^{p-1}\left\{M_1b_1\sum_j(A_j+|r'_j|)
 +(M_1+M_2)b_2\left[\sum_jA_j^2+\|r_G\|_2^2\right]\right\}
 +H(1-a^p).
\]

This exact formula identifies the only basis-dependent factors in the direct certificate. The geometric squared norm and tail are invariant. The sums of coordinate maxima and geometric absolute loadings need not be. Reducing a coordinate-envelope sum is neither reducing the intrinsic covariance risk nor proving a quantum cost reduction.

For every orthogonal basis,

\[
\max_i\|R_i\|_2^2\le\sum_j\max_i|(R_iO)_j|^2
 \le\sum_i\|R_i\|_2^2.
\]

The lower bound follows by fixing any asset row before summing coordinate maxima. The upper bound follows by bounding the maximum by the sum over asset squares in each coordinate. Both outer bounds are basis-independent. A finite search does not prove attainment of the lower bound or global optimality. In one outer dimension the only orthogonal transformations are sign changes. Symmetric finite laws, box domains and the direct loading envelopes are unchanged, so the two-asset case is a useful negative control.

## 4. Financial and resource recomposition

The positive-R1 implemented probability is bounded by a basis-specific arithmetic certificate and financial transport. The invariant continuous second moment and \(c_1\) may tighten its prior, but \(\beta\), the normaliser \(B\), earlier finite mean bound and arithmetic enlargement \(\eta\) are new. Reusing the old implemented flag interval after rotating the finite law is unjustified. The existing prior formula properly divides by the external mass and expands by arithmetic error; its complete inputs must be new.

A reduced deterministic price error can increase the remaining statistical allowance and reduce the required phase width, but widths change discretely and full circuit charges may instead rise. Rotation changes rational coefficient numerators and denominators, rounded fixed-point constants, addition chains, nonzero loading counts, scalar table spans, headroom, actual word widths and register layout. Exact orthogonality of the financial input does not mean rounded arithmetic preserves covariance; arithmetic error is handled separately by the certificate. Resource conclusions require emitted complete financial preparations, their actual phase adjoints, controlled iterates, zero and good reflections, global-minus phases, Fourier tapes, repetitions and any failures. Replacing them with loading-sum or derivative-bound proxies is insufficient.

The direct certificate is conservative. The minimum of independently valid direct and residual-decomposition certificates remains valid and incurs no statistical selection penalty if selection is deterministic from shared mathematical information. However, complete older certificates for every new basis must actually be recomputed if this minimum is claimed. A cheap phase restricted to the direct certificate is defensible only if that scope is explicit and the identity comparison uses the same certificate family.

## 5. Classical comparisons and selection

For a fixed orthogonal basis and continuous independent standard Gaussian factors, the conditional payoff has exactly the same distribution. Therefore ordinary iid Monte Carlo using this same conditional payoff has identical ideal sampling variance before implementation effects. This statement does not extend to QMC, finite product quadrature, coordinate-conditioned controls or finite-grid score laws, all of which can exploit basis selection. QMC linear-transform and PCA methods are especially relevant classical opponents.

Every classical opponent must receive each permitted basis, exact moment/coefficient information and derivative certificates. Charging quantum basis search whilst denying it to a classical comparator, or showing a quantum change relative to a deliberately weak unrotated classical integration rule, would confound the result. If a classical interval already prices the invariant continuous target, improving quantum resource caps does not create a practical quantum advantage.

The finite rotation pool, objective, tie-breaking, fallback, candidate failure retention and basis-specific recomputation need to be registered before fresh outcomes. Source hashes and candidate exact matrices must be frozen. A selection criterion based solely on mathematical certificates and deterministic complete charges can be evaluated without using financial price answers and creates no stochastic multiplicity penalty. Selection using sampled quantum histories or classical price estimates requires separate simultaneous confidence accounting or fresh validation. Computational cost and timing of basis search and its table/certificate construction must be retained separately, with amortisation assumptions stated before treating repeated-use T-gate reductions as an operational improvement. Classical preprocessing is not automatically represented by logical T counts.

## 6. Minimum audit before an affirmative result

The independent auditor is staged as `review/audit_factor_basis.py`, with callable `audit_pair(original_model, transformed_model, O)`. It deliberately does not import the producer's transformation routine. Its retained self-test covers 17 pairs across all four originals, both signs of rational plane rotations in selected early and late coordinate planes, the identity and one-dimensional sign reversal. It checks each continuous multinomial exponent through degree six. It also checks stale derived covariance/variance/dispersion fields where present. It is not yet an audit of the final registered producer schema or complete quantum resources.

For final evidence, I require every declared matrix and transformed model to pass `audit_pair`, all original/rotated invariant fields and every selected global direct allowance to be replayed independently, a full candidate and failure count, a pointwise quantum-cost comparison to the identity under the same candidate methods/architectures/precision/confidence rule, complete arithmetic-profile and probability-prior reconstruction for costed candidates, and honest accounting of classical preprocessing. Any shortlist based on a cheaper proxy requires an explicit shortlist limitation and retained proxy outcomes, rather than a global optimum claim.

## 7. Literature leads and novelty limits

The strongest obvious antecedent is the linear-transform approach to derivative-pricing dimension reduction associated with Imai and Tan, including the title *A general dimension reduction technique for derivative pricing* and the earlier *Minimizing effective dimension using linear transformation*. The first is believed to have a published Journal of Computational Finance version; its exact volume, issue, pages and source status need verification. These methods optimise orthogonal Gaussian representations for QMC. PCA factor generation and Brownian bridge constructions are also established Gaussian reparameterisations. The project is a terminal multiasset factor problem, so Brownian path-dimension methods are related antecedents rather than automatically identical problems.

The exact rotation, Gaussian invariance, geometric coefficient invariance, derivative chain rule, finite-pool selection and coordinate-envelope bounds do not establish a new quantum transformation. Potential value lies in a checked project-specific optimisation of complete certified quantum pricing plans under unchanged classical risk-neutral targets, with rounding and financial-error composition accounted for and realistic classical comparisons. Whether a published paper has already combined orthogonal basis optimisation with complete fault-tolerant arithmetic resource/error budgets is an open literature question. Quantum finance state preparation, numerical integration, importance sampling, control variates and QMC comparisons should also be searched before claiming even a specialised methodological priority. This reviewer explicitly does not certify originality or positive contribution before those comparisons and fresh resource outcomes.
