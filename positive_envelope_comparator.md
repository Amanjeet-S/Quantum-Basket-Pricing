# 1. Setting

This is an original project-specific derivation for comparison, based on established Gaussian tilting and classical/quantum rejection. I make no theoretical-priority claim. Let \(Y\sim N(0,I_d)\), \(p(y)=(2\pi)^{-d/2}\exp(-\|y\|^2/2)\), and
\[
T(y)=\sum_{i=1}^{n}w_i\exp(d_i+v_i\cdot y),\qquad T(y)\geq1.
\]
Fix \(s\geq1\) and write
\[
M_-=E[(T-1)^{2s}],\qquad M_+=E[(T+1)^{2s}].
\]
If \(M_-=0\), then \(T=1\) almost surely and the controlled remainder is zero almost surely. No importance state or quantum correction is required. I assume \(M_->0\) below. Both moments are finite because a finite sum of lognormal variables has every finite integer moment.

# 2. A positive Gaussian-mixture envelope

For a multi-index \(\alpha\in\mathbb N_0^n\), put \(j=|\alpha|\), \(v_\alpha=\sum_i\alpha_i v_i\), \(d_\alpha=\sum_i\alpha_i d_i\), and \(w^\alpha=\prod_iw_i^{\alpha_i}\). The binomial and multinomial identities give
\[
(T+1)^{2s}=\sum_{j=0}^{2s}\sum_{|\alpha|=j}
\binom{2s}{j}\frac{j!}{\alpha!}w^\alpha
\exp(d_\alpha+v_\alpha\cdot y).
\]
Define
\[
a_\alpha=\binom{2s}{j}\frac{j!}{\alpha!}w^\alpha
\exp(d_\alpha+\|v_\alpha\|^2/2).
\]
All coefficients are nonnegative. Completing the square yields
\[
e^{v_\alpha\cdot y}p(y)
=e^{\|v_\alpha\|^2/2}p(y-v_\alpha).
\]
Consequently \(M_+=\sum_\alpha a_\alpha\) and
\[
q_+(y)=\frac{(T(y)+1)^{2s}p(y)}{M_+}
=\sum_\alpha\frac{a_\alpha}{M_+}p(y-v_\alpha).
\]
This is a genuine positive mixture of shifted standard Gaussians. It can be sampled classically by choosing \(\alpha\) with probability \(a_\alpha/M_+\) and drawing \(N(v_\alpha,I_d)\). There are at most \(\binom{n+2s}{2s}\) multi-indices, before identical shifts are merged. Parameter-table construction and selection costs belong in the benchmark.

# 3. Classical acceptance and a quantum purification

Define the bounded acceptance probability
\[
h(y)=\left(\frac{T(y)-1}{T(y)+1}\right)^{2s}\in[0,1].
\]
It satisfies
\[
q_+(y)h(y)=\frac{(T(y)-1)^{2s}p(y)}{M_+}.
\]
The classical acceptance rate is therefore \(\eta=M_-/M_+\), and the conditional accepted density is
\[
q_-(y)=\frac{(T(y)-1)^{2s}p(y)}{M_-}.
\]
Ordinary repeated rejection has expected \(1/\eta=M_+/M_-\) trials per accepted draw.

The formal continuous quantum proposal is a purification
\[
|\Omega_+\rangle=\sum_\alpha\sqrt{a_\alpha/M_+}\,|\alpha\rangle
\int_{\mathbb R^d}\sqrt{p(y-v_\alpha)}\,|y\rangle\,dy.
\]
The index states are orthogonal, so measuring only \(y\) produces \(q_+\). A controlled flag rotation with success amplitude \(\sqrt{h(y)}=((T-1)/(T+1))^s\) has success probability \(\eta\). The normalised success state retains \(\alpha\), but its \(y\)-marginal is exactly \(q_-\). Any diagonal financial observable depending only on \(y\) has the desired expectation in that state. Erasing the mixture index is unnecessary for this calculation. A coherent preparation and its adjoint must preserve the full index register for the amplification and amplitude-estimation reflections.

Standard state-specific amplitude amplification changes the trial-call scale to \(O(\eta^{-1/2})\), subject to its usual known-success or robust approximate-success formulation. This is an established quantum rejection mechanism. The statement is about calls to a realised proposal/flag circuit; it does not price its gates or claim deterministic exact loading when only an interval for \(\eta\) is available.

# 4. Comparison with the signed Gaussian LCU

Put \(\psi_0(y)=\sqrt{p(y)}\). Completing the square at the amplitude level gives
\[
e^{v_\alpha\cdot y}\psi_0(y)
=e^{\|v_\alpha\|^2}\psi_{2v_\alpha}(y),\qquad
\psi_a(y)=\sqrt{p(y-a)}.
\]
Expanding \((T-1)^s\psi_0\) therefore yields coefficients
\[
\beta_\alpha=(-1)^{s-j}\binom{s}{j}\frac{j!}{\alpha!}w^\alpha
e^{d_\alpha+\|v_\alpha\|^2},\qquad j\leq s.
\]
Its squared norm is \(M_-\), and the coefficient one-norm is \(L_s=\sum_{|\alpha|\leq s}|\beta_\alpha|\). The usual signed LCU succeeds with probability \(M_-/L_s^2\).

Since replacing the signs gives exactly \((T+1)^s\psi_0\),
\[
M_+=\left\|\sum_\alpha|\beta_\alpha|\psi_{2v_\alpha}\right\|_2^2
\leq\left(\sum_\alpha|\beta_\alpha|\|\psi_{2v_\alpha}\|_2\right)^2=L_s^2.
\]
Thus \(M_-/M_+\geq M_-/L_s^2\), and the positive-mixture rejection scheme has no larger amplification overhead in the abstract trial-call comparison. This does not establish lower gate cost: its mixture contains degrees up to \(2s\), whereas the signed LCU contains degrees up to \(s\), and the acceptance arithmetic also has a cost. It is nevertheless a necessary competitor to a proposal whose claimed benefit depends on signed LCU.

# 5. Monetary estimator and the competing costs

Let the order \(2s-1\) Taylor residual satisfy
\[
|R(y)|\leq b_s(T(y)-1)^{2s},\qquad
b_s=\frac{H D_{2s}}{(2s)!}.
\]
Define \(F(y)=R(y)/(b_s(T(y)-1)^{2s})\) where \(T>1\), and set \(F=0\) where \(T=1\). The assumed bound forces \(R=0\) on the latter set, so the definition is consistent and \(|F|\leq1\). Then
\[
E_p[R]=b_s M_- E_{q_-}[F].
\]
Writing \(B_s=b_sM_-\), the correction equals \(B_s(2a-1)\), with \(a=E_{q_-}[(1+F)/2]\). A monetary error allowance \(\varepsilon\) therefore requires amplitude error at most \(\varepsilon/(2B_s)\). State loading is multiplied by this repeated-estimation effort, not paid only once.

At the level of abstract trial calls, signed LCU plus amplitude estimation has scale \(b_sL_s\sqrt{M_-}/\varepsilon\), while positive-mixture rejection plus amplitude estimation has scale \(b_s\sqrt{M_-M_+}/\varepsilon\), up to confidence and amplification precision factors. An equally informed classical rejection sampler followed by bounded Monte Carlo has trial-call scale \(b_s^2M_-M_+/\varepsilon^2\). These are not hardware or wall-clock comparisons. A proposal sampled without rejection, an optimised drift and existing deterministic price certificates may be stronger classical baselines.

If the existing deterministic interval already resolves the requested price tolerance, the optimal quantum correction effort is zero. Any superiority claim must therefore test tolerances or cases where that interval leaves a material unresolved correction.

# 6. Finite-register obligations

The continuous purification is not a finite-qubit implementation. I must choose a common finite register and retain certified errors for tail truncation, grid/cell approximation, component preparation, mixture weights, flag rotations and residual arithmetic. Truncated shifted Gaussians have component-dependent retained mass; using continuous mixture weights without charging that difference changes the density. A modular shift of a discrete Gaussian is not automatically the sampled translated continuous Gaussian. Both discrepancies need explicit bounds.

For the signed LCU, if every normalised component has vector error at most \(\delta\), then the unnormalised vector error is at most \(L_s\delta\). Provided the perturbed vector is nonzero, the normalised vector error is at most \(2L_s\delta/\sqrt{M_-}\). Severe cancellation thus amplifies component error as well as preparation cost. Declaring the components accurate at the baseline tolerance is insufficient.

For either method, evaluating \(F\) by subtracting nearly equal prices and dividing by a small \((T-1)^{2s}\) is numerically hazardous. An integral Taylor-remainder representation can provide a bounded observable without this subtraction/division, but it needs its own reversible implementation and error certificate. The current mathematical identity is a design specification until those circuits, synthesis errors and repeated costs are retained.

# 7. A bounded integral observable and a finite auxiliary register

Write \(r=2s\), \(U=T-1\), and keep \(m\) fixed when differentiating \(g_m\). Taylor's integral remainder gives
\[
R(y)=\frac{U(y)^r}{(r-1)!}\int_0^1(1-z)^{r-1}
g_m^{(r)}(1+zU(y))\,dz.
\]
If \(Z\) has density \(r(1-z)^{r-1}\) on \([0,1]\), it follows that
\[
F(y)=E_Z\left[\frac{g_m^{(r)}(1+ZU(y))}{H D_r}\right].
\]
The integrand has absolute value at most \((1+ZU)^{-r}\leq1\). It never divides by \(U\). At \(U=0\), its value may be defined by this same expression; the \(q_-\) density is zero there, so the pricing identity remains unaffected.

A finite implementation can obtain a discretised \(Z\) from the minimum of \(r\) independent uniform registers. Indeed, if \(V_1,\ldots,V_r\) are independent uniforms on \([0,1]\), then
\[
\Pr\{\min_i V_i>z\}=(1-z)^r.
\]
Use midpoint values \(V_i^{(h)}=(j_i+1/2)h\), with \(h=2^{-m}\) and independent uniform integer registers \(j_i\in\{0,\ldots,2^m-1\}\). Coupling each continuous uniform with its containing grid cell gives \(|V_i-V_i^{(h)}|\leq h/2\), and hence \(|Z-Z_h|\leq h/2\). A reversible comparator network computes their minimum; the uniform registers and work bits must be retained or coherently uncomputed inside the full preparation unitary.

The derivative of the bounded observable in \(z\) has absolute value at most \(U D_{r+1}/D_r\). Therefore the expected correction bias from this auxiliary-grid substitution, under the ideal \(q_-\), is bounded by
\[
B_s E_{q_-}[U]\frac{D_{r+1}}{D_r}\frac h2
=\frac{H D_{r+1}M_{r+1}}{r!}\frac h2,
\qquad M_{r+1}=E[U^{r+1}].
\]
All quantities are finite and can be bounded by the existing Gaussian moment engine. A finite truncated financial register may instead use its certified maximum \(U\). This proof charges only auxiliary-grid bias. Gaussian discretisation, derivative evaluation, synthesized flag rotation, amplification error and coefficient rounding remain separate charges.
