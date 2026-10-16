# Independent shared-threshold derivation and exclusions

## 1. Exact bounded positive identity

Let K,H,lambda>0 and let N be a standard normal variable. Define

\[
P(x)=\mathbb E\left[\min\{(K\exp(x+\lambda N)-K)^+,H\}\right].
\]

For every realisation of N, the layer-cake identity for the bounded nonnegative payoff gives

\[
\min\{(K e^{x+\lambda N}-K)^+,H\}
 =\int_0^H 1\{K e^{x+\lambda N}>K+v\}\,dv.
\]

Tonelli's theorem applies to the nonnegative integrand. Substituting v=Hu and taking an independent uniform U on (0,1) gives

\[
P(x)=H\Pr\{k(U)-\lambda N<x\},\qquad
k(u)=\log(1+(H/K)u).
\]

The law of k(U)-lambda N is continuous because lambda>0 and N is independent and has a continuous density. Replacing the strict inequality by a non-strict inequality leaves this continuous price unchanged. For any mA>=mG,

\[
P(m_A)-P(m_G)
 =H\Pr\{m_G<k(U)-\lambda N\le m_A\}.
\]

The same U,N must be used in both comparisons. Independent flags combined by multiplication would compute the product of marginal probabilities, which is not this price difference. The event representation has a genuine pointwise normaliser H: its integrand is an indicator on every label. It is not obtained by clipping a derivative integrand of range exceeding H.

For each financial outer coordinate y, weighted arithmetic-geometric mean yields mA(y)>=mG(y). The event's N-coordinate form is

\[
\frac{k(U)-m_A}{\lambda}\le N<\frac{k(U)-m_G}{\lambda}.
\]

A circuit can equivalently compare mG+lambda N<k(U)<=mA+lambda N. Comparator equality conventions must be retained in the execution record. For the finite N law equality atoms need not be absent; their effect is covered by the CDF approximation below.

## 2. Uniform strike quadrature

Let c=H/K and let

\[
\psi(u)=\Phi((k(u)-m_G)/\lambda)
       -\Phi((k(u)-m_A)/\lambda).
\]

Then F/H=integral_0^1 psi(u) du. Set xG=(k-mG)/lambda and xA=(k-mA)/lambda. The chain rule gives

\[
\psi''(u)=\frac{k'(u)^2}{\lambda^2}
  [\phi'(x_G)-\phi'(x_A)]
 +\frac{k''(u)}{\lambda}[\phi(x_G)-\phi(x_A)].
\]

Since k'=c/(1+cu), k''=-c^2/(1+cu)^2, u in [0,1], and

\[
\sup_z\phi(z)=1/\sqrt{2\pi},\qquad
\sup_z|\phi'(z)|=1/\sqrt{2\pi e},
\]

one obtains the uniform, outer-coordinate-independent bound

\[
\|\psi''\|_\infty
 \le 2c^2\left(\frac{1}{\lambda^2\sqrt{2\pi e}}
              +\frac{1}{\lambda\sqrt{2\pi}}\right)=M_U.
\]

The density maximum follows by differentiating exp(-z^2/2); the first density-derivative maximum follows by maximising |z| exp(-z^2/2), whose stationary positive point is z=1. The composite midpoint rule with nU equal-width cells therefore has probability bias at most M_U/(24 nU^2). The monetary bias is H times this quantity. This statement concerns exact midpoint strikes k((j+1/2)/nU); strike-table encoding error is charged separately.

## 3. A sharp finite common-Gaussian CDF bound

Let L>0, a=Phi(L)-Phi(-L), t=(1-a)/2, and let nN midpoint quantile nodes satisfy

\[
\Phi(z_j)=t+a(j+1/2)/n_N,\qquad j=0,\ldots,n_N-1.
\]

The finite law assigns probability 1/nN to each node. It is normalised even though its quantiles lie inside a truncated Gaussian interval. At z_j, the left and right finite-CDF values are j/nN and (j+1)/nN. Their discrepancies from Phi(z_j) are

\[
\frac{(1-a)j}{n_N}-t-\frac{a}{2n_N},\qquad
\frac{(1-a)j}{n_N}-t+\frac{1}{n_N}-\frac{a}{2n_N}.
\]

The first is never below -t-a/(2nN), and the second is never above t+a/(2nN). Their other extremes have smaller absolute values. Between nodes, the finite CDF is constant whilst Phi is increasing, so the same endpoint bound holds. The tails also obey it. Thus

\[
\sup_x|F_{n_N}(x)-\Phi(x)|\le t+\frac{a}{2n_N}.
\]

The same inequality holds for left CDFs by taking left limits and using continuity of Phi. If the actual rounded nodes obey |zHat_j-z_j|<=deltaN, coupling the corresponding labels gives F_node(x-deltaN)<=F_rounded(x)<=F_node(x+deltaN). The normal CDF is globally 1/sqrt(2pi)-Lipschitz. Hence

\[
D_N:=\sup_x|F_{\rm rounded}(x)-\Phi(x)|
 \le\frac{1-a}{2}+\frac{a}{2n_N}
      +\frac{\delta_N}{\sqrt{2\pi}}.
\]

For any interval with any combination of closed or strict endpoints, its finite probability differs from the Gaussian interval probability by at most 2D_N, since both endpoint CDF errors are bounded by D_N. Gaussian endpoint probabilities are zero. In particular, there is no claim that the finite atomic CDF is Lipschitz.

The common Gaussian finite label state must remain normalised. Its excluded-tail effect is already included in D_N. Multiplying the decoded signal by an additional common-factor mass a would target a different price. The outer retained mass remains a_y^p from the separate outer-grid construction.

## 4. Threshold arithmetic

Fix an ideal strike midpoint and a financial outer grid point. Let the circuit's two endpoint residual calculations approximate the ideal values mG+lambda n-k and mA+lambda n-k with uniform absolute errors eG,eA for every actual common-node label. These errors include the financial log values, exact strike-table approximation, encoded lambda times n, signed-affine shifts, and all fixed-point operations involved. Uniform nonnegative eG,eA may be used even when the actual error varies with n.

The computed interval event is contained in the expanded ideal interval

\[
(k-m_A-e_A)/\lambda\le n<(k-m_G+e_G)/\lambda,
\]

and contains the contracted ideal interval

\[
(k-m_A+e_A)/\lambda\le n<(k-m_G-e_G)/\lambda
\]

whenever that contracted interval is nonempty. If it is empty, its probability is defined as zero. This follows directly from the signed error inequalities and uses no smoothness of a binary flag or monotonicity of the computed arithmetic in n.

Each ideal endpoint displacement changes the corresponding continuous Gaussian CDF by at most its absolute displacement divided by sqrt(2pi). If the contracted interval is empty, the original interval length is at most (eG+eA)/lambda, so the same bound still follows. Applying the finite-law 2D_N bound directly to the expanded/contracted intervals yields

\[
|\Pr\{\hbox{computed event}\}-\psi(u)|
 \le 2D_N+\frac{e_G+e_A}{\lambda\sqrt{2\pi}}.
\]

A second 2D_N is unnecessary when this sandwich is composed directly. A looser bound can remain valid, but must be identified as such. The theorem requires valid endpoint error certificates over every finite input label and register headroom avoiding modular aliasing. Exact rational financial covariance preservation does not imply the rounded arithmetic meets these conditions.

## 5. Complete price composition

Let outer-grid transfer beta compare the exact conditional F against the full continuous outer expectation and the unnormalised retained node rule of mass a_y^p. Let etaU=M_U/(24nU^2), and let etaEvent=etaU+2D_N+(eG+eA)/(lambda sqrt(2pi)). The exact implemented binary event probability is q, its unnormalised residual is a_y^p Hq, and its discrepancy from the ideal outer-grid residual is at most a_y^p H etaEvent.

The complete undiscounted monetary error must include beta, coefficient rounding, external-mass interval uncertainty, a_y^p H etaEvent, and the mean-estimation error multiplied by a_y^p H. Distinct new approximations must not be charged twice or omitted. The new probability prior must enlarge the ideal finite residual interval by etaEvent, and must use B=H rather than an older B based on U_star unless the two are exactly identical and this is proved. The quantum failure probability applies to the implemented q after all deterministic allowances, with the monetary tolerance left after those allowances.

The generic layer-cake identity, threshold representation, midpoint theorem and CDF arguments are established probability and numerical-analysis mechanisms. The project-specific contribution, if a complete coherent resource comparison succeeds, concerns this finite financial-oracle configuration and its certified composition. It would not establish a new financial asset model, a new general quantum state-preparation theorem or empirical market evidence.
