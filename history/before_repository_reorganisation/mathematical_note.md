# Monetary certificates for a conditional digital evaluator

## 1. Target and scope

I investigate a conditional payoff evaluator whose accuracy is certified in currency units. I hold the financial model and contract fixed before comparing arithmetic costs. The candidate replaces two conditional call calculations with a random strike, one exercise threshold and one normal tail probability. I also consider its nonnegative geometric-control residual. This note establishes the identities and error bounds needed for that feasibility test. It does not establish an original theorem, a compiled quantum circuit or a quantum advantage.

Let the positive monitored prices be

\[
S_j(y,z)=\exp(\ell_j(y)+b_jz),\qquad
A_y(z)=\sum_{j=1}^n w_jS_j(y,z),
\]

where \(w_j>0\), \(\sum_jw_j=1\), \(Z\sim N(0,1)\), and \(Y\) is an independent standard Gaussian vector. These are rotated coordinates of the continuous Gaussian log-price model. Initially I require \(b_j>0\), so \(A_y\) is strictly increasing from zero to infinity. The strike \(K\) and cap \(H\) are positive and have the same currency units as the average. The undiscounted payoff and discounted price are

\[
F(a)=\min\{(a-K)^+,H\},\qquad
V=d_T\,\mathbb E[F(A_Y(Z))],\qquad d_T=e^{-rT}>0.
\]

Every accuracy allowance below is undiscounted unless stated otherwise. Its contribution to price error is multiplied by \(d_T\).

## 2. One threshold and one normal probability

The layer-cake identity gives

\[
F(a)=\int_K^{K+H}\mathbf1_{\{a>k\}}\,dk
=H\,\mathbb E_U[\mathbf1_{\{a>K+HU\}}],
\qquad U\sim\mathrm{Unif}[0,1].
\]

For each \(y,k\), let \(z_A(y,k)\) be the unique root of \(A_y(z)=k\), and write \(\Phi\) for the standard normal CDF. Integrating the Gaussian coordinate gives

\[
h_A(y,u)=H\Phi[-z_A(y,K+Hu)],\qquad
\mathbb E_{Y,U}[h_A(Y,U)]=\mathbb E[F(A_Y(Z))].
\]

This output lies in \([0,H]\). Its formula uses one root and one normal probability, rather than two roots and a tilted-CDF sum over every monitored term. The additional uniform coordinate is part of the estimator input and its cost must be charged.

If \(g_A(y)=\mathbb E_Z[F(A_y(Z))]\), then \(\mathbb E_U[h_A(y,U)]=g_A(y)\). Consequently

\[
\mathrm{Var}(h_A(Y,U))
=\mathrm{Var}(g_A(Y))
+\mathbb E_Y[\mathrm{Var}_U(h_A(Y,U))].
\]

Randomising the strike can increase variance relative to full conditioning. There is no general variance ordering between \(h_A\) and the original payoff. Ordinary amplitude estimation with a fixed normalisation has the same mean amplitude for both. A resource claim based on standard deviation requires an implemented variance-sensitive mean estimator; it cannot be inferred from a conventional amplitude-estimation demonstration.

## 3. A nonnegative geometric-control interval

I use the weighted geometric average

\[
G_y(z)=\exp\left(\sum_jw_j\ell_j(y)+b_Gz\right),
\qquad b_G=\sum_jw_jb_j>0.
\]

The weighted arithmetic–geometric mean inequality gives \(A_y(z)\geq G_y(z)\) for every input. Define

\[
R(y,z)=F(A_y(z))-F(G_y(z)).
\]

Since \(F\) is increasing and one-Lipschitz,

\[
0\leq R\leq\min\{H,A-G\}.
\]

The geometric threshold is explicit:

\[
z_G(y,k)=\frac{\log k-\sum_jw_j\ell_j(y)}{b_G}.
\]

Because \(z_A(y,k)\leq z_G(y,k)\), the controlled digital output is the positive interval probability

\[
h_R(y,u)
=H\{\Phi[-z_A(y,K+Hu)]-\Phi[-z_G(y,K+Hu)]\}
=H\Pr\{z_A(y,K+Hu)<Z\leq z_G(y,K+Hu)\}.
\]

Thus \(0\leq h_R\leq H\), \(\mathbb E[h_R]=\mathbb E[R]\), and

\[
V=d_T\{\mathbb E[F(G)]+\mathbb E[h_R]\}.
\]

The geometric correction is lognormal and can be calculated separately. Its numerical error belongs in the complete price allowance. The interval probability should be evaluated directly or through a stable normal-tail difference, rather than by subtracting two rounded probabilities close to one. This changes the arithmetic design, not the exact identity.

For common loading \(b_j=\lambda>0\), write

\[
A_y(z)=B_A(y)e^{\lambda z},\qquad
G_y(z)=B_G(y)e^{\lambda z},
\]

where \(B_A(y)=\sum_jw_je^{\ell_j(y)}\) and \(B_G(y)=\exp(\sum_jw_j\ell_j(y))\). Both roots are logarithms:

\[
z_A(y,k)=\frac{\log(k/B_A(y))}{\lambda},\qquad
z_G(y,k)=\frac{\log(k/B_G(y))}{\lambda}.
\]

No iterative root solver is required in this case. I retain this established analytical construction as a mandatory comparator. Its statistical effectiveness can weaken with dense Asian monitoring, so computational convenience does not establish an overall resource benefit.

## 4. Residual certificates give a monetary mean bound

Let \(q_y(k)=\Pr\{A_y(Z)>k\}\), extended by \(q_y(k)=1\) for \(k\leq0\). Suppose an approximate threshold \(t_y(k)\) satisfies the true residual certificate

\[
|A_y(t_y(k))-k|\leq\rho
\]

for every retained \(y\) and every \(k\in[K,K+H]\). Monotonicity implies

\[
q_y(k+\rho)\leq\Phi[-t_y(k)]\leq q_y(k-\rho).
\]

Integrating over the strike interval bounds the approximate mean between the exact capped prices with strikes \(K+\rho\) and \(K-\rho\). A capped payoff is one-Lipschitz in its strike. Therefore

\[
\left|\mathbb E_U[H\Phi(-t_y(K+HU))]-g_A(y)\right|\leq\rho.
\]

This is a mean certificate after averaging over the uniform strike. It is not a pointwise bound \(|h_A-\widetilde h_A|\leq\rho\). Near a flat exercise root, the digital probability can be sensitive even when the integrated capped price is stable. The approximate threshold map need not itself be monotone; the residual sandwich is sufficient.

For midpoint strikes \(k_j=K+H(j+1/2)/N\), \(j=0,\ldots,N-1\), monotonicity of \(q_y\) bounds the midpoint quadrature error by \(H/(2N)\). Applying the residual sandwich to each midpoint yields

\[
\left|\frac{H}{N}\sum_{j=0}^{N-1}\Phi[-t_y(k_j)]-g_A(y)\right|
\leq\rho+\frac{H}{2N}.
\]

This bound does not require a total variation comparison between a continuous uniform law and its discrete grid.

For the controlled output, suppose the arithmetic and geometric thresholds have true residual bounds \(\rho_A,\rho_G\). Form the exact-CDF difference and project it to \([0,H]\) if the approximate thresholds reverse order. The upper sandwich

\[
q_A(k-\rho_A)-q_G(k+\rho_G)
\]

is nonnegative because \(A\geq G\). The lower sandwich is \(q_A(k+\rho_A)-q_G(k-\rho_G)\); replacing a negative lower value by zero only strengthens the lower bound. Integrating, or summing midpoint strikes, gives respectively

\[
|\mathbb E_{Y,U}[\widetilde h_R]-\mathbb E[R]|
\leq\rho_A+\rho_G,
\]

\[
|\mathbb E_Y[N^{-1}\sum_j\widetilde h_R(Y,(j+1/2)/N)]-\mathbb E[R]|
\leq\rho_A+\rho_G+\frac{H}{N}.
\]

The two grid allowances are charged separately. If each normal probability has uniform absolute error at most \(\eta_\Phi\), projection to the payoff range adds at most \(H\eta_\Phi\) for the single digital, or \(2H\eta_\Phi\) for the controlled difference. Coefficient, exponential, logarithm and comparison errors must be included when certifying the true residual. A small rounded residual alone does not establish either certificate.

## 5. Boundary handling and reversible implementation

I restrict the first prototype to positive stochastic loadings. Zero-loading terms contribute an offset. When the threshold is at or below that offset, the digital probability is one; an entirely deterministic average is evaluated directly. Mixed-sign loadings require a different crossing analysis and are outside this initial construction.

A finite inner search domain \([-R_Z,R_Z]\) needs certified endpoint comparisons. If a root lies above the domain, replacing its normal probability by zero incurs at most \(\Phi(-R_Z)\); if it lies below, replacing it by one incurs the same bound. Thus the integrated payoff allowance is at most \(H\Phi(-R_Z)\) for one digital, with separate allowances for two controlled probabilities. Ambiguous endpoint comparisons must either supply a true residual certificate or retain an uncertainty allowance.

The initial constructive target is the root-free common-loading interval evaluator, including its logarithms, stable normal-probability difference and uncomputation. General-direction root solving remains an extension. If I introduce it, a fixed reversible schedule needs guarded exponential-sum evaluation, certified comparisons and a bounded root bracket. Early-convergence flags do not remove gates from an ordinary fixed-depth circuit. Classical average iteration counts therefore cannot certify its logical resource cost.

## 6. Continuous and finite distributions

A continuous Gaussian law and a finitely supported register law are mutually singular, with total variation distance one. I therefore use coupling or quadrature to transfer a continuous price to a finite numerical target. Each rotation and discretisation has its own target-bias bound; finite product grids are not exactly rotationally invariant.

For the fully averaged conditional payoff \(g_A(y)\), a useful coupling bound follows from

\[
\|\nabla_yF(A_y(z))\|
\leq(K+H)\max_j\|r_j\|
\]

almost everywhere, where \(\ell_j(y)=\ell_{j,0}+r_j^\mathsf Ty\). Hence \(g_A\) has the same Lipschitz bound. Coupling \(Y\) to its discretised version \(Y_d\) gives a mean error at most

\[
(K+H)\max_j\|r_j\|\,\mathbb E\|Y-Y_d\|.
\]

For \(g_R=g_A-g_G\), the resulting global Lipschitz constant is

\[
L_R=(K+H)\left(\max_j\|r_j\|+
\left\|\sum_jw_jr_j\right\|\right).
\]

The residual and finite-strike certificates hold uniformly in retained \(y\), so this transfer can be applied to the averaged payoff without asserting uniform smoothness of an individual digital probability. In particular, an individual common-loading digital can involve a derivative proportional to \(1/\lambda\), while this fully strike-integrated residual bound does not.

I can specify a finite Gaussian target without a separate outer truncation. For one coordinate, let \(N_Y=2^{b_Y}\) and

\[
Y_{b_Y}=\Phi^{-1}\left(
\frac{\lfloor N_Y\Phi(Y)\rfloor+1/2}{N_Y}\right).
\]

This target has uniform mass at the midpoint quantiles
\(z_j=\Phi^{-1}((j+1/2)/N_Y)\). Put
\(a_j=\Phi^{-1}(j/N_Y)\), with \(a_0=-\infty\),
\(a_{N_Y}=+\infty\), and \(\phi(\pm\infty)=0\). Integrating the absolute coupling error within a probability cell gives

\[
\int_{a_j}^{a_{j+1}}|y-z_j|\phi(y)\,dy
=2\phi(z_j)-\phi(a_j)-\phi(a_{j+1}).
\]

The terms proportional to \(z_j\) cancel because the midpoint quantile splits the cell's probability equally. Therefore

\[
D(N_Y):=\mathbb E|Y-Y_{b_Y}|
=2\sum_{j=0}^{N_Y-1}\phi(z_j)
-2\sum_{j=1}^{N_Y-1}\phi(a_j).
\]

For \(p=\dim Y\), coordinatewise coupling and
\(\|x\|_2\leq\|x\|_1\) yield

\[
|\mathbb E[g_R(Y_{b_Y})]-\mathbb E[g_R(Y)]|
\leq L_RpD(N_Y).
\]

The extreme probability cells already include the Gaussian tails. An additional outer truncation allowance is required only if the implementation changes this stated discretisation. Uniform error at most \(\eta_Y\) in each implemented inverse-normal coordinate adds at most \(L_R\sqrt p\,\eta_Y\). The formula for \(D\) should be evaluated with a certified upper enclosure; subtracting large aggregate sums without error control is insufficient. The retained [certificate diagnostics](results/certificate_checks.json) compare it with direct absolute-error integration for \(N_Y=1,2,4,16,256\); the largest observed difference is about \(2.0\times10^{-15}\). This numerical agreement is an implementation diagnostic rather than a rounding certificate.

For a controlled residual with \(N_U\) midpoint strikes, the directly identified mean-bias components consequently include

\[
\rho_A+\rho_G+\frac{H}{N_U}
+2H\eta_\Phi+L_RpD(N_Y)+L_R\sqrt p\,\eta_Y.
\]

This expression still needs the declared inner-search boundary allowance, any remaining arithmetic allowance and the finite-preparation term below. I do not replace these terms with a continuous-to-discrete total variation argument.

Total variation is useful between an intended finite law \(P_d\) and a prepared finite law \(Q_d\). If their distance is at most \(\tau\), an output in \([0,H]\) has mean perturbation at most \(H\tau\) and variance perturbation at most \(H^2\tau\). A consistently implemented approximate sampler with its exact adjoint defines a changed finite target; I charge that target bias once. Inconsistent repeated operations, approximate reflections and physical errors require a separate implementation analysis.

## 7. Conservative variance bounds and a checked adverse case

The positive geometric residual supplies an analytical bound without a variance pilot:

\[
\mathbb E[R^2]\leq H\mathbb E[R]\leq H\mathbb E[A-G],
\qquad
\mathbb E[h_R^2]\leq H\mathbb E[A-G].
\]

If the log prices are \(\mu_j+a_j^\mathsf TX\), their moment difference is

\[
\mathbb E[A-G]
=\sum_jw_je^{\mu_j+\|a_j\|^2/2}
-\exp\left(\sum_jw_j\mu_j+
\frac12\left\|\sum_jw_ja_j\right\|^2\right).
\]

For a clipped implemented residual with complete mean-bias allowance \(\beta_R\), a safe standard-deviation bound is

\[
\overline\sigma_R
=\min\left\{\frac H2,
\sqrt{H[\mathbb E(A-G)+\beta_R]}\right\}.
\]

The moment calculation itself must have an upper numerical enclosure. This bound is conservative and applies to the direct geometric residual too. It cannot establish a variance advantage for the conditional digital. Its role is to permit an initial guaranteed resource comparison without treating an empirical sample variance as known.

I checked a scalar adverse case with \(A(Z)=105e^{bZ}\), \(K=100\), \(H=10\), and no outer random coordinate. One-dimensional positive-payoff quadrature and normal-tail calculations give the following undiscounted variances, in squared payoff currency units.

| Loading \(b\) | Direct capped-payoff variance | Random-strike digital variance |
| ---: | ---: | ---: |
| 0.0001 | 0.00011025 | 24.94076 |
| 0.1 | 18.80507 | 1.15150 |
| 0.8 | 24.20837 | 0.01880 |

Full conditioning has zero variance here because there is no outer coordinate. These are mathematical diagnostics, not market observations. They demonstrate both a favourable digital-conditioning regime and a severe adverse regime. A separate check of 45 basket, midpoint-grid and alternating residual-shift combinations satisfied the stated residual-plus-grid bound. These numerical checks support implementation review; the preceding arguments establish the bounds.

## 8. Literature position and decisive comparison

I treat geometric controls and analytical Gaussian smoothing as established methods. Related primary sources include Kemna and Vorst, *A pricing method for options based on average asset values*, Journal of Banking & Finance (1990), [DOI](https://doi.org/10.1016/0378-4266(90)90039-5), and Bayer, Siebenmorgen and Tempone, *Smoothing the payoff for efficient computation of basket option prices*, [author manuscript, Section 3](https://arxiv.org/abs/1607.05572). Montanaro's [quantum mean-estimation manuscript](https://arxiv.org/abs/1504.06987) provides relevant estimator context. In this step I consulted the supplied critiques and checked the displayed calculations directly; I have not represented these links as final published texts read during this subtask. A source-version and section record from primary-paper reading remains necessary before a manuscript contribution claim.

The [source record](SOURCES.md) retains the primary-source versions and sections consulted during the wider redesign. I have not established theoretical novelty for these calculations.

The proposed contribution is the prospective implementation and comparison of a certified conditional interval evaluator. The first decisive test is whether its root-free common-loading form reduces complete reversible resources at a fixed monetary accuracy after accounting for strike randomisation, finite distributions and estimator variance. The comparison includes stable fully conditioned common-loading pricing and the direct geometric residual. General-direction root methods remain a later comparison. Selecting an existing cheapest method, reproducing a known identity or reporting a favourable classical runtime would not establish that contribution.

## 9 Constructive representation update

I now use the checked bounded-variation certificate in the [finite-grid supplement](finite_grid_supplement.md) as an additional, often sharper transfer bound. For positive baskets, each outer coordinate contributes at most 3H/(2N), while the direct common Gaussian coordinate and the controlled strike coordinate contribute at most H/N. The supplement proves these constants, handles the exact comparison-register probability and gives the primitive generator enclosures. I retain the earlier coupling derivations as valid alternative bounds rather than interpret their loose constants as inherent discretisation costs.

The [reversible study](reversible_encoding_study.md) constructs all three evaluators for a stipulated two-asset instance and reports their matched bias budgets and emitted logical resources. These constructive results do not certify a statistical pricing interval or a literature novelty claim.

## 10 Specialised conditioning update

I now implement the scalar normalised conditional payoff \(P_h(m)=h^{-1}\int_1^{1+h}\Phi((m-\log k)/\lambda)\,dk\). Its log-level Lipschitz constant is at most \(1/(\lambda\sqrt{2\pi})\), and its second derivative is bounded by \(1/(\lambda^2\sqrt{2\pi e})\). The [certificate review](checks/kernel_certificate_review.md) checks differentiation, the interpolation remainder, financial input error, realised domains and amplified offline ordinate enclosures. The [stronger benchmark](conditioning_benchmark.md) applies analogous scalar specialisation to the interval method. The revised full comparator reverses the initial resource ranking in all six cases. These derivations use established mathematical ingredients and do not establish theoretical novelty or a complete statistical pricing guarantee.

## 11 Keywords

Basket call spreads; Gaussian conditioning; geometric controls; residual certificates; reversible arithmetic; quantum mean estimation.
