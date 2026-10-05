# Independent certificate review of the normalised capped-lognormal kernel

## 1. Keywords and scope

Keywords: capped lognormal payoff; conditional expectation; interpolation certificate; common Gaussian factor; fixed-point arithmetic; reversible evaluation.

I review a scalar analytical conditioning kernel as a stronger fully conditioned comparator for the stipulated two-asset instance. I treat the layer-cake identity, Gaussian conditioning and elementary interpolation inequalities as established mathematics. I make no theoretical novelty claim. This memorandum concerns the numerical certificate and the implemented resource boundary; it does not establish end-to-end quantum advantage.

## 2. Exact financial target and derivatives

I use the price scale of 100, strike 1, relative cap \(h=H/100>0\) and common loading \(\lambda=0.2\). For an unshifted conditional log level \(m\), I define

\[
P_h(m)=\frac{1}{h}\mathbb E_Z\min\{(e^{m+\lambda Z}-1)^+,h\}
=\frac1h\int_1^{1+h}\Phi\!\left(\frac{m-\log k}{\lambda}\right)\,dk.
\]

The second equality follows by integrating the exercise indicator over the strike interval. It proves \(0\leq P_h\leq1\) and monotonicity directly. Differentiation under this finite integral is justified by the bounded Gaussian derivatives. Therefore

\[
P_h'(m)=\frac1{h\lambda}\int_1^{1+h}\phi\!\left(\frac{m-\log k}{\lambda}\right)\,dk,
\qquad
P_h''(m)=\frac1{h\lambda^2}\int_1^{1+h}\phi'\!\left(\frac{m-\log k}{\lambda}\right)\,dk.
\]

I obtain the global bounds

\[
0\leq P_h'(m)\leq L:=\frac1{\lambda\sqrt{2\pi}},
\qquad |P_h''(m)|\leq C:=\frac1{\lambda^2\sqrt{2\pi e}}.
\tag{1}
\]

Here \(L\approx1.994711402\) and \(C\approx6.049268113\). These bounds are independent of the positive cap. The corresponding unnormalised monetary conditional value has log-level Lipschitz constant \(HL\), rather than a coefficient proportional to the strike. This distinction prevents artificial amplification by \(1/h\) when I perturb one coherent financial input.

The exact residual is

\[
g_R(y)=H\{P_h(m_A(y))-P_h(m_G(y))\},
\qquad
m_A=\log B_A,\quad m_G=\log B_G.
\]

If the implemented log inputs obey \(|\widehat m_A-m_A|\leq e_A\) and \(|\widehat m_G-m_G|\leq e_G\), I can charge their financial effect by \(HL(e_A+e_G)\). This bound treats the exact financial function as a whole. It does not permit me to ignore unrelated approximation errors in independently implemented CDF calls in the old eight-CDF formula.

## 3. Kernel interpolation and the complete output certificate

I tabulate the scalar function on \([-1,1]\), with \(N=2^b\) cells and spacing \(\Delta=2/N\). The circuit can address this table using \(m+1\). For the existing log output \(\log B+2\), subtracting the exactly representable constant 1 supplies this biased kernel input without an additional rounding error.

Suppose each high-precision kernel ordinate has certified absolute generation error at most \(\epsilon_K\), and is rounded to the nearest \(q\)-fractional-bit value, with \(\delta=2^{-q}\). Endpoint linear interpolation, followed by the existing integer interpolation product, obeys

\[
|\widehat P_h(x)-P_h(x)|
\leq \frac{C\Delta^2}{8}+\frac32\delta+\epsilon_K
=:E_K
\tag{2}
\]

at every represented input in the table domain. The first term is the linear interpolation remainder, the second contains endpoint rounding and the product's floor error, and the last is the ordinate enclosure. Argument error is charged separately through (1). Thus a sufficient monetary evaluator allowance for two kernel calls, exact equal-scale subtraction and final projection to \([0,1]\) is

\[
\beta_{\rm arithmetic}
\leq H\{L(e_A+e_G)+2E_K\}.
\tag{3}
\]

Projection is non-expansive because the exact normalised residual belongs to \([0,1]\). Its endpoints 0 and 1 are exactly representable, so that final projection has no cap-rounding allowance. The final comparison-register encoding introduces no new floor if it compares directly with the already represented integer threshold.

I can make a minor, optional refinement when both kernel calls use the same fractional step. The interpolation-product floor residual lies in \([-\delta,0]\); the difference of two such residuals has magnitude at most \(\delta\). Each rounded linear interpolant contributes ordinate error at most \(\delta/2+\epsilon_K\), so their difference is bounded by \(\delta+2\epsilon_K\). A valid direct bound on the pair's kernel approximation errors is therefore

\[
\frac{C\Delta^2}{4}+2\delta+2\epsilon_K,
\tag{4}
\]

instead of \(2E_K=C\Delta^2/4+3\delta+2\epsilon_K\). The simpler nodewise use of (2) is already valid. Equation (4) requires the stated common-step floor construction and should not be transferred to arbitrary rounded arithmetic.

## 4. Realised input range

For the stipulated loadings and drifts,

\[
m_A(y)=\log\!\left\{\frac{e^{-.015-.1y}+e^{-.03+.2y}}2\right\},
\qquad m_G(y)=-.0225+.05y.
\]

I restrict the represented Gaussian coordinate only after the retained inverse-quantile certificate has established \(|y|<4.75\). The arithmetic average is convex and has its unique minimum at

\[
y_*=(.015-\log2)/.3\approx-2.260490602.
\]

Direct evaluation gives \(m_A(y_*)\approx-0.07663301227\), \(m_A(-4.75)\approx-0.02051648927\) and \(m_A(4.75)\approx0.4452999865\). Hence its exact range on this interval is contained in \([-0.076634,0.445301]\). The geometric log lies in \([-0.26,0.215]\). In particular, the arithmetic lower limit is not \(-0.03\).

The chosen \([-1,1]\) table domain is amply wide. Nevertheless, I must check the implemented input ranges including quantile representation, coefficient rounding, exponential interpolation, weighted averaging and logarithm interpolation. A nominal mathematical range does not establish safe table addressing after perturbation. A narrower \([-0.5,0.5]\) table is also possible for this instance only after these perturbed ranges are certified; it would reduce the squared-spacing term by a factor of four at the same address width. I do not use that narrower domain in (2).

## 5. High-precision ordinate generation

For generation only, I may use the exact four-CDF closed form. With

\[
d_0=m/\lambda,\qquad
d_1=(m-\log(1+h))/\lambda,
\]

it is

\[
P_h(m)=\frac{e^{m+\lambda^2/2}
\{\Phi(d_0+\lambda)-\Phi(d_1+\lambda)\}
-\{\Phi(d_0)-\Phi(d_1)\}+h\Phi(d_1)}h.
\tag{5}
\]

For \(m\in[-1,1]\) and the retained caps \(h\in\{0.01,0.1,0.25\}\), all four CDF arguments lie in \([-6.115718,5.2]\), so the existing \([-8,8]\) decimal-CDF generator certificate applies. A \([-2,2]\) kernel table would not satisfy that generator-domain requirement.

I charge cancellation in this offline formula explicitly. If each CDF has absolute allowance \(\epsilon_\Phi\), the exponential coefficient has allowance \(\epsilon_{\exp}\), and the remaining generator arithmetic has an absolute numerator allowance \(\epsilon_{\rm arith}\), then

\[
\epsilon_K\leq
\frac{(2e^{1.02}+2+h)\epsilon_\Phi+
\epsilon_{\exp}+\epsilon_{\rm arith}}h.
\tag{6}
\]

Coefficient multiplication introduces a negligible cross term if the exponential and CDF enclosures are applied separately; I include it in \(\epsilon_{\rm arith}\). With the existing per-CDF allowance \(10^{-38}\), the CDF part of (6) is below \(7.8\,10^{-36}\) across the three caps. It is therefore incorrect to reuse \(10^{-38}\) unmodified as the kernel ordinate allowance. For a seventy-digit implementation with at most 100 elementary operations outside the separately enclosed CDF calls, correctly rounded logarithm and exponential centres, and magnitudes below 32, I can conservatively retain \(\epsilon_{\exp}=10^{-68}\) and \(\epsilon_{\rm arith}=10^{-64}\). The latter also covers the Gaussian derivative bound applied to argument-rounding errors. These allowances make \(10^{-34}\) a sufficient kernel ordinate enclosure. The actual generator must satisfy the operation-count, rounding and range assumptions; the allowance is negligible beside the searched fixed-point steps.

I generate the kernel at the exact declared relative cap, not at a rounded circuit copy of the cap. The high-precision closed form can be checked against the positive integral as a numerical diagnostic, but agreement alone does not replace (6). I retain every endpoint enclosure and check that rounded ordinates are in \([0,2^q]\) and nondecreasing. The represented table then supplies a monotone integer interpolator.

## 6. Common-coordinate sensitivity and the old closed form

For the new scalar kernel, common coordinate errors should not be needlessly charged through large independent intermediate coefficients. The exact derivative is

\[
\frac{d}{dy}\{P_h(m_A(y))-P_h(m_G(y))\}
=P_h'(m_A)m_A'-P_h'(m_G)m_G'.
\]

Here \(m_A'\in[-0.1,0.2]\), \(m_G'=0.05\) and both kernel derivatives lie in \([0,L]\). Consequently the residual probability derivative is bounded in magnitude by \(0.2L\), or approximately \(0.3989422804\). The monetary derivative is bounded by \(0.2HL\). This is a direct bound for a coherent perturbation of the same coordinate. It does not absorb separate errors in the exponential and logarithm implementations.

The previous four-CDF value formula has an additional stationarity property which explains why independent nodewise argument propagation can be unfair. Write its exact unnormalised numerator as

\[
J(b,t_0,t_1)=b e^{\lambda^2/2}
\{\Phi(t_0+\lambda)-\Phi(t_1+\lambda)\}
-\Phi(t_0)+(1+h)\Phi(t_1).
\]

At \(t_0=\log b/\lambda\) and \(t_1=\log\{b/(1+h)\}/\lambda\), both partial threshold derivatives vanish. If only these thresholds are perturbed by at most \(\eta_0,\eta_1\), while the coefficient, tilted offset and CDF remain exact, Taylor's theorem and the Gaussian density give

\[
|\Delta J|\leq
\frac{\lambda}{2\sqrt{2\pi}}
\{e^{\lambda\eta_0}\eta_0^2+
(1+h)e^{\lambda\eta_1}\eta_1^2\}.
\tag{7}
\]

This is a threshold-only second-order certificate. It does not make a genuine log-level perturbation second order, and does not remove first-order error from inconsistent exponential coefficients, rounded tilted offsets or CDF interpolation. The scalar kernel avoids implementing those separately perturbed pieces in the reversible evaluator. It therefore supplies a stronger opponent without requiring me to claim that the old cancellation-prone arithmetic is free of error.

## 7. Resource and interpretation boundary

I retain the two exponential evaluations, arithmetic average, logarithm and geometric log formation that construct the scalar kernel inputs. I count both one-dimensional cap-specific kernel lookups, their interpolation arithmetic, subtraction, projection, comparison and the exact inverse. I also retain cap-specific table identity, word counts, program bits, generator precision and measured classical setup work. This kernel table is a function of a realised scalar conditional log level, not a joint table over the market-input indices or a shortcut using the Gaussian index as a payoff address.

The resulting comparison remains a realised architecture comparison. A favourable result against the old eight-CDF formula would establish that the revised full-conditioning benchmark is stronger, not that scalar-kernel tabulation is new or that either quantum construction has an end-to-end advantage. Independent arithmetic precision selection, complete distribution-transfer allowances and a separately certified geometric correction remain necessary.
