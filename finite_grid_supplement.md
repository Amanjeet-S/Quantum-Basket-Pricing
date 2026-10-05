# Finite-grid and coherent-evaluator certificates for a capped basket residual

## 1. Keywords

Keywords: capped basket option; geometric control variate; bounded variation; midpoint quadrature; conditional expectation; fixed-point arithmetic; coherent comparison; quantum resource accounting.

## 2. Target and status of the argument

I compare three constructive evaluators for the same continuous stipulated lognormal model. I write the arithmetic and geometric averages as

\[
A(y,z)=\sum_{j=1}^{n}w_j\exp(\ell_j(y)+\lambda z),
\qquad
G(y,z)=\exp\!\left(\sum_{j=1}^{n}w_j\ell_j(y)+\lambda z\right),
\]

where \(w_j>0\), \(\sum_jw_j=1\), \(\ell_j(y)=\mu_j+r_j^{\mathsf T}y\), \(\lambda>0\), and the coordinates of \(Y\) and \(Z\) are independent standard normal variables. The contract is

\[
F(a)=\min\{(a-K)^+,H\},\qquad
R(y,z)=F(A(y,z))-F(G(y,z)).
\]

The arithmetic–geometric mean inequality and monotonicity of \(F\) give \(0\leq R\leq H\). I recover the continuous undiscounted basket value from \(\mathbb E F(G)+\mathbb E R\), and multiply the complete sum and its monetary error allowance by the discount factor.

The direct evaluator uses \(R(Y,Z)\). The fully conditioned evaluator uses \(g_R(Y)=\mathbb E_ZR(Y,Z)\). The conditional-strike evaluator uses

\[
h_R(y,u)=H\left[
\Phi\!\left(\frac{\log B_A(y)-\log(K+Hu)}{\lambda}\right)
-
\Phi\!\left(\frac{\log B_G(y)-\log(K+Hu)}{\lambda}\right)
\right],
\]

where \(B_A(y)=\sum_jw_je^{\ell_j(y)}\) and \(B_G(y)=\exp(\sum_jw_j\ell_j(y))\). The layer-cake identity implies \(\mathbb E_Uh_R(y,U)=g_R(y)\) for a uniform \(U\). Common-loading analytical conditioning and its computational advantages are established antecedents, rather than proposed theoretical contributions; the retained source record identifies the author manuscript actually consulted for [Bayer, Siebenmorgen and Tempone, Quantitative Finance 18, 491–505 (2018)](https://doi.org/10.1080/14697688.2017.1308003).

The following proofs are elementary project-specific derivations used to certify these constructors. I distinguish them from a claim to a new quadrature theorem, a new conditioning identity or a quantum advantage. Total variation below means variation of a real-valued function along a coordinate. It does not mean total variation distance between a continuous law and a finite register law.

## 3. Midpoint quadrature for a function of bounded variation

For a function \(f:[0,1]\to\mathbb R\), let \(\operatorname{Var}(f)\) denote its total variation. With \(N\) equal cells and midpoints \(u_j=(j+1/2)/N\), I use

\[
Q_Nf=\frac1N\sum_{j=0}^{N-1}f(u_j).
\]

Put \(a_j=j/N\) and \(b_j=(j+1)/N\). On the left half of a cell, \(|f(u)-f(u_j)|\leq\operatorname{Var}_{[a_j,u_j]}(f)\); on its right half the analogous bound is \(\operatorname{Var}_{[u_j,b_j]}(f)\). Integrating these inequalities gives

\[
\left|\int_{a_j}^{b_j}f(u)\,du-\frac{f(u_j)}N\right|
\leq
\frac{1}{2N}\operatorname{Var}_{[a_j,b_j]}(f).
\]

Adding over cells proves

\[
\left|\int_0^1f(u)\,du-Q_Nf\right|
\leq\frac{\operatorname{Var}(f)}{2N}.
\tag{1}
\]

Endpoint limits suffice when a Gaussian quantile transformation initially defines \(f\) only on \((0,1)\). Bounded monotone branches have such limits. A monotone change of argument cannot increase total variation, so composition with \(\Phi^{-1}\) preserves the bounds used below. Equation (1) is the familiar bounded-variation midpoint argument, included here as a complete proof of the exact constant needed for the finite register target.

## 4. Outer Gaussian coordinates

I fix every variable except one outer coordinate \(y_i\). The arithmetic average then has the form

\[
A(y_i)=\sum_j c_j e^{r_{ji}y_i},\qquad c_j>0,
\]

and

\[
A''(y_i)=\sum_jc_jr_{ji}^{\,2}e^{r_{ji}y_i}\geq0.
\]

Its derivative is nondecreasing. Therefore \(A\) is nonincreasing and then nondecreasing, allowing either branch to be empty or the whole function to be constant. If \(\psi\) is nondecreasing and takes values in \([0,H]\), the composition \(\psi(A(y_i))\) has variation at most \(H\) on each branch, hence at most \(2H\) overall. The geometric average is exponential affine in \(y_i\); it is monotone, so \(\psi(G(y_i))\) has variation at most \(H\). Subadditivity gives

\[
\operatorname{Var}_{y_i}\{\psi(A)-\psi(G)\}\leq3H.
\tag{2}
\]

For the direct evaluator I take \(\psi=F\) and hold \(z\) fixed. For the conditional-strike evaluator I hold \(u\) fixed and take

\[
\psi(a)=H\Phi\!\left(\frac{\log a-\log(K+Hu)}{\lambda}\right).
\]

This function is nondecreasing and bounded by zero and \(H\), so the same argument applies. For the fully conditioned evaluator I can either use the nondecreasing bounded function \(\psi(a)=\mathbb E_ZF(ae^{\lambda Z})\), or apply

\[
\operatorname{Var}(\mathbb E_Zf_Z)
\leq\mathbb E_Z\operatorname{Var}(f_Z),
\]

which follows by taking any finite partition, applying the triangle inequality to its increments and then taking the supremum.

Let coordinate \(i\) use \(N_i\) midpoint Gaussian quantiles. The tensor integration and midpoint operators are averages and have supremum-norm operator norm at most one. Telescoping the product of these operators, and applying (1) and (2) to every coordinate slice, yields

\[
\beta_{\mathrm{outer}}
\leq\sum_{i=1}^{p}\frac{3H}{2N_i}.
\tag{3}
\]

For equal grids this is \(3pH/(2N_Y)\). It is uniform in any retained inner Gaussian or strike coordinate. It therefore applies to all three evaluator targets before their additional inner discretisation allowances are added.

The direct outer proof needs positive basket weights, but does not require common loading. Holding all other Gaussian factors fixed still leaves a positive sum of exponentials in one coordinate. Conditioning also preserves the variation bound through the expectation inequality just proved. Common loading is used to construct the inexpensive root-free conditional evaluator; it is not needed for this outer quadrature certificate. I retain this broader scope as a mathematical qualification, rather than a theoretical novelty claim.

## 5. Inner Gaussian and strike coordinates

For fixed \(y\), both \(F(B_A(y)e^{\lambda z})\) and \(F(B_G(y)e^{\lambda z})\) are nondecreasing, bounded functions of \(z\), each with variation at most \(H\). Their difference therefore has variation at most \(2H\). Composition with the normal quantile and (1) give the direct evaluator's additional inner allowance

\[
\beta_{\mathrm{direct},Z}\leq\frac{H}{N_Z}.
\tag{4}
\]

For fixed \(y\), each conditional digital component in \(h_R(y,u)\) is nonincreasing in \(u\), bounded by zero and \(H\). The residual has variation at most \(2H\), giving

\[
\beta_{\mathrm{strike}}\leq\frac{H}{N_U}.
\tag{5}
\]

The single arithmetic digital component alone has the sharper bound \(H/(2N_U)\). I charge both components in the controlled residual, rather than silently transferring the single-component constant.

For equal outer grids, complete distribution-transfer allowances are consequently \(3pH/(2N_Y)+H/N_Z\) for direct evaluation, \(3pH/(2N_Y)\) for full conditioning, and \(3pH/(2N_Y)+H/N_U\) for conditional-strike evaluation. These bounds refer to exact midpoint quantiles and exact functions on the finite grid. Quantile-coordinate representation and function arithmetic are separate errors; I charge them through the evaluator certificate rather than count them a second time here.

## 6. Relation to the coupling certificate

I retain coupling as an independently valid alternative where its Lipschitz constants give a smaller result. Let \(q(u)=\Phi^{-1}(u)\), let \(q_R\) be its clipping to \([-R,R]\), and assume every midpoint quantile lies in this interval. The clipped quantile is monotone with variation \(2R\). The half-cell argument used in section 3, now applied to absolute deviations, gives

\[
\int_0^1|q_R(u)-q_R(u_{\lfloor Nu\rfloor})|\,du
\leq\frac{R}{N}.
\]

The remaining clipping error is

\[
\mathbb E(|Z|-R)^+
=2\{\phi(R)-R[1-\Phi(R)]\}.
\]

Integration by parts gives

\[
1-\Phi(R)=\frac{\phi(R)}R-
\int_R^\infty\frac{\phi(t)}{t^2}\,dt
\geq\frac{\phi(R)}R-\frac{1-\Phi(R)}{R^2},
\]

so \(1-\Phi(R)\geq R\phi(R)/(R^2+1)\). Therefore

\[
\mathbb E|Z-Z_N|
\leq D_R(N):=\frac RN+\frac{2\phi(R)}{R^2+1}.
\tag{6}
\]

I use \(R=5\) only after establishing that all midpoint quantiles in the searched grid widths lie strictly inside that domain. In the two-asset instance, the separately checked \(\Phi(-4.75)<1/(2\,2^{18})\) establishes \(|z_j|<4.75\) for all Gaussian index widths up to eighteen bits.

For coordinate \(i\), an outer residual Lipschitz constant is

\[
L_i=(K+H)\left(\max_j|r_{ji}|+
\left|\sum_jw_jr_{ji}\right|\right).
\]

Thus \(\sum_iL_iD_5(N_i)\) is a valid outer coupling allowance for the fully averaged residual. Direct inner coupling gives \(2(K+H)\lambda D_5(N_Z)\). I may use the minimum of complete applicable coupling and variation bounds. For a finite strike target I first transfer the continuously strike-averaged residual in \(Y\), then add (5); I do not assert that the same Lipschitz constant holds pointwise for an individual digital probability.

Equation (6) has a positive fixed tail term and can be substantially looser than (3) to (5). The variation argument does not discard Gaussian tails: it transforms the entire continuous Gaussian integral to \((0,1)\) and bounds the midpoint quadrature directly.

## 7. Fixed-point arithmetic, projection and normalisation

I represent fixed-point values by integers with fractional step \(\delta=2^{-q}\). The implementation uses nonnegative biased coordinates for elementary functions and modular integer registers for arithmetic. The six integer bits give a signed interpretation in \([-32,32)\), although many intermediate registers are deliberately interpreted as unsigned. Every signed intermediate and every nonnegative function or multiplication input requires an explicit range certificate. A small floating-point error allowance cannot repair modular overflow.

Integer addition and subtraction at a common scale are exact. A nonnegative product followed by removal of \(q\) low bits implements floor truncation and contributes less than \(\delta\) of additional absolute error. A right shift used for a weighted average has its corresponding floor allowance. If exact operands \(a,b\) have approximation errors \(\epsilon_a,\epsilon_b\), a conservative nonnegative-product allowance is

\[
b_{\max}\epsilon_a+a_{\max}\epsilon_b+
\epsilon_a\epsilon_b+\delta.
\]

Coefficient approximation is charged separately. In the instantiated two-asset model, I scale currency values by \(100\), use \(h=H/100\), and convert the residual to probability through division by the fixed exact contract cap. Approximating the constants \(h\) or \(1/h\) does not change the declared contract: their representation errors are propagated to the final normalised output.

For a twice differentiable primitive, linear interpolation over a cell of width \(\Delta\) contributes at most \(\Delta^2\sup|f''|/8\). An ordinate rounded to the nearest fixed-point value contributes at most \(\delta/2\), in addition to its high-precision generator enclosure. The implemented interpolation product adds a floor allowance below \(\delta\). I therefore include the interpolation term, the generator allowance and \(3\delta/2\), together with the error in the primitive's argument. For the normal CDF, the global derivative and curvature bounds are \(1/\sqrt{2\pi}\) and \(1/\sqrt{2\pi e}\).

I use the finite domains \([-2,2]\) for the exponential, \([0.25,4.25]\) for the logarithm, and \([-8,8]\) for the CDF. The logarithm's stored output is \(\log x+2\). Range checks apply to the perturbed implemented arguments, not just nominal values. The interpolation table must be generated and enclosed over all its ordinates; it is a one-dimensional primitive table, not a joint market-input payoff lookup.

Projection to an exact interval is non-expansive: for \(x\) in that interval, \(|P(\widehat x)-x|\leq|\widehat x-x|\). A rounded cap contributes its additional threshold error, at most \(\delta/2\). I use this property for individual capped payoffs, analytical conditional values and the final residual, while constructing and counting their comparators and selectors. I do not use projection to justify an uncharged domain error.

In the full-conditioning constructor, subtraction of equal-scale fixed-point integers preserves the represented fractional bits. I nevertheless propagate the errors in the ordinary and tilted CDF differences and their multiplying coefficients, because cancellation in the exact formula can amplify those input errors. Ordered CDF arguments and a monotone integer interpolator guarantee nonnegative probability differences before unsigned multiplication. This structural guarantee must accompany the range checks; a nominal nonnegative real bound alone is insufficient.

## 8. Normal-CDF generator and inverse-quantile enclosure

The generator uses seventy significant decimal digits and

\[
\Phi(x)=\frac12+
\frac1{\sqrt{2\pi}}\sum_{n=0}^{\infty}
\frac{(-1)^n x^{2n+1}}{2^n n!(2n+1)},
\qquad |x|\leq8.
\tag{7}
\]

Writing the summand as \(t_n\), its recurrence is

\[
t_0=x,\qquad
t_n=-t_{n-1}\,
\frac{x^2(2n-1)}{2n(2n+1)}.
\]

Once \(n>x^2/2\), the magnitudes decrease. The implementation waits beyond this point and stops only when the last included term has magnitude below \(10^{-50}\). The alternating-series remainder is bounded by the first omitted term, which is smaller still. At \(|x|=8\), the implementation requires at most 167 recurrence iterations; the last exact term has magnitude below \(3.64\,10^{-51}\). Smaller absolute arguments are dominated term by term by the \(|x|=8\) series. This gives ample margin for the rounding perturbation in the stopping test.

The sum of absolute terms is bounded by

\[
\sum_{n\geq0}|t_n|
\leq |x|\exp(x^2/2)
\leq8e^{32}<6.32\,10^{14}.
\tag{8}
\]

I use a generous upper bound of 4,000 decimal arithmetic operations for this CDF computation and its normalising constant. This covers the recurrence, accumulation and loop-test arithmetic, together with the Machin-formula calculation of \(\pi\), its arctangent series and the square root. It replaces the insufficient literal bound of 500 operations.

For rounded seventy-digit decimal operations, I can use relative unit roundoff \(u=5\,10^{-70}\), with no underflow or overflow in this domain. The ordinary forward-product bound

\[
\gamma_{4000}=\frac{4000u}{1-4000u}<2.01\,10^{-66}
\]

controls recurrence coefficient products. Accumulation errors obey the same conservative order after multiplication by the absolute-term bound (8). The resulting series-rounding error is below \(1.3\,10^{-51}\). The Machin formula

\[
\pi=16\arctan(1/5)-4\arctan(1/239)
\]

uses decreasing alternating series whose stopping and first-omitted-term bounds are below \(10^{-68}\) before the fixed coefficients are applied. Including seventy-digit input division, accumulation and square-root errors gives a conservative absolute enclosure below \(10^{-64}\) for the computed \(1/\sqrt{2\pi}\). Multiplying this allowance by (8), and adding the alternating remainder and final operations, leaves the absolute CDF error below \(10^{-48}\). The retained \(10^{-38}\) enclosure is therefore deliberately wider. Agreement with a double-precision error function and CDF symmetry are useful diagnostics, but do not replace this bound.

For a midpoint probability \(p_j=(j+1/2)/2^b\), let \(\widehat z_j\) be the rounded quantile coordinate before adding its storage bias. I require the exhaustive enclosure check

\[
\widehat\Phi(\widehat z_j-\delta)+10^{-38}
\leq p_j
\leq
\widehat\Phi(\widehat z_j+\delta)-10^{-38}.
\tag{9}
\]

The CDF enclosure and strict monotonicity imply
\(|\widehat z_j-\Phi^{-1}(p_j)|\leq\delta\). I check every nonnegative-half quantile and verify exact integer reflection of the negative half; CDF symmetry then proves the remaining enclosures. Each distinct index-width and fractional-width table used by a retained constructor must satisfy (9). A floating-point inverse-CDF seed, an unchecked table or a precision-search-only run does not establish this coordinate certificate.

## 9. Exact comparison flag and clean uncomputation

Let the compiled arithmetic return an integer \(T(x)\in\{0,\ldots,2^q\}\) for the normalised residual on model input \(x\). I prepare a \(q\)-bit comparison register uniformly over \(r=0,\ldots,2^q-1\), and XOR the predicate \(r<T(x)\) into a flag. I append a fixed zero high bit to \(r\) and compare against the \(q+1\) low bits of the threshold. This represents \(T=2^q\) correctly rather than wrapping probability one to zero. Exactly \(T(x)\) comparison values satisfy the predicate, hence

\[
\Pr(\mathrm{flag}=1\mid x)=\frac{T(x)}{2^q}.
\tag{10}
\]

The comparison itself introduces no additional rounding. Financial approximation has already been charged through the difference between the exact normalised residual and \(T(x)/2^q\). If I chose instead to floor a separately evaluated real probability, that floor would require its own allowance; I do not charge the same output quantisation twice.

I compute the threshold reversibly, perform a clean comparison that preserves the threshold, and apply the exact inverse of the entire arithmetic body. This leaves the model indices, comparison integer and flag, with all evaluator work restored to zero. This comparison-register encoding is established in the quantum-finance literature; the final published paper actually consulted is [Stamatopoulos and Zeng, Quantum 8, 1322 (2024)](https://quantum-journal.org/papers/q-2024-04-30-1322/), sections 4, 5.1, 6, 7 and Appendix A. The project's requirement is to instantiate, certify and count its whole evaluator rather than treat a financial function as a free oracle.

The selected arithmetic library emits X, CNOT and exact Toffoli gates. I declare the seven-T, T-depth-three Toffoli implementation as a resource schedule, not an optimal decomposition. The supporting source is [Selinger, Physical Review A 87, 042302 (2013)](https://doi.org/10.1103/PhysRevA.87.042302). Detailed material actually consulted was the [author version, arXiv:1210.0974v2, dated 3 April 2013](https://www.mathstat.dal.ca/~selinger/papers/toffoli.pdf), sections 1 to 3; I do not imply that the final journal text was read.

Every quantile table, primitive table, coefficient transform, product, cap, flag comparison and inverse is part of the constructor. The current exact prefix-tree table implementation is an explicit selected architecture, with all its gates counted. Published unary-iteration QROM offers a stronger alternative: [Babbush et al., Physical Review X 8, 041015 (2018)](https://doi.org/10.1103/PhysRevX.8.041015), sections III.A and III.C and figures 3 to 10, gives the relevant four-T logical-AND and \(4L-4\) T-count construction for an \(L\)-row lookup, under its stated ancillary and uncomputation assumptions. The final journal PDF was actually consulted through this [accessible retrieval copy](https://inspirehep.net/files/8e4f34cd299ada4bdc0b8441dbfb73af). I do not transfer that count to a different circuit without implementing its assumptions.

## 10. Reflections on the clean-work invariant subspace

Let \(X\) contain all finite model and strike indices, \(Q\) the comparison register, \(F\) the flag, and \(W\) the evaluator work. If \(E\) is the arithmetic body and \(C\) the clean comparison, the flag evaluator is \(B=E^\dagger C E\). For every basis input and either initial flag value, \(B\) preserves \(X,Q\), toggles \(F\) by the declared predicate, and returns \(W\) to zero. By linearity it preserves the invariant subspace \(\mathcal H_{XQF}\otimes|0_W\rangle\).

The preparation unitary is \(A=B(H_X\otimes H_Q)\), acting on an initially zero flag and work. Its inverse also preserves the clean-work subspace. The good-state reflection is the Pauli Z on the flag. On this invariant subspace, the initial-state reflection may act only on \(X,Q,F\):

\[
S_0=
\left(I-2|0_{XQF}\rangle\langle0_{XQF}|\right)\otimes I_W.
\]

It agrees there with a reflection about the full all-zero state. Thus a conventional Grover iterate formed from \(A\), \(A^\dagger\), \(S_0\) and the flag reflection is valid without controlling the zero test on every arithmetic-work qubit. This statement relies on exact work reset; a dirty evaluator does not justify the reduced reflection.

The reflection is still an additional circuit. For \(m\geq3\) reflected register qubits, a straightforward all-zero phase construction uses X gates around an all-one test, an AND ladder on the first \(m-1\) bits, a controlled Z with the remaining bit, and the inverse ladder. It uses \(2(m-2)\) exact Toffolis, \(m-2\) clean auxiliary bits, \(2m\) X gates and one controlled Z. This is another explicit upper-bound construction, not a minimum. Initial uniform preparation contributes one parallel Hadamard layer. Controlled iterates, confidence procedures and all estimator repetitions require their own accounting; the present evaluator resource totals do not include them.

Because the compiled Boolean circuit is exact for its declared fixed tables and integers, the deterministic financial error is a bias in the represented mean. It is charged once in the price certificate. I do not multiply that bias by the number of ideal estimator calls. Physical gate error and approximate unitary synthesis would create different accumulation questions and are outside these logical constructor bounds.

## 11. Contribution and comparison limits

I retain a positive constructor when every table ordinate and inverse quantile is enclosed, every realised arithmetic range is checked, the entire flag evaluator and inverse are specified, and its resources refer to actual emitted blocks or a proved serial composition. I retain source hashes for both study scripts and their arithmetic dependency, and round monetary upper bounds outwards when exporting them to binary floating point. Numerical basis replay and sampled financial errors support implementation checks, but do not replace the uniform analytical certificate.

This construction does not establish the strongest classical or quantum baseline, optimal precision allocation or an end-to-end advantage. Stable analytical conditioning, positive numerical integration, range reduction, alternative reversible arithmetic and more efficient QROM can improve the opponents. The conditional-strike construction can also increase estimator variance relative to full conditioning. Complete mean-estimation resources still require preparation and reflection costs, repeated and controlled evaluator calls where used, a valid accuracy and confidence rule, and a separately certified geometric correction. I therefore interpret any favourable result as a comparison within the realised constructor family and stipulated synthetic model. I treat the checked implementation and its retained evidence as a research deliverable. A contribution to the literature remains contingent on a literature comparison and an informative, reproducible project-specific result.
