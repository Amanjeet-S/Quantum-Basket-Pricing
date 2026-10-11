# Independent review of the broader financial applications

## 1 Scope and verdict

I reviewed financial_applications_proof.tex and general_dispersion_study.py in the root staging directory. The risk-neutral drift algebra, common-factor covariance criterion, gamma-normal mixture moment domain and dispersion constants are mathematically correct under the stated hypotheses. The code evaluates the same symmetric financial residual and retains its continuous quadrature status. It does not establish certified numerical integration, a finite-grid limit or empirical validation.

## 2 Risk-neutral shift and kernel hypotheses

Writing ψ(v)=log E exp(v·Y), the shift μ_i=r−d−λ²/2−ψ(R₀+εb_i) gives the required terminal mean because the independent normal shock contributes exp(λ²/2). Taylor's integral remainder gives exactly ψ(R₀+εb_i)=ψ(R₀)+ε∇ψ(R₀)·b_i+ε²e_i(ε). Subtracting the weighted drift leaves D_i=εb_i·(Y−∇ψ(R₀))−ε²(e_i−Σ_k w_k e_k). These intercepts are weighted centred and uniformly of order ε². The geometric base differs from m₀ by −ε²Σ_iw_ie_i, so it converges without needing individual limits for e_i. The general dispersion theorem therefore applies to V=Y−∇ψ(R₀) under its explicit norm-exponential moment hypothesis.

The proposition should replace “a suitable bounded log-kernel” by an explicit reference to the hypotheses of Theorem general-weak-dispersion: P maps the real line into [0,1], is increasing and continuously differentiable, and has globally bounded derivative. Positive weights summing to one should also be stated or explicitly inherited. For the existing capped-lognormal P with λ>0 those hypotheses have already been proved. The correct terminal normalisation does not itself specify or calibrate a time-dependent financial process; the current text appropriately presents a stipulated terminal pricing family.

## 3 Common-factor covariance

For a non-singular covariance matrix C, positive semidefiniteness implies positive definiteness. The displayed condition C−λ²11ᵀ≥0 is equivalent to λ²≤1/(1ᵀC⁻¹1). Cauchy–Schwarz after inserting C^(1/2) and C^(−1/2) proves sufficiency; x=C⁻¹1 proves necessity. A residual square root and independent standard normals then give exactly the required covariance and common loading. At equality the residual can be singular, which is permitted. The proof's separate warning about singular original C is necessary and correct. Arbitrary deterministic monitoring means belong in the affine constants; this application does not incorrectly impose a one-year drift on every monitored price.

## 4 Gamma-normal mixture and limiting constants

For G with shape 1/ν and scale ν, conditioning on G and rescaling the gamma integral gives M_Y(t)=(1−νt²/2)^(−1/ν) exactly for |t|<sqrt(2/ν). The integral diverges at and outside the boundary. Gamma recursion gives E G=1 and E G²=1+ν, hence E Y²=1 and E Y⁴=3(1+ν). At ν=0.4 the fourth moment is 4.2 and Var(Y²)=3.2. The stated bound E exp(c|Y|)≤M_Y(c)+M_Y(−c), followed by a strictly positive exponential margin absorbing |Y|⁴, proves the actual polynomial-times-exponential hypothesis at c=0.4. It is stronger than merely checking the two moment arguments used in the finite second-moment formula.

For equal weights and opposite loadings ±0.2ε, the two drifts are equal for every symmetric latent law. Consequently T=cosh(0.2εY), J=log T and the base is deterministic. The coefficients E H₀=0.02P′(m₀) and Var(H₀)=0.0004P′(m₀)²Var(Y²) are correct. A standardised uniform has fourth moment 9/5, giving variance factor 0.8; a standard normal gives 2; the mixture gives 3.2; a symmetric two-point law gives zero. In the last law the complete fractional residual is constant, since both the base and cosh(0.2εY) are constant. Its Bernoulli flag can nevertheless have positive variance.

## 5 Diagnostic implementation

The four implemented MGFs match the analytical laws. The second-moment expression M_Y(2t)/2−2M_Y(t)+3/2 equals E(cosh(tY)−1)². Its ε⁴ limit is b⁴E Y⁴/4. The financial drift uses each law's own MGF and the decimal stock-mean check verifies the terminal normalisation. The Laguerre Jacobi coefficients correspond to gamma shape 5/2; scaling its nodes by 0.4 and squaring the first eigenvector components gives the normalised positive gamma quadrature. Combining it with the normal quadrature matches the gamma-normal mixture.

The positive density formula integrates the derivative along the log-gap and the uniform strike coordinate. The factors gap/λ, the two half-interval weights, λ=0.2 and h=0.1 are correct. The formula log1p(2sinh(z/2)²)=log cosh(z) avoids subtraction of nearly equal payoffs at small dispersion.

One diagnostic limitation remains: the reported coarse/fine comparison refines the latent-law and strike orders, but keeps the along-gap order at 12. A shared error in that unrefined dimension would not appear in the refinement difference. I recommend exposing its order and using 12/24 or 16/32 in the corresponding coarse/fine runs, while continuing to label all quadrature output diagnostic. This does not affect the analytical theorem or the separate 273 finite-output confidence certificates.
