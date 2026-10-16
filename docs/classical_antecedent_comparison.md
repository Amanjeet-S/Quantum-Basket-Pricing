# Residual certificate comparison with classical pricing antecedents

## 1. Scope and result

I compare the exact second-correction grid certificate in `tensor_transfer_proof.tex` with four established basket-pricing sources. The result is a narrowing of the claim. Gaussian conditioning, conditional-price expansion, control subtraction and finite Gaussian moments are established. The checked project statement concerns the actual rounded-grid transfer error of a specific residual. The source versions below do not establish its theoretical priority or superiority over the cited pricing algorithms.

## 2. Published polynomial-pricing antecedent

Olivares, P. and Alvarez, A. (2016), *Pricing Basket Options by Polynomial Approximations*, Journal of Applied Mathematics, article 9747394, DOI [10.1155/2016/9747394](https://onlinelibrary.wiley.com/doi/10.1155/2016/9747394). I read final publisher HTML, Sections 1–5, Proposition 1, equations 6–18 and Theorem 6, equations 39–54. It strips many formula bodies. The [Wiley-provided final journal body on ResearchGate](https://www.researchgate.net/publication/261597975_Pricing_of_Basket_Options_Using_Polynomial_Approximations) supplies mathematical OCR, with damaged glyphs. I additionally read the [2014 author preprint](https://arxiv.org/pdf/1404.3229), Section 2, equations 5–16 and Proposition 1. Its numbering and content are distinct from the final article.

The final article conditions on the last asset log-returns under constant-covariance Black–Scholes dynamics, represents the remaining conditional contract analytically, expands the conditional price, and integrates Gaussian exponential-power moments. It also uses truncated moments in its Chebyshev construction. Spread examples admit negative weights. This is a direct antecedent for the project's broad architecture. It does not justify a claim that conditional polynomial pricing or truncated Gaussian moments are new. The final Chebyshev formula coefficients have not all been independently transcribed and checked here.

## 3. Further published antecedents and access limits

Ju, N. (2002), *Pricing Asian and basket options via Taylor expansion*, Journal of Computational Finance, 5(3), 79–103, DOI 10.21314/JCF.2002.088. I read the [final publisher abstract and issue record](https://www.risk.net/journal-of-computational-finance/2160548/pricing-asian-and-basket-options-via-taylor-expansion). The DOI is recorded in the final Olivares–Alvarez reference 12. Ju uses a moment-matched lognormal reference and a characteristic-function ratio expanded around zero volatility through order six. This differs from the project's exact conditional multiplicative residual identity. The final body remained inaccessible, so I cannot exclude overlapping coefficient, remainder or error arguments.

Curran, M. (1994), *Valuing Asian and Portfolio Options by Conditioning on the Geometric Mean Price*, Management Science, 40(12), 1705–1711, DOI [10.1287/mnsc.40.12.1705](https://pubsonline.informs.org/doi/10.1287/mnsc.40.12.1705). I first read final metadata and abstract. Safari subsequently obtained the [author-uploaded final journal scan](https://www.researchgate.net/publication/333906878_Valuing_Asian_and_portfolio_options_by_conditioning_on_the_g). I visually read all seven printed pages, 1705–1711, including Sections 2.1–2.3 and equations 1–4. Gaussian conditioning gives conditional covariance \(C-Cww^{\mathsf T}C/(w^{\mathsf T}Cw)\). The component where the geometric average satisfies \(G\geq K\) is handled exactly. Below the strike, Curran approximates the conditional arithmetic-to-geometric gap by a lognormal law with moments frozen at \(G=K\); the conditional-mean construction gives a Jensen lower bound. The project instead conditions on a common independent Gaussian factor and retains an exact whole-domain integral remainder \(R_2\). This distinguishes targets and exactness, while leaving conditioning and computable basket moments as established mechanisms.

Dingeç, K.D. and Hörmann, W. (2013), *Control variates and conditional Monte Carlo for basket and Asian options*, Insurance: Mathematics and Economics, 52(3), 421–434, DOI [10.1016/j.insmatheco.2013.03.002](https://www.sciencedirect.com/science/article/pii/S0167668713000309). I read final abstract, introduction and snippets for Sections 2–5 and conclusions. These establish a GBM geometric-exercise-region control combined with conditional Monte Carlo and quadratic controls. The numerical comparisons concern variance reduction. Formula bodies were mostly stripped and the full final text was not read. Neither generic conditional variance reduction nor quadratic financial controls are original here. I have not performed a matched implementation comparison with this algorithm.

## 4. Exact project distinction

In the project's notation, the weighted arithmetic-to-geometric ratio is (T\geq1), (U=T-1), and the conditional geometric log-price is (m). With (t=1+zU), the exact residual is

\[
R_2(y)=U(y)^2\int_0^1(1-z)t^{-2}C_2(m(y)+\log t)\,dz.
\]

The project's derivative calculation bounds its first coordinate derivative by terms (U|T_j|,U^2,U^2|T_j|), and its second derivative by terms (T_j^2,UT_{jj},U|T_j|,UT_j^2,U^2,U^2|T_j|,U^2T_{jj},U^2T_j^2). These follow by product and chain rules, (t\geq1), bounded (C_2,C_3,C_4) and their relevant derivatives. Cauchy's inequality under each finite positive transport measure replaces the absolute-value terms by nonnegative polynomial moments through degree four. This is a residual-specific certified calculation, rather than a new conditional-price Taylor principle.

For (X_i=\exp(d_i+s_i\cdot y)), any resulting monomial has the elementary identity

\[
\int\prod_iX_i^{\alpha_i}\,d\Bigl(\bigotimes_j\eta_j\Bigr)
=\exp\Bigl(\sum_i\alpha_i d_i\Bigr)
 \prod_j\mathcal M_{\eta_j}\Bigl(\sum_i\alpha_i s_{ij}\Bigr).
\]

Finite sums and integrability justify Fubini. The project uses actual truncated Gaussian, discrete rounded-node and positive transport factors, rather than only a multivariate Gaussian probability law. This change supplies the grid-transfer certificate. It does not make the product-MGF identity itself novel. The signed centred terms and unnormalised masses must remain intact, because the target is the difference between the continuous expectation and the finite-grid expectation.

## 5. Contribution judgement

The extra checked conclusion is a price-independent, explicit error allowance for that exact residual and rounded grid, obtained without Cartesian enumeration through degree-four exponential moments, retained tail terms and outward arithmetic. This is useful certified numerical work. Its degree bound is explained by the particular residual's derivative algebra and standard moment inequalities. The present comparisons do not establish a fundamentally new pricing model, a general Gaussian integration theory, literature priority, or a physical quantum advantage. The final Ju and Dingeç–Hörmann mathematical bodies remain needed for a complete proposition-level exclusion. Curran’s final body has now been read, and its exercise-region partition and below-strike approximation have been compared with the exact integral residual.

## 6. Retrieval record

On 8 October 2026, the paper-index search failed with a DNS error; I used primary web retrieval. Wiley/Hindawi final PDF and XML retrievals failed, while final HTML and the Wiley-provided indexed journal body were accessible. Local arXiv download failed through DNS; the 2014 PDF text remained accessible through a separate retrieval route. I subsequently obtained Curran’s seven-page author-uploaded final scan in Safari and read it visually. The separate access ledger retains its file hash. No copyrighted PDF or copied OCR is installed with this comparison. No subscription, paid content or additional entitlement was acquired. This record preserves remaining access limits rather than treating them as evidence of non-overlap.
