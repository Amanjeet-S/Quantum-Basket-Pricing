# 1. Scope and access record

I checked the closest accessible polynomial and bounded basket-pricing predecessors against the proposed cancellation-preserving Gaussian derivative norm. The comparison identifies established components and a narrower testable extension. It does not establish theoretical priority or identify a universally strongest classical pricer. All access used free sources or my existing Safari access at no additional cost. I attempted a paper-index search first; it failed with a DNS error, and primary retrieval and Safari supplied the accessible bodies. Checked on 7 October 2026.

# 2. Source versions actually consulted

## 2.1 Olivares and Alvarez (2016)

Olivares, P. and Alvarez, A. (2016). Pricing basket options by polynomial approximations. *Journal of Applied Mathematics*, 2016, article 9747394, 12 pages. [DOI](https://doi.org/10.1155/2016/9747394). I read the final publisher HTML and, through Safari, the final publisher ePDF, particularly Section 4, equations (39) to (54), Theorem 6 and its proof, and Section 5, pp. 9 to 10. The publisher records received 2 May 2016, revised 21 July 2016, accepted 26 July 2016 and first published 5 October 2016. The final paper supersedes the closely related 2014 Taylor and polynomial preprints for the present comparison.

The method approximates a conditional Black–Scholes price, including a truncated Chebyshev construction and Gaussian exponential-power moments. Theorem 6 evaluates the polynomial expectation using a truncated Gaussian moment generating function. Section 5 chooses the finite interval heuristically and suggests endpoint continuation. These passages do not provide a retained, outward-rounded joint certificate for tail omission, polynomial approximation and coefficient arithmetic. I therefore attribute conditional Gaussian polynomial pricing and truncated Gaussian moments to this published antecedent, while treating the proposed norm certificate as a different object.

Primary consulted body: [publisher ePDF](https://onlinelibrary.wiley.com/doi/epdf/10.1155/2016/9747394).

## 2.2 Caldana et al. (2016)

Caldana, R., Fusai, G., Gnoatto, A. and Grasselli, M. (2016). General closed-form basket option pricing bounds. *Quantitative Finance*, 16(4), 535–554. [DOI](https://doi.org/10.1080/14697688.2015.1073854). I read the City repository accepted manuscript through primary retrieval and Safari, particularly Section 2, Proposition 2.1 and Corollary 2.2, Section 3, and Appendices B and D. The repository explicitly calls it the accepted version; its running heading is dated 20 November 2014. I did not read the final publisher body during this check.

Corollary 2.2 gives the positive-basket call interval from the geometric-basket call and the arithmetic-geometric expectation gap. Appendix B proves the underlying pathwise inequalities; Section 3 supplies lognormal formulas. Section 1 gives an optimised event-set lower bound and a control variate. These are established bounded-pricing mechanisms. Independently subtracting two call intervals for a capped call loses information shared by the two strikes, so a direct capped-payoff certificate remains worth comparing. The cited text is not evidence that the new derivative norm or its implementation is original.

Primary consulted body: [accepted manuscript](https://openaccess.city.ac.uk/id/eprint/15729/7/basket_option_qf.pdf).

## 2.3 Carr (2001): direct derivative-operator antecedent

Carr, P. (2001). Deriving derivatives of derivative securities. *Journal of Computational Finance*, 4(2), Winter 2000/2001. [DOI](https://doi.org/10.21314/JCF.2001.075). I read the publicly available author-uploaded journal text via primary retrieval, especially Section 3, equation (9), Section 6 and Appendix 4. The consulted body has final journal headers and pagination starting at p. 5. Its ResearchGate landing-page metadata mixes the CIFEr 2000 proceedings DOI and pp. 101–128 with the journal record, so those landing-page fields must not be copied into the bibliography. Final page extent requires publisher metadata verification before a new reference is inserted.

Equation (9) explicitly relates stock-price derivatives to log-price derivatives through signed Stirling numbers. Appendix 4 records their recurrence. Section 6 discusses limits of Taylor expansions. Thus the signed-Stirling identity in the proposed polynomial has a direct finance antecedent and cannot support a novelty claim. The proposed supremum norm, its certified global maximisation and finite basket remainder are separate comparison targets.

Consulted author body: [public author-uploaded text](https://www.researchgate.net/publication/3847827_Deriving_Derivatives_of_Derivative_Securities). No reliance on unattributed third-party copies was needed for the above algebraic claim.

## 2.4 Gaß et al. (2018)

Gaß, M., Glau, K., Mahlstedt, M. and Mair, M. (2018). Chebyshev interpolation for parametric option pricing. *Finance and Stochastics*, 22(3), 701–731. [DOI](https://doi.org/10.1007/s00780-018-0361-y). I read the final journal PDF deposited at the German National Library through primary retrieval: Section 2.1, Proposition 2.1, Section 2.3, Theorem 2.5, and Section 3, Condition 3.1 and Theorem 3.2 with its proof. I did not read this paper in Safari during the present check.

The results bound tensor Chebyshev interpolation on a parameter box using a complex analyticity domain, and account for distorted nodal values. They establish that rigorous polynomial pricing and propagation of input approximation errors are existing methods. Their principal approximation variable is the financial parameter, whereas the proposed certificate expands the arithmetic-geometric ratio within one expectation. The new norm is not an improvement over their theorem without checking its different assumptions and costs.

Consulted final body: [final journal PDF](https://d-nb.info/1165853183/34).

## 2.5 Ju (2002): full-text gap remains

Ju, N. (2002). Pricing Asian and basket options via Taylor expansion. *Journal of Computational Finance*, 5(3), 79–103. [DOI](https://doi.org/10.21314/JCF.2002.088). I verified the publisher metadata and abstract in Safari. It records first publication on 1 April 2002. The abstract describes a zero-volatility expansion of a ratio of characteristic functions. The publisher displayed a subscription barrier, with no full text or established institutional entitlement. I did not purchase access or submit a subscription enquiry. The inaccessible final body remains a source-level priority gap; it is not evidence of originality.

Consulted metadata and abstract: [publisher record](https://www.risk.net/journal-of-computational-finance/2160548/pricing-asian-and-basket-options-via-taylor-expansion).

# 3. Mathematical comparison of the proposed norm

I distinguish the operator identity, the norm statement and the pricing consequence. Carr's signed-Stirling identity and the standard Hermite representation of Gaussian derivatives produce the polynomial

\[
A_{r,\lambda}(z)=\sum_{\ell=1}^{r}s(r,\ell)\lambda^{-\ell}(-1)^{\ell-1}\operatorname{He}_{\ell-1}(z).
\]

For a finite positive threshold measure \(\nu\) of mass \(H\), Gaussian convolution gives \(t^r g_m^{(r)}(t)=\int A_{r,\lambda}(z)\phi(z)\nu(dx)\), with \(z=(m+\log t-x)/\lambda\). The immediate bound is \(H D_r\), where \(D_r=\sup_z|A_{r,\lambda}(z)\phi(z)|\). Positivity and total mass are sufficient for this step; no basket dimension occurs in the constant. A Dirac measure attains the operator bound after the location parameter is selected appropriately. That is sharpness for the stated measure class, not sharpness for a particular capped-call measure or for the expected Taylor remainder.

The comparison against the existing certificate must retain both sources of looseness. Replacing the sum of separate Gaussian derivative envelopes by the norm of their signed combination preserves cancellation. It may also sharpen each Fourier-derived envelope. An intermediate comparator

\[
\sum_{\ell=1}^{r}|s(r,\ell)|\lambda^{-\ell}\sup_z|\operatorname{He}_{\ell-1}(z)\phi(z)|
\]

would distinguish these two improvements. Comparing only the old Fourier triangle constant and the new norm does not isolate cancellation as the entire cause of a numerical saving.

The stationary polynomial is \(A'_{r,\lambda}-zA_{r,\lambda}\). Its rational coefficients permit exact real-root counting and isolation when \(\lambda\) is rational. Outward-rounded evaluation at every isolated stationary root, together with the zero limit at both infinities, certifies a global maximum. A floating-point optimiser, root plot or dense grid alone cannot provide this certificate. The mathematical facts about Hermite polynomials, polynomial root isolation and interval arithmetic are established; project-specific retained evidence is needed to show what their combination buys for basket pricing.

The finite remainder can use \(H D_{k+1}E[(T-1)^{k+1}]/(k+1)!\), since all points between 1 and \(T\) are at least 1. Gaussian tilting still evaluates the same finite coefficients. This is compatible with the current model and directly extends to positive threshold mixtures, including cash digitals and capped spreads. It does not automatically handle negative payoff measures, unbounded payoffs, negative basket weights or arbitrary stochastic volatility. The original raw-remainder divergence theorem concerns the old triangle constant and must not silently be relabelled as a theorem for \(D_r\).

# 4. Defensible next contribution test

I would present the extension as a cancellation-preserving, rigorously computed constant for the project's finite certificate, with a sharp operator statement over positive threshold mixtures. I would attribute the established differentiation and Gaussian-moment mechanisms at their use. A successful numerical test must show whether the constant changes a required degree, resolves a previously uncertified case or reduces total certified work at the same tolerance. Replacing a constant without changing any decision is a mathematical refinement with limited numerical value. Reduced work is an implementation finding until a wider matched comparison substantiates an algorithmic claim.

A compact-domain polynomial plus analytic Gaussian-tail certificate remains a viable alternative, but it has clear published basket-pricing ancestors in Olivares and Alvarez and a substantial general polynomial-error literature. The inspected sources do not establish the exact tail certificate proposed here, yet that limited observation does not establish its absence from the literature. I would not introduce a second polynomial method merely to make the manuscript longer. The sharper global norm is a cleaner immediate test because it preserves the existing coefficient engine, isolates a concrete source of conservatism and provides a reusable statement beyond one payoff.

# 5. Safari preservation

In Safari I consulted the Ju publisher page, the Caldana accepted manuscript and the Olivares publisher ePDF. A later Carr DOI navigation in Safari did not proceed, so the Carr reading above rests on primary retrieval. No account settings, subscriptions, purchases or messages were changed.
