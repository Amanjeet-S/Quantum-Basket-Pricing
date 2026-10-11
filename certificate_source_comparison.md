# 1. Scope and finding

I compare a deterministic certificate for the same retained Gaussian measure, equal-CDF cells, actual rounded nodes, conditioned basket residual and monetary accuracy. This comparison concerns financial discretisation error and certificate construction. It does not compare full pricing algorithms or imply practical quantum advantage.

The strongest new source finding is adverse to a broad novelty claim: the positive second-order cell representation in `tensor_transfer_proof.tex` is an exact specialisation of a published weighted identity. Consequently, positive transport, tensor telescoping and a transform of that transport must be treated as established mechanisms. The remaining contribution candidate is the residual-specific derivative envelope and its verified finite-moment closure, including its concrete construction cost. A matched numerical improvement can support that narrower claim, but cannot establish priority by itself.

# 2. Published weighted identity with an exact proposition mapping

Cerone, P. (2002), “Weighted Three Point Identities and their Bounds”, *SUT Journal of Mathematics*, 38(1), 17–37. I read the final [publisher PDF](https://www.rs.tus.ac.jp/sutjmath/_userdata/38-1/02-cerone.pdf), particularly Theorem 3, equations (2.13)–(2.15), and Theorem 5, equation (3.21), pages 23–30. Theorem 3 permits the endpoint limits below; its remarks explicitly discuss them. I set \(n=1,\alpha=a,\beta=b,x=c,w=\phi\). The middle derivative coefficient becomes
\[
-U_1(a,c;\phi)+L_1(c,b;\phi)=\int_a^b(y-c)\phi(y)\,dy=\mu_1.
\]
Its kernel is exactly
\[
\kappa_2(z)=
\begin{cases}
\int_a^z(z-y)\phi(y)\,dy,&z\leq c,\\
\int_z^b(y-z)\phi(y)\,dy,&z\geq c.
\end{cases}
\]
Thus the project's \(\mu_1 f'(c)+\int f''\,d\rho_2\) formula is already covered, including arbitrary nonnegative weights and interior nodes. Theorem 5 supplies its norm bound. The printed \(n=0\) kernel branch appears inconsistent; I do not rely on it. The second-order equivalence is unaffected.

# 3. The quantum-finance supplement and the grid distinction

Herman et al. (2026), [arXiv:2602.03725v1](https://arxiv.org/html/2602.03725v1), 3 February 2026. I read sections 4.2–4.3, Theorem 4.1 and Appendix A.1. The current [arXiv record](https://arxiv.org/abs/2602.03725) lists only v1 and no journal reference; exact-title searches identified no published successor as of 8 October 2026. I retain this as an unpublished supplement, not a peer-reviewed theorem.

Equation (4.13) concerns uniform physical-coordinate cells, left-endpoint values and density-weighted integration. Theorem 4.1 gives an asymptotic GBM discretisation guarantee with logarithmic bits per Gaussian under a bounded-slope piecewise-linear payoff assumption. Appendix A bounds cellwise derivatives using exponential/Gaussian relative variation. It does not state a numerical allowance for rounded equal-CDF midpoint nodes. Its unspecified asymptotic constants cannot simply become a certified monetary number.

I therefore implement a first-order weighted FTC adaptation on the project's actual grid. I do not label it a verbatim implementation of Theorem 4.1 or interpret its failure to meet an allocation as refuting that theorem. This preserves the derivative mechanism while changing the measure and quadrature rule transparently.

# 4. Independent common-grid first-order construction

For one cell \([a,b]\) and its actual node \(c\in[a,b]\), put
\[
k_c(z)=
\begin{cases}
-\{\Phi(z)-\Phi(a)\},&a\leq z<c,\\
\Phi(b)-\Phi(z),&c<z\leq b.
\end{cases}
\]
Directly applying the fundamental theorem of calculus to \(f(y)-f(c)\), and then Fubini on the compact cell, gives
\[
\int_a^b\{f(y)-f(c)\}\phi(y)\,dy
=\int_a^b k_c(z)f'(z)\,dz.
\]
The positive competitor measure \(\rho_0\) has density \(|k_c|\) summed over the disjoint cells. Its exact mass is
\[
\mathcal M_{\rho_0}(0)=\sum_c\int_a^b |y-c|\phi(y)\,dy.
\]
For a nonzero rate \(v\), one cell contributes
\[
\frac{ e^{vc}\{2\Phi(c)-\Phi(a)-\Phi(b)\}
+J_v(c,b)-J_v(a,c)}{v},
\quad
J_v(l,u)=e^{v^2/2}\{\Phi(u-v)-\Phi(l-v)\}.
\]
This follows by integration by parts on each segment; it retains the correct opposite signs of the two weighted Gaussian integrals. At zero the cell mass is
\[
2\phi(c)-\phi(a)-\phi(b)
+c\{2\Phi(c)-\Phi(a)-\Phi(b)\}.
\]
Near zero the certificate must use a rigorous limiting series or a support bound, rather than an unstable floating-point quotient.

With the unnormalised retained continuous measure \(\nu\) and rounded measure \(\tau\), both of mass \(A\), product telescoping gives the common-grid allowance
\[
\left|\int R\,d\nu^{\otimes p}-\int R\,d\tau^{\otimes p}\right|
\leq\sum_j\int|\partial_jR|\,
d\{\nu^{\otimes(j-1)}\otimes\rho_0\otimes\tau^{\otimes(p-j)}\}.
\]
The other-coordinate retained masses are already inside these measures. Their transforms factor independently. Thus a first-order competitor can use the same residual dispersion envelope, Cauchy bounds, signed exponential-polynomial coefficients and certified finite moments as the second-order construction. It need not enumerate Cartesian cells or replace \(U\) by a box maximum. I share that information in the comparison, so any resulting improvement concerns the transport order and resulting allowance, rather than an artificial information advantage.

This derivation is included to specify the implemented adaptation. I do not claim the weighted FTC mechanism, its Gaussian antiderivatives or product-measure factorisation as new theory.

# 5. The published transform-kernel antecedent

Blair, J., Luttman, A. and Machorro, E. (2016), “A generalized Peano Kernel Theorem for distributions of exponential decay”, *Journal of Mathematical Analysis and Applications*, 433(1), 622–641, [DOI 10.1016/j.jmaa.2015.07.072](https://www.sciencedirect.com/science/article/pii/S0022247X15007271). I verified the publisher abstract and metadata. It covers Laplace-transformable tempered distributions and a Fourier-domain kernel formula beyond compact support. A transform-kernel representation is therefore not a defensible broad novelty claim. Publisher retrieval returned HTTP 403; no final body was obtained. The advertised open-manuscript link could not be resolved. The earlier authorised Safari attempt was blocked by ongoing user interaction and yielded no final body. I do not infer that its full theorem lacks the project's specialisation.

The simple project formula
\[
v^2\mathcal M_{\rho_2}(v)
=\mathcal M_\nu(v)-\mathcal M_\tau(v)
-v\sum_c\mu_{1,c}e^{vc_c}
\]
also follows immediately by applying the already established cell identity to \(f(y)=e^{vy}\). It is a useful computable Gaussian specialisation, not evidence of a new transform theorem.

# 6. Executed matched benchmark and its interpretation

I execute first-order positive FTC transport and centred second-order Peano transport on the same 36 frozen second-correction payoff cases. Both receive identical rounded equal-CDF nodes, retained mass, \(L=6,N=1024,q=28\), boundary norms, first- and second-derivative dispersion envelopes \(E_1,E_2\), signed finite moments and tail allowances. The first-order row is a competent independently derived adaptation of the derivative mechanism described in section 3. It is not a verbatim implementation of Herman et al.'s equation (4.13), which uses a different grid and weighting rule.

I refine each transported-coordinate positive measure into \(B=1,8,32,128,1024\) nested blocks and apply Cauchy's inequality separately within each block. At \(B=1024\), every original cell in the transported coordinate has its own block. The other coordinates remain factorised under their continuous or rounded measures. Neither comparator enumerates a Cartesian grid. This is a cellwise refinement of the same shared derivative envelope, rather than a comparison against a deliberately expensive full tensor enumeration or direct residual-derivative sampling.

At the financial-transfer allocation \(0.005\), the FTC construction passes 32 of 36 cases and the Peano construction passes 33 of 36 cases at every tested block level. At the finest level, the Peano monetary allowances are 25.7–72.3 per cent smaller than the matched FTC allowances. This gap compares different established transport orders. Refining the second-order construction from one block to 1024 blocks reduces its own allowances by only 0.0702–0.3480 per cent and rescues no allocation failures. The larger first-order comparison must not be attributed to the block refinement, nor described as evidence that a new generic Peano theorem has been discovered.

The retained timings are diagnostic joint construction times for both transport orders and all twelve payoffs within each asset-count group. The \(B=1\) runs begin with fresh grouping caches, whereas the higher block levels reuse those caches. They exclude the common 1.276-second preparation of the original rounded nodes, prefix moments and scalar transform queries. Consequently, these timings do not support an order-specific speed claim or a comparison of cold-run cost at different block levels. The retained archive identifies the common preparation, actual scalar queries, source hashes and outward precision; it does not turn elapsed setup time into quantum gate savings or practical pricing advantage.

The remaining contribution candidate is narrower: for a positive-weight conditioned lognormal basket residual with bounded payoff-boundary derivatives, I derive explicit dispersion envelopes whose hybrid-transport integrals close in scalar transforms and monomials of degree at most four. I provide a reproducible, price-independent certificate construction with a transparent refinement parameter. I have not found that residual-specific statement in the sources actually read. This limited comparison supports its usability and the stated numerical findings, but does not settle wider literature priority. Since both admitted competitors use the same finite closure, this benchmark cannot support a claim that only the second-order method avoids Cartesian enumeration.

# 7. Access and retained evidence

The initial paper-index search failed with a DNS error. Primary web retrieval supplied the exact arXiv sections and publisher metadata. I downloaded the freely available final Cerone PDF to a temporary path and extracted and checked its relevant equations there. The source checksum and sections consulted are retained in `certificate_source_access.json`. No new subscription, purchase, paid dataset or separate document compilation was obtained. The related 2001 multivariate weighted Ostrowski record, DOI 10.1216/rmjm/1020171572, was discovered but its publisher route was blocked; I do not use it as a body-level comparison.
