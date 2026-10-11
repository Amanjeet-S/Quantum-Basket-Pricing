# 1. Independent review of the clipped-tail moment problem

I independently checked the optimisation problem

\[
B(r;u,v,\rho)=\sup\{\mathbb E(X-r)_+:0\le X\le\rho,
\ \mathbb EX\le u,\ \mathbb EX^2\le v\}.
\]

The proposed formula is correct for nonnegative finite parameters. I assume \(r>0\) in the nontrivial regimes below. The case \(r=0\) requires a separate statement, namely \(B(0;u,v,\rho)=\min\{u,\sqrt v,\rho\}\). A constant random variable attains this value. If \(\rho\le r\), \(u=0\), or \(v=0\), the value is zero. Negative upper moment parameters describe an infeasible problem and must be rejected by the implementation.

Write \(U=\min\{u,\rho\}\) and \(V=\min\{v,\rho U,\rho^2\}\). These reductions preserve the feasible set because \(X^2\le\rho X\), \(\mathbb EX\le\rho\), and \(X^2\le\rho^2\). The last entry in the minimum defining \(V\) is redundant once \(U\le\rho\), but harmless. With \(\rho>r>0\) and \(U,V>0\), the value is

\[
B=\begin{cases}
(\rho-r)V/\rho^2,&\rho<2r,\\
V/(4r),&\rho\ge2r,\ V\le4r^2,\ U\ge V/(2r),\\
\sqrt V-r,&\rho\ge2r,\ V>4r^2,\ U^2\ge V,\\
U-rU^2/V,&\text{otherwise}.
\end{cases}
\]

# 2. Proof by feasible quadratic majorants and attaining laws

I verify each upper bound using an explicit polynomial \(Q\) satisfying \(Q(x)\ge(x-r)_+\) on \([0,\rho]\). Whenever a coefficient multiplies an upper moment bound, that coefficient is nonnegative. I then exhibit a feasible law attaining the value. This proves the result directly without requiring an unproved strong-duality assertion or interior-moment assumption.

If \(r<\rho<2r\), take \(Q(x)=(\rho-r)x^2/\rho^2\). For \(x\le r\), the comparison is immediate. For \(r\le x\le\rho\), the function \((x-r)/x^2\) has derivative \((2r-x)/x^3\ge0\), and therefore \((x-r)/x^2\le(\rho-r)/\rho^2\). Thus \(\mathbb E(X-r)_+\le(\rho-r)V/\rho^2\). The law with mass \(V/\rho^2\) at \(\rho\) and the remaining mass at zero is feasible, since \(V/\rho^2\le1\) and its mean \(V/\rho\le U\), and attains the value.

If \(\rho\ge2r\), \(V\le4r^2\), and \(U\ge V/(2r)\), take \(Q(x)=x^2/(4r)\). For \(x\ge r\), the difference is \((x-2r)^2/(4r)\ge0\), and for \(x\le r\) the polynomial is nonnegative. The bound is \(V/(4r)\). A law with mass \(V/(4r^2)\) at \(2r\) and the remaining mass at zero is feasible and attains it. Both the mass constraint and the upper mean constraint follow from the regime assumptions.

If \(V>4r^2\) and \(U^2\ge V\), set \(s=\sqrt V\). The reduction \(V\le\rho^2\) gives \(s\le\rho\), and the regime gives \(s>2r\). Take \(Q(x)=x^2/(2s)+s/2-r\). Its constant term is nonnegative. For \(x\ge r\), the difference from \(x-r\) is \((x-s)^2/(2s)\). Hence \(\mathbb E(X-r)_+\le V/(2s)+s/2-r=s-r\). The constant law \(X=s\) is feasible because \(s\le U\), and attains the bound.

In the remaining regime, I first check \(V>2rU\) and \(U^2<V\). If \(V\le4r^2\), failure of the second regime gives \(U<V/(2r)\); this implies \(U^2<V\) because \(V\le4r^2\). If \(V>4r^2\), failure of the third regime gives \(U<\sqrt V<V/(2r)\). Define \(a=V/U\) and \(p=U^2/V\). Then \(2r<a\le\rho\), where the upper bound follows from \(V\le\rho U\), and \(0<p<1\). Take

\[
Q(x)=\left(1-\frac{2rU}{V}\right)x+\frac{rU^2}{V^2}x^2.
\]

Both coefficients are nonnegative. For \(0\le x\le r\), this proves \(Q(x)\ge0\). For \(x\ge r\),

\[
Q(x)-(x-r)=r\left(1-\frac{Ux}{V}\right)^2\ge0.
\]

Taking expectations gives \(B\le U-rU^2/V\). The law with mass \(p\) at \(a\) and the remaining mass at zero has mean \(pa=U\) and second moment \(pa^2=V\), and its payoff expectation is \(p(a-r)=U-rU^2/V\). This proves attainment.

# 3. Boundaries, numerical checks and grid scope

The regimes meet continuously. At \(\rho=2r\), the endpoint formula and the second-moment formula both give \(V/(4r)\). At \(U=V/(2r)\), the final formula equals \(V/(4r)\). At \(V=4r^2\) with \(U\ge2r\), the second-moment and constant-law formulas both give \(r\). At \(U^2=V>4r^2\), the final and constant-law formulas both give \(\sqrt V-r\). At \(V=\rho U\), the final attaining law has upper support point exactly \(\rho\). The implementation should select branches by exact inequalities or use outward arithmetic at boundary comparisons.

I additionally executed an independent standard-library rational enumeration. It verified feasibility, exact attainment, nonnegative dual coefficients, the dual objective and a 17-point majorant check for 481,241 rational configurations: 157,456 zero cases, 22,138 zero-strike mean cases, 114,688 endpoint cases, 92,128 second-moment cases, and 94,831 both-active cases. Irrational constant-law cases were excluded from this enumeration and are covered by the explicit analytic square identity above. This is a diagnostic supplement to the proof, not an LP-based proof or an empirical validation claim.

The result is sharp over all probability distributions on the continuous interval \([0,\rho]\), including atomic distributions. It is a valid upper bound for a fixed-point numerical grid, but can be strictly conservative there. For example, with \(\rho=1\), \(r=1/4\), \(U=2/5\) and \(V=1/5\), the continuous-support optimum is \(1/5\), attained by mass \(4/5\) at \(1/2\). On the restricted support \(\{0,1\}\), the maximum is \(3/20\), since the second-moment constraint restricts the mass at one to \(1/5\). A finite-grid sharpness claim therefore requires a separate discrete-support argument.

# 4. Published antecedents and originality assessment

Bertsimas and Popescu (2002) already formulate sharp option-payoff bounds as moment optimisation and prove a closed-form call-payoff bound by quadratic majorisation in Section 2.2, Theorem 3, page 362. With fixed moments \(m=\mathbb EX\) and \(s=\mathbb EX^2\), the branch \(r<s/(2m)\) is exactly \(m-rm^2/s\). Their Section 2.1, Theorem 2, addresses piecewise polynomial payoffs. The paper attributes the classical call bound to Scarf (1958) and its application to option pricing to Lo (1987). I consulted the final typeset journal PDF hosted on the author's MIT site, including the introductory attribution on page 359 and Sections 2.1–2.2 on pages 360–362; I did not consult the original Scarf chapter or Lo paper. The metadata was checked against the publisher page. [Final journal PDF](https://web.mit.edu/dbertsim/www/papers/Finance/On%20the%20relation%20between%20option%20and%20stock%20prices-%20a%20convex%20optimization%20approach.pdf), [publisher record](https://pubsonline.informs.org/doi/abs/10.1287/opre.50.2.358.424).

Bertsimas and Popescu (2005), Section 2, pages 784–785, formulate polynomial majorants and explicitly state that upper moment constraints require nonnegative corresponding dual coefficients. Their Theorems 2.1–2.2 distinguish weak from strong duality and describe the interior-moment condition for the latter. I consulted the final typeset SIAM journal PDF hosted on the author's MIT site, not an earlier INSEAD working paper. [Final journal PDF](https://www.mit.edu/~dbertsim/papers/MomentProblems/Optimal-inequalities-in-probability-theory-A-convex-optimization-approach-SIAM15.pdf), [publisher record](https://epubs.siam.org/doi/abs/10.1137/S1052623401399903).

Liu and Li (2009) are an especially close bounded-support antecedent. The final publisher HTML, abstract, introduction and Section 2 text identify sharp upper call-payoff bounds under fixed first and second moments and bounded support in Proposition 2, equations (2.2)–(2.4), with quadratic majorisation and attaining two-point distributions. The mathematical formulae were stripped from the retrieved HTML, and a separately indexed final PDF could not be opened. I therefore have not verified the exact Proposition 2 equations and must not claim to have read its mathematical proof in full. This is a material outstanding comparison for any novelty claim about bounded-support refinements. [Publisher article](https://www.sciencedirect.com/science/article/abs/pii/S0167715209002120), [DOI](https://doi.org/10.1016/j.spl.2009.06.002).

The audit used free primary publication copies and publisher records. The paper-index search failed with a DNS error, so I used primary web retrieval. I made no Safari navigation and did not change publisher entitlements. No purchase or new subscription was attempted.

The derived upper-moment envelope can be presented as a checked specialisation of established moment-problem methods. The exact expression, its attainability and its use in certified quantum clipping are useful for this project. They do not establish a new general moment theory. A claim of theoretical originality would require completing the Liu–Li comparison and a wider search on bounded-support and moment-inequality stop-loss problems. Any strengthened whole-mean or clipped-tail prior should be applied to eligible classical and single-flag quantum comparators as well.

# 5. Bibliographic metadata

Bertsimas, D. and Popescu, I. (2002). 'On the relation between option and stock prices: A convex optimization approach'. Operations Research, 50(2), pp. 358–374. DOI: [10.1287/opre.50.2.358.424](https://doi.org/10.1287/opre.50.2.358.424).

Bertsimas, D. and Popescu, I. (2005). 'Optimal inequalities in probability theory: A convex optimization approach'. SIAM Journal on Optimization, 15(3), pp. 780–804. DOI: [10.1137/S1052623401399903](https://doi.org/10.1137/S1052623401399903).

Liu, G. and Li, W. V. (2009). 'Moment bounds for truncated random variables'. Statistics & Probability Letters, 79(18), pp. 1951–1956. DOI: [10.1016/j.spl.2009.06.002](https://doi.org/10.1016/j.spl.2009.06.002).

The historical attribution to Scarf, H. (1958), 'A min-max solution of an inventory problem', in Arrow, K. J., Karlin, S. and Scarf, H. (eds), Studies in the Mathematical Theory of Inventory and Production, Stanford University Press, pp. 201–209, was checked in Bertsimas and Popescu's final bibliography. I have not independently read that chapter and recommend citing the consulted Bertsimas–Popescu material for the attribution unless the original source is verified.

# 6. Retained audit and review of the project proof

The original rational enumeration is retained as `checks/audit_clipped_tail_independent.py`, with default machine-readable output `results/clipped_tail_independent_audit.json`. It imports only the Python standard library and does not import `clipping_certificate.py` or another project certificate. The exact parameter ranges, effective-moment reductions, branch decisions and checks are preserved. The machine-readable report records the skipped square-root branch counts, the exact rational branch counts and the executed source's SHA-256 hash. The 17-point quadratic checks remain diagnostics; the universal quadratic comparison is proved in Section 2 above.

I read `clipping_proof.tex` independently. I found no inconsistency in its effective upper moments, exhaustive cases, feasible attaining laws, upward square-root enclosure charge or clipping-confidence interval combination. Its use of a positive discount factor is also correct, including the express requirement to adjust the undiscounted tolerance for negative rates. The only minor mathematical wording repair is to define the ratio \(g(x)/x^2\) for \(x>0\): it vanishes on \(0<x\le r\), while \(x=0\) is handled directly. The current phrase saying the ratio vanishes for \(x\le r\) includes the undefined quotient at zero. I recommend declaring \(u,v,\rho,r\ge0\) explicitly as theorem parameters, although feasible nonnegative \(X\) already forces nonnegative moment bounds. Neither issue affects the derived formula.
