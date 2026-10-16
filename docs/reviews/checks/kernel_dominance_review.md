# Independent review of the matched-parameter kernel comparison

## 1. Scope

I check a dominance argument for the realised `full_kernel` and `interval_kernel` constructors in this repository. I compare the same cap, fractional width, Gaussian index width, exponential and logarithm prefix precisions, and probability-kernel address width. I restrict the conclusion to the emitted prefix-tree lookup, the selected interpolation arithmetic and the declared seven-T exact Toffoli assignment. I make no claim of optimal arithmetic, superiority over every quantum architecture, or end-to-end advantage.

## 2. Ordering of the actual rounded kernel ordinates

I write \(\Psi(m)=\Phi(5m)\) and

\[
P_h(m)=\frac1h\int_1^{1+h}\Phi(5[m-\log k])\,dk.
\]

Because \(\log k\geq0\), I have \(P_h(m)\leq\Psi(m)\). Both functions are nonnegative and increasing. Nearest rounding at a common fractional width is a nondecreasing map, including ties under the same rounding rule. It therefore preserves ordering of exact ordinates.

The implementation rounds high-precision centres, so I also check their generation allowances. For \(m\in[-1,1]\) and \(h\in[0.01,0.25]\), restricting the integral to \(k\in[1+h/2,1+h]\) gives

\[
\Psi(m)-P_h(m)
\geq\frac52\log(1.005)\,
\phi\{5+5\log(1.25)\}
>3.7583\,10^{-11}.
\tag{1}
\]

The last lower bound exceeds the combined kernel and scaled-CDF centre allowances \(10^{-34}+10^{-38}\). Consequently the actual centres remain ordered throughout this domain for the retained caps, and their same-width rounded integer ordinates are ordered. This step would require a new allowance or a retained endpoint check for caps outside the stated range. A claim based only on exact-function ordering would otherwise omit generator error.

## 3. Prefix-tree selector counts

Each packed table row stores its nonnegative left ordinate \(v_j\) and nonnegative difference \(v_{j+1}-v_j\). The row is zero exactly when both endpoint ordinates are zero. Endpoint ordering implies

\[
\{j:\text{scaled-CDF row }j=0\}
\subseteq
\{j:\text{capped-kernel row }j=0\}.
\tag{2}
\]

Equivalently, every nonzero capped-kernel row is a nonzero scaled-CDF row. The implemented tree visits an internal node precisely when its address subtree contains a nonzero row. It emits two exact Toffolis at each visited internal node and no Toffolis for leaf word copies. Thus the capped-kernel selector count \(T_P\) is no larger than the scaled-CDF selector count \(T_\Psi\).

I do not infer CNOT dominance. Integer ordering does not order binary popcounts, so a smaller stored word can require more CNOTs. Nor does (2) prove dominance for alternative QROM schemes, table decompositions, compiler exploitation of reachable financial addresses or physical spacetime cost.

## 4. Complete flag evaluator count

The two graphs have the same financial prefix in the retained comparisons. Both probability primitives have domain span 2, identical packed word width, identical output width and identical fractional-remainder width. Their generic interpolation multipliers and additions therefore emit the same Toffoli counts. The argument subtractions, probability subtraction, final clamp and clean comparison also have the same register widths and costs. None of their operand registers is aliased.

The full graph prepares its exact bias constant 1 using X gates. The interval graph additionally prepares a uniform-strike coordinate and evaluates the `strike_log` interpolation block. Uniform preparation and the midpoint representation emit no Toffolis, but the strike-log block has Toffoli count \(T_{\log K}\). The arithmetic body is emitted once and reversed once; the clean comparator is counted only once. Hence the exact complete-flag difference is

\[
T_{\rm interval}-T_{\rm full}
=2T_{\log K}+4(T_\Psi-T_P)\geq0.
\tag{3}
\]

This equation includes the arithmetic inverse. Multiplication by seven gives the assigned logical T-count difference under the stated decomposition. A different constant-aware interpolation or table architecture requires its own comparison.

## 5. Certificate comparison

Both probability kernels use the same global derivative bounds \(L=5/\sqrt{2\pi}\) and \(C=25/\sqrt{2\pi e}\), the same span-2 cell spacing, and the same endpoint and interpolation-product rounding allowances. Under the current conservative certificate, the full graph charges its bias constant at most \(\delta/2\). The interval strike-log output error contains at least \(3\delta/2\), even before its positive interpolation term is added. The same financial-prefix input error and subtraction allowance are present in both graphs.

Per probability call, the interval certificate therefore exceeds the full certificate by at least

\[
L\delta-(10^{-34}-10^{-38}).
\tag{4}
\]

For the searched fractional widths \(18\leq q\leq30\), this is greater than \(1.8577\,10^{-9}\). The final subtraction and projection allowances coincide, so the full arithmetic monetary certificate is no larger. This compares analytical upper bounds and does not assert an ordering of observed pointwise approximation errors.

At the same outer Gaussian width, both use the same applicable outer transfer allowance. Full conditioning has no retained strike coordinate, whereas interval evaluation adds the nonnegative strike allowance \(H/2^{b_U}\). Thus full conditioning has no larger distribution allowance. I require valid perturbed addressing domains for both mapped configurations; interval feasibility alone does not abstractly prove full-domain feasibility for an arbitrary model or prefix enclosure.

## 6. Independent emitted-code checks

I regenerated the matched full and interval graphs at all six retained interval-kernel configurations. I checked every actual table-row endpoint ordering, zero-row inclusion, selector lookup count and arithmetic certificate. I also checked the selected full configurations use the same matched parameters. The complete emitted Toffoli gaps were 6,472; 6,356; 6,380; 10,676; 7,600; and 15,740. Every gap equals (3) exactly. Both addressing certificates passed, and every full arithmetic and distribution allowance was lower.

The retained [evidence](../results/kernel_dominance_audit.json) is generated by [audit_kernel_dominance.py](audit_kernel_dominance.py), with hashes of the constructor, resource certificate, gate library and independently enclosed table record. These are actual emitted-block counts, not asymptotic gate formulas. Equation (3) explains those counts structurally; the six comparisons alone would not prove an unrestricted architecture theorem.

## 7. Interpretation

I can retain this as a checked argument about the selected constructor family: a full scalar conditional-payoff kernel supplies a no-more-expensive assigned T-count evaluator at matching precision and no larger retained error allowance. It invalidates a savings claim based solely on the weaker eight-CDF full-conditioning opponent in this instance.

I do not claim cheaper classical setup, fewer programme bits over a collection of caps, fewer CNOTs or universally lower depth. The cap-specific kernel needs its own classical generation and compiled programme, while the scaled-CDF table can be shared across caps. Controlled estimator operations, confidence rules, fault tolerance, hardware and alternative compiler families are outside this proof. The underlying probability identity and elementary monotonicity are established mathematics; the project-specific argument is the checked mapping and cost identity for the retained executable constructors. Its theoretical or literature novelty remains unestablished.

## 8. Keywords

Conditional basket pricing; scalar functions; certified endpoint rounding; coherent lookup; prefix-tree arithmetic; matched precision; Toffoli resources.
