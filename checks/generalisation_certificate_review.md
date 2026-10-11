# Independent review of the generalised conditioning certificates

## 1. Keywords, scope and conclusion

Keywords: capped basket option; coherent evaluator; ordinate enclosure; fixed-point interpolation; directed reporting; precision selection; constructive resource analysis.

I reviewed the staged `general_kernels.py`, `generalisation.py` and `resource_study.py`, the underlying Decimal normal-series implementation, the prespecified model construction and the independent audit. I did not change their source or perform Git operations. I checked the 108 retained resource rows and the corresponding independent audit, whose configuration hash matches that resource file. I found no material violation of the declared probability ordinate allowance or retained monetary bias budget. I did find an independently important logical weakness in the endpoint audit assertion and several reporting qualifications. The already retained audit metrics satisfy the stronger assertion, so the weak assertion does not presently imply a failed constant or a need to regenerate circuit counts.

## 2. Probability ordinate allowance

For the declared kernel arguments \(|m|\leq2\), loadings \(0.05\leq\lambda\leq0.4\) and relative caps \(0.01\leq h\leq0.25\), I bound the growth coefficient by \(B=e^{2.08}\). Each clipped CDF has absolute error at most \(6.4\times10^{-16}\), including saturation at arguments of magnitude at least eight. The four-CDF expression therefore has CDF contribution no greater than

\[
\frac{(2B+2+h)\,6.4\times10^{-16}}h
\leq\frac{(2e^{2.08}+2.25)\,6.4\times10^{-16}}{0.01}
<1.169\times10^{-12}.
\]

The difference from the retained \(2\times10^{-12}\) allowance exceeds \(8.31\times10^{-13}\). The seventy-digit exponential, logarithm, argument and final arithmetic errors are negligible relative to this margin, conditional on the retained normal-series enclosure on \([-8,8]\). I checked that the retained supplement uses the appropriate absolute-term bound \(8e^{32}<6.32\times10^{14}\) and a generous 4,000-operation allowance, rather than the insufficient literal 500-operation count. The normal-series implementation stops after the terms are decreasing and charges a first-omitted-term remainder.

Projection to the approximate upper CDF does not invalidate the allowance. If the raw cap centre is within \(\epsilon_K\) of the exact \(P\), the approximate upper centre is within \(\epsilon_\Psi\) of the exact \(\Psi\), and \(0\leq P\leq\Psi\), then clipping the cap centre to \([0,\widetilde\Psi]\) has error at most \(\max(\epsilon_K,\epsilon_\Psi)\), provided the upper centre is nonnegative. Here \(\epsilon_\Psi\leq\epsilon_K\). Taking the cumulative maximum also preserves the uniform enclosure of an increasing exact function. Pointwise projection before this operation and identical address grids preserve cap-centre ordering relative to the cumulative-max CDF centres. Common nearest rounding preserves that order.

The numerical guarantee is contextual. `conditional_centre` validates loading and cap but does not itself assert \(|m|\leq2\). The current table generator and function-address check impose that domain, including its table endpoints. A future caller outside that domain cannot reuse the same constant without a new proof. Adding an explicit generator range assertion would prevent accidental extrapolation, but I found no such extrapolation among the retained rows.

## 3. Endpoint audit assertion

The independent audit's `certify_table` routine currently asserts

\[
|T-c|\leq\delta/2+\epsilon_K+e,
\]

where \(c\pm e\) encloses the exact ordinate. This is only an interval compatibility test. It permits some exact ordinates in the reference interval to exceed the claimed allowance. A sufficient full-enclosure test is

\[
|T-c|+e\leq\delta/2+\epsilon_K.
\]

The routine already stores the maximum of the left side as `maximum_stored_vs_independent_enclosure`. I independently checked that recorded maximum against the claimed budget at 100 decimal digits for all 121 retained function-table records. Every record passes. The smallest margin is approximately \(1.3868\times10^{-15}\), for the general exponential at \(q=29\). The audit includes 348,537 function endpoints. Thus the retained numbers support the stronger statement, but the executable assertion and its explanation should be strengthened before describing the routine itself as independently certifying that statement. The apparent subtraction of reference uncertainty in the comment should not be retained as the justification.

This correction is logically separate from the analytical generator argument. A comparison with an uncertain reference cannot prove a uniform generator bound by merely allowing its uncertainty on the permissive side of the inequality. Conversely, the present analytic allowance remains defensible even before an independent audit is considered.

## 4. Derivative and interpolation bounds

I checked the positive-integral derivative representation and the integration-by-parts identity in the companion mathematical memorandum. They justify the implemented constants

\[
L=\min\{\phi(0)/\lambda,(1+h)/h\},\qquad
C=\min\{\phi(1)/\lambda^2,(1+h)\phi(0)/(h\lambda)\}.
\]

The curvature constant controls interpolation of the exact smooth function; the clipped numerical generator need not itself be twice differentiable because its endpoint discrepancy is separately charged. The table spacing is correctly based on the declared span, and the endpoint rounding plus fixed-point interpolation floor is charged by \(3\delta/2\). For two probability calls the error propagation charges both calls before multiplication by the monetary cap. Non-expansive clipping of the exact nonnegative residual is legitimate because the weighted arithmetic basket dominates its weighted geometric counterpart.

The cap-dependent constants improve curvature only in the prescribed weak-loading cases \((\lambda,h)=(0.05,0.1)\) and \((0.05,0.25)\), with the slope improvement at the latter. I checked all nine scalar combinations against the implementation. No reduction occurs for the original \(\lambda=0.2\) cases. This is a regime certificate, not a general improvement over the generic Gaussian bound.

## 5. Domains, arithmetic and outward reporting

The propagation checks perturbed function inputs, not only nominal financial inputs, against their half-open table address domains. It separately checks signed intermediate magnitude against 64 and uses physical clamp bounds for a clipped value followed by constant multiplication. The table endpoints include the right endpoint for interpolation while circuit input addresses remain below it. Biased quantiles, affine drifts, weights and common loading are handled consistently with the stipulated frozen-decimal models. I rebuilt all 108 selected graphs with a constant-alias diagnostic and found no collision between distinct ideal constant values. This matters because the underlying graph otherwise deduplicates constants by their represented integer.

The internal bounds use nearest-rounded Decimal arithmetic with manually conservative allowances. They are not a general directed-rounding interval library. In particular, the logarithm derivative and curvature centres are returned without separate outward increments, although the existing generator allowance and the final outward conversion have ample margin relative to seventy-digit arithmetic in these finite cases. I recommend describing this as analytical bounds with explicit conservative arithmetic allowances, rather than implying that every primitive operation is directionally rounded. If a formal local enclosure at each operation is intended, its rounding allowance should be explicit at those operations too.

The top-level arithmetic and distribution price bounds are converted with `nextafter` towards positive infinity. For conversion of a positive Decimal bound, one nextafter suffices to move past any nearest binary conversion error. I additionally checked all 108 combined output bounds against the exact Decimal sum of their two stored binary input bounds and the exact correction allowance \(10^{-12}\). Every output is outward, with smallest margin approximately \(2.13\times10^{-19}\). This retained-case check does not establish that one nextafter is sufficient for an arbitrary longer floating-point reduction. The largest retained arithmetic allowance is approximately 0.000997474, below its 0.001 budget, and the largest combined deterministic bias allowance is approximately 0.00488087, below 0.005.

Several `domain_checks` fields are nevertheless exported using ordinary `float`, including argument minima, maxima, errors and generator allowances. These are numerical summaries of the internally checked enclosures, not independently outward-rounded enclosure endpoints. For example, converting the exact Decimal \(2\times10^{-12}\) to binary float rounds down by approximately \(4.02\times10^{-29}\). This discrepancy has no material effect on the actual price budget, but the literal exported value should not be asserted to enclose the declared exact allowance. Retain Decimal strings for certificate constants and interval endpoints, or round lower endpoints down and upper/error values up. Otherwise label the exported domain fields as diagnostic summaries.

## 6. Distribution transfer and research claims

The shared outer distribution bound is appropriately applied to the continuously inner-averaged financial function for every opponent. For direct evaluation, I can first replace the outer continuous variables while the common Gaussian remains continuous, then charge the inner Gaussian grid separately. For interval evaluation, I can similarly retain the uniform strike integral during the outer replacement and charge its grid afterwards. Thus the cap-aware outer sensitivity is not unfairly reserved for the full method. The separate bounded-variation alternatives and the loose quantile coupling are consistent with this telescoping order.

The resource selector emits only two proxy-shortlisted feasible configurations per method, with the full-cap case additionally retaining the full-global winner. It does not identify the minimum over all feasible configurations or all quantum implementations. The full-cap winner cannot cost more than that included global configuration by construction. The meaningful numerical result is strict saving, not merely weak dominance. I checked three strict paired savings among the 27 full-kernel cases, all at \(\lambda=0.05\), cap 25: approximately 13.11 per cent for two assets, 3.87 per cent for four assets and 1.03 per cent for eight assets. Each reduces probability address width from 13 to 12; the four-asset winner also changes its fixed-point width. The other 24 cases must remain visible.

I would describe these as contract-dependent certificate consequences for the best retained emitted configuration within a prespecified constructor and search, at a deterministic monetary bias target. I would not call them a new conditioning method, globally optimal mixed precision, a general dimension-scaling law, empirical validation, an end-to-end statistical pricing guarantee or quantum advantage. The cap-aware and generic methods evaluate the same scalar financial function; they differ in a checked certificate and the finite implementation it permits. That is the appropriate project-specific interpretation, subject to the separate literature comparison.

## 7 Retained correction

I strengthened the endpoint audit to require the entire independent enclosure to fit inside the charged ordinate allowance, and reran it successfully. All 121 distinct function-table records satisfy this stronger condition. The ordinary-float domain fields remain labelled as diagnostic summaries in the report and proof. The portable 100-digit certificate recurrence separately checks the original 108 budgets and the 54 post-study matched-prefix certificates.
