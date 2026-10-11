# Critical structure review, 10 October 2026

## 1. Scope and evidence

This read-only review concerns the current generalisation_proofs.tex. The abstract has 234 word tokens under a conservative hyphen-preserving count, so it exceeds the requested 200-word ceiling. Seven complete proof environments occur before the appendix. Thirteen unnumbered display-math delimiters and several align* blocks remain in the main argument. The keywords contain manual line breaks and ragged-right formatting. These are objective presentation defects, independent of the numerical findings.

## 2. Essential narrative

The main paper needs a coherent sequence: the risk-neutral basket law and ordered geometric control; the smoothed conditional payoff and controlled corrections; the derivative/transport theorem explaining how finite preparation can be certified without a Cartesian price answer; the exact implemented encoding and full monetary contract; the matched estimation/resource assumptions; and the positive result together with classical competitors. Moving every formula to the appendix would damage that chain. In particular, the main text must retain the exact financial target, the relation A=GT and T≥1, the smoothed kernel, the R1 and R2 identities, the relevant theorem hypotheses and conclusions, the actual Bernoulli decoder, the aggregate monetary budget and the stronger rival comparison. The main text also needs the complete logical-T qualification and the fact that strong classical methods price all studied contracts. These are conclusions, not optional asides.

## 3. Proof relocations

The physical-to-risk-neutral exponential density and admissible-strategy no-arbitrage derivation at the start of Section 2 belongs in a financial foundations appendix. Keep the stipulated equivalent measure, martingale pricing equation and admissibility/completeness qualification in the main model. The separate common-factor covariance proof can move to the same appendix, while its existence statement and positive-definite covariance hypothesis remain in the main model.

All seven full proof environments before the appendix can move without losing the argument, provided each theorem/proposition ends with an explicit proof-appendix reference and one sentence of structural intuition. The joint norm proof should preserve the intuition that the falling log-derivative operator combines terms before taking its norm. The finite-price proof should preserve the finite-order character of the Taylor remainder. The boundary proof should preserve the fact that affine payoff terms disappear and only kink curvature remains. The convex-screen proof needs only a one-line sketch after the statement. The variance result can become a short numbered identity and one sentence, with its square expansion in the appendix. The complete-error proof is standard triangle composition, so its full derivation belongs in the appendix while the decoder-dependent factor d remains explicit. The shared-threshold proof is short but can still move, leaving the nested-event mechanism and shared auxiliary draw immediately after the proposition.

## 4. Algebra relocations

The Hermite/Stirling recurrences and Fourier-norm integration recurrence are background details suitable for the joint norm appendix. Keep definitions of A_r and D_r with the theorem, so a reader can identify exactly what is being certified. The completion-of-square steps beneath the Gaussian tilt can move while the tilt and coefficient/moment finite sums remain accessible in the main paper. The explicit closed capped-kernel formula and its boundary cancellation differentiation can move to an evaluation appendix; the main text needs a reference to the validated scalar evaluator, not the full Black-Scholes call-spread calculation.

Finite MGF factorisation is central to the absence of Cartesian enumeration. Keep its finite-sum identity or an equivalent precise statement in the main text; move the multinomial expansion and cancellation arithmetic discussion. The step-by-step signed prior interval propagation and mean-informed prior recombination can move to an information-priors appendix. Keep the final prior formulas and statement that the same analytic controls are supplied to every rival. Otherwise the main comparison could appear to receive unequal financial information.

## 5. Critical derivations to retain

The raw-payoff smoothing-first step is central to the strengthened rival, rather than routine background algebra. Keep its bound and the order of operations: charge intrinsic scalar/inner arithmetic pointwise, transfer the inner law in distribution, then integrate the continuous normal before charging outer shifts. This ordering avoids applying a continuous derivative to a finite atomic law and materially changes estimator resolution. Move the detailed clamp identity, power-series headroom proof, register-offset arithmetic and exact constants to the raw-payoff appendix.

The shared-threshold preparation, the actual ordered endpoints with a common offset, and the complete event error need to remain stated in the main paper. The midpoint and CDF coupling derivations can move. Reader understanding requires that the binary flag retains normaliser H, while inner-tail error is charged inside D_N and only the outer mass is decoded. These are design assumptions which cannot be hidden in a footnote.

## 6. Footnotes and supplementary discussions

Good short footnotes include: the exact source version consulted when it differs from the final published citation; unavailability of Ju's full article and the resulting incomplete exclusion; full-row-rank sufficiency rather than necessity; the distinction between a single model price certificate and separately rounded cross-strike quotes; future LSEG entitlement/redistribution detail; put-call parity being an applicability consequence rather than a novelty claim; and the distinction between recorded construction times and per-price execution times. Avoid multi-paragraph footnotes that merely relocate clutter to the page bottom.

The main introduction should retain one concise literature-scope paragraph. Repeated source-access logistics are better collected in a source-access appendix or ledger reference with brief footnotes at affected claims. Moving access limitations must not make the prose imply full published text was consulted. The detailed 64-bit Fourier tape diagnostic, historical 32-bit count, minor constant-simplifier small-fixture outcome and transport-block audit counts can move to implementation/replay appendices. Keep the actual estimator success and confidence conditions and resource components in the main account.

The factor-basis trial and failed bounded Chebyshev adaptation are adjacent investigations, not the main pricing-method gain. Their substantial protocols can move to a supplementary comparisons appendix, with one main paragraph recording limited/no gains and the continuing modern-estimator gap. The same is true of long setup timing and block-refinement tables if the main narrative is shortened further. By contrast, the dispersion-stress figure and stronger classical pricing results explain why transfer certification is not classical hardness and must remain prominent.

## 7. Presentation specifics

Replace manual keyword line breaks with a normal full-width paragraph, conventional spacing and no forced ragged-right. JEL classification C63/C65/G13 is appropriate to numerical/computational methods, mathematical methods and derivative pricing, subject to a final classification check. The abstract should be about 160–190 words, state the changed pricing preparation, how the error is certified, the strongest current quantitative result, and the limit that classical pricing succeeds. Excluding historical comparator counts and side investigations improves clarity.

Number displayed equations with equation or align. Use aligned inside equation for identities that form one logical object. Figure label definitions and mathematical labels must not be duplicated when moving source blocks. Keep every cross-reference attached to the same equation or result; add separate proof labels rather than changing existing theorem labels. If renumbering appendix equations, distinguish appendix labels visually and avoid automatic renumbering being represented as changed mathematical results.

## 8. Open review conditions

The current Table four-complete-prices is the preceding baseline and needs its historical role stated, because later stronger rivals have lower caps. The later shared-threshold table correctly compares the cheaper admitted rival per case. Do not let shortened narrative falsely call the preceding baseline the final lowest cost. The possible future wider precision pool would require updated table, abstract and figure cross-references together. The finite policy boundary, strongest classical conclusion and lack of established literature priority must survive editorial shortening, although repeated cautions can be consolidated.

A source-text audit cannot verify PDF pixels or page count. Native compiler success establishes compilation only. The keyword line break can be repaired structurally, and balanced source layout can be checked, but final PDF line wrapping remains subject to the native preview.

## 9. Definition repairs found during close review

The continuous moment symbol M_r is defined as E U^r in Section 4 and reused as a derivative-norm bound M_r≥||C_r||∞ later. The ordinary-call remainder already switches meaning immediately after the capped theorem. Rename the continuous moment consistently, for example as mathcal M_r, leaving M_r for payoff-boundary derivative constants. The finite moments M_{r,N} can retain their present explicit distinction.

The symbols omega and C_L first appear in the signed-prior argument without a preceding definition. For the actual retained truncated-quantile law, the score weight is one and the external mass is alpha. If the weighted uniform-box variant remains, its weight and separate monetary factor need a precise first-use definition or a supplementary-reference sentence. Its zeroth raw moment is not one. Moving its algebra must not silently erase the mass factor.

The complete monetary theorem should explicitly say that its bound holds on the jointly certified estimation event, rather than relying on the later proof to supply this qualifier. A shortened main statement otherwise looks unconditional.
