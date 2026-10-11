# Literature audit of cap-aware coherent conditioning

## 1. Keywords and assessment

Keywords: capped basket option; analytical conditioning; call spread; quantum arithmetic; finite precision; interpolation certificate; constructive resource analysis.

I regard contract-dependent sensitivity bounds and component-dependent precision selection as established methodological ingredients. I have not established a new conditioning identity, a new general quantum approximation method or a new theory of precision allocation. A defensible next question is narrower: do explicit cap and common-loading certificates change the certified, fully emitted cost of a declared reversible evaluator, after the strongest fully conditioned scalar kernel is independently tuned under the same monetary bias budget? This question is falsifiable and permits unfavourable results. Its originality would lie in a checked financial and architectural consequence, if one survives, rather than in the sophistication of its presentation.

## 2. Closest published quantum predecessors

Chakrabarti et al., [A Threshold for Quantum Advantage in Derivative Pricing](https://doi.org/10.22331/q-2021-06-01-463), Quantum 5, 463 (2021), is the strongest direct obstacle to a broad precision-allocation claim. I consulted the final publisher PDF through Safari, particularly sections 3.3.2 and 4.1.1 and Appendix C.1–C.2. Section 3.3.2 uses second-derivative bounds to choose discretisation precision. Appendix C.2 propagates arithmetic errors using represented domains, function slopes and multiplication errors. Appendix C.1 explicitly identifies component-dependent integer and fractional precision as a possible cost reduction, while retaining uniform precisions for its study. Thus a proposal to combine financial derivatives, fixed-point error budgets and cheaper circuits substantially follows an already identified direction. Some displayed propagation steps discard higher-order terms, so I distinguish their resource methodology from a strict enclosure retaining every remainder.

Stamatopoulos and Zeng, [Derivative Pricing using Quantum Signal Processing](https://quantum-journal.org/papers/q-2024-04-30-1322/), Quantum 8, 1322 (2024), supplies the published payoff-encoding comparison. The retained record documents final publisher PDF reading in Safari, sections 4, 5.1, 6 and 7 and Appendix A. Its QSP and comparison-register constructions prevent treating coherent payoff evaluation or an explicit precision-resource trade-off as new. My exact-Toffoli arithmetic family is one constructive implementation family; its counts do not determine the best achievable QSP or quantum oracle cost.

Bhaskar, Hadfield, Papageorgiou and Petras, [Quantum Algorithms and Circuits for Scientific Computing](https://doi.org/10.26421/QIC16.3-4-2), Quantum Information and Computation 16(3–4), 197–236 (2016), provides published finite-precision scientific arithmetic. I read the final publisher PDF through retrieval, including section 2 and the logarithm construction in section 2.4. Section 2 already permits operands with different widths when their magnitudes differ; section 2.4 explicitly truncates intermediate computations and states an error bound. I did not read the detailed appendix proof of that logarithm bound. A Safari reading attempt was interrupted, so I do not record the final text as read there.

O'Brien and Sünderhauf, [Quantum state preparation via piecewise QSVT](https://doi.org/10.22331/q-2025-07-03-1786), Quantum 9, 1786 (2025), is a stronger published known-function baseline than an unrestricted appeal to arithmetic tables. I read the final publisher PDF through retrieval, sections 2.3–2.4, 3.3 and 4.1. Section 3.3 optimises dyadic segmentation for a fixed polynomial degree and error target, while section 2.4 includes precision, segments and amplification in its cost analysis. Its restricted segmentation optimum is not a global quantum circuit optimum. My lookup interpolation does not supersede this method merely because it supplies a rigorous scalar error bound.

McArdle, Gilyén and Berta, [Quantum State Preparation without Coherent Arithmetic](https://doi.org/10.1103/ntvs-c48s), Physical Review Letters 136, 240603 (2026), is the published successor of arXiv:2210.14892. I verified the final publisher metadata and abstract through retrieval; I did not access the detailed final text. Its stated arithmetic-free function-state construction is a relevant competing direction, but I do not assign its theorem or cost assumptions from the abstract alone.

## 3. Finite-representation pricing and classical conditioning

Chen, Li and Neufeld, [Quantum Monte Carlo Algorithm for Option Pricing and Its Complexity Analysis](https://doi.org/10.1111/mafi.70036), Mathematical Finance, Early View, 1–42, Version of Record 10 June 2026, is the published successor of arXiv:2301.09241. I read the final metadata and abstract in Safari. The final PDF required access that was not available in that session; I did not purchase or upgrade access. For details I read author version v6, 3 May 2026, sections 2.4.1–2.4.2 and section 5.3, including Theorem 1, Lemma 5.7 and Proposition 5.8. It separately represents the probability grid and payoff coefficients and gives explicit approximation bounds for continuous piecewise affine payoffs, subject to distribution-loading assumptions. Consequently finite monetary error certificates for quantum pricing are not a gap in the literature. I do not imply that the final journal theorem text was consulted.

Bayer, Siebenmorgen and Tempone, [Smoothing the payoff for efficient computation of basket option prices](https://doi.org/10.1080/14697688.2017.1308003), Quantitative Finance 18, 491–505 (2018), establishes the analytical basket-conditioning comparator. The retained source record documents the accessible author manuscript, including section 3, and verification of the published successor. For this project I therefore require the fully conditioned scalar kernel as a baseline. Root-based smoothing or random-strike conditioning cannot be presented as improvements by comparing only with an unnecessarily expensive unspecialised call-spread formula.

Jäckel, [Let's Be Rational](https://doi.org/10.1002/wilm.10395), Wilmott 2015(75), 40–53, addresses stable normalised Black evaluation. I verified the publisher metadata and read section 6 of the author's postpublication version dated 25 March 2016 through retrieval. I did not read the final journal text. That section explicitly treats cancellation in weighted normal-CDF differences with scaled functions and specialised expansions. I therefore treat numerical stability of Black-type formulas as a classical predecessor, rather than a new quantum-finance discovery. A coherent fixed-point implementation still needs its own absolute error certificate; this version does not itself supply the particular cap-dependent coherent certificate considered here.

## 4. Precise project question and analytical scope

I use the established conditional-payoff function

\[
P_{h,\lambda}(m)=\frac1h\int_1^{1+h}\Phi\!\left(\frac{m-\log k}{\lambda}\right)\,dk.
\]

Differentiation under this finite positive integral and integration by parts give elementary global constants

\[
L_{h,\lambda}=\min\left\{\frac{\phi(0)}\lambda,\frac{1+h}{h}\right\},\qquad
C_{h,\lambda}=\min\left\{\frac{\phi(1)}{\lambda^2},\frac{(1+h)\phi(0)}{h\lambda}\right\}.
\]

I treat these as checked elementary project derivations, not as novel theorems. The cap-dependent term improves the generic curvature constant precisely when \(\lambda<e^{-1/2}h/(1+h)\); the slope improvement requires \(\lambda<\phi(0)h/(1+h)\). Neither changes the current \(\lambda=0.2\), \(h\in\{0.01,0.1,0.25\}\) global certificate. In the prespecified extension \(\lambda\in\{0.05,0.2,0.4\}\), curvature improves at \((0.05,0.1)\) and \((0.05,0.25)\), with slope improvement at the latter. These thresholds offer a genuine analytical regime statement, but an interpolation-bound improvement is not yet a circuit-cost finding.

My proposed experiment compares cap-aware and generic bounds using identical models, domains, gate primitives and monetary budgets, with each configuration independently selected from the same declared search. I report the best retained feasible configuration, not a global optimum if the selector is greedy. I retain table generation, coefficient enclosures, compiled selectors, arithmetic, the comparison flag, exact inverse and reflections. An actual precision change must alter the emitted circuit; changing an allowance in a report is insufficient. The global certificate is independent of asset count, but realised prefix domains and prefix costs are not. A two-, four- and eight-asset family can therefore test whether a scalar interpolation saving survives complete prefix costs without asserting that its synthetic models reproduce equal empirical difficulty.

## 5. Blockers and stop rules

I stop an efficiency claim if independent retuning leaves emitted costs unchanged, if savings arise only from an inferior conditioning comparator, or if precision changes are not realised in the constructor. I reject certificates based only on sampled derivative maxima, first-order error propagation with omitted remainders, or unconstrained independently perturbed table ordinates whose claimed order is essential to the selector proof. Weak common loading requires explicit normal-tail clipping and amplified generator errors; a previous small CDF allowance cannot silently be reused outside its certified argument range.

The inequality \(P_{h,\lambda}(m)\leq\Phi(m/\lambda)\) can support an architecture-specific selector argument after a retained, error-preserving monotone generator constraint. It does not order CNOT counts, depth or every quantum lookup design. Exact seven-T Toffoli assignment is a declared constructive upper bound, not the strongest published QROM implementation. Unless QSP, unary-iteration QROM and estimator implementation are actually costed, I make no claim of globally optimal evaluation or end-to-end quantum advantage.

The closest published predecessors establish conditioning, derivative-based discretisation, mixed arithmetic precision, stable Black evaluation and approximation optimisation separately, with several expressly anticipating their combination. I have not demonstrated that the particular contract-dependent circuit consequence is absent from every predecessor. A retained strict improvement or a carefully scoped dominance counterexample could become positive project-specific knowledge. Before claiming a literature contribution, I require the checked mathematical regime, reproducible emitted evidence, and a source-level explanation of precisely which antecedent assumptions or conclusions do not contain that consequence.

## 6. Access and reproducibility record

I used my existing Safari access for academic reading. Some academic navigation attempts were interrupted, so I distinguish successful final-PDF Safari reading, retrieved final text, author versions and abstract-only access above. No unsuccessful navigation is recorded as a reading, and no new subscription, purchase or entitlement was sought. This record was prepared on 5 October 2026 from the retained project source ledger and the additional version checks described here.
