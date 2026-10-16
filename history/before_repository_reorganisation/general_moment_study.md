# Financial moment certificates, dispersion and estimation after encoding

## 1 Research question and keywords

I investigate which financial variance reductions remain useful after numerical representation and choice of estimator. The question concerns positive-asset payoffs under general pricing laws with the required moments, rather than only the original GBM grid or one reversible compiler. I derive conditional moment and implementation-error certificates, identify an encoding-dependent dispersion consequence, and test a complete finite pricing workflow. This is a candidate contribution about a transferable certificate and its computational consequences. I do not claim to have invented geometric controls, conditional Monte Carlo, transform-based financial frameworks or variance-aware quantum mean estimation.

Keywords: financial mathematics; geometric controls; conditional expectation; joint moment transforms; weak dispersion; quantum mean estimation; finite precision; computational finance.

## 2 A theorem that does not require the original model

Let positive asset values S_i have finite second moments, let positive weights sum to one, and let A and G be their weighted arithmetic and geometric averages. For an increasing bounded payoff F with Lipschitz constant c_F and upper bound H, define X_0=(F(A)−F(G))/H. Under any chosen conditioning information, X=E[X_0|𝓗] lies in [0,1]. If Δ=E[A−G|𝓗], the proved certificate is

$$
E[X^2]\leq\min\{1,c_F^2E[\Delta^2]/H^2\}
\leq\min\{1,c_F^2E[(A-G)^2]/H^2\}.
$$

This statement requires neither Gaussian returns, affine logs, independent coordinates nor a particular quantum method. Its practical use requires a computable conservative moment bound. For a finite-law approximation with pointwise error at most η, the variance is at most min(1/4,(√E[X²]+η)²), using the same underlying finite law. I prove the conditioning and numerical-transfer steps explicitly in the [proof supplement](generalisation_proofs.tex). A continuous mean-bias allowance cannot silently replace that finite-law second moment.

An affine-log-price model supplies an additional, computable route. Write ℓ_i=a_i+b_i·Y, let bars denote weighted averages, and define the joint transform Λ(v)=E[exp(v·Y)]. If P is an increasing [0,1]-valued log-kernel with global Lipschitz constant L, then X=P(logΣ_i w_i expℓ_i)−P(Σ_i w_iℓ_i) satisfies E[X²]≤min(1,L²M), where

$$
M=\sum_{i,k}w_iw_k e^{d_i+d_k}\Lambda(s_i+s_k)
-2\sum_iw_i e^{d_i}\Lambda(s_i)+1,
\qquad d_i=a_i-\bar a,\quad s_i=b_i-\bar b.
$$

The proof uses T=A/G≥1 and logT≤T−1, then evaluates E[(T−1)²]. It allows dependent latent variables through their joint transform. Required transform arguments must lie in the moment domain. Heavy-tailed models without those moments are not covered by this formula, although the positive-asset theorem can still apply when its different second-moment assumptions hold.

General financial transform controls and conditional Monte Carlo already appear in Dingeç et al. (2015). The distinction I investigate is a checked finite-output certificate, its explicit numerical-error transfer and the estimation resources that follow. The published abstract alone does not establish whether the precise inequality or dispersion theorem is absent from that paper.

## 3 A distribution-general weak-dispersion argument

I consider centred relative log-prices D_i,ε=εb_i·V+a_i(ε), where Σ_iw_i b_i=0, Σ_iw_i a_i(ε)=0 and |a_i(ε)|≤Cε². Let the geometric base m_ε converge almost surely to m_0. For an increasing continuously differentiable bounded kernel with globally bounded derivative, and under the explicitly stated polynomial-times-exponential integrability condition, I prove

$$
\frac{X_\varepsilon}{\varepsilon^2}
\longrightarrow H_0
=\tfrac12P'(m_0)V^{\mathsf T}\Sigma_bV
\quad\text{in }L^2,
\qquad\Sigma_b=\sum_iw_i b_ib_i^{\mathsf T}.
$$

Consequently the mean is ε²E[H_0]+o(ε²) and the fractional variance is ε⁴Var(H_0)+o(ε⁴). Strict Θ orders require positive leading coefficients. The proof includes the weighted cancellation, Taylor remainder, integrable square dominator and first- and second-moment limits. The latent law need not be Gaussian or independent.

For a canonical flag B=1{U<X}, with U independent and uniform, the exact identity is Var(B)=Var(X)+E[X(1−X)]. In the nondegenerate regime, single-shot flag variance has order ε² while fractional-output variance has order ε⁴. This concerns observing the flag. It is neither fundamental information loss in the entire coherent state nor a lower bound on all quantum algorithms. Montanaro (2015) already supplies a variance-aware quantum mean mechanism using a numerical value oracle; that method remains a valid further opponent.

At fixed residual-relative mean-square accuracy, the fractional sample mean can require a dispersion-independent statistical count, whereas unamplified flag averaging requires order ε⁻² draws. This does not establish constant arithmetic runtime. Maintaining the continuous limit in a finite implementation requires total L² representation error o(ε²). The existing bounded-range Bernstein formula also has a linear range term, so its relative-error sample count does not automatically share that rate. At fixed absolute price tolerance the residual eventually becomes small enough to omit, which is why I distinguish relative-residual and absolute-price experiments.

Korn and Zeytun (2013) provide a close arithmetic/geometric limiting-control antecedent. Their accessible theorem discussion uses an additive shift tending to infinity, whereas this experiment scales centred loadings towards zero. Some publisher formulae were unavailable in the retrieved rendering. I retain that distinction and the access limitation rather than assert a complete theorem-level novelty exclusion.

## 4 Financial applications and model consistency

If ψ(v)=logΛ(v), the appropriate risk-neutral drift for a rebased terminal asset with loading R_i is r−d−λ²/2−ψ(R_i), with an independent common normal shock λZ. Reusing the Gaussian correction −‖R_i‖²/2 after changing the latent law would alter the model. I prove the drift normalisation and show that the general weak-dispersion shift is Y−∇ψ(R_0), with an explicit bounded-Hessian Taylor remainder.

The mathematical class includes Gaussian mixtures, finite scenario laws and suitable non-Gaussian return models. I demonstrate a symmetric variance-gamma latent factor, an established financial model class from Madan et al. (1998). Its gamma-normal representation gives Λ(t)=(1−νt²/2)^(-1/ν) only inside its stated domain. The proof checks that domain, the required exponential integrability and the risk-neutral correction. These are stipulated model experiments, with no market calibration or empirical pricing validation.

The affine-log-price theorem also applies to positively weighted monitored GBM prices, including arithmetic Asian baskets, once their moments and conditional kernel are supplied. For a non-singular log-covariance matrix C, I prove that a common Gaussian factor can be extracted exactly when λ²≤(1ᵀC⁻¹1)^(-1). Singular covariances require a separate condition. Geometric-average controls and Asian conditioning have established antecedents in Kemna and Vorst (1990) and Curran (1994). This mathematical application does not imply that new Asian circuits have been executed.

## 5 Continuous diagnostics beyond Gaussian returns

The [four-law experiment](results/general_dispersion_results.json) retains 28 cases across seven dispersion levels. Equal-weight loadings ±0.2, λ=0.2, cap 10 and law-specific risk-neutral drifts give a common limiting derivative a=P′(−0.01). The exact mean coefficient is 0.02a and the variance coefficient is 0.0004a²Var(Y²).

| Latent law | Var(Y²) | Fractional mean / ε² at ε=1/64 | Fractional variance / ε⁴ at ε=1/64 | Limiting variance coefficient |
| --- | ---: | ---: | ---: | ---: |
| Gaussian | 2.0 | 0.03789984 | 0.00287283 | 0.00287280 |
| Standardised uniform | 0.8 | 0.03789975 | 0.00114912 | 0.00114912 |
| Symmetric variance gamma, ν=0.4 | 3.2 | 0.03789992 | 0.00459660 | 0.00459649 |
| Symmetric two-point | 0.0 | 0.03789970 | 0.00000000 | 0.00000000 |

The two-point example is retained as a degenerate boundary: the entire conditioned fractional output is constant while the canonical flag can still be random. The continuous calculations use positive density integration, refining latent quadrature 64→128, strike quadrature 48→96 and along-gap quadrature 12→24. They are numerical diagnostics, distinct from the proved theorem, the high-precision analytical moment identities and the original finite-grid confidence certificates.

## 6 Complete classical prices and charged setup

I constructed an exact vectorised integer evaluator for the same finite financial graphs. Proved limb arithmetic handles products exceeding64 bits, and exact integer reductions retain sample sums and squared sums. An independent interpreter replayed all 108 complete seeded streams, all 49,788,913 outputs, all monetary intervals and all 273 candidate-selection records. Selection uses certified sample counts times conservative scalar-lane work, with no sample outcomes used to choose a precision or stopping count. It is a finite-pool operation objective, rather than a claim of minimum CPU time.

Each price receives absolute error 0.01 and 99 per cent confidence under the declared ideal independent finite-uniform premise. The executed streams use reproducible PCG64DXSM pseudorandom numbers; their reconstruction verifies computation, not mathematical independence. The guarantee is per price, not a joint 99 per cent guarantee for all 108 outputs. Exact finite-sample means, separately certified bias, geometric correction and outward discounting supply the reported intervals. No unchecked floating finance calculation enters the pricing evaluator.

| Method | Executed finite evaluations across 27 prices | Integer evaluation seconds | Complete sampling and interval execution seconds |
| --- | ---: | ---: | ---: |
| Direct geometric residual | 18,766,096 | 1.368 | 1.959 |
| Randomised-strike interval | 17,041,317 | 1.300 | 1.854 |
| Generic full conditioning | 6,991,442 | 0.477 | 0.660 |
| Cap-aware full conditioning | 6,990,058 | 0.454 | 0.626 |

These are measurements of this NumPy implementation on the recorded host. Total sampling and interval execution is 5.099 seconds. Cold variance-certificate generation is separately recorded at 53.708 seconds, candidate table generation and selection at 57.118 seconds, including 55.049 seconds of unique-table generation, and the earlier financial-prefix certification at 2.098 seconds. The setup records do not include every research audit, code import or source-reading action. Chosen standalone table-generation costs are reconstructed from first-generation records and are labelled accordingly. None of these CPU quantities is commensurate with logical T counts.

Across matching finite prefixes, the stronger full-conditioning variance bound reduces the previous mean-only certified sample budget by factors 3.312–46.071. The relevant [certificates](results/classical_variance_certificates.json), [independent finite-law audit](results/classical_variance_audit.json), [pricing protocol](results/classical_pricing_protocol.json), [prices](results/classical_pricing_results.json) and [complete replay](results/classical_pricing_audit.json) retain the evidence. This improvement benefits fractional-output sampling; it is not assigned to the Bernoulli flag.

## 7 Stronger quantum comparison

I retain the published journal IQAE theorem, with its actual physical-angle stopping rule and a proved prior-intersection consequence, separately from the executable MIQAE-informed variant. The latter follows published confidence-allocation and greedy-power ideas, but its rational radius, outward inverses, prior bootstrap and integer chain caps have their own proof. Neither the published coefficient fifty nor sixty-two is transferred to the modified implementation. The actual MIQAE text consulted is the author version of the published proceedings chapter; the final chapter was inaccessible at no additional cost.

All 819 candidate prefixes are retained across the three architectures. An independent audit checked 13,104 complete charge profiles and 486 minima/tie sets, including preparations at zero Grover power, lookup measurements, feedback and resets. The ratio of the earlier controller's minimum guaranteed T bound to the new variant's bound ranges from 0.601 to 4.277 across case/architecture selections. Values below one mean that the new bound is worse. These are differences between proved upper budgets, not measured hardware costs or a proof of algorithmic dominance from loose bounds.

The equally supplied financial mean prior is used by both quantum controllers. The new fractional-output second moment is not assigned to their Bernoulli measurements. The [quantum comparison](results/published_comparison.json), [independent charge audit](results/published_comparison_audit.json) and [source ledger](estimator_literature_audit.md) preserve the complete record. The source IQAE envelope assumes its published real postprocessing; the finite numerical variant is verified separately.

## 8 Research contribution and the next discriminating test

The retained argument is a distribution-general financial moment and numerical-transfer certificate, a precise dispersion comparison between fractional observations and flag observations, and a complete matched implementation experiment. It is applicable beyond the original GBM cases and compiler because its hypotheses concern positive asset values, conditional information, joint transforms, kernel regularity and implementation error. The finite experiments are one checked instantiation, not the definition of its entire scope.

The literature already contains broad financial controls, transform methods, low-volatility approximations and variance-aware quantum means. I therefore retain the candidate distinction at the certificate and encoding/error interface. The available source comparison does not yet establish theoretical priority. A positive literature claim needs the precise distinction to survive full mathematical comparison with those predecessors, together with a useful complete-cost consequence. Publication at a top journal cannot be inferred from presentation or implementation sophistication.

The next discriminating implementation is a variance-aware quantum value estimator with a fully charged numerical-output oracle, alongside a stronger classical continuous conditional-pricing or quadrature routine. That comparison can determine whether the proved fractional moment improves a quantum budget once value encoding and representation costs are included. I will preserve a negative result if existing methods already supply the same argument or eliminate its computational benefit.

## 9 Reproduction and references

The new workflow uses the existing Python and NumPy environment, with no paid data or new entitlement:

```sh
python variance_certificate.py
python checks/audit_variance_certificate.py
python classical_pricing_study.py
python checks/audit_classical_pricing.py
python general_dispersion_study.py
python build_comparison.py
python checks/audit_published_comparison.py
python -m unittest discover -s checks -p 'test_*.py' -v
```

Source versions and material actually consulted are recorded in [the new literature ledger](estimator_literature_audit.md) and [SOURCES.md](SOURCES.md). References follow the established author–date convention and alphabetical order.

Curran, M. (1994). Valuing Asian and portfolio options by conditioning on the geometric mean price. *Management Science*, 40(12), 1705–1711. [DOI](https://doi.org/10.1287/mnsc.40.12.1705).

Dingeç, K. D., Sak, H. and Hörmann, W. (2015). Variance reduction for Asian options under a general model framework. *Review of Finance*, 19(2), 907–949. [DOI](https://doi.org/10.1093/rof/rfu005).

Fukuzawa, S., Ho, C., Irani, S. and Zion, J. (2023). Modified iterative quantum amplitude estimation is asymptotically optimal. *2023 Proceedings of the Symposium on Algorithm Engineering and Experiments*, 135–147. SIAM. [DOI](https://doi.org/10.1137/1.9781611977561.ch12).

Grinko, D., Gacon, J., Zoufal, C. and Woerner, S. (2021). Iterative quantum amplitude estimation. *npj Quantum Information*, 7, 52. [DOI](https://doi.org/10.1038/s41534-021-00379-1).

Kemna, A. G. Z. and Vorst, A. C. F. (1990). A pricing method for options based on average asset values. *Journal of Banking & Finance*, 14(1), 113–129. [DOI](https://doi.org/10.1016/0378-4266(90)90039-5).

Korn, R. and Zeytun, S. (2013). Efficient basket Monte Carlo option pricing via a simple analytical approximation. *Journal of Computational and Applied Mathematics*, 243(1), 48–59. [DOI](https://doi.org/10.1016/j.cam.2012.10.035).

Madan, D. B., Carr, P. P. and Chang, E. C. (1998). The variance gamma process and option pricing. *European Finance Review*, 2(1), 79–105. [DOI](https://doi.org/10.1023/A:1009703431535).

Montanaro, A. (2015). Quantum speedup of Monte Carlo methods. *Proceedings of the Royal Society A*, 471(2181), 20150301. [DOI](https://doi.org/10.1098/rspa.2015.0301).
