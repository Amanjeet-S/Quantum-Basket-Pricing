# Independent review of the deterministic comparators

## 1. Finding

I found no material mathematical error in the retained deterministic pricing proof or a discrepancy in the numerical evidence. Exact rational algebra and file provenance pass for all 47 finite-grid price records and nine continuous price records. The inherited 273-prefix tensor inventory, prespecified enumeration scope and selected-prefix minimisation rule agree with the input protocol. The retained full Python-integer replay audit covers all 330,752 enumerated outputs. I checked its hashes and counts rather than repeating that historical computation.

The maximum finite-price radius is approximately 0.004679193, below the reserved 0.005 deterministic allowance. Every finite-grid interval contains the corresponding narrower continuous interval. Across all 47 matched prefixes, the enumerated-output to certified-Monte-Carlo sample ratio is between 4096/1053823 and 256/441. This compares financial evaluator invocations with a particular conservative sample budget, not execution time or the best possible classical estimator.

## 2. Continuous enclosure and derivation

I checked the positive strike representation, integration-by-parts expression, weighted arithmetic–geometric residual sign, geometric Gaussian correction and one-sided tail proof. Their use is mathematically consistent for the stipulated positive weights, Gaussian common factor and exact rational risk-neutral parameters. The retained code verifies these parameter identities.

The primary residual prices use 128 bits and L=8; the direct conditional-reference prices use 176 bits and L=10. The primary continuous interval widths are approximately 1.2442×10^−15 to 3.1313×10^−14, computed from the retained exact rational endpoints. All nine representations overlap. The primary calculation has 2,112 integrand calls. Its roughly 0.185-second combined primary and reference timing is a host diagnostic that excludes process start and imports.

The callback forwards analytic=True into the complex logarithm. Recorded non-finite evaluations trigger the integrator's neighbourhood handling. A regular real path would not justify removing this guard. No midpoint extraction, requested-tolerance substitution or output clipping constructs the certificate. Exporting exact outward endpoints avoids reliance on rounded diagnostic ball strings. The 126 positive strike-integral kernel checks are separate representations but share the Arb library with the outer price calculation.

## 3. Correction and limitations

I corrected the source ledger's misidentified Johansson journal reference. The Arb paper was published in IEEE Transactions on Computers 66(8), 1281–1292, with DOI 10.1109/TC.2017.2690633. The review verifies primary-source bibliographic metadata and actual python-flint 0.9.0 API contracts. It does not claim that the inaccessible final journal text was read.

These are constructive competitors for nine synthetic two-asset prices. They falsify an argument that improving upon this Monte Carlo budget alone demonstrates quantum pricing advantage. They do not establish an empirical pricing result, general superiority, a high-dimensional infeasibility result, a classical complexity lower bound, a new quadrature theorem or a physical quantum runtime advantage. Confidence one remains conditional on correct validated-library enclosures, exact evaluator execution and the stipulated model and bias certificates.

## 4. Retained evidence

The machine-readable deterministic_evidence_review.json contains the exact endpoint and ratio fractions, each checked configuration, protocol hashes, source hashes and result hashes. The read-only verification source is retained alongside it. I made no Git changes and did not edit the main project repository.

## 5. Keywords

Validated quadrature; deterministic pricing; numerical certification; basket options; comparator fairness; source provenance.

## 6. Subsequent access update

After the independent review, I accessed the author-hosted typeset Arb copy through Safari and consulted the sections recorded in continuous_quadrature_sources.md. The historical review JSON records the proof source before its bibliographic acknowledgement was added. Its numerical evidence hashes remain unchanged.
