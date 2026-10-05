# Certified quantum encoding of bounded basket derivatives

## 1 Research question and status

I investigate how certified scalar-function specialisation changes the complete coherent cost of direct sampling, analytical conditioning and randomised-strike conditioning of a geometrically controlled basket call spread. The proposed interval advantage has failed the stronger-comparator test. The probability identities, geometric control, Gaussian conditioning and interpolation are established ingredients. An original positive literature contribution and complete quantum resource savings remain unverified.

I replaced the earlier direction and precision selector with one concrete evaluator. For positive normalised weights, arithmetic basket A and geometric basket G satisfy A >= G. I write the capped payoff as F(x) = min(max(x - K, 0), H), retain the analytical correction E[F(G)], and estimate the nonnegative residual E[F(A) - F(G)]. A uniformly sampled strike turns that residual into a conditional Gaussian interval probability. Under the established common-loading factor, both interval boundaries are logarithms and no iterative root solver is required.

The initial financial object is a terminal equity basket call spread under stipulated correlated GBM. REITs remain a possible later LSEG parameter scenario. Asian monitoring is an extension only after the core evaluator survives its benchmarks. Neither asset substitution nor increased nominal dimension is a contribution.

## 2 Completed feasibility work

The [stronger conditioning benchmark](conditioning_benchmark.md) reports five complete constructed evaluators for three caps at two monetary targets. Specialised full conditioning has the lowest assigned T count in all six cases. At cap £10 and target £0.01, it uses 597,478 T gates against 672,210 for the comparably specialised interval evaluator and 1,207,388 for direct evaluation. These are constructor-specific logical costs with half the tolerance reserved for unimplemented statistical estimation. The [current resource record](results/conditioning_resources.json) and [independent audit](results/conditioning_independent_audit.json) preserve the evidence. The [initial study](reversible_encoding_study.md) and its results remain as the superseded comparison.

The [research design](research_design.md) sets the comparison and decision rules. The [mathematical note](mathematical_note.md) checks the identities and monetary bias certificates. The [source record](SOURCES.md) distinguishes published versions from the versions actually consulted. The [pilot report](results/pilot_report.md) describes twelve fixed synthetic scenarios and an adverse variance example. The [machine-readable results](results/pilot.json) retain inputs, seeds, diagnostics and a script hash.

The small state-vector check verifies finite-output Born probabilities and Grover interference. Its host-evaluated probability truth table is explicitly recorded. It does not implement scalable reversible arithmetic, demonstrate a financial accuracy target or establish a speedup. The pilot standard deviations and batch standard errors are screening diagnostics, not certified variance bounds or pricing confidence intervals.

## 3 Reproduction

Use Python 3.11 or later with NumPy. In a virtual environment, install [requirements.txt](requirements.txt), then run:

```sh
python pilot.py
python audit_certificates.py
python -m unittest discover -s checks -p "test_*.py" -v
python resource_study.py --certify-quantiles --small-checks --output results/conditioning_resources.json
python checks/reference.py
python reflections.py
python checks/audit_composition.py
python checks/shadow_audit.py --config results/conditioning_resources.json --output results/conditioning_independent_audit.json --samples 262144 --certify-quantiles --new-variant-gates
python checks/audit_kernel_dominance.py
```

The new resource study uses only the Python standard library. Its independent reference and replay audit require NumPy. It emits complete primitive-block gates, checks function domains and price-bias budgets, and independently encloses every retained quantile and function constant. The two-asset resource instance is distinct from the earlier four-asset pilot.

The default pilot run uses eight independent batches of 4,096 paths for each of twelve scenarios. It writes results/pilot.json. Runtime is recorded separately from the numerical outputs and is not a quantum cost estimate. The exact script and NumPy version are retained in the output.

## 4 Next decision

I have constructed and checked stronger conditioned arithmetic. The measured interval/full gap reverses, so I withdraw the interval improvement claim. I will next assess whether a checked precision-allocation or scalar-kernel argument generalises across a prespecified financial family and supplies an original, informative contribution against equally specialised opponents. Selecting the cheapest established method is insufficient. A large mean-estimation implementation should follow a surviving contribution argument, with a valid confidence and monetary-error guarantee. Observed pilot variance cannot supply its stopping certificate.

I will stop the claimed improvement if its savings depend on tuning withheld from competitors, if the finite-output monetary certificate cannot be met at useful precision, or if stronger constructed comparators eliminate the measured benefit. This repository is an initial feasibility study, not a completed paper or an assured positive contribution.

## 5 Data and access

All current experiments use explicitly stipulated synthetic parameters. No LSEG observations, credentials, purchased data or publisher PDFs are included. Existing LSEG access may later support market-informed scenarios after a field and adjustment audit. Raw and reconstructible provider observations remain local and excluded from Git, including when the remote repository is private.

## 6 Keywords

Quantum computational finance; basket derivatives; geometric control variates; conditional expectation; reversible arithmetic; monetary error certification.
