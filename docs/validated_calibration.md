# Validated parameter-to-price calculation

## 1. Financial target

I turn a scalar Gaussian-return volatility interval into an outward-enclosed capped-GBM price interval, with known spot, contractual inputs and rates. The [companion paper](../paper/capped_gbm_companion.pdf), Appendix E's validated-bridge subsection, derives the observation-rounding allowance, quantile inversion, global price extrema and probability statement. The [implementation](../src/gbm_price_certificate.py) extends the existing inference study without changing its frozen source or numerical evidence. This is a checked project-specific implementation of established probability and interval-arithmetic principles; it is not a new quantum theory.

## 2. Exact input and numerical conventions

Decimal strings and fractions denote exact rational inputs. A Python float denotes its stored binary64 value. The engine forms sample scatter exactly, propagates supplied observation-error radii through the Euclidean centring projection, and brackets chi-squared quantiles using strictly proved CDF inequalities. Arb encloses the endpoint prices and the clamped stationary maximum. Requested numerical widths must be attained; unresolved comparisons, nonpositive scatter lower bounds and exhausted budgets produce explicit abstentions. Monetary price tolerances and dimensionless quantile tolerances are separate.

The rounding helper bounds one final nearest-binary64 rounding of an exact return. It does not automatically cover returns computed from rounded prices, earlier arithmetic or provider transformations. Such observations require their own error propagation. The synthetic study deliberately rounds generated returns to six decimal places, with an exact radius of half a decimal unit in the final place. Approximate pseudorandom Gaussian generation is a diagnostic, not a certified implementation of continuous Gaussian sampling.

## 3. Reproduction

Use Python 3.13 with the [declared environment](../configs/environments/requirements-price-certificate.txt). From the organised repository root:

```sh
python -m pip install -r configs/environments/requirements-price-certificate.txt
python -m unittest discover -s tests -p 'test_gbm_price_certificate.py'
python scripts/run_gbm_certificate_study.py --output /tmp/qbp-validated-calibration
```

The output must be a new or empty directory outside the repository. The runner rejects evidence overwrites and records source hashes before and after execution. The [independent numerical representations](../src/gbm_certificate_references.py) supply high-precision call-spread and positive digital-integral diagnostics; they are not the logical basis of interval validity. The [amendment record](../configs/gbm_certificate_amendments.json) retains the corrected stationary-point refinement floor, tiny-tail integration diagnostic and separation of tolerance units. This is informed numerical development, not an untouched prospective experiment.

## 4. Evidence and interpretation

The [retained validation](../results/capped_gbm/validated_calibration/validation.json) passes 13 adverse image checks, 12 quantile checks and numerical enclosure checks for 180 synthetic datasets. All 180 calculations return certificates. The 23 tests include mathematical boundaries, abstentions, exact CDF witnesses and evidence-protection checks. Near a CDF equality, the quantile witness increases precision to 1,024 bits before retaining strict rational inequalities.

Observed true-price containment is 165/180. The SciPy full-image comparator has identical containment counts. In the near-stationary case the full image covers 54/60 targets and endpoint-only propagation covers 0/60. These small-sample counts are retained without being described as an empirical 95 per cent guarantee. Verified outward arithmetic is the improvement over the floating-point comparator; higher observed coverage is not established. The [source record](../results/capped_gbm/validated_calibration/source_record.json) distinguishes reused published foundations from the local code/API review actually performed.

If A denotes successful return, the proof establishes P(A and exclusion of the true model price) at most alpha under the stated Gaussian, constant-diffusion and rounding premises. It does not establish coverage conditional on A. The interval concerns the model price under parameter uncertainty, rather than a future payoff. It does not certify market-model adequacy, Greeks, basket calibration or practical quantum advantage. The [previous companion](../replication/history/before_validated_calibration/) remains available for historical reproduction; the quantum working paper is unchanged.
