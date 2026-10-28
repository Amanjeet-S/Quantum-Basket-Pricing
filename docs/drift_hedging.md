# Physical-drift uncertainty in static hedging

## 1. Financial question

I extend the [hedging companion](../paper/constrained_hedging_companion.pdf) from a stipulated physical drift to a compact drift interval. The objective minimises the largest discounted mean-square loss over that interval, subject to the existing position and funding limits. One known volatility is common to the physical and pricing measures. Risk-neutral premiums, instrument costs and admissible holdings therefore stay fixed as physical drift varies.


## 2. Mathematical and statistical safeguards

The [engine](../src/gbm_drift_hedging.py) evaluates exact rational portfolios with outward arithmetic. Gaussian density interpolation and Hölder's inequality bound every between-node loss by the largest endpoint loss multiplied by an explicit factor. The factor tends to one quadratically with grid spacing. Exact probability-mixture weights and non-negative funding/box multipliers provide an independently checked dual lower bound. The resulting gap certifies the candidate against the continuous uncertainty interval, rather than assuming its worst drift is an endpoint.

Independent historical Gaussian logreturns give a known-volatility drift interval, enlarged for declared observation-rounding radii, a validated chi-square quantile and a rational square-root bound. Under the stated historical/future model, the probability of both successful return and a reported upper bound below the chosen hedge's true future population mean-square loss is at most the declared significance level. This is neither conditional-on-success coverage nor a bound on each realised squared error.

The argument demonstrates conditional probability, Gaussian measures, estimation, compactness, strong convexity and primal-dual optimisation. It is a checked application of established mathematics, without a theoretical-priority claim. Unknown volatility, changing drift and misspecified returns remain outside the guarantee.

## 3. Retained evidence

The [audit](../results/capped_gbm/drift_uncertainty/validation.json) passes 12 finite-grid optimisation comparisons, 204 node intervals, 12 mixture duals and 18,504 dense baseline-risk diagnostics. All 21 new tests pass. Worst-grid risk improves in 10/12 cases, but true-drift risk worsens in those same cases. The confidence interval contains truth in 27/30 synthetic datasets. Three predetermined inference-to-decision examples have no true-risk understatement; each robust hedge nevertheless has greater true-drift risk than its plug-in comparator. Original studies and source hashes are preserved. The [alignment amendment](../configs/drift_hedge_protocol.json) makes the historical and future dividend inputs exactly identical. The [earlier execution and paper](../replication/history/drift_dividend_alignment/) remain archived; all reported comparisons were rerun after the correction.

The [protocol](../configs/drift_hedge_protocol.json) separates designed drift intervals, confidence diagnostics and the inference-to-decision examples. High-precision optimisation and dense drift evaluations check implementations; the analytic envelope establishes continuum validity. Synthetic containment counts do not demonstrate nominal confidence coverage.

## 4. Reproduction

Use the existing [certificate environment](../configs/environments/requirements-price-certificate.txt). From the repository root:

```sh
python -m unittest discover -s tests -p 'test_gbm_drift_hedging.py'
python scripts/run_drift_hedge_study.py --output /tmp/qbp-drift-hedging
```

The [runner](../scripts/run_drift_hedge_study.py) requires a fresh external destination and retains failures and source hashes. [Tests](../tests/test_gbm_drift_hedging.py) cover arithmetic, feasibility, statistical inputs and evidence protection. Unmet accuracy targets return unresolved bounds; unresolved arithmetic raises an explicit failure. Neither is silently reported as successful certification.
