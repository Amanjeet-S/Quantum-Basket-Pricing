# Validated constrained-hedging calculation

## 1. Financial target

I certify the funded mean-square hedging loss for the retained capped-GBM contracts, using stock and the lower-strike call under box and premium-budget constraints. The [companion paper](../paper/constrained_hedging_companion.pdf) derives the objective, feasible set and numerical certificate. Each original Python float denotes its exact stored binary64 rational value. New calculations use the corresponding exact model fair quotes. Archived losses used rounded binary64 quotes, so differences can reflect both financing prices and arithmetic.

I distinguish original holdings, their robust feasible repairs and refined holdings. Repair clips holdings to the box and scales them using a call-price upper bound and budget lower bound; it need not reduce risk. Refinement solves a conservative rational proposal problem and then certifies feasibility and optimality against the true model.

## 2. Mathematical certificate

The [engine](../src/gbm_hedge_certificate.py) encloses a three-gain Gram representation of the residual, avoiding unstable subtraction of expanded losses. An outward-evaluated feasible risk supplies the upper bound; convex Lagrange duality supplies a lower bound. A determinant/trace eigenvalue bound converts the primal-dual gap into a holding-distance bound. Exact fractions are authoritative; exhausted numerical budgets produce explicit failures.

GBM moments, convex duality and interval arithmetic are established methods. The contribution is the checked project-specific certificate and retained numerical comparison. I make no theoretical-priority claim. The auxiliary upper-strike call is an algebraic representation, not an additional permitted trading instrument.

## 3. Evidence and interpretation

The [retained reconstruction](../results/capped_gbm/validated_hedging/validation.json) checks 168 fixed portfolios, 126 refined portfolios against independent high-precision quadratic-programme references, 126 SLSQP candidates and 60 positive-residual integrals. All 21 tests pass. At 256 bits, the refined rare-tail losses remain strictly positive, approximately $1.61\times10^{-43}$ to $1.53\times10^{-41}$.

Among the 126 original constrained portfolios, 75 satisfy both box and budget constraints and 51 fail at least one. Budget-only classification gives 85 feasible and 41 infeasible. Tiny negative stock holdings explain the distinction. These are numerical feasibility findings, not evidence of economically material constraint violations.

The references evaluate the same stipulated synthetic law. They establish neither empirical hedge performance nor practical quantum advantage. Earlier outcomes were known during reconstruction; I do not describe it as an untouched prospective experiment.

## 4. Reproduction

Use Python 3.13 and the existing [price-certificate requirements](../configs/environments/requirements-price-certificate.txt). From the repository root:

```sh
python -m pip install -r configs/environments/requirements-price-certificate.txt
python -m unittest discover -s tests -p 'test_gbm_hedge_certificate.py'
python scripts/run_hedge_certificate_study.py --output /tmp/qbp-validated-hedging
```

Choose a fresh or empty external output directory. The [runner](../scripts/run_hedge_certificate_study.py) rejects internal destinations and evidence overwrites, records source hashes and fails if required comparisons or counts are missing. Independent references are diagnostics; interval validity follows from the checked mathematical construction.
