# Initial feasibility results

## 1 Interpretation

I ran twelve prespecified synthetic four-asset cases, each with eight independent batches of 4,096 paths. The figures below are sample diagnostics. They do not establish a gate-count advantage, certified variance bounds, reliable pricing confidence coverage or empirical market validity.

The comparison uses the nonnegative geometric-control residual in both the candidate and its direct opponent. If a variance-sensitive algorithm had the measured standard deviations and all other costs were equal, the candidate's evaluation overhead would need to remain below the stated ratio. Complete circuit and estimator costs are not yet available.

| Correlation | Strike | Cap | Maximum overhead versus direct control | Maximum overhead versus full conditioning with control |
| --- | --- | --- | --- | --- |
| 0.1 | 100 | 1 | 4.723 | 1.000 |
| 0.1 | 100 | 10 | 2.623 | 0.990 |
| 0.1 | 100 | 25 | 1.643 | 0.915 |
| 0.1 | 115 | 10 | 2.962 | 0.980 |
| 0.6 | 100 | 1 | 8.116 | 1.000 |
| 0.6 | 100 | 10 | 3.480 | 0.992 |
| 0.6 | 100 | 25 | 2.152 | 0.959 |
| 0.6 | 115 | 10 | 3.442 | 0.993 |
| 0.9 | 100 | 1 | 12.226 | 1.000 |
| 0.9 | 100 | 10 | 3.832 | 0.982 |
| 0.9 | 100 | 25 | 2.133 | 0.933 |
| 0.9 | 115 | 10 | 3.185 | 0.993 |

The fully conditioned controlled opponent has essentially the same standard deviation for narrow caps and lower standard deviation for wider caps. The proposed representation must therefore provide a cheaper complete oracle to beat that opponent. I do not describe the direct-opponent screen as evidence of superiority over the best quantum construction.

## 2 Adverse variance diagnostic

For the scalar conditional basket 105 exp(lambda Z), strike 100 and cap 10, I integrated the direct payoff and strike-randomised estimator separately with transition-aware positive quadrature. At lambda=0.0001, the direct payoff variance is about 0.00011025, while strike randomisation produces variance about 24.94076. At lambda=0.1 the corresponding variances are about 18.80507 and 1.15150. These are synthetic numerical findings, not observations. They rule out a universal variance-reduction claim.

## 3 Representation checks

The 13-qubit state-vector check uses a host-evaluated finite probability table. Its success probability is 0.00634765625. Across eight Grover iteration counts, the maximum discrepancy from the theoretical probability is 2.78e-15. This verifies finite-output quantum mechanics only. It does not execute the proposed reversible arithmetic or meet a financially useful error tolerance.

For sixteen selected conditioning inputs in each scenario, doubling positive quadrature order from 48 to 96 changed the result by at most 1.15e-14 basket points. Midpoint strike integration with 4,096 strikes satisfied its conservative H/(2N) diagnostic bound in all those checks. Quadrature convergence is numerical evidence, not a universal rigorous certificate.

## 4 Decision

I retain the conditional interval representation as a candidate for a bounded implementation investigation. I do not promote it to a completed contribution. The next decisive result is a stable reversible circuit with a monetary certificate and complete resource counts against direct geometric control and stable fully conditioned geometric control. Published comparison-register encoding and digital preintegration remain close antecedents.

## 5 Reproduction record

The machine-readable file retains all seeds, model parameters, statistics, checks, NumPy version and script SHA-256. The generation runtime measures this classical pilot workflow only. No provider observations or private data were used.
