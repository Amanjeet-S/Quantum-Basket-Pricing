# 1 Research objective and financial target

I develop a quantum pricing method for derivatives under an established classical financial model. Development of an existing quantum method is within scope. I retain positive-weight capped basket payoffs under constant-coefficient risk-neutral GBM as the first explicit implementation target. The quantum component is a coherent pricing oracle and estimator; it does not supply asset dynamics derived from quantum mechanics. The paper now gives the equivalent-measure, discounted-gain and admissible-strategy argument that fixes its no-arbitrage price. V denotes the undiscounted mean and Π its discounted time-zero value.

The revised question is: under what certified representation, arithmetic and estimation budgets does an exact integral-remainder oracle improve complete quantum basket pricing relative to equally informed direct correction evaluation? The improvement could concern a resource reduction for an established quantum algorithm. It need not demonstrate physical quantum advantage, but it must establish a substantive result beyond a renamed encoding or a cheap auxiliary circuit.

# 2 The construction selected for development

I subtract the finite analytic Taylor coefficients from the financial target and express its exact order-r correction using a Beta(1,r) variable and a uniformly sampled strike. The integrand is U^r t^(-r) A_r(ξ) φ(ξ), with t=1+ZU and U equal to the arithmetic/geometric ratio minus one. This representation replaces a conditional-price subtraction by a polynomial and Gaussian density. It does not establish that its reversible evaluator is cheaper. Its auxiliary randomness adds nonnegative conditional variance against exact residual evaluation, as proved in the manuscript.

On a retained Gaussian box I use the certified bound B_r=H Dbar_r U_*^r/r! and encode the signed integrand by a finite comparison flag. The proof treats saturated endpoints, rounding, original input labels and clean uncomputation. The positive first-order case has a separate flag without a signed offset. The full financial transfer includes tilted tail bounds, Beta and strike-grid bias, quantile-grid coupling, coefficient and normaliser arithmetic, value precision and statistical estimation. Approximate state preparation would require an additional certified charge. The existing financial flag's error allocation cannot be transferred without checking this new residual.

# 3 Completed mathematics and executed components

The cancellation-preserving derivative theorem is sharp over nonnegative threshold measures of fixed mass. It is not necessarily sharp for a capped-call measure or its expected remainder. Every required stationary point is enclosed by an exact rational root-counting certificate and outward interval arithmetic. Independent checks cover 27 joint profiles,9 individual Hermite profiles and 180 real stationary roots.

The price refinement checks 243 finite-degree enclosures. The same 25 of 27 prices certify, with three first-degree reductions. Executed stopping uses 793 Gaussian kernels and 760 geometry exponentials, compared with 1169 and 1552. Norm generation takes approximately 4.52 seconds of additional shared setup; hierarchy execution takes approximately 0.156 seconds. The reduction in kernel work is approximately 32.16 per cent, and a faster cold total run is not established. The frozen financial cases remain synthetic.

The minimum circuit realises the finite Beta law, retaining the original uniform labels and clearing private work. The production audit checks 10,076 dirty-output basis cases; the independent audit checks 38,828,40 emitted profiles and three complex Clifford and T replays with a reference system. At order 2 and width 8 it costs 784 T gates for the auxiliary minimum alone. The signed-threshold interface separately passes 64 endpoint cases,84 rounding checks,596 midpoint-coupling endpoints and 24 polynomial recurrences. These checks do not emit the financial value calculation.

# 4 Rejected loading candidate and the QSP boundary

The retained signed Gaussian-tilt toy is an actual small compiled quantum signal. Its signal costs 4416 T gates before amplification and 173288 after its specified coherent amplification, with ideal success approximately 0.99836. The equally informed direct finite-support signal costs 736 T gates and wins the signal-only correction-weighted comparison by approximately 18.44. The financial observable and full price estimator were not compiled in this toy. Original evidence and a regenerated portable version are both retained. Tape replay and its independent certificate do not require the optional synthesis generator.

The signed loader uses established LCU and rejection mechanisms. Its weak-dispersion asymptotics require explicit centring, a nondegenerate quadratic limit, a matched fixed base law and the active estimation regime. Exact merging of identical components can change its coefficient one-norm. A positive Gaussian-mixture competitor and coherent amplitude-normalisation requirements are derived separately; repeated measured acceptance alone does not give the clean preparation and adjoint required by estimation.

The QSP audit identifies the precise comparator singular value, the identity case, endpoint convention, projector-phase costs and a higher-dimensional obstruction to using geometric log price alone. It is not a nontrivial implemented QSP baseline. A complete QSP competitor must supply its actual financial scalar, phase sequence, admissibility and monetary-error certificate. The source record distinguishes QSP payoff encoding from QSP amplitude estimation.

# 5 The next decisive experiment

I will first emit the second-order financial value oracle, using the linear A_2 polynomial, component exponentials, log t, the strike logarithm and the Gaussian density. A certified budget must propagate every table and arithmetic error to the signed output and then to price. The finite financial support, input preparation, clean inverse, success reflection and estimator must be included. The auxiliary circuit is not a cost proxy for these operations.

I will compare it with direct evaluation of the same second-order correction and with the established positive first-order residual, giving each the same analytic coefficients, support, prior, monetary tolerance and confidence. Classical Taylor selection and continuous integration receive the same information. No quantum calls are required if the whole correction is already certified within tolerance. All 27 original synthetic cases are retained, including cases that defeat the method. A prospective finer-tolerance study must be specified before running it and must retain setup and failure costs; it cannot be selected retrospectively from favourable results.

I retain an improvement claim only after the full matched result survives this comparison. The checked derivative theorem, its pricing consequences and the oracle-error derivation are project-specific arguments; the consulted sources do not yet establish world-first theoretical priority. The full Ju article and compatible modern quantum constructions remain material comparison gaps. A top-journal contribution is an objective, not an achievement inferred from mathematical notation or manuscript length.

# 6 Data and repository state

The completed work uses synthetic parameters and public mathematics. A later REIT illustration may use existing LSEG Workspace access, with corporate actions and dividends checked; historical estimates do not establish market-consistent exotic prices. I will acquire no new subscription or data entitlement, and OptionMetrics is unavailable.

I removed the unpublished latest commit f6f74da by a soft reset to 1a25d67, which matched remote main. All 316 pre-existing research files were verified unchanged by that rollback. The changes in this phase remain uncommitted. Earlier source and result records retain their historical hashes and commit context.
