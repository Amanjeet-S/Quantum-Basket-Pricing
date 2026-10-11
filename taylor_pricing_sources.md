# 1. Source scope

I use the unchanged frozen synthetic-model protocol and the existing outward interval utilities in `continuous_quadrature.py`. I parse the retained decimal model values as exact rational numbers; no market or provider data enters this calculation. The installed runtime is python-flint 0.9.0, which exposes Arb, Acb and formal power-series arithmetic. I consulted the local `arb_series` documentation and available methods to implement the independent formal differentiation check. The existing project's Arb/Acb source ledger covers the previously consulted software documentation and publication metadata. I do not imply that a previously inaccessible final publication was read afresh during this calculation.

# 2. Mathematical status

Taylor's theorem, the multinomial theorem, Gaussian change of measure, Gaussian convolution and the Fourier representation of the normal density are established mathematical ingredients. The proof snippet derives their project-specific combination, every coefficient, and the finite-order remainder explicitly. I derive the project-specific certificate directly; this is not a priority claim for a new basket-pricing expansion. The normal Fourier constants are conservative uniform bounds, and I do not represent them as sharp derivative suprema.

Published basket moment/expansion and approximation results need to be compared before any literature contribution claim is made. The main paper's literature review should state the actual version and sections consulted for those comparisons. Numerical success of this certificate does not establish absence of an earlier equivalent result.
