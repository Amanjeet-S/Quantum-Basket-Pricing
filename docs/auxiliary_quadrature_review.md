# Independent auxiliary quadrature review

## 1. Beta register

The minimum of two independent uniform midpoint labels has exactly the continuous Beta(1,2) mass of each midpoint cell. Its cell-centred first moment is -h^3/6, rather than zero. Taylor's formula therefore gives h^2(L1/6+L2/24). The affine example f(z)=z attains h^2/6 exactly. Omitting the first-derivative term would incorrectly assign zero error to this example.

For the financial second-order remainder, the prefactor is H U^2/2. Differentiating its Beta integrand once and twice contributes U D3 and U^2 D4, since 1+zU>=1. The undiscounted bound is consequently H h^2(D3 MU3/12+D4 MU4/48). The moments refer to the ideal financial function on the actual outer grid. The error of its fixed-point implementation is charged separately.

## 2. Strike register

The uniform midpoint rule has bound h^2 times the supremum of the second derivative divided by 24. The quadratic example f(s)=s^2 attains h^2/12, which verifies this coefficient. Differentiating the financial strike argument gives a first derivative of magnitude at most H/(K lambda) and a second derivative of magnitude at most H^2/(K^2 lambda). The resulting monetary bound is H^3 MU2 h^2(F2/lambda^2+E2/lambda)/(48 K^2). The Hermite envelopes for E2 and F2 are valid triangle bounds on the first two derivatives of A2 phi.

The Beta and strike discretisation errors can be compared sequentially. Their uniform derivative bounds hold for either continuous or discretised values of the other auxiliary variable, so a triangle inequality combines them. No independence between statistical estimators enters this deterministic argument.

## 3. Support and likelihood weights

On the finite outer support, 0<=U<=U_star implies U^4<=U_star U^3 pointwise. This proves MU4<=U_star MU3. The conclusion remains valid when the moments contain a fixed nonnegative Gaussian likelihood weight. An external retained Gaussian mass or uniform-box volume factor remains outside the helper and must be applied by the monetary reconstruction. Neither that factor nor the oracle arithmetic allowance can be omitted.

## 4. Retained verification

The independent audit checks eight exact minimum-label laws, including non-dyadic diagnostic cell counts, all cell masses and both centred moments. It checks 64 signed-polynomial quadrature comparisons, 108 financial norm and scaling cases with ordinary or likelihood-weighted moments, and nine invalid-input rejections. It reconstructs the selected minimum of the older linear coupling cap and the new quadratic cap using exact rational arithmetic and retained norm certificates.

The proof and helper are correct under their stated range, smoothness, moment and normalisation hypotheses. The quadrature mechanism is established mathematics. The project-specific application supplies a stronger error certificate for the existing emitted auxiliary laws, rather than a new generic midpoint theorem or a changed quantum circuit.
