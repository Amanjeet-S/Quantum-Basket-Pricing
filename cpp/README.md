# Native capped-GBM pricing engine

## 1. Build

I implement the established capped-GBM valuation and five Greeks in C++17, with no third-party C++ dependency. The reusable interface is [capped_gbm.hpp](include/qbp/capped_gbm.hpp). The [mathematical and numerical note](../docs/cpp_pricing.md) gives the derivations, discrete residual proof and comparison limits.

From the repository root, build into a fresh directory outside the repository:

```sh
mkdir /tmp/qbp-native-build
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Icpp/include \
  cpp/src/capped_gbm.cpp cpp/src/main.cpp -o /tmp/qbp-native-build/qbp_gbm
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -Icpp/include \
  cpp/src/capped_gbm.cpp cpp/tests/native_tests.cpp -o /tmp/qbp-native-build/qbp_native_tests
/tmp/qbp-native-build/qbp_native_tests
```

Alternatively, use `cmake -S cpp -B /tmp/qbp-native-build -DCMAKE_BUILD_TYPE=Release`, followed by `cmake --build /tmp/qbp-native-build` and `ctest --test-dir /tmp/qbp-native-build --output-on-failure`. The retained run verified direct compiler builds; CMake/CTest were unavailable.

## 2. Interface and numerical conventions

```sh
/tmp/qbp-native-build/qbp_gbm analytic --spot 110 --strike 100 --cap 20 \
  --maturity 1 --rate 0.03 --yield 0.01 --volatility 0.2
/tmp/qbp-native-build/qbp_gbm pde --spot 110 --lower -2 --upper 2 \
  --intervals 256 --time-steps 2048
```

Each successful call emits JSON with price, delta, gamma, vega, rho and theta. Vega/rho use unit parameter changes; calendar theta is minus the remaining-maturity derivative. Yield is continuous. Negative gamma/vega are possible for the capped payoff. `floating_point` records actual type precision; `long double` has only 53 mantissa bits on the tested arm64 system.

Analytic mode uses tail-aware normal intervals and density differences, retaining the cap width even when `K+H` rounds to `K`. Severe price cancellation invokes a disclosed positive digital-strike quadrature. `price_method`, `cancellation_ratio`, `quadrature_estimated_error` and `quadrature_evaluations` report the path taken. The error estimate is not a rigorous enclosure. Underflow remains a limitation; an unresolved fallback transformation is explicitly rejected. Prices are not clipped.

The analytic interface supports zero cap, strike, maturity or volatility with explicit limiting conventions. Undefined derivatives at deterministic payoff kinks are `null`; vega at zero volatility is a right derivative. PDE mode requires strictly positive parameters, both payoff kinks strictly inside its fixed heat domain, and all bumped queries inside it. Invalid or unresolved grids, bumps and non-finite arithmetic intermediates are rejected.

PDE Greeks are central differences on the same fixed domain/grid. Supply all four optional sizes together using `--spot-bump`, `--volatility-bump`, `--rate-bump` and `--maturity-bump`; otherwise the companion's shrinking defaults apply. `analytic_same_bumps` separates analytic bump bias from the PDE discrepancy. `domain_price_bound` is only a domain component, not a total price/Greek certificate.

`--dump-grid` emits all final nodes and values. `--dump-trace` additionally emits initial and subsequent grids, limited to one million entries. The [exact auditor](../scripts/certify_cpp_heat_trace.py) certifies the exported-coefficient discrete target. Its bound excludes continuum and original-input errors. The [reproduction guide](../docs/reproducing_paper.md) supplies the comparison/audit commands and retained evidence links. There is no language-speedup, new-model or quantum-advantage claim.
