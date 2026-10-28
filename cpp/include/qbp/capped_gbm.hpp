#ifndef QBP_CAPPED_GBM_HPP
#define QBP_CAPPED_GBM_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace qbp {

using Real = long double;

struct Contract {
    Real spot = 100;
    Real strike = 100;
    Real cap = 20;
    Real maturity = 1;
    Real rate = 0.03L;
    Real dividend = 0.01L;
    Real volatility = 0.2L;
    void validate() const;
};

// Calendar theta is -dV/d(maturity); vega/rho use unit parameter changes.
// A missing value denotes a derivative undefined at a degenerate payoff kink.
struct Greeks {
    Real price = 0;
    std::optional<Real> delta;
    std::optional<Real> gamma;
    std::optional<Real> vega;
    std::optional<Real> rho;
    std::optional<Real> theta;
};

struct AnalyticResult {
    Greeks values;
    std::string price_method;
    std::string derivative_convention;
    Real cancellation_ratio = 1;
    Real quadrature_estimated_error = 0;
    std::size_t quadrature_evaluations = 0;
};

struct HeatGrid {
    Real lower = -3;
    Real upper = 3;
    std::size_t intervals = 512;
    std::size_t time_steps = 2048;
    void validate(const Contract& c) const;
    Real spacing() const;
};

struct HeatSolution {
    std::vector<Real> nodes;
    std::vector<Real> values;
    Real lambda = 0;
    std::vector<std::vector<Real>> trace;
};

struct Bumps {
    Real spot = 0;
    Real volatility = 0;
    Real rate = 0;
    Real maturity = 0;
    void validate(const Contract& c) const;
};

struct PDEResult {
    Greeks values;
    Bumps bumps;
    Real base_domain_bound = 0;
    Real lambda = 0;
    Real heat_min = 0;
    Real heat_max = 0;
    // These are central differences of the exact price at identical bumps.
    // They expose bump bias separately from the PDE discrepancy.
    Greeks analytic_same_bumps;
};

Real normal_cdf(Real x);
Real normal_pdf(Real x);
Real payoff(Real terminal_spot, const Contract& c);
AnalyticResult analytic(const Contract& c);
Bumps default_bumps(const Contract& c, const HeatGrid& grid);
HeatSolution heat_backward_euler(const Contract& c, const HeatGrid& grid,
                                 bool capture_trace = false);
Real pde_price(const Contract& c, const HeatGrid& grid);
Real domain_price_bound(const Contract& c, const HeatGrid& grid);
PDEResult pde(const Contract& c, const HeatGrid& grid,
              const std::optional<Bumps>& bumps = std::nullopt);

} // namespace qbp

#endif
