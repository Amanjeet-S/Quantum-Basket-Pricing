#include "qbp/capped_gbm.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string escape_json(const std::string& text) {
    std::string out;
    for (char c : text) {
        if (c == '"' || c == '\\') { out += '\\'; out += c; }
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else out += c;
    }
    return out;
}

void number(qbp::Real value) {
    if (!std::isfinite(value)) throw std::runtime_error("A non-finite number cannot be emitted as JSON.");
    std::cout << value;
}

void optional_number(const std::optional<qbp::Real>& value) {
    if (value) number(*value);
    else std::cout << "null";
}

void greeks(const qbp::Greeks& g) {
    std::cout << "{\"price\":"; number(g.price);
    std::cout << ",\"delta\":"; optional_number(g.delta);
    std::cout << ",\"gamma\":"; optional_number(g.gamma);
    std::cout << ",\"vega\":"; optional_number(g.vega);
    std::cout << ",\"rho\":"; optional_number(g.rho);
    std::cout << ",\"theta\":"; optional_number(g.theta);
    std::cout << '}';
}

void array(const std::vector<qbp::Real>& values) {
    std::cout << '[';
    for (std::size_t i = 0; i < values.size(); ++i) { if (i) std::cout << ','; number(values[i]); }
    std::cout << ']';
}

qbp::Real real_value(const std::map<std::string, std::string>& options, const std::string& name, qbp::Real fallback) {
    auto found = options.find(name);
    if (found == options.end()) return fallback;
    std::size_t used = 0;
    qbp::Real value;
    try { value = std::stold(found->second, &used); }
    catch (...) { throw std::invalid_argument("Invalid numeric value for --" + name + '.'); }
    if (used != found->second.size() || !std::isfinite(value))
        throw std::invalid_argument("--" + name + " must be a finite number.");
    return value;
}

std::size_t count_value(const std::map<std::string, std::string>& options, const std::string& name, std::size_t fallback) {
    const auto found = options.find(name);
    if (found == options.end()) return fallback;
    if (found->second.empty()) throw std::invalid_argument("Empty integer for --" + name + '.');
    for (char c : found->second) if (c < '0' || c > '9')
        throw std::invalid_argument("--" + name + " requires an unsigned integer.");
    std::size_t used = 0;
    const auto value = std::stoull(found->second, &used);
    if (used != found->second.size() || value > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Integer overflow for --" + name + '.');
    return static_cast<std::size_t>(value);
}

void usage() {
    std::cout << "Usage: qbp_gbm analytic|pde [options]\n"
                 "Contract: --spot S --strike K --cap H --maturity T --rate r --yield q --volatility sigma\n"
                 "PDE: --lower L --upper U --intervals M --time-steps N [--dump-grid] [--dump-trace]\n"
                 "Bumps: --spot-bump dS --volatility-bump ds --rate-bump dr --maturity-bump dT\n"
                 "All four bump options must be supplied together. PDE heat boundaries remain fixed under bumps.\n"
                 "Output: one JSON record. Calendar theta; vega/rho per unit change. No external dependencies.\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
            usage(); return argc < 2 ? 2 : 0;
        }
        const std::string mode = argv[1];
        if (mode != "analytic" && mode != "pde") throw std::invalid_argument("The first argument must be analytic or pde.");
        std::map<std::string, std::string> options;
        bool dump_grid = false, dump_trace = false;
        const std::vector<std::string> allowed = {"spot", "strike", "cap", "maturity", "rate", "yield", "volatility",
            "lower", "upper", "intervals", "time-steps", "spot-bump", "volatility-bump", "rate-bump", "maturity-bump"};
        for (int i = 2; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--dump-grid") { dump_grid = true; continue; }
            if (option == "--dump-trace") { dump_trace = dump_grid = true; continue; }
            if (option == "--help" || option == "-h") { usage(); return 0; }
            if (option.rfind("--", 0) != 0) throw std::invalid_argument("Unexpected argument: " + option);
            const std::string name = option.substr(2);
            if (std::find(allowed.begin(), allowed.end(), name) == allowed.end())
                throw std::invalid_argument("Unknown option: " + option);
            if (i + 1 >= argc) throw std::invalid_argument("Missing value for " + option);
            if (!options.emplace(name, argv[++i]).second) throw std::invalid_argument("Duplicate option: " + option);
        }
        if (mode == "analytic") {
            if (dump_grid || dump_trace) throw std::invalid_argument("Grid/trace output requires pde mode.");
            for (const std::string& option : {"lower", "upper", "intervals", "time-steps", "spot-bump", "volatility-bump", "rate-bump", "maturity-bump"})
                if (options.count(option)) throw std::invalid_argument("--" + option + " requires pde mode.");
        }
        qbp::Contract c;
        c.spot = real_value(options, "spot", c.spot);
        c.strike = real_value(options, "strike", c.strike);
        c.cap = real_value(options, "cap", c.cap);
        c.maturity = real_value(options, "maturity", c.maturity);
        c.rate = real_value(options, "rate", c.rate);
        c.dividend = real_value(options, "yield", c.dividend);
        c.volatility = real_value(options, "volatility", c.volatility);
        c.validate();
        std::cout << std::setprecision(std::numeric_limits<qbp::Real>::max_digits10);
        const auto start = std::chrono::steady_clock::now();
        if (mode == "analytic") {
            const auto result = qbp::analytic(c);
            const auto elapsed = std::chrono::duration<qbp::Real>(std::chrono::steady_clock::now() - start).count();
            std::cout << "{\"mode\":\"analytic\",\"values\":"; greeks(result.values);
            std::cout << ",\"price_method\":\"" << result.price_method << "\",\"derivative_convention\":\"" << escape_json(result.derivative_convention) << '"';
            std::cout << ",\"cancellation_ratio\":"; number(result.cancellation_ratio);
            std::cout << ",\"quadrature_estimated_error\":"; number(result.quadrature_estimated_error);
            std::cout << ",\"quadrature_evaluations\":" << result.quadrature_evaluations;
            std::cout << ",\"elapsed_seconds\":"; number(elapsed);
        } else {
            qbp::HeatGrid grid;
            grid.lower = real_value(options, "lower", grid.lower);
            grid.upper = real_value(options, "upper", grid.upper);
            grid.intervals = count_value(options, "intervals", grid.intervals);
            grid.time_steps = count_value(options, "time-steps", grid.time_steps);
            grid.validate(c);
            if (dump_trace && grid.time_steps >= 1000000 / (grid.intervals + 1))
                throw std::invalid_argument("Trace capture is limited to one million grid entries.");
            unsigned bump_count = 0;
            for (const std::string& name : {"spot-bump", "volatility-bump", "rate-bump", "maturity-bump"}) bump_count += options.count(name) != 0;
            if (bump_count != 0 && bump_count != 4) throw std::invalid_argument("Supply all four bump sizes or none.");
            std::optional<qbp::Bumps> bumps;
            if (bump_count == 4) bumps = qbp::Bumps{real_value(options, "spot-bump", 0), real_value(options, "volatility-bump", 0),
                                                  real_value(options, "rate-bump", 0), real_value(options, "maturity-bump", 0)};
            const auto result = qbp::pde(c, grid, bumps);
            const auto elapsed = std::chrono::duration<qbp::Real>(std::chrono::steady_clock::now() - start).count();
            std::cout << "{\"mode\":\"pde\",\"values\":"; greeks(result.values);
            std::cout << ",\"analytic_same_bumps\":"; greeks(result.analytic_same_bumps);
            std::cout << ",\"bumps\":{\"spot\":"; number(result.bumps.spot);
            std::cout << ",\"volatility\":"; number(result.bumps.volatility);
            std::cout << ",\"rate\":"; number(result.bumps.rate);
            std::cout << ",\"maturity\":"; number(result.bumps.maturity); std::cout << '}';
            std::cout << ",\"domain_price_bound\":"; number(result.base_domain_bound);
            std::cout << ",\"heat_min\":"; number(result.heat_min);
            std::cout << ",\"heat_max\":"; number(result.heat_max);
            std::cout << ",\"grid\":{\"lower\":"; number(grid.lower);
            std::cout << ",\"upper\":"; number(grid.upper);
            std::cout << ",\"intervals\":" << grid.intervals << ",\"time_steps\":" << grid.time_steps;
            std::cout << ",\"spacing\":"; number(grid.spacing());
            std::cout << ",\"lambda\":"; number(result.lambda);
            std::optional<qbp::HeatSolution> solution;
            if (dump_grid) {
                solution = qbp::heat_backward_euler(c, grid, dump_trace);
                std::cout << ",\"nodes\":"; array(solution->nodes);
                std::cout << ",\"values\":"; array(solution->values);
            }
            std::cout << '}';
            if (dump_trace) {
                std::cout << ",\"heat_trace\":[";
                for (std::size_t j = 0; j < solution->trace.size(); ++j) { if (j) std::cout << ','; array(solution->trace[j]); }
                std::cout << ']';
            }
            std::cout << ",\"elapsed_seconds\":"; number(elapsed);
            std::cout << ",\"derivative_convention\":\"central bumped prices on a fixed heat domain; calendar theta\"";
        }
        std::cout << ",\"floating_point\":{\"long_double_digits\":" << std::numeric_limits<qbp::Real>::digits
                  << ",\"double_digits\":" << std::numeric_limits<double>::digits
                  << ",\"long_double_bytes\":" << sizeof(qbp::Real)
                  << "},\"numerical_error_certificate\":false}" << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "{\"error\":\"" << escape_json(error.what()) << "\"}" << '\n';
        return 2;
    }
}
