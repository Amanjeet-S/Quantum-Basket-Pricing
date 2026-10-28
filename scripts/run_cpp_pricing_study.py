"""Verify the C++ capped-GBM engine against separate numerical representations.

All inputs are synthetic. Cross-language agreement is an implementation check;
neither that agreement nor high-precision quadrature validates the GBM model
against market data. This study contains no cross-language speed comparison.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad


QUANTITIES = ("price", "delta", "gamma", "vega", "rho", "theta")
PDE_MESHES = ((64, 512), (128, 1024), (256, 2048))
PROTOCOL = {
    "version": 2,
    "design_status": "informed amendment after the initial C++/Python subtraction study; not an untouched preregistration",
    "amendment": {
        "reason": "The first study exposed cancellation in the existing subtractive Python call-spread baseline. A binary64 Python/SciPy digital-strike integral was therefore added to compare against the stable mathematical representation in the same ordinary precision.",
        "changes": "add six Python binary64 quadrature quantities and their high-precision errors; retain original contracts, meshes, thresholds and subtraction comparison; add two disclosed implementation stress contracts covering the otherwise unexecuted positive-integral C++ branch and an extremely broad strike interval",
        "retained_record": "amendment_history.json contains the initial protocol, its hash and the observed outcomes motivating each amendment",
        "binary64_implementation_development": "The first stronger-baseline run retained a small-volatility Gamma cancellation failure. A subsequently inspected actual-fallback contract exposed precision loss in log(S/K); the standard log1p((S-K)/K) expression for nearby spot and strike was therefore used. Earlier results remain retained. No tolerance was relaxed.",
        "broad_interval_addition": "A numerical review then found loss of the upper-strike normal coordinate for K=1e-300, cap=spot=100, volatility=1e-8, rate=yield=0. The native endpoint evaluation was corrected and this additional informed stress case retained without changing any earlier tolerance. The direct binary64 integral is evaluated as written, with any resulting cancellation failures retained.",
        "prior_record_status": "self-contained machine-readable history supplied with this study; initial outcomes informed this amendment",
        "threshold_status": "original 2e-11 scaled absolute plus 2e-11 relative analytic tolerance unchanged for both C++ and stronger Python baseline",
    },
    "data_status": "entirely synthetic; no market observations or licensed data",
    "model": "risk-neutral GBM with continuous dividend yield; capped European call",
    "quantities": list(QUANTITIES),
    "greek_units": {
        "delta": "per unit spot", "gamma": "per unit spot squared",
        "vega": "per unit volatility", "rho": "per unit continuously compounded rate",
        "theta": "per year of calendar time; minus maturity derivative",
    },
    "comparison_roles": {
        "cpp_python_analytic": "cross-language consistency; shared formula is not an independent financial baseline",
        "python_binary64_integral": "ordinary binary64 SciPy quadrature of the same derived digital-strike representation; stronger stable Python baseline, not independent empirical evidence",
        "cpp_python_pde": "matched backward-Euler discretisation; distinct Thomas/DST implementations",
        "analytic_reference": "closed-form dominant baseline for this specified model",
        "mp_integral": "80-digit integral of digital calls over strikes, differentiated under the integral; separate numerical representation under the same GBM",
        "mp_call_spread": "100-digit call-spread check on the quadrature reference",
    },
    "pde_meshes": [{"intervals": n, "time_steps": nt} for n, nt in PDE_MESHES],
    "pde_scope": "fixed nine-sigma domain, ordinary kinked payoff, identical finite bumps; numerical refinement is not a convergence proof or a complete error certificate",
    "relative_error_rule": "report only when abs(reference)>1e-12*quantity_scale; otherwise null",
    "quantity_scale_rule": "positive dimensional scales derived from discounted cap, spot, maturity and volatility, including the narrow-spread scale",
    "analytic_cpp_tolerance": {"absolute_scaled": 2e-11, "relative": 2e-11},
    "analytic_python_integral_tolerance": {"absolute_scaled": 2e-11, "relative": 2e-11},
    "fallback_branch_price_relative_tolerance": 1e-8,
    "fallback_branch_check_scope": "additional informed price-only relative check for the new positive_integral_fallback contract; standard near-zero threshold would otherwise make the price check uninformative; no threshold changes to the original eighteen contracts",
    "failure_roles": "Native/reference/PDE failures determine implementation status. Stronger Python integral failures remain explicit baseline diagnostics. The comparison does not claim that every method passed when baseline diagnostics fail.",
    "analytic_reference_tolerance": {"absolute_scaled": 1e-25, "relative": 1e-25},
    "pde_implementation_tolerance": {"absolute_scaled": 2e-9, "relative": 2e-10},
    "pde_grid_tolerance": {"absolute_cap_scaled": 2e-10, "relative": 2e-11},
    "financial_identity_tolerance_scaled": 2e-10,
    "limitations": [
        "High-precision arithmetic and quadrature are numerical references, not outward-rounded certificates.",
        "No empirical model validation or quantum advantage follows from these comparisons.",
        "All failure rows are retained; a stress failure cannot be removed from the aggregate status.",
        "Fixed-domain, finite-bump errors are separated from analytic bump bias.",
        "These timing-free comparisons do not identify a language speedup.",
        "The stronger Python baseline was added after observing the original subtraction failures; it is an informed amendment.",
        "QUADPACK estimates are diagnostic error estimates and omit a proof of binary64 rounding error; they are not certificates.",
    ],
}


def load_engine(root: Path):
    path = root / "src" / "capped_gbm.py"
    spec = importlib.util.spec_from_file_location("cpp_study_python_reference", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scales(c) -> dict[str, float]:
    # The cap-based scale detects errors hidden by an enormous underlying spot.
    cap_value = c.cap * math.exp(-c.rate * c.maturity)
    scale_price = min(cap_value, c.spot * math.exp(-c.dividend * c.maturity))
    return {
        "price": scale_price,
        "delta": scale_price / c.spot,
        "gamma": scale_price / (c.spot * c.spot),
        "vega": scale_price / c.volatility,
        "rho": c.maturity * scale_price,
        "theta": scale_price / c.maturity,
    }


def errors(estimate: float, reference: float, scale: float) -> dict:
    error = estimate - reference
    return {
        "signed_error": error,
        "absolute_error": abs(error),
        "quantity_scale": scale,
        "scaled_absolute_error": abs(error) / scale,
        "relative_error": abs(error) / abs(reference) if abs(reference) > 1e-12 * scale else None,
        "reference_near_zero": abs(reference) <= 1e-12 * scale,
    }


def within(estimate: float, reference: float, scale: float, tolerance: dict) -> bool:
    return (math.isfinite(estimate) and math.isfinite(reference)
            and abs(estimate - reference) <= tolerance["absolute_scaled"] * scale
            + tolerance["relative"] * abs(reference))


def tolerance_ratio(estimate: float, reference: float, scale: float, tolerance: dict) -> float:
    return abs(estimate-reference)/(tolerance["absolute_scaled"]*scale
                                    + tolerance["relative"]*abs(reference))


def stress_contracts(engine):
    C = engine.Contract
    return [
        ("deep_cap_large_spot", C(spot=1e12)),
        ("deep_cap_extreme_spot", C(spot=1e16)),
        ("narrow_cap", C(cap=1e-8)),
        ("cap_below_strike_ulp", C(cap=1e-15)),
        ("wide_cap", C(cap=1e8)),
        ("deep_out_of_money", C(spot=1.0)),
        ("near_expiry_lower_kink", C(maturity=1e-8)),
        ("near_expiry_upper_kink", C(spot=120.0, maturity=1e-8)),
        ("small_positive_volatility", C(volatility=1e-4)),
        ("negative_rate_high_dividend", C(rate=-0.05, dividend=0.08)),
        ("very_high_volatility", C(volatility=2.0, maturity=5.0)),
        ("positive_integral_fallback", C(spot=99.999995, rate=0.0, dividend=0.0,
                                          volatility=1e-8)),
        ("extremely_broad_strike_interval", C(strike=1e-300, cap=100.0,
                                                rate=0.0, dividend=0.0,
                                                volatility=1e-8)),
    ]


def binary64_positive_integral(c) -> tuple[dict[str, float], dict]:
    """Evaluate the differentiated digital-strike representation in binary64.

    Price and Delta have non-negative integrands. The remaining derivatives
    may have signed integrands and cancellation. Neither the transformation
    nor the QUADPACK error estimates provide an arithmetic certificate.
    """
    S, K, H, T, r, q, sig = (c.spot, c.strike, c.cap, c.maturity,
                            c.rate, c.dividend, c.volatility)
    root_t = math.sqrt(T)
    discount = math.exp(-r*T)
    width = math.log1p(H/K)
    v = sig*root_t
    drift = r-q-0.5*sig*sig
    log_moneyness = (math.log1p((S-K)/K) if 0.5*K < S < 2.*K
                     else math.log(S/K))
    a = (log_moneyness+drift*T)/v
    cuts = [0.0, 1.0]
    for threshold in (-12., -6., 0., 6., 12.):
        t = (a-threshold)*v/width
        if 0.0 < t < 1.0:
            cuts.append(t)
    cuts = sorted(set(cuts))
    norm = math.sqrt(2.*math.pi)
    quantity_scales = scales(c)

    def integrand(t, quantity):
        d = a-width*t/v
        weight = discount*K*width*math.exp(width*t)
        density = math.exp(-0.5*d*d)/norm
        if quantity == "price":
            z = 0.5*math.erfc(-d/math.sqrt(2.))
        elif quantity == "delta":
            z = density/(S*v)
        elif quantity == "gamma":
            z = -density*(1.+d/v)/(S*S*v)
        elif quantity == "vega":
            z = density*(-d/sig-root_t)
        elif quantity == "rho":
            z = density*root_t/sig
        else:
            z = density*(drift/v-d/(2.*T))
        return weight*z

    values, estimated_errors, counts, messages = {}, {}, {}, {}
    for quantity in QUANTITIES:
        pieces, error_pieces, evaluations, diagnostics = [], [], [], []
        for left, right in zip(cuts, cuts[1:]):
            result = quad(lambda t: integrand(t, quantity), left, right,
                          epsabs=1e-13*quantity_scales[quantity]/(len(cuts)-1),
                          epsrel=1e-13, limit=300, full_output=1)
            value, estimated_error, info = result[:3]
            pieces.append(value)
            error_pieces.append(estimated_error)
            evaluations.append(info["neval"])
            diagnostics.extend(result[3:])
        values[quantity] = math.fsum(pieces)
        estimated_errors[quantity] = math.fsum(error_pieces)
        counts[quantity] = sum(evaluations)
        messages[quantity] = diagnostics
    values["rho"] -= T*values["price"]
    estimated_errors["rho"] += T*estimated_errors["price"]
    values["theta"] = r*values["price"]-values["theta"]
    estimated_errors["theta"] += abs(r)*estimated_errors["price"]
    return values, {
        "method": "binary64 log1p normalised log-strike integral; split at normal coordinates -12,-6,0,6,12; SciPy QUADPACK",
        "arithmetic": "ordinary Python float and SciPy binary64; no higher precision in this comparator",
        "estimated_quadrature_errors": estimated_errors,
        "integrand_evaluations": counts,
        "quadpack_messages": messages,
        "normalised_log_strike_cuts": cuts,
        "error_estimate_scope": "diagnostic only, including the rho/theta combination of estimated integral errors; no directed rounding or rigorous enclosure",
    }


def high_precision_reference(c) -> tuple[dict[str, float], dict[str, float], dict]:
    """Digital-strike integral, plus a higher-precision call-spread check.

    For d(x)=[log(S/x)+(r-q-sigma^2/2)T]/(sigma*sqrt(T)),
    V=exp(-rT)*integral_K^(K+H) Phi(d(x))dx. On the compact
    strike interval all five derivatives may pass under the integral for
    S,K,H,T,sigma>0. The integrands below are those exact derivatives.
    The normalised log-strike transformation avoids a short interval formed
    by subtracting two rounded binary64 strike endpoints.
    """
    with mp.workdps(80):
        S, K, H, T, r, q, sig = (mp.mpf(value) for value in
            (c.spot, c.strike, c.cap, c.maturity, c.rate, c.dividend, c.volatility))
        root_t, discount = mp.sqrt(T), mp.exp(-r * T)
        width = mp.log1p(H / K)
        v = sig * root_t
        drift = r - q - sig * sig / 2
        a = (mp.log(S / K) + drift * T) / v
        cuts = [mp.mpf(0), mp.mpf(1)]
        for threshold in (-12, -6, 0, 6, 12):
            t = (a - threshold) * v / width
            if 0 < t < 1:
                cuts.append(t)
        cuts = sorted(set(cuts))
        norm = mp.sqrt(2 * mp.pi)

        def integrand(t, quantity):
            d = a - width * t / v
            weight = K * width * mp.exp(width * t)
            density = mp.exp(-d * d / 2) / norm
            if quantity == "price":
                z = mp.erfc(-d / mp.sqrt(2)) / 2
            elif quantity == "delta":
                z = density / (S * v)
            elif quantity == "gamma":
                z = -density * (1 + d / v) / (S * S * v)
            elif quantity == "vega":
                z = density * (-d / sig - root_t)
            elif quantity == "rho":
                z = density * root_t / sig
            else:
                z = density * (drift / v - d / (2 * T))
            return weight * z

        integral = {key: discount * mp.quad(lambda t: integrand(t, key), cuts)
                    for key in QUANTITIES}
        integral["rho"] -= T * integral["price"]
        integral["theta"] = r * integral["price"] - integral["theta"]
        integral_strings = {key: mp.nstr(value, 75) for key, value in integral.items()}

    # The exact binary64 input values are reused, without a rounded K+H.
    with mp.workdps(100):
        S, K, H, T, r, q, sig = (mp.mpf(value) for value in
            (c.spot, c.strike, c.cap, c.maturity, c.rate, c.dividend, c.volatility))

        def call(strike):
            rt = mp.sqrt(T)
            d1 = (mp.log(S / strike) + (r - q + sig * sig / 2) * T) / (sig * rt)
            d2 = d1 - sig * rt
            n1 = mp.erfc(-d1 / mp.sqrt(2)) / 2
            n2 = mp.erfc(-d2 / mp.sqrt(2)) / 2
            density = mp.exp(-d1 * d1 / 2) / mp.sqrt(2 * mp.pi)
            dq, dr = mp.exp(-q * T), mp.exp(-r * T)
            return {
                "price": S * dq * n1 - strike * dr * n2,
                "delta": dq * n1,
                "gamma": dq * density / (S * sig * rt),
                "vega": S * dq * density * rt,
                "rho": strike * T * dr * n2,
                "theta": -S * dq * density * sig / (2 * rt) + q * S * dq * n1 - r * strike * dr * n2,
            }

        low, high = call(K), call(K + H)
        spread = {key: low[key] - high[key] for key in QUANTITIES}
        spread_strings = {key: mp.nstr(value, 95) for key, value in spread.items()}
        discrepancy = {key: mp.nstr(abs(mp.mpf(integral_strings[key]) - spread[key]), 25)
                       for key in QUANTITIES}
    detail = {
        "input_interpretation": "exact binary64 values supplied to both engines; upper strike summed at reference precision",
        "digital_integral_decimal_digits": 80,
        "call_spread_decimal_digits": 100,
        "digital_integral": integral_strings,
        "call_spread": spread_strings,
        "absolute_reference_discrepancy": discrepancy,
        "quadrature_status": "high precision numerical reference; no interval enclosure or directed rounding",
    }
    return ({key: float(value) for key, value in integral_strings.items()},
            {key: float(value) for key, value in spread_strings.items()}, detail)


def call_cpp(executable: Path, mode: str, c, grid=None, bumps=None, dump_grid=False) -> dict:
    args = [str(executable), mode]
    for name, value in asdict(c).items():
        flag = "yield" if name == "dividend" else name
        args.extend(("--" + flag, repr(value)))
    if grid is not None:
        for name in ("lower", "upper", "intervals", "time_steps"):
            args.extend(("--" + name.replace("_", "-"), str(getattr(grid, name))))
    if bumps is not None:
        for name, value in bumps.items():
            args.extend(("--" + name + "-bump", repr(value)))
    if dump_grid:
        args.append("--dump-grid")
    result = subprocess.run(args, check=False, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(json.dumps({"arguments": args[1:], "exit_code": result.returncode,
                                      "stderr": result.stderr, "stdout": result.stdout}))
    payload = json.loads(result.stdout)
    values = payload.get("values", payload.get("greeks", payload))
    if not all(name in values for name in QUANTITIES):
        raise ValueError("C++ output lacks the six price/Greek quantities: " + repr(payload))
    payload["values"] = {name: float(values[name]) for name in QUANTITIES}
    return payload


def run(root: Path, executable: Path, output: Path, pde: bool = True,
        amendment_history: Path | None = None) -> dict:
    root, output = root.resolve(), output.resolve()
    if output.is_relative_to(root):
        raise ValueError("Use a new external output directory; the study must not write over repository evidence.")
    executable_hash = digest(executable)
    runner_hash = digest(Path(__file__))
    python_source_hash = digest(root/"src"/"capped_gbm.py")
    amendment_history = amendment_history or root/"configs"/"cpp_pricing_amendment_history.json"
    history_hash = digest(amendment_history)
    history_bytes = amendment_history.read_bytes()
    json.loads(history_bytes)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("The output directory must be new or empty; retained evidence cannot be overwritten.")
    engine = load_engine(root)
    output.mkdir(parents=True, exist_ok=True)
    (output/"amendment_history.json").write_bytes(history_bytes)
    protocol = {**PROTOCOL,
                "contracts": [{"case": name, "group": group, "parameters": asdict(c)}
                              for group, cases in (("existing", list(engine.synthetic_contracts())),
                                                    ("stress", stress_contracts(engine)))
                              for name, c in cases]}
    # Persist this informed amendment before this run. Earlier outcomes are
    # known and retained, so this is not a claim of untouched preregistration.
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    analytic_rows, references, financial_rows, python_integral_details = [], [], [], []
    failures, baseline_failures, branch_rows = [], [], []
    cpp_floating_point = None
    for group, cases in (("existing", list(engine.synthetic_contracts())),
                         ("stress", stress_contracts(engine))):
        for name, c in cases:
            ref, ref_spread, detail = high_precision_reference(c)
            references.append({"case": name, "group": group, **detail})
            cpp = call_cpp(executable, "analytic", c)
            if cpp_floating_point is None:
                cpp_floating_point = cpp.get("floating_point", {})
            elif cpp_floating_point != cpp.get("floating_point", {}):
                raise RuntimeError("The C++ floating-point environment changed during the study.")
            py = engine.analytic(c)
            py_integral, py_integral_detail = binary64_positive_integral(c)
            python_integral_details.append({"case": name, "group": group,
                                            "values": py_integral, **py_integral_detail})
            scale = scales(c)
            for quantity in QUANTITIES:
                cpp_ok = within(cpp["values"][quantity], ref[quantity], scale[quantity],
                                PROTOCOL["analytic_cpp_tolerance"])
                python_integral_ok = within(py_integral[quantity], ref[quantity], scale[quantity],
                                            PROTOCOL["analytic_python_integral_tolerance"])
                with mp.workdps(100):
                    integral_mp = mp.mpf(detail["digital_integral"][quantity])
                    spread_mp = mp.mpf(detail["call_spread"][quantity])
                    reference_gap = abs(integral_mp-spread_mp)
                    ref_tolerance = PROTOCOL["analytic_reference_tolerance"]
                    ref_ok = bool(reference_gap <= mp.mpf(ref_tolerance["absolute_scaled"])*scale[quantity]
                                  + mp.mpf(ref_tolerance["relative"])*abs(spread_mp))
                    reference_gap_string = mp.nstr(reference_gap, 25)
                if not cpp_ok:
                    failures.append(f"analytic:{name}:{quantity}")
                if not ref_ok:
                    failures.append(f"reference:{name}:{quantity}")
                if not python_integral_ok:
                    baseline_failures.append(f"python_integral:{name}:{quantity}")
                analytic_rows.append({
                    "case": name, "group": group, "quantity": quantity,
                    "cpp": cpp["values"][quantity], "python_call_spread": py[quantity],
                    "python_positive_integral": py_integral[quantity],
                    "reference_digital_integral": ref[quantity],
                    "reference_high_precision_call_spread": ref_spread[quantity],
                    "reference_absolute_discrepancy_high_precision": reference_gap_string,
                    **errors(cpp["values"][quantity], ref[quantity], scale[quantity]),
                    "python_absolute_error": abs(py[quantity] - ref[quantity]),
                    "python_scaled_absolute_error": abs(py[quantity] - ref[quantity]) / scale[quantity],
                    "python_positive_integral_absolute_error": abs(py_integral[quantity]-ref[quantity]),
                    "python_positive_integral_scaled_absolute_error": abs(py_integral[quantity]-ref[quantity])/scale[quantity],
                    "python_positive_integral_relative_error": errors(py_integral[quantity], ref[quantity], scale[quantity])["relative_error"],
                    "python_positive_integral_pass": python_integral_ok,
                    "cpp_tolerance_ratio": tolerance_ratio(cpp["values"][quantity], ref[quantity], scale[quantity], PROTOCOL["analytic_cpp_tolerance"]),
                    "python_positive_integral_tolerance_ratio": tolerance_ratio(py_integral[quantity], ref[quantity], scale[quantity], PROTOCOL["analytic_python_integral_tolerance"]),
                    "python_positive_integral_estimated_quadrature_error": py_integral_detail["estimated_quadrature_errors"][quantity],
                    "python_positive_integral_quadpack_messages": json.dumps(py_integral_detail["quadpack_messages"][quantity]),
                    "cpp_python_absolute_difference": abs(cpp["values"][quantity] - py[quantity]),
                    "cpp_pass": cpp_ok, "reference_check_pass": ref_ok,
                    "cpp_backend": cpp.get("price_method", "not_reported"),
                    "cpp_quadrature_estimated_error": cpp.get("quadrature_estimated_error"),
                    "cpp_quadrature_evaluations": cpp.get("quadrature_evaluations"),
                    "cpp_cancellation_ratio": cpp.get("cancellation_ratio"),
                    "status": "numerical validation, not a mathematical certificate",
                })
            if name == "positive_integral_fallback":
                branch_executed = cpp.get("price_method") == "positive_digital_quadrature"
                if not branch_executed:
                    failures.append("fallback_branch_not_executed")
                for method, estimate in (("cpp", cpp["values"]["price"]),
                                         ("python_positive_integral", py_integral["price"]),
                                         ("python_original_call_spread", py["price"])):
                    relative = abs(estimate-ref["price"])/abs(ref["price"])
                    passed = relative <= PROTOCOL["fallback_branch_price_relative_tolerance"]
                    branch_rows.append({
                        "case": name, "quantity": "price", "method": method,
                        "estimate": estimate, "reference": ref["price"],
                        "relative_error_without_near_zero_suppression": relative,
                        "relative_tolerance": PROTOCOL["fallback_branch_price_relative_tolerance"],
                        "relative_check_pass": passed,
                        "cpp_price_method": cpp.get("price_method"),
                        "required_cpp_branch_executed": branch_executed,
                        "scope": "additional informed check only for this disclosed branch-coverage contract",
                    })
                    if not passed and method == "cpp":
                        failures.append("fallback_branch_cpp_price_relative")
                    if not passed and method == "python_positive_integral":
                        baseline_failures.append("fallback_branch_python_integral_price_relative")
            a = cpp["values"]
            residual = (a["theta"] + 0.5 * c.volatility**2 * c.spot**2 * a["gamma"]
                        + (c.rate-c.dividend)*c.spot*a["delta"] - c.rate*a["price"])
            cap_discount = c.cap * math.exp(-c.rate*c.maturity)
            tolerance = PROTOCOL["financial_identity_tolerance_scaled"]
            checks = {
                "non_negative_price": a["price"] >= -tolerance*scale["price"],
                "discounted_cap_upper_bound": a["price"] <= cap_discount + tolerance*scale["price"],
                "non_negative_delta": a["delta"] >= -tolerance*scale["delta"],
                "delta_upper_bound": a["delta"] <= math.exp(-c.dividend*c.maturity) + tolerance,
                "black_scholes_pde_identity": abs(residual) <= tolerance*scale["theta"] + 1e-12*abs(a["theta"]),
            }
            for check, passed in checks.items():
                financial_rows.append({"case": name, "group": group, "check": check,
                                       "pass": passed, "pde_residual": residual,
                                       "pde_residual_scale": scale["theta"]})
                if not passed:
                    failures.append(f"financial_identity:{name}:{check}")
    (output / "high_precision_references.json").write_text(json.dumps(references, indent=2) + "\n")
    (output / "python_binary64_integral_diagnostics.json").write_text(json.dumps(python_integral_details, indent=2) + "\n")
    write_csv(output / "analytic_comparison.csv", analytic_rows)
    write_csv(output / "financial_identities.csv", financial_rows)
    write_csv(output / "fallback_branch_comparison.csv", branch_rows)

    pde_rows, grid_rows, refinement_rows = [], [], []
    if pde:
        for name, c in engine.synthetic_contracts():
            lower, upper = engine.study_domain(c)
            exact, scale = engine.analytic(c), scales(c)
            for intervals, nt in PDE_MESHES:
                grid = engine.HeatGrid(lower, upper, intervals, nt)
                bumps = engine.default_bumps(c, grid)
                cpp = call_cpp(executable, "pde", c, grid, bumps)
                py = engine.bumped_greeks(c, grid, bumps)
                exact_bump = engine.analytic_bumped_greeks(c, grid, bumps)
                bounds = engine.pde_greek_domain_bounds(c, grid, bumps)
                for quantity in QUANTITIES:
                    passed = within(cpp["values"][quantity], py[quantity], scale[quantity],
                                    PROTOCOL["pde_implementation_tolerance"])
                    if not passed:
                        failures.append(f"pde_implementation:{name}:{intervals}:{quantity}")
                    pde_rows.append({
                        "case": name, "quantity": quantity, "space_intervals": intervals,
                        "time_steps": nt, "heat_lower": lower, "heat_upper": upper,
                        **{key+"_bump": val for key, val in bumps.items()},
                        "cpp": cpp["values"][quantity], "python_dst": py[quantity],
                        "analytic": exact[quantity], "analytic_bumped": exact_bump[quantity],
                        **errors(cpp["values"][quantity], exact[quantity], scale[quantity]),
                        "pde_bumped_minus_analytic_bumped": cpp["values"][quantity]-exact_bump[quantity],
                        "analytic_bump_bias": exact_bump[quantity]-exact[quantity],
                        "cpp_python_absolute_difference": abs(cpp["values"][quantity]-py[quantity]),
                        "domain_component_only_upper_bound": bounds[quantity],
                        "implementation_pass": passed,
                    })
            for nt in (1, 3, 64):
                grid = engine.HeatGrid(lower, upper, 64, nt)
                cpp = call_cpp(executable, "pde", c, grid, engine.default_bumps(c, grid), True)
                cpp_nodes = np.asarray(cpp["grid"]["nodes"], dtype=float)
                cpp_values = np.asarray(cpp["grid"]["values"], dtype=float)
                nodes, thomas = engine.heat_backward_euler_thomas(c, grid)
                _, dst = engine.heat_backward_euler(c, grid)
                cpp_thomas = float(np.max(np.abs(cpp_values-thomas)))
                cpp_dst = float(np.max(np.abs(cpp_values-dst)))
                tol = (PROTOCOL["pde_grid_tolerance"]["absolute_cap_scaled"]*c.cap
                       + PROTOCOL["pde_grid_tolerance"]["relative"]*float(np.max(np.abs(thomas))))
                passed = bool(cpp_values.shape == thomas.shape and cpp_nodes.shape == nodes.shape
                              and np.max(np.abs(cpp_nodes-nodes)) <= 1e-12*max(1., upper-lower)
                              and cpp_thomas <= tol and cpp_dst <= tol
                              and cpp_values.min() >= -tol and cpp_values.max() <= c.cap+tol
                              and cpp_values[0] == 0.0 and cpp_values[-1] == c.cap)
                grid_rows.append({"case": name, "intervals": 64, "time_steps": nt,
                                  "lambda": 0.5*c.volatility**2*(c.maturity/nt)/grid.spacing**2,
                                  "cpp_python_thomas_linf": cpp_thomas,
                                  "cpp_python_dst_linf": cpp_dst,
                                  "minimum_heat_value": float(cpp_values.min()),
                                  "maximum_heat_value": float(cpp_values.max()), "cap": c.cap,
                                  "pass": passed})
                if not passed:
                    failures.append(f"grid:{name}:{nt}")
        for intervals, nt in PDE_MESHES:
            for quantity in QUANTITIES:
                rows = [row for row in pde_rows if row["quantity"] == quantity
                        and row["space_intervals"] == intervals]
                refinement_rows.append({
                    "quantity": quantity, "space_intervals": intervals, "time_steps": nt,
                    "case_count": len(rows),
                    "maximum_absolute_total_error": max(row["absolute_error"] for row in rows),
                    "maximum_scaled_absolute_total_error": max(row["scaled_absolute_error"] for row in rows),
                    "rms_scaled_total_error": math.sqrt(sum(row["scaled_absolute_error"]**2 for row in rows)/len(rows)),
                    "maximum_absolute_pde_bump_error": max(abs(row["pde_bumped_minus_analytic_bumped"]) for row in rows),
                    "maximum_absolute_analytic_bump_bias": max(abs(row["analytic_bump_bias"]) for row in rows),
                    "status": "descriptive coupled refinement; no assumed monotonicity or fitted convergence order",
                })
        write_csv(output / "pde_comparison.csv", pde_rows)
        write_csv(output / "solver_crosscheck.csv", grid_rows)
        write_csv(output / "pde_refinement_summary.csv", refinement_rows)

    if digest(executable) != executable_hash:
        raise RuntimeError("The executable changed during the study; retain diagnostic files and rerun with a fixed build.")
    if digest(Path(__file__)) != runner_hash:
        raise RuntimeError("The study runner changed during execution; retain diagnostic files and rerun with a fixed source.")
    if digest(root/"src"/"capped_gbm.py") != python_source_hash:
        raise RuntimeError("The Python reference source changed during execution; retain diagnostic files and rerun with a fixed source.")
    if digest(amendment_history) != history_hash:
        raise RuntimeError("The amendment history changed during execution; retain diagnostic files and rerun with a fixed record.")
    manifest = {
        "status": "pass" if not failures else "fail",
        "failures": failures,
        "stronger_python_baseline_status": "pass" if not baseline_failures else "diagnostic_failures_retained",
        "baseline_diagnostic_failures": baseline_failures,
        "all_compared_stable_methods_pass": not failures and not baseline_failures,
        "analytic_output_rows": len(analytic_rows),
        "pde_output_rows": len(pde_rows), "solver_output_rows": len(grid_rows),
        "maximum_cpp_scaled_analytic_error": max(row["scaled_absolute_error"] for row in analytic_rows),
        "maximum_python_scaled_analytic_error": max(row["python_scaled_absolute_error"] for row in analytic_rows),
        "maximum_python_positive_integral_scaled_analytic_error": max(row["python_positive_integral_scaled_absolute_error"] for row in analytic_rows),
        "python_positive_integral_failed_output_rows": sum(not row["python_positive_integral_pass"] for row in analytic_rows),
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "mpmath": mp.__version__, "scipy": scipy.__version__, "platform": platform.platform(),
                        "cpp_floating_point": cpp_floating_point},
        "executable_sha256": executable_hash,
        "amendment_history_sha256": history_hash,
        "code_sha256": {"run_cpp_pricing_study.py": runner_hash,
                        "src/capped_gbm.py": python_source_hash},
        "files": {path.name: {"sha256": digest(path), "bytes": path.stat().st_size}
                  for path in sorted(output.iterdir()) if path.is_file() and path.name != "study_manifest.json"},
        "scope": "implementation verification and deterministic numerical stress study under specified GBM; no empirical validation or new pricing theory",
    }
    (output / "study_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="new or empty directory outside the project; retained evidence is never overwritten")
    parser.add_argument("--skip-pde", action="store_true")
    parser.add_argument("--amendment-history", type=Path,
                        help="defaults to configs/cpp_pricing_amendment_history.json in the project")
    args = parser.parse_args()
    output = args.output
    result = run(args.project_root.resolve(), args.executable.resolve(), output.resolve(), not args.skip_pde,
                 args.amendment_history.resolve() if args.amendment_history else None)
    print(json.dumps({key: result[key] for key in ("status", "failures", "stronger_python_baseline_status",
                                                 "baseline_diagnostic_failures", "all_compared_stable_methods_pass",
                                                 "analytic_output_rows", "pde_output_rows")}, indent=2))
    sys.exit(0 if result["status"] == "pass" else 1)
