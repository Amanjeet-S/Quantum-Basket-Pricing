"""Exact residual bounds for published binary64 backward-Euler heat traces.

The target uses the published coefficient, boundary and initial values.
It is not a certificate for the whole-line financial price or input transforms.
"""

from __future__ import annotations

from fractions import Fraction
import math
from numbers import Real
from typing import Sequence


def binary64(value: Real) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError("A finite real binary64 value is required.")
    try:
        rounded = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError("A finite real binary64 value is required.") from error
    if not math.isfinite(rounded):
        raise ValueError("A finite real binary64 value is required.")
    return Fraction.from_float(rounded)


def upper_binary64(value: Fraction) -> float | None:
    """Return an upward bound, or None if no finite binary64 bound exists."""
    if value < 0:
        raise ValueError("An error bound cannot be negative.")
    try:
        result = float(value)
    except OverflowError:
        return None
    if not math.isfinite(result):
        return None
    if Fraction.from_float(result) < value:
        result = math.nextafter(result, math.inf)
    return result if math.isfinite(result) else None


def record(value: Fraction) -> dict:
    result = upper_binary64(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator),
            "upper_binary64": result, "finite_binary64_upper_bound": result is not None}


def certify_trace(trace: Sequence[Sequence[Real]], lambda_value: Real, cap: Real,
                  *, price: Real | None = None, left_index: int | None = None,
                  right_weight: Real | None = None, discount: Real | None = None) -> dict:
    """Bound final grid error using exact rational residual propagation.

    A=I+lambda*L has diagonal 1+2*lambda and off-diagonal -lambda.
    The exact target starts at the published initial interior values and has
    fixed boundaries 0/cap. All trace values are interpreted as binary64.
    A supplied scalar price is checked against the exact interpolation target
    with the separately supplied binary64 weight and discount.
    """
    lam, height = binary64(lambda_value), binary64(cap)
    if lam < 0 or height <= 0:
        raise ValueError("lambda must be nonnegative and cap positive.")
    if len(trace) < 2 or len(trace[0]) < 3:
        raise ValueError("At least one step and one interior node are required.")
    width = len(trace[0])
    vectors = []
    for row in trace:
        if len(row) != width:
            raise ValueError("Every trace row must have the same dimension.")
        vector = [binary64(value) for value in row]
        if vector[0] != 0 or vector[-1] != height:
            raise ValueError("Every trace boundary must equal 0/cap exactly.")
        vectors.append(vector)
    total = Fraction(0)
    step_bounds = []
    for previous, current in zip(vectors, vectors[1:]):
        # Both boundary contributions are included when there is one interior node.
        residuals = [(1 + 2 * lam) * current[j] - lam * current[j - 1]
                     - lam * current[j + 1] - previous[j]
                     for j in range(1, width - 1)]
        bound = max(map(abs, residuals))
        step_bounds.append(record(bound))
        total += bound
    range_excess = max(max(Fraction(0), -value, value-height) for row in vectors for value in row)
    initial_in_range = all(0 <= value <= height for value in vectors[0])
    result = {
        "schema": "qbp-exact-heat-trace-certificate-v1",
        "target": "Exact backward-Euler recurrence defined by published binary64 lambda, cap and initial vector",
        "steps": len(vectors) - 1, "interior_nodes": width - 2,
        "initial_error": record(Fraction(0)), "step_residual_norms": step_bounds,
        "final_grid_error_bound": record(total),
        "observed_range_violation": range_excess > 0,
        "initial_range_valid": initial_in_range,
        "max_observed_range_excess": record(range_excess),
        "range_excess_within_residual_bound": range_excess <= total if initial_in_range else None,
        "arithmetic": "Exact integer/Fraction operations on binary64 input values",
        "limitations": "Excludes original-parameter/coordinate/payoff errors, continuum discretisation, finite-domain error and market-model adequacy.",
    }
    scalar = (price, left_index, right_weight, discount)
    if any(value is not None for value in scalar):
        if not all(value is not None for value in scalar):
            raise ValueError("All scalar interpolation arguments are required together.")
        if isinstance(left_index, bool) or not isinstance(left_index, int) or not 0 <= left_index < width - 1:
            raise ValueError("Interpolation index must identify an adjacent node pair.")
        weight, factor, published = binary64(right_weight), binary64(discount), binary64(price)
        if not 0 <= weight <= 1 or factor < 0:
            raise ValueError("Interpolation weight must be in [0,1] and discount nonnegative.")
        interpolated = factor * ((1 - weight) * vectors[-1][left_index]
                                 + weight * vectors[-1][left_index + 1])
        mismatch = abs(published - interpolated)
        result["scalar_price_error_bound"] = record(mismatch + factor * total)
        result["scalar_evaluation_mismatch"] = record(mismatch)
        result["scalar_target"] = "Discounted exact-scheme interpolation with the supplied binary64 index, weight and discount"
    return result
