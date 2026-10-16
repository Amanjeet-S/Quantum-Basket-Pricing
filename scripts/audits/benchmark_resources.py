"""Emit exact large primitive gate counts; no whole large tape is retained."""

import json
from math import erf, sqrt
from pathlib import Path
from statistics import NormalDist
from time import perf_counter

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from reversible import Circuit


def count_primitive(name, build):
    c = Circuit(store_gates=False)
    start = perf_counter()
    build(c)
    resources = c.resources()
    resources["compilation_seconds"] = perf_counter() - start
    resources["name"] = name
    assert resources["retained_gate_tape_bytes"] == 0
    return resources


def arithmetic(method, n=36):
    def build(c):
        a, b = c.alloc(n, "a"), c.alloc(n, "b")
        if method in ("multiply", "multiply_low_space"):
            getattr(c, method)(a, b, output_bits=n, drop_lsb=24)
        else:
            getattr(c, method)(a, b)
    return build


def main():
    results = [count_primitive("unsigned_36bit_" + method, arithmetic(method))
               for method in ("unsigned_add", "unsigned_subtract", "unsigned_less_than", "multiply", "multiply_low_space")]
    for label, integer, shift in (
            ('fraction_0p1_q30', round(.1 * (1 << 30)), 30),
            ('fraction_0p05_q30', round(.05 * (1 << 30)), 30),
            ('integer_5', 5, 0), ('integer_10', 10, 0), ('integer_100', 100, 0),
            ('integer_100_times_2pow30_drop30', 100 << 30, 30)):
        results.append(count_primitive('unsigned_36bit_constant_' + label,
            lambda c, integer=integer, shift=shift:
                c.multiply_constant(c.alloc(36), integer, output_bits=36, drop_lsb=shift)))
    normal = NormalDist()
    quantile_table = [int(round(normal.inv_cdf((j + 0.5) / 65536) * (1 << 24))) & ((1 << 32) - 1)
                      for j in range(65536)]
    cdf_table = [int(round(0.5 * (1 + erf((-8 + 16 * j / 4095) / sqrt(2))) * ((1 << 32) - 1)))
                 for j in range(4096)]
    results.append(count_primitive("mid_quantile_16address_32word_24fraction_tree",
                                   lambda c: c.lookup(c.alloc(16), quantile_table, 32)))
    results.append(count_primitive("cdf_12address_32word_tree",
                                   lambda c: c.lookup(c.alloc(12), cdf_table, 32)))
    output = {
        "status": "expanded streaming resource counts, not large basis replay",
        "tables": {"quantile": "Standard-normal midpoint quantiles; 24 fractional bits; signed two's-complement 32-bit words.",
                   "cdf": "Standard-normal CDF on 4096 equally spaced points in [-8,8], mapped to unsigned 32-bit words."},
        "limitations": "Host-generated table contents and traversal; no table approximation or interpolation error certificate implied. No complete pricing architecture or quantum algorithm overhead included.",
        "primitives": results,
    }
    path = Path(__file__).resolve().parents[1]/"results"/"arithmetic_resource_counts.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"output": str(path), "primitives": [
        {k: r[k] for k in ("name", "expanded_gates", "toffoli", "wires", "toffoli_layer_depth", "t_count_7t_decomposition", "compilation_seconds")}
        for r in results]}, indent=2))


if __name__ == "__main__":
    main()
