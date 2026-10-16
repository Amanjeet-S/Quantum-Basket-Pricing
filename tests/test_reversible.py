"""Independent integer references and basis-state replay of expanded gates."""

import itertools
import random
import unittest

from reversible import Circuit


class ArithmeticTests(unittest.TestCase):
    def verify(self, circuit, initial, output, expected, other_outputs=()):
        state = circuit.run(initial)
        self.assertEqual(circuit.value(state, output), expected)
        for reg, integer in initial.items():
            self.assertEqual(circuit.value(state, reg), integer)
        retained = set(output)
        for reg in initial:
            retained.update(reg)
        for reg in other_outputs:
            retained.update(reg)
        self.assertTrue(all(state[w] == 0 for w in range(circuit.wire_count) if w not in retained))
        return state

    def test_primitive_truth_tables(self):
        for kind in ("x", "cnot", "toffoli"):
            c = Circuit()
            reg = c.alloc(3)
            if kind == "x":
                c.x(reg[2])
            elif kind == "cnot":
                c.cnot(reg[0], reg[2])
            else:
                c.toffoli(reg[0], reg[1], reg[2])
            for integer in range(8):
                trigger = kind == "x" or (kind == "cnot" and integer & 1) or (kind == "toffoli" and integer & 3 == 3)
                state = c.run({reg: integer})
                self.assertEqual(c.value(state, reg), integer ^ (4 if trigger else 0))

    def test_exhaustive_add_subtract_compare(self):
        for n in range(1, 5):
            for operation in ("unsigned_add", "unsigned_subtract", "unsigned_less_than"):
                c = Circuit()
                a, b = c.alloc(n), c.alloc(n)
                out = getattr(c, operation)(a, b)
                for x, y in itertools.product(range(1 << n), repeat=2):
                    expected = ((x + y) if operation == "unsigned_add" else
                                ((x - y) if operation == "unsigned_subtract" else int(x < y)))
                    if operation != "unsigned_less_than":
                        expected %= 1 << n
                    self.verify(c, {a: x, b: y}, out, expected)

    def test_exhaustive_product_and_truncation(self):
        for n, m, shift, out_bits in ((1, 1, 0, 2), (2, 3, 0, 5), (3, 3, 2, 3),
                                      (4, 4, 0, 8), (3, 2, 4, 4), (2, 2, 8, 2)):
            c = Circuit()
            a, b = c.alloc(n), c.alloc(m)
            out = c.multiply(a, b, out_bits, shift)
            for x, y in itertools.product(range(1 << n), range(1 << m)):
                self.verify(c, {a: x, b: y}, out, ((x * y) >> shift) % (1 << out_bits))

    def test_exhaustive_inplace_add(self):
        for n in range(1, 5):
            c = Circuit()
            a, b = c.alloc(n), c.alloc(n)
            c.add_inplace(a, b)
            for x, y in itertools.product(range(1 << n), repeat=2):
                state = c.run({a: x, b: y})
                self.assertEqual(c.value(state, a), x)
                self.assertEqual(c.value(state, b), (x + y) % (1 << n))
                self.assertTrue(all(state[w] == 0 for w in range(c.wire_count) if w not in a + b))

    def test_low_space_multiply(self):
        for n, m, shift, out_bits in ((1, 1, 0, 2), (2, 3, 0, 5), (3, 3, 2, 3),
                                      (4, 4, 0, 8), (3, 2, 4, 4), (2, 2, 8, 2)):
            c = Circuit()
            a, b = c.alloc(n), c.alloc(m)
            out = c.multiply_low_space(a, b, out_bits, shift)
            for x, y in itertools.product(range(1 << n), range(1 << m)):
                self.verify(c, {a: x, b: y}, out, ((x * y) >> shift) % (1 << out_bits))

    def test_aliased_product_inputs(self):
        for method in ("multiply", "multiply_low_space"):
            c = Circuit()
            a = c.alloc(4)
            out = getattr(c, method)(a, a, 8)
            for x in range(16):
                self.verify(c, {a: x}, out, x * x)

    def test_exhaustive_constant_product(self):
        for n, integer, shift, out_bits in (
                (4, 0, 0, 6), (4, 1, 0, 4), (4, 2, 0, 6), (4, 8, 4, 4),
                (4, 5, 0, 8), (4, 10, 0, 8), (4, 100, 0, 12),
                (4, 100 << 8, 8, 12), (4, 13, 7, 4), (4, 2, 20, 4),
                (4, 8, 0, 2), (4, 205, 11, 7), (4, 7, 0, 2)):
            c = Circuit()
            a = c.alloc(n)
            out = c.multiply_constant(a, integer, out_bits, shift)
            for x in range(1 << n):
                self.verify(c, {a: x}, out, ((x * integer) >> shift) % (1 << out_bits))

    def test_fixed_point_constant_and_large_boundary_cases(self):
        rng = random.Random(9104)
        for n, fractional in ((12, 6), (36, 30)):
            mask = (1 << n) - 1
            for integer, shift in ((round(.1 * (1 << fractional)), fractional),
                                   (round(.05 * (1 << fractional)), fractional),
                                   (5, 0), (10, 0), (100, 0), (100 << fractional, fractional)):
                c = Circuit()
                a = c.alloc(n)
                out = c.multiply_constant(a, integer, n, shift)
                values = [0, 1, mask, mask - 1, 1 << (n - 1)]
                values += [rng.randrange(mask + 1) for _ in range(8)]
                for x in values:
                    self.verify(c, {a: x}, out, ((x * integer) >> shift) & mask)

    def test_constant_multiplier_streaming_and_pipeline(self):
        reports = []
        for stored in (True, False):
            c = Circuit(store_gates=stored)
            a, flag = c.alloc(8), c.alloc(1)
            snap = c.snapshot()
            out = c.multiply_constant(a, 205, 8, 11)
            body_end = c.gate_count
            c.compare_into(a, out, flag[0])
            c.uncompute(snap, end_gate_count=body_end)
            reports.append(c.resources())
            if stored:
                for x in range(256):
                    state = c.run({a: x})
                    self.assertEqual(c.value(state, a), x)
                    self.assertEqual(c.value(state, flag), int(x < ((x * 205) >> 11)))
                    self.assertTrue(all(state[w] == 0 for w in range(c.wire_count) if w not in a + flag))
        for key in ('x', 'cnot', 'toffoli', 'wires', 'all_gate_dependency_depth', 'toffoli_layer_depth'):
            self.assertEqual(reports[0][key], reports[1][key])

    def test_larger_boundary_and_seeded_cases(self):
        rng = random.Random(613)
        for n in (8, 12):
            limit = (1 << n) - 1
            pairs = [(0, 0), (0, limit), (limit, 0), (limit, limit), (limit, 1),
                     (1, limit), (1 << (n - 1), 1 << (n - 1))]
            pairs += [(rng.randrange(limit + 1), rng.randrange(limit + 1)) for _ in range(12)]
            for operation in ("unsigned_add", "unsigned_subtract", "unsigned_less_than", "multiply"):
                c = Circuit()
                a, b = c.alloc(n), c.alloc(n)
                out = c.multiply(a, b, n, n // 2) if operation == "multiply" else getattr(c, operation)(a, b)
                for x, y in pairs:
                    if operation == "unsigned_add":
                        expected = (x + y) & limit
                    elif operation == "unsigned_subtract":
                        expected = (x - y) & limit
                    elif operation == "unsigned_less_than":
                        expected = int(x < y)
                    else:
                        expected = (x * y >> (n // 2)) & limit
                    self.verify(c, {a: x, b: y}, out, expected)

    def test_compare_into_nonzero_flag(self):
        c = Circuit()
        a, b, flag = c.alloc(3), c.alloc(3), c.alloc(1)
        c.compare_into(a, b, flag[0])
        for x, y, f in itertools.product(range(8), range(8), range(2)):
            state = c.run({a: x, b: y, flag: f})
            self.assertEqual(c.value(state, flag), f ^ (x < y))
            self.assertEqual(c.value(state, a), x)
            self.assertEqual(c.value(state, b), y)
            self.assertTrue(all(state[w] == 0 for w in range(c.wire_count) if w not in a + b + flag))

    def test_lookup_tree_against_minterm(self):
        tables = [([], 0, 0), ([3], 0, 2), ([0, 1], 1, 1),
                  ([0, 7, 2, 5, 4, 3, 6, 1], 3, 3), ([1, 0, 4], 3, 3),
                  ([0] * 16, 4, 3), ([(j * 31 + 11) % 64 for j in range(32)], 5, 6)]
        for table, n, bits in tables:
            for method in ("lookup_tree", "lookup_minterm"):
                c = Circuit()
                a = c.alloc(n)
                out = getattr(c, method)(a, table, bits)
                for i in range(1 << n):
                    self.verify(c, {a: i}, out, table[i] if i < len(table) else 0)

    def test_clean_wire_reuse_and_whole_pipeline_inverse(self):
        c = Circuit()
        a, b = c.alloc(4), c.alloc(4)
        snap = c.snapshot()
        summation = c.unsigned_add(a, b)
        product = c.multiply(summation, b, 4, 2)
        difference = c.unsigned_subtract(product, a)
        out = c.copy_uncompute(snap, difference)
        self.assertEqual(c.active_wires, frozenset(a + b + out))
        for x, y in itertools.product(range(16), repeat=2):
            expected = (((((x + y) & 15) * y) >> 2) - x) & 15
            self.verify(c, {a: x, b: y}, out, expected)
        # Reuse cleaned wire numbers in a subsequent preserved-output operation.
        second = c.unsigned_add(a, b)
        for x, y in ((0, 15), (15, 15), (8, 9)):
            state = self.verify(c, {a: x, b: y}, second, (x + y) & 15, (out,))
            self.assertEqual(c.value(state, out), (((((x + y) & 15) * y) >> 2) - x) & 15)

    def test_uncompute_pipeline(self):
        c = Circuit()
        a, b = c.alloc(3), c.alloc(3)
        snap = c.snapshot()
        s = c.unsigned_add(a, b)
        c.lookup(s, [1, 2, 3, 4, 5, 6, 7, 0], 3)
        c.uncompute(snap)
        for x, y in itertools.product(range(8), repeat=2):
            state = c.run({a: x, b: y})
            self.assertEqual(c.value(state, a), x)
            self.assertEqual(c.value(state, b), y)
            self.assertTrue(all(state[w] == 0 for w in range(c.wire_count) if w not in a + b))

    def test_clean_compare_retained_after_body_inverse(self):
        for stored in (True, False):
            c = Circuit(store_gates=stored)
            a, b, flag = c.alloc(3), c.alloc(3), c.alloc(1)
            snap = c.snapshot()
            s = c.unsigned_add(a, b)
            body_end = c.gate_count
            c.compare_into(s, a, flag[0])
            c.uncompute(snap, end_gate_count=body_end)
            self.assertEqual(c.active_wires, frozenset(a + b + flag))
            if stored:
                for x, y, f in itertools.product(range(8), range(8), range(2)):
                    state = c.run({a: x, b: y, flag: f})
                    self.assertEqual(c.value(state, flag), f ^ (((x + y) % 8) < x))
                    self.assertEqual(c.value(state, a), x)
                    self.assertEqual(c.value(state, b), y)
                    self.assertTrue(all(state[w] == 0 for w in range(c.wire_count) if w not in a + b + flag))

    def test_streaming_resource_agreement(self):
        results = []
        for stored in (True, False):
            c = Circuit(store_gates=stored)
            a, b = c.alloc(5), c.alloc(5)
            snap = c.snapshot()
            summation = c.unsigned_add(a, b)
            product = c.multiply(summation, b, 5, 2)
            c.lookup(product, [(j * 7) % 32 for j in range(32)], 5)
            c.uncompute(snap)
            results.append(c.resources())
        for key in ("x", "cnot", "toffoli", "wires", "expanded_gates", "all_gate_dependency_depth",
                    "toffoli_layer_depth", "t_count_7t_decomposition", "t_depth_3layer_schedule"):
            self.assertEqual(results[0][key], results[1][key])
        self.assertEqual(results[1]["retained_gate_tape_bytes"], 0)

    def test_scheduling_and_t_count(self):
        c = Circuit()
        a = c.alloc(6)
        c.toffoli(a[0], a[1], a[2])
        c.toffoli(a[3], a[4], a[5])
        self.assertEqual(c.resources()["toffoli_layer_depth"], 1)
        c.cnot(a[2], a[3])
        c.toffoli(a[3], a[4], a[5])
        r = c.resources()
        self.assertEqual(r["all_gate_dependency_depth"], 3)
        self.assertEqual(r["toffoli_layer_depth"], 2)
        self.assertEqual(r["t_count_7t_decomposition"], 21)
        self.assertEqual(r["t_depth_3layer_schedule"], 6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
