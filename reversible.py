"""Executable reversible Boolean arithmetic, using LSB-first registers.

All operations expand to X, CNOT and exact Toffoli gates.  Arithmetic is a
deliberately conservative Bennett construction, not an efficient adder design.
``store_gates=False`` streams resource counts; temporary compute histories are
retained only while a snapshot requires their subsequent inverse.  The declared
Clifford+T costing is seven T gates and three T layers per exact Toffoli.
This is a specified decomposition schedule, not an optimality claim.
"""

from array import array
from collections import Counter
from dataclasses import dataclass
from typing import Callable, Iterable, Mapping, Sequence

Register = tuple[int, ...]
X, CNOT, TOFFOLI = 1, 2, 3


@dataclass(frozen=True)
class Snapshot:
    owner: int
    gate_start: int
    emitted_start: int
    active: frozenset[int]
    serial: int


class Circuit:
    """A Boolean permutation circuit with allocation and exact basis replay."""

    def __init__(self, store_gates: bool = True):
        self.store_gates = store_gates
        self._tape = array("i")  # four integers per expanded gate
        self._counts = Counter()
        self._next_wire = 0
        self._active: set[int] = set()
        self._free: list[int] = []
        self._depth: list[int] = []
        self._toffoli_depth: list[int] = []
        self._max_active = 0
        self._snapshots: list[Snapshot] = []
        self._serial = 0
        self.labels: dict[Register, str] = {}

    @property
    def wire_count(self) -> int:
        return self._next_wire

    @property
    def active_wires(self) -> frozenset[int]:
        return frozenset(self._active)

    @property
    def gate_count(self) -> int:
        return sum(self._counts.values())

    @property
    def gates(self):
        """Iterate tuples (operation, control1, control2, target)."""
        if not self.store_gates:
            raise RuntimeError("The streaming circuit has no complete gate tape")
        for i in range(0, len(self._tape), 4):
            yield tuple(self._tape[i:i + 4])

    def alloc(self, n: int, label: str = "", *, fresh: bool = False) -> Register:
        if n < 0:
            raise ValueError("Negative register width")
        reg = []
        for _ in range(n):
            if self._free and not fresh:
                w = self._free.pop()
            else:
                w = self._next_wire
                self._next_wire += 1
                self._depth.append(0)
                self._toffoli_depth.append(0)
            self._active.add(w)
            reg.append(w)
        result = tuple(reg)
        self._max_active = max(self._max_active, len(self._active))
        if label:
            self.labels[result] = label
        return result

    def release(self, reg: Iterable[int]):
        """Release known-zero wires.  Zero is a caller precondition, not inferred."""
        wires = tuple(reg)
        if len(set(wires)) != len(wires) or not set(wires) <= self._active:
            raise ValueError("Release must name distinct allocated wires")
        for w in wires:
            self._active.remove(w)
            self._free.append(w)

    def _emit(self, op: int, c1: int, c2: int, target: int):
        touched = (target,) if op == X else ((c1, target) if op == CNOT else (c1, c2, target))
        if len(set(touched)) != len(touched) or not set(touched) <= self._active:
            raise ValueError("Gate wires must be distinct and allocated")
        self._counts[op] += 1
        if self.store_gates or self._snapshots:
            self._tape.extend((op, c1, c2, target))
        layer = 1 + max(self._depth[w] for w in touched)
        tlayer = max(self._toffoli_depth[w] for w in touched) + (op == TOFFOLI)
        for w in touched:
            self._depth[w] = layer
            self._toffoli_depth[w] = tlayer

    def x(self, target: int):
        self._emit(X, -1, -1, target)

    def cnot(self, control: int, target: int):
        self._emit(CNOT, control, -1, target)

    def toffoli(self, control1: int, control2: int, target: int):
        if control1 == control2:
            self.cnot(control1, target)
        else:
            self._emit(TOFFOLI, control1, control2, target)

    def snapshot(self) -> Snapshot:
        self._serial += 1
        snap = Snapshot(id(self), len(self._tape), self.gate_count, frozenset(self._active), self._serial)
        self._snapshots.append(snap)
        return snap

    def _check_snapshot(self, snap: Snapshot):
        if snap.owner != id(self) or not self._snapshots or self._snapshots[-1] != snap:
            raise ValueError("Snapshots must be closed in last-in first-out order")

    def _reverse_segment(self, start: int, end: int):
        # The fixed end permits the emitted inverse to append to the same tape.
        for i in range(end - 4, start - 1, -4):
            op, c1, c2, target = self._tape[i:i + 4]
            # Released wires in historical gates may have been recycled. They
            # are valid physical wires and must be live for inverse scheduling.
            touched = (target,) if op == X else ((c1, target) if op == CNOT else (c1, c2, target))
            for w in touched:
                if w not in self._active:
                    self._active.add(w)
                    self._free.remove(w)
            self._max_active = max(self._max_active, len(self._active))
            self._emit(op, c1, c2, target)

    def _finish_snapshot(self, snap: Snapshot, retained: Iterable[int] = ()):
        keep = snap.active | set(retained)
        self.release(sorted(self._active - keep))
        self._snapshots.pop()
        if not self.store_gates and not self._snapshots:
            self._tape = array("i")

    def uncompute(self, snap: Snapshot, end_gate_count: int | None = None):
        """Reverse snapshot's body, optionally excluding a later clean action.

        To retain a clean comparison into a flag allocated before the snapshot,
        save ``end=c.gate_count`` before ``compare_into``, then specify
        ``end_gate_count=end`` here. The excluded suffix must preserve every
        body wire and leave its work zero; this is a caller precondition.
        """
        self._check_snapshot(snap)
        if end_gate_count is None:
            end = len(self._tape)
        else:
            if not snap.emitted_start <= end_gate_count <= self.gate_count:
                raise ValueError("Body end is outside this snapshot's emitted gates")
            end = snap.gate_start + 4 * (end_gate_count - snap.emitted_start)
        self._reverse_segment(snap.gate_start, end)
        self._finish_snapshot(snap)

    def copy_uncompute(self, snap: Snapshot, result: Sequence[int]) -> Register:
        """Copy result to virgin wires, then reverse the whole compute history."""
        self._check_snapshot(snap)
        end = len(self._tape)
        # Virgin wires are essential: a released wire could occur in the inverse.
        output = self.alloc(len(result), "copied result", fresh=True)
        for src, dst in zip(result, output):
            self.cnot(src, dst)
        self._reverse_segment(snap.gate_start, end)
        self._finish_snapshot(snap, output)
        return output

    def _clean(self, builder: Callable[[], Sequence[int]], n: int) -> Register:
        # Allocate before snapshot, preventing any inverse gate from touching it.
        output = self.alloc(n, "arithmetic output")
        snap = self.snapshot()
        result = tuple(builder())
        if len(result) != n:
            raise ValueError("Unexpected builder output width")
        end = len(self._tape)
        for src, dst in zip(result, output):
            self.cnot(src, dst)
        self._reverse_segment(snap.gate_start, end)
        self._finish_snapshot(snap)
        return output

    def _majority(self, a: int, b: int, c: int) -> int:
        out = self.alloc(1)[0]
        self.toffoli(a, b, out)
        self.toffoli(a, c, out)
        self.toffoli(b, c, out)
        return out

    def _ripple(self, a: Sequence[int], b: Sequence[int], subtract: bool = False):
        if len(a) != len(b) or not a:
            raise ValueError("Arithmetic requires equal positive widths")
        carry = self.alloc(1)[0]
        bits = []
        for ai, bi in zip(a, b):
            si = self.alloc(1)[0]
            for src in (ai, bi, carry):
                self.cnot(src, si)
            bits.append(si)
            if subtract:
                not_a = self.alloc(1)[0]
                self.x(not_a)
                self.cnot(ai, not_a)
                carry = self._majority(not_a, bi, carry)
            else:
                carry = self._majority(ai, bi, carry)
        return tuple(bits), carry

    def unsigned_add(self, a: Sequence[int], b: Sequence[int]) -> Register:
        """Clean out-of-place (a+b) modulo 2**width, preserving both inputs."""
        return self._clean(lambda: self._ripple(a, b)[0], len(a))

    clean_add = unsigned_add

    def add_inplace(self, a: Sequence[int], b: Sequence[int]):
        """Reversibly replace b by a+b modulo 2**n, preserving disjoint a.

        Carries are computed before any sum bit changes.  To clear each carry
        afterwards, temporarily restore its original b bit, reverse that carry
        computation, and restore its sum bit.  Work returns to zero.  This uses
        n-1 carry wires, with more gates than an optimised Cuccaro construction.
        """
        n = len(a)
        if not n or len(b) != n or set(a).intersection(b):
            raise ValueError("In-place addition requires disjoint equal positive widths")
        if n == 1:
            self.cnot(a[0], b[0])
            return
        carries = self.alloc(n - 1, "in-place carry work")
        self.toffoli(a[0], b[0], carries[0])
        for i in range(1, n - 1):
            self.toffoli(a[i], b[i], carries[i])
            self.toffoli(a[i], carries[i - 1], carries[i])
            self.toffoli(b[i], carries[i - 1], carries[i])
        for i in range(n):
            self.cnot(a[i], b[i])
            if i:
                self.cnot(carries[i - 1], b[i])
        for i in range(n - 2, -1, -1):
            self.cnot(a[i], b[i])
            if i:
                self.cnot(carries[i - 1], b[i])
                self.toffoli(b[i], carries[i - 1], carries[i])
                self.toffoli(a[i], carries[i - 1], carries[i])
                self.toffoli(a[i], b[i], carries[i])
                self.cnot(carries[i - 1], b[i])
            else:
                self.toffoli(a[0], b[0], carries[0])
            self.cnot(a[i], b[i])
        self.release(carries)

    def unsigned_subtract(self, a: Sequence[int], b: Sequence[int]) -> Register:
        """Clean out-of-place (a-b) modulo 2**width, preserving inputs."""
        return self._clean(lambda: self._ripple(a, b, True)[0], len(a))

    clean_subtract = unsigned_subtract

    def unsigned_less_than(self, a: Sequence[int], b: Sequence[int]) -> Register:
        """A clean one-bit output encoding unsigned a<b, preserving inputs."""
        return self._clean(lambda: (self._ripple(a, b, True)[1],), 1)

    def compare_into(self, a: Sequence[int], b: Sequence[int], flag: int):
        """XOR a<b into an existing flag; all work returns to zero."""
        snap = self.snapshot()
        result = self._ripple(a, b, True)[1]
        end = len(self._tape)
        self.cnot(result, flag)
        self._reverse_segment(snap.gate_start, end)
        self._finish_snapshot(snap)

    def multiply(self, a: Sequence[int], b: Sequence[int], output_bits: int | None = None,
                 drop_lsb: int = 0) -> Register:
        """Clean unsigned shift-add product, then exact binary truncation.

        Returns (a*b >> drop_lsb) modulo 2**output_bits. For two fixed-point
        inputs with f fractional bits, drop_lsb=f retains f fractional bits.
        This is floor truncation, not rounding to nearest or saturation.
        """
        if not a or not b or drop_lsb < 0:
            raise ValueError("Positive input widths and nonnegative shift required")
        width = len(a) + len(b)
        if output_bits is None:
            output_bits = max(0, width - drop_lsb)
        if output_bits < 0:
            raise ValueError("Negative output width")

        def build():
            accum = self.alloc(width)
            for j, bj in enumerate(b):
                row = self.alloc(width)
                for i, ai in enumerate(a):
                    self.toffoli(ai, bj, row[i + j])
                accum = self._ripple(accum, row)[0]
            padding = self.alloc(max(0, drop_lsb + output_bits - width))
            extended = tuple(accum) + padding
            return extended[drop_lsb:drop_lsb + output_bits]

        return self._clean(build, output_bits)

    def multiply_low_space(self, a: Sequence[int], b: Sequence[int], output_bits: int | None = None,
                           drop_lsb: int = 0) -> Register:
        """Shift-add product with row and carry work cleared between iterations.

        This exchanges more gates for substantially less peak workspace than
        ``multiply``. It has the same floor-truncation output contract.
        """
        if not a or not b or drop_lsb < 0:
            raise ValueError("Positive input widths and nonnegative shift required")
        width = len(a) + len(b)
        if output_bits is None:
            output_bits = max(0, width - drop_lsb)
        if output_bits < 0:
            raise ValueError("Negative output width")

        def build():
            accum = self.alloc(width, "product accumulator")
            row = self.alloc(width, "product row")
            for j, bj in enumerate(b):
                for i, ai in enumerate(a):
                    self.toffoli(ai, bj, row[i + j])
                self.add_inplace(row, accum)
                for i in range(len(a) - 1, -1, -1):
                    self.toffoli(a[i], bj, row[i + j])
            self.release(row)
            padding = self.alloc(max(0, drop_lsb + output_bits - width))
            return (tuple(accum) + padding)[drop_lsb:drop_lsb + output_bits]

        return self._clean(build, output_bits)

    def multiply_constant(self, a: Sequence[int], integer: int, output_bits: int | None = None,
                          drop_lsb: int = 0) -> Register:
        """Clean multiplication by a known nonnegative classical integer.

        Returns ((a*integer)>>drop_lsb) modulo 2**output_bits. The constant
        need not fit the input/output register. Set-bit partial rows are copied
        with CNOTs. Trailing zero bits are implemented by the selected shift,
        and the first set-bit row requires no addition to a zero accumulator.
        Remaining additions use the same conservative clean-history ripple
        construction. No quantum constant register or partial-product Toffoli
        gates are required.
        """
        if not a or not isinstance(integer, int) or integer < 0 or drop_lsb < 0:
            raise ValueError("Positive input width and nonnegative integer/shift required")
        if output_bits is None:
            output_bits = max(0, len(a) + integer.bit_length() - drop_lsb)
        if output_bits < 0:
            raise ValueError("Negative output width")
        if integer == 0 or output_bits == 0:
            return self.alloc(output_bits, "zero constant product")
        trailing = (integer & -integer).bit_length() - 1
        odd = integer >> trailing
        shift = drop_lsb - trailing
        width = len(a) + odd.bit_length()
        if shift >= width or -shift >= output_bits:
            return self.alloc(output_bits, "zero shifted constant product")
        if odd == 1:
            output = self.alloc(output_bits, "power-of-two constant product")
            for i, dst in enumerate(output):
                src_index = i + shift
                if 0 <= src_index < len(a):
                    self.cnot(a[src_index], dst)
            return output

        def build():
            accum = self.alloc(width, "constant product accumulator")
            for i, ai in enumerate(a):
                self.cnot(ai, accum[i])
            for j in range(1, odd.bit_length()):
                if not (odd >> j) & 1:
                    continue
                row = self.alloc(width, "constant partial row")
                for i, ai in enumerate(a):
                    self.cnot(ai, row[i + j])
                accum = self._ripple(accum, row)[0]
            if shift >= 0:
                padding = self.alloc(max(0, shift + output_bits - width))
                return (tuple(accum) + padding)[shift:shift + output_bits]
            low_zeros = self.alloc(-shift)
            padding = self.alloc(max(0, output_bits - (width - shift)))
            return (low_zeros + tuple(accum) + padding)[:output_bits]

        return self._clean(build, output_bits)

    def lookup_minterm(self, address: Sequence[int], table: Sequence[int], out_bits: int) -> Register:
        """Exact truth-table XOR lookup using fully emitted minterm ladders.

        Missing rows are zero.  A streaming circuit iterates every emitted gate
        and schedules it, but retains no lookup tape unless a snapshot is open.
        The exponential table traversal is an explicit cost, never free QRAM.
        """
        n = len(address)
        if out_bits < 0 or len(table) > 1 << n:
            raise ValueError("Invalid table or output width")
        if any(v < 0 or v >= 1 << out_bits for v in table):
            raise ValueError("Table value does not fit output width")
        output = self.alloc(out_bits, "lookup output")
        ladder = self.alloc(max(0, n - 1), "lookup work")
        for index, value in enumerate(table):
            if value == 0:
                continue
            negatives = [address[i] for i in range(n) if not (index >> i) & 1]
            for w in negatives:
                self.x(w)
            if n == 0:
                for j in range(out_bits):
                    if (value >> j) & 1:
                        self.x(output[j])
            else:
                if n == 1:
                    control = address[0]
                else:
                    self.toffoli(address[0], address[1], ladder[0])
                    for i in range(2, n):
                        self.toffoli(ladder[i - 2], address[i], ladder[i - 1])
                    control = ladder[n - 2]
                for j in range(out_bits):
                    if (value >> j) & 1:
                        self.cnot(control, output[j])
                if n > 1:
                    for i in range(n - 1, 1, -1):
                        self.toffoli(ladder[i - 2], address[i], ladder[i - 1])
                    self.toffoli(address[0], address[1], ladder[0])
            for w in reversed(negatives):
                self.x(w)
        self.release(ladder)
        return output

    def lookup_tree(self, address: Sequence[int], table: Sequence[int], out_bits: int) -> Register:
        """Exact prefix-tree QROM with one reusable selector per address level.

        For a prefix p and address bit b, child starts as p AND b.  A CNOT
        from p changes it to p AND NOT b, permitting both children to share
        a wire.  Each visited internal node uses two exact Toffolis and two
        CNOTs; leaf word writes are explicit CNOTs.  Empty subtrees are pruned.
        This is not the relative-phase four-T unary QROM implementation.
        """
        n = len(address)
        if out_bits < 0 or len(table) > 1 << n:
            raise ValueError("Invalid table or output width")
        if any(v < 0 or v >= 1 << out_bits for v in table):
            raise ValueError("Table value does not fit output width")
        output = self.alloc(out_bits, "tree lookup output")
        if n == 0:
            v = table[0] if table else 0
            for j in range(out_bits):
                if (v >> j) & 1:
                    self.x(output[j])
            return output
        # Prefix counts identify all-zero subtrees without emitting circuits.
        nonzero_prefix = [0]
        for v in table:
            nonzero_prefix.append(nonzero_prefix[-1] + (v != 0))

        def nonempty(lo, hi):
            lo = min(lo, len(table))
            hi = min(hi, len(table))
            return nonzero_prefix[hi] != nonzero_prefix[lo]

        if not nonempty(0, 1 << n):
            return output
        one = self.alloc(1, "tree root selector")[0]
        self.x(one)
        selectors = self.alloc(n, "tree selectors")

        def visit(prefix, bit, lo, hi):
            if not nonempty(lo, hi):
                return
            if bit < 0:
                v = table[lo]
                for j in range(out_bits):
                    if (v >> j) & 1:
                        self.cnot(prefix, output[j])
                return
            child = selectors[bit]
            middle = (lo + hi) // 2
            self.toffoli(prefix, address[bit], child)
            visit(child, bit - 1, middle, hi)
            self.cnot(prefix, child)
            visit(child, bit - 1, lo, middle)
            self.cnot(prefix, child)
            self.toffoli(prefix, address[bit], child)

        visit(one, n - 1, 0, 1 << n)
        self.x(one)
        self.release((one,) + selectors)
        return output

    lookup = lookup_tree

    def run(self, initial: Mapping[Register, int] | None = None) -> list[int]:
        if not self.store_gates:
            raise RuntimeError("Basis replay requires store_gates=True")
        state = [0] * self.wire_count
        assigned: set[int] = set()
        for reg, value in (initial or {}).items():
            if value < 0 or value >= 1 << len(reg) or assigned.intersection(reg):
                raise ValueError("Initial registers overlap or value is out of range")
            assigned.update(reg)
            for i, wire in enumerate(reg):
                state[wire] = (value >> i) & 1
        for op, c1, c2, target in self.gates:
            if op == X or (op == CNOT and state[c1]) or (op == TOFFOLI and state[c1] and state[c2]):
                state[target] ^= 1
        return state

    @staticmethod
    def value(state: Sequence[int], reg: Sequence[int]) -> int:
        return sum(state[w] << i for i, w in enumerate(reg))

    def resources(self) -> dict:
        layers = max(self._toffoli_depth, default=0)
        return {
            "x": self._counts[X], "cnot": self._counts[CNOT],
            "toffoli": self._counts[TOFFOLI], "expanded_gates": self.gate_count,
            "wires": self.wire_count, "peak_active_wires": self._max_active,
            "active_wires_at_end": len(self._active),
            "all_gate_dependency_depth": max(self._depth, default=0),
            "toffoli_layer_depth": layers,
            "t_count_7t_decomposition": 7 * self._counts[TOFFOLI],
            "t_depth_3layer_schedule": 3 * layers,
            "complete_gate_tape_retained": self.store_gates,
            "retained_gate_tape_bytes": self._tape.buffer_info()[1] * self._tape.itemsize,
            "costing_note": "Actual emitted gates; qubit-disjoint dependency schedule; exact 7-T, 3-T-layer Toffoli decomposition specified, not optimised; no noise or physical runtime model.",
        }


def value(state: Sequence[int], register: Sequence[int]) -> int:
    return Circuit.value(state, register)
