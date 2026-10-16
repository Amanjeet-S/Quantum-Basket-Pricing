"""Independent finite-unitary witness for deferred financial cleanup.

This imports no financial producer, gate emitter or reflection-controller code.
The proof is algebraic; these literal small-state checks catch parity/control
sign and work-zero invariant failures. This is not a hardware experiment.
"""
from __future__ import annotations
import json
import math
from pathlib import Path

LABELS = (0, 1)
FLAG = 2
WORK = (3, 4)
DIM = 32


def basis(i, size=DIM):
    x = [0j] * size
    x[i] = 1 + 0j
    return x


def permute(v, fn):
    out = [0j] * len(v)
    for i, z in enumerate(v):
        out[fn(i)] += z
    return out


def body(v):
    return permute(v, lambda i: i ^ (((i >> 0) & 1) << 3)
                   ^ (((i >> 1) & 1) << 4))


def predicate(v, k):
    return permute(v, lambda i: i ^ ((1 << FLAG)
                   if (((i >> 3) & 1) + 2 * ((i >> 4) & 1)) < k
                   else 0))


def label_hadamards(v):
    out = list(v)
    for bit in LABELS:
        nxt = list(out)
        for i in range(len(v)):
            if not (i & (1 << bit)):
                j = i | (1 << bit)
                nxt[i] = (out[i] + out[j]) / math.sqrt(2)
                nxt[j] = (out[i] - out[j]) / math.sqrt(2)
        out = nxt
    return out


def a(v, k):
    return body(predicate(body(label_hadamards(v)), k))


def a_inverse(v, k):
    return label_hadamards(body(predicate(body(v), k)))


def b(v, k):
    return predicate(body(label_hadamards(v)), k)


def b_inverse(v, k):
    return label_hadamards(body(predicate(v, k)))


def good_reflection(v):
    return [(-z if i & (1 << FLAG) else z) for i, z in enumerate(v)]


def core_reflection(v):
    return [(-z if (i & 7) == 0 else z) for i, z in enumerate(v)]


def rank_one_reflection(v):
    return [(-z if i == 0 else z) for i, z in enumerate(v)]


def borrowed_reflection(v):
    # A literal clean-scratch three-core-bit zero-reflection tape:
    # X on all core controls, CCX into borrowed work bit3, CZ between
    # bit3 and flag, CCX inverse, then restore all core controls.
    out = permute(v, lambda i: i ^ 7)
    ccx = lambda x: x ^ ((1 << 3) if (x & 3) == 3 else 0)
    out = permute(out, ccx)
    out = [(-z if (i & (1 << 3)) and (i & (1 << FLAG)) else z)
           for i, z in enumerate(out)]
    return permute(permute(out, ccx), lambda i: i ^ 7)


def controlled_borrowed_reflection(v):
    # Actual retained control-last ladder, with an extra borrowed wire.
    out = permute(v, lambda i: i ^ 7)
    ccx1 = lambda i: i ^ ((1 << 3) if (i & 3) == 3 else 0)
    ccx2 = lambda i: i ^ ((1 << 4) if (i & (1 << 3))
                         and (i & (1 << FLAG)) else 0)
    out = permute(permute(out, ccx1), ccx2)
    out = [(-z if (i & (1 << 4)) and (i & (1 << 5)) else z)
           for i, z in enumerate(out)]
    return permute(permute(permute(out, ccx2), ccx1), lambda i: i ^ 7)


def actual_controlled_iterate(v, k, deferred):
    prep, inv = (b, b_inverse) if deferred else (a, a_inverse)
    out = [(-z if (i & (1 << FLAG)) and (i & (1 << 5)) else z)
           for i, z in enumerate(v)]
    out = prep(controlled_borrowed_reflection(inv(out, k)), k)
    return [(-z if i & (1 << 5) else z) for i, z in enumerate(out)]


def q(v, k, deferred, rank_one=False, retain_minus=True,
      literal_borrowing=True):
    prep, inv = (b, b_inverse) if deferred else (a, a_inverse)
    zero = rank_one_reflection if rank_one else (
        borrowed_reflection if literal_borrowing else core_reflection)
    out = prep(zero(inv(good_reflection(v), k)), k)
    return [-z for z in out] if retain_minus else out


def assert_close(x, y):
    assert max(abs(a - b) for a, b in zip(x, y)) < 2e-12


def zero_work_mass(v):
    return sum(abs(z) ** 2 for i, z in enumerate(v) if (i >> 3) == 0)


def check():
    clean_basis_cases = 0
    parity_cases = 0
    invariant_cases = 0
    controlled_sign_cases = 0
    full_operator_cases = 0
    actual_controlled_operator_cases = 0
    controlled_clean_domain_cases = 0
    for k in range(5):
        for i in range(8):
            v = basis(i)
            assert abs(zero_work_mass(a(v, k)) - 1) < 2e-12
            assert abs(zero_work_mass(a_inverse(v, k)) - 1) < 2e-12
            assert_close(b(v, k), body(a(v, k)))
            assert_close(b_inverse(b(v, k), k), v)
            assert_close(borrowed_reflection(v), core_reflection(v))
            clean_basis_cases += 1
        for i in range(DIM):
            v = basis(i)
            assert_close(q(v, k, True), body(q(body(v), k, False)))
            full_operator_cases += 1
        initial = a(basis(0), k)
        dirty_initial = b(basis(0), k)
        p = k / 4
        theta = math.asin(math.sqrt(p))
        for degree in range(1, 13):
            r = degree // 2
            clean = list(initial)
            dirty = list(dirty_initial)
            for _ in range(r):
                clean = q(clean, k, False)
                dirty = q(dirty, k, True)
            assert_close(dirty, body(clean))
            assert abs(zero_work_mass(b_inverse(dirty, k)) - 1) < 2e-12
            assert_close(core_reflection(b_inverse(dirty, k)),
                         rank_one_reflection(b_inverse(dirty, k)))
            assert_close(borrowed_reflection(b_inverse(dirty, k)),
                         core_reflection(b_inverse(dirty, k)))
            assert_close(q(dirty, k, True),
                         body(q(body(dirty), k, False)))
            invariant_cases += 1
            if degree % 2:
                observed = sum(abs(z) ** 2 for i, z in enumerate(dirty)
                               if i & (1 << FLAG))
                expected = math.sin(degree * theta) ** 2
            else:
                decoded = b_inverse(dirty, k)
                observed = sum(abs(z) ** 2 for i, z in enumerate(decoded)
                               if (i & 7) == 0)
                expected = math.cos(degree * theta) ** 2
            assert abs(observed - expected) < 2e-12
            parity_cases += 1
        # A controlled operation retains the global minus on branch one.
        # An omitted minus differs in vector norm even though a single
        # uncontrolled terminal good observation cannot see that phase.
        branch_one = q(dirty_initial, k, True)
        broken_branch = q(dirty_initial, k, True, retain_minus=False)
        assert_close(branch_one, [-z for z in broken_branch])
        assert sum(abs(x - y) ** 2 for x, y in
                   zip(branch_one, broken_branch)) > 3.999999999
        # Explicit control-zero and control-one branches of the conjugacy.
        control_zero = list(dirty_initial)
        control_one = list(branch_one)
        assert_close(control_zero, body(initial))
        assert_close(control_one, body(q(initial, k, False)))
        controlled_sign_cases += 1
        for i in range(64):
            v = basis(i, 64)
            assert_close(actual_controlled_iterate(v, k, True),
                         body(actual_controlled_iterate(body(v), k, False)))
            actual_controlled_operator_cases += 1
        for i in range(8):
            for c in (0, 1):
                v = basis(i | (c << 5), 64)
                expected = [(-z if c and (j & 7) == 0 else z)
                            for j, z in enumerate(v)]
                assert_close(controlled_borrowed_reflection(v), expected)
                controlled_clean_domain_cases += 1
    # The small-core reflection is NOT the global rank-one reflection.
    # Its use is justified by work-zero at the reflection, not everywhere.
    outside = basis(8)
    assert core_reflection(outside) != rank_one_reflection(outside)
    assert borrowed_reflection(outside) != core_reflection(outside)
    # A dirty EXTRA conjunction wire breaks off-domain equality between
    # the actual controlled ladder and controlled(uncontrolled ladder).
    dirty_extra = basis((1 << 4) | (1 << 5), 64)
    assert controlled_borrowed_reflection(dirty_extra) != borrowed_reflection(dirty_extra)
    return {
        'status': 'passed',
        'clean_basis_cases': clean_basis_cases,
        'parity_cases': parity_cases,
        'invariant_cases': invariant_cases,
        'controlled_sign_cases': controlled_sign_cases,
        'full_actual_operator_conjugacy_cases': full_operator_cases,
        'actual_controlled_operator_conjugacy_cases': actual_controlled_operator_cases,
        'actual_controlled_clean_domain_cases': controlled_clean_domain_cases,
        'probabilities': [0, .25, .5, .75, 1],
        'degrees': list(range(1, 13)),
        'state_dimension': DIM,
        'controlled_state_dimension': 64,
        'global_rank_one_equivalence_outside_invariant': False,
        'borrowed_tape_equals_symbolic_core_reflection_outside_invariant': False,
        'borrowed_tape_equals_symbolic_core_reflection_on_M': True,
        'actual_controlled_tape_equals_controlled_uncontrolled_tape_off_M': False,
        'controlled_global_minus_retained': True,
        'numerical_tolerance': '2e-12',
        'scope': 'Independent literal small-unitary witness only; algebraic proof and actual financial opcode/resource invariants remain required.',
    }


if __name__ == '__main__':
    out = check()
    path = Path(__file__).with_name('deferred_invariant_independent_checks.json')
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out))
