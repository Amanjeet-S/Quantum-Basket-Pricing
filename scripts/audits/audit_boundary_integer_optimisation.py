"""Portable independent integer, clean value/flag and accounting fixtures."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json, sys, time
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import boundary_value_oracle as financial
import boundary_integer_optimisation as optimisation
import integral_value_oracle as original
from integral_audit_support import packed_replay, register_values, independent_gate_accounting


def nearest(value):
    value = F(value) + F(1, 2)
    return value.numerator // value.denominator


def interpret(graph, sample, tables):
    """Independent faithful masked semantics; no graph evaluate call."""
    values = dict(sample)
    scale = 1 << graph.q
    modulus = 1 << graph.width
    for name, node in graph.nodes.items():
        a = [values[x] for x in node.args]
        op = node.kind
        if op == 'constant':
            value = nearest(F(node.meta) * scale)
        elif op == 'add':
            value = a[0] + a[1]
        elif op in ('sub', 'signed_sub'):
            value = a[0] - a[1]
        elif op == 'positive_sub':
            value = max(0, a[0] - a[1])
        elif op == 'mul':
            value = a[0] * a[1] // scale
        elif op == 'scale':
            coefficient = F(node.meta)
            dyadic = coefficient.denominator & (coefficient.denominator - 1) == 0
            integer = coefficient.numerator if dyadic else nearest(coefficient * scale)
            drop = coefficient.denominator.bit_length() - 1 if dyadic else graph.q
            value = a[0] * integer // (1 << drop)
        elif op == 'uniform':
            value = (2 * a[0] + 1) * (scale // (1 << (node.meta + 1)))
        elif op == 'minimum':
            value = min(a)
        elif op in ('quantile', 'function'):
            if op == 'quantile':
                key = ('quantile', node.bits, str(node.meta))
                value = tables[key][a[0]]
            else:
                key = ('scalar', json.dumps(node.meta.as_dict(), sort_keys=True), node.bits)
                span = F(node.meta.span)
                exponent = span.numerator.bit_length() - span.denominator.bit_length()
                fractional = graph.q + exponent - node.bits
                address, remainder = divmod(a[0], 1 << fractional)
                packed = tables[key][address]
                base = packed % modulus
                magnitude = (packed // modulus) % modulus
                sign = packed // (modulus * modulus) % 2
                correction = magnitude * remainder // (1 << fractional)
                value = base - correction if sign else base + correction
        elif op in ('clamp', 'clamp_signed'):
            source = a[0] - modulus if op == 'clamp_signed' and a[0] >= modulus // 2 else a[0]
            value = max(0, min(source, nearest(F(node.meta) * scale)))
        else:
            raise ValueError(op)
        values[name] = value % (1 << node.output_bits)
    return values[graph.output], values


def normalise(counts):
    result = Counter()
    for gate, count in counts.items():
        result[{'TD': 'Tdg', 'SD': 'Sdg'}.get(gate, gate)] += count
    return result


def check_clean(graph, tables, mode, samples, outputs, counts_cache):
    value_circuit, inputs, output = original.compile_value(graph, tables, mode)
    assignments, expected = [], []
    for sample, value in zip(samples, outputs):
        for dirty in range(1 << len(output)):
            lane = {inputs[n]: sample[n] for n in inputs}
            lane[output] = dirty
            assignments.append(lane)
            expected.append(value ^ dirty)
    state, _ = packed_replay(value_circuit, assignments)
    if register_values(state, output, len(assignments)) != expected:
        raise AssertionError('Optimised value XOR differs from the independent integer answer')
    keep = set(output + tuple(w for register in inputs.values() for w in register))
    if any(word != 0 for wire, word in enumerate(state) if wire not in keep):
        raise AssertionError('Optimised value work did not return to zero')
    for register in inputs.values():
        if register_values(state, register, len(assignments)) != [lane[register] for lane in assignments]:
            raise AssertionError('An original value input register changed')
    restored, _ = packed_replay(value_circuit, assignments, reverse=True, initial_state=state)
    if register_values(restored, output, len(assignments)) != [lane[output] for lane in assignments]:
        raise AssertionError('Optimised value adjoint failed dirty-output restoration')
    value_lanes = len(assignments)
    del value_circuit, assignments, expected, state, restored
    flag_circuit, inputs, comparison, flag = original.compile_flag(graph, tables, mode)
    assignments, expected = [], []
    for sample, value in zip(samples, outputs):
        for comparator in range(1 << graph.q):
            for dirty in (0, 1):
                lane = {inputs[n]: sample[n] for n in inputs}
                lane[comparison], lane[(flag,)] = comparator, dirty
                assignments.append(lane)
                expected.append(dirty ^ int(comparator < value))
    state, _ = packed_replay(flag_circuit, assignments)
    if register_values(state, (flag,), len(assignments)) != expected:
        raise AssertionError('Optimised flag XOR differs from the independent threshold answer')
    keep = set(comparison + (flag,) + tuple(w for register in inputs.values() for w in register))
    if any(word != 0 for wire, word in enumerate(state) if wire not in keep):
        raise AssertionError('Optimised flag work did not return to zero')
    for register in list(inputs.values()) + [comparison]:
        if register_values(state, register, len(assignments)) != [lane[register] for lane in assignments]:
            raise AssertionError('An original flag/comparison input register changed')
    restored, _ = packed_replay(flag_circuit, assignments, reverse=True, initial_state=state)
    if register_values(restored, (flag,), len(assignments)) != [lane[(flag,)] for lane in assignments]:
        raise AssertionError('Optimised flag adjoint failed dirty-target restoration')
    accounting = independent_gate_accounting(flag_circuit)
    actual_counts = normalise(accounting['expanded_literal_gate_counts'])
    actual_counts['H'] += len(flag_circuit.uniform_preparation_gates)
    profile = original.counted_composition(graph, tables, counts_cache, mode)
    if actual_counts != Counter(profile['primitive_counts']) or accounting['logical_T_count'] != profile['logical_T_count']:
        raise AssertionError('Optimised counted profile differs from independently expanded gates')
    if profile['logical_qubits_upper'] < flag_circuit.wire_count:
        raise AssertionError('Optimised wire bound is below the emitted circuit')
    if tuple(flag_circuit.reflection_core) != tuple(w for register in inputs.values() for w in register) + comparison + (flag,):
        raise AssertionError('The original input/comparison/flag reflection core changed')
    return {'dirty_value_basis_lanes': value_lanes, 'dirty_flag_basis_lanes': len(assignments),
        'T_count': profile['logical_T_count'], 'emitted_logical_wires': flag_circuit.wire_count,
        'reflection_core_qubits': len(flag_circuit.reflection_core),
        'literal_tape_sha256': accounting['expanded_unitary_tape_sha256']}


def arithmetic_fixture(kind):
    graph = original.IntegralGraph(3, 1)
    graph.width = 8
    graph.inputs = {'unused_or_variable': 3}
    graph.bounds = {'unused_or_variable': (F(0), F(7))}
    graph.configuration = {'kind': kind, 'fixture': 'masked arithmetic rule coverage'}
    def node(name, op, args=(), meta=None, output_bits=8):
        graph.nodes[name] = original.IntegralNode(op, args, meta=meta, output_bits=output_bits)
        graph.bounds[name] = (F(0), F(255, 8))
        return name
    if kind == 'zero_rounded_scale':
        result = node('scale', 'scale', ('unused_or_variable',), F(1, 10000))
    elif kind == 'masked_negative_signed_clamp':
        constant = node('negative_bits', 'constant', meta=F(255, 8))
        result = node('signed_clip', 'clamp_signed', (constant,), F(1))
    elif kind == 'masked_constant_overflow':
        a = node('a', 'constant', meta=F(250, 8))
        b = node('b', 'constant', meta=F(20, 8))
        result = node('sum', 'add', (a, b))
    elif kind == 'narrow_known_minimum':
        a = node('a', 'constant', meta=F(6, 8))
        b = node('b', 'constant', meta=F(3, 8))
        low = node('narrow', 'minimum', (a, b), output_bits=3)
        result = node('lift', 'scale', (low,), F(1))
    elif kind == 'masked_subtraction':
        a = node('a', 'constant', meta=F(2, 8))
        b = node('b', 'constant', meta=F(5, 8))
        result = node('difference', 'signed_sub', (a, b))
    elif kind == 'zero_multiplicand':
        zero = node('zero', 'constant', meta=F(0))
        result = node('product', 'mul', ('unused_or_variable', zero))
    else:
        raise ValueError(kind)
    graph.output = result
    return graph


def main():
    start = time.perf_counter()
    source_hash = optimisation.digest(optimisation.graph_description(arithmetic_fixture('zero_rounded_scale')))
    protocol = json.loads((ROOT / 'results/generalisation_protocol.json').read_text())
    models = protocol['models']
    table_cache, counts_cache = {}, {}
    records = []
    input_states = value_lanes = flag_lanes = 0
    for dimension in (3, 7):
        model = next(m for m in models if m['outer_dimension'] == dimension and F(m['lambda']) == F(1, 5))
        for cap in (25, None):
            for method in ('boundary_r2', 'stable_gauss_r2', 'raw_direct_r2', 'positive_r1'):
                encodings = ('centred', 'positive_part', 'negative_part') if cap is not None and method != 'positive_r1' else ('centred',)
                for encoding in encodings:
                    graph = financial.build_graph(model, cap, q=5, b_y=1, beta_bits=1, table_bits=2,
                        outer_bound=2, U_star=F(1, 8), kind=method, output_encoding=encoding)
                    tables = financial.make_tables(graph, table_cache)
                    optimised, retained, report = optimisation.optimise(graph, tables)
                    names = tuple(graph.inputs)
                    samples = [dict(zip(names, values)) for values in product(*(range(1 << graph.inputs[n]) for n in names))]
                    original_outputs = [interpret(graph, sample, tables)[0] for sample in samples]
                    new_outputs = [interpret(optimised, sample, retained)[0] for sample in samples]
                    if original_outputs != new_outputs:
                        raise AssertionError('Optimisation changed a financial integer output')
                    if graph.inputs != optimised.inputs or graph.width != optimised.width or graph.q != optimised.q:
                        raise AssertionError('Optimisation changed the original registers or widths')
                    input_states += len(samples)
                    for mode in ('exact7t', 'unitary8t'):
                        checks = check_clean(optimised, retained, mode, samples, new_outputs, counts_cache)
                        baseline = original.counted_composition(graph, tables, counts_cache, mode)
                        if baseline['reflection_core_qubits'] != checks['reflection_core_qubits']:
                            raise AssertionError('Optimised and baseline reflection cores differ')
                        value_lanes += checks['dirty_value_basis_lanes']
                        flag_lanes += checks['dirty_flag_basis_lanes']
                        records.append({'model_id': model['model_id'], 'outer_dimension': dimension,
                            'payoff': 'vanilla' if cap is None else 'capped', 'cap': cap, 'method': method,
                            'encoding': encoding, 'mode': mode, 'finite_input_states': len(samples),
                            'original_T_count': baseline['logical_T_count'], 'optimised_T_count': checks['T_count'],
                            'T_reduction': baseline['logical_T_count'] - checks['T_count'],
                            'independent_circuit_checks': checks, 'optimisation_report': report,
                            'scope': 'Coarse complete-circuit fixture, not an attained financial monetary tolerance.'})
                    print(dimension, cap, method, encoding, 'nodes', len(graph.nodes), len(optimised.nodes), flush=True)
    rule_records = []
    for kind in ('zero_rounded_scale', 'masked_negative_signed_clamp', 'masked_constant_overflow',
                 'narrow_known_minimum', 'masked_subtraction', 'zero_multiplicand'):
        graph = arithmetic_fixture(kind)
        optimised, retained, report = optimisation.optimise(graph, {})
        samples = [{'unused_or_variable': j} for j in range(8)]
        outputs = [interpret(graph, s, {})[0] for s in samples]
        if outputs != [interpret(optimised, s, retained)[0] for s in samples]:
            raise AssertionError('A masked arithmetic rule fixture changed output')
        # Comparator/value targets expose q+1 bits. This fixture family must
        # clamp larger modular results before those interfaces are exercised.
        if max(outputs) <= 1 << graph.q:
            for mode in ('exact7t', 'unitary8t'):
                check_clean(optimised, retained, mode, samples, outputs, counts_cache)
        rule_records.append({'kind': kind, 'integer_outputs': outputs, 'rewrite_report': report})
    geometry = []
    for model in models:
        if model['outer_dimension'] not in (3, 7) or F(model['lambda']) != F(1, 20):
            continue
        rg = original.model_parameters(model)[4]
        geometry.append({'model_id': model['model_id'], 'R_G_exact': list(map(str, rg)),
            'R_G_decimal_diagnostic': list(map(float, rg)),
            'absolute_R_G_q28_scale_integers': [original.scale_integer(abs(v), 28)[0] for v in rg],
            'all_q28_scale_integers_nonzero': all(original.scale_integer(abs(v), 28)[0] != 0 for v in rg),
            'correction': 'The frozen R_G is not zero or nearly zero at q28. No constant-PG or constant-PGprime assumption is made. Any future exact simplification must follow actual integer semantics.'})
    result = {'schema': 'qbp-boundary-integer-optimisation-independent-audit-v1', 'status': 'passed',
        'financial_graphs': 28, 'financial_graph_architecture_configurations': len(records),
        'exhaustive_financial_input_states': input_states,
        'dirty_value_basis_lanes': value_lanes, 'dirty_flag_basis_lanes': flag_lanes,
        'arithmetic_rule_fixtures': rule_records, 'records': records,
        'frozen_geometric_loading_correction': geometry,
        'optimiser_source_sha256': hashlib.sha256((ROOT / 'boundary_integer_optimisation.py').read_bytes()).hexdigest(),
        'audit_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in
            ['boundary_value_oracle.py', 'integral_value_oracle.py', 'checks/integral_audit_support.py',
             'results/generalisation_protocol.json', 'results/boundary_derivative_norms.json']},
        'generation_seconds_diagnostic': time.perf_counter() - start,
        'scope': 'Exact integer equivalence, preserved input/reflection registers, clean dirty-target value and flag, literal adjoint and independently counted emitted primitives. Primitive phase correctness is inherited from the unchanged previously checked compiler. Neither quantum execution nor practical quantum advantage is claimed.'}
    (ROOT / 'results/boundary_integer_optimisation_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k in ('status', 'financial_graphs',
        'financial_graph_architecture_configurations', 'exhaustive_financial_input_states',
        'dirty_value_basis_lanes', 'dirty_flag_basis_lanes', 'generation_seconds_diagnostic')}))


if __name__ == '__main__':
    main()
