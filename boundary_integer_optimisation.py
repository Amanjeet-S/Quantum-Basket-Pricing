"""Exact masked constant propagation for the retained financial integer DAG.

This pass changes neither the source graph's certificate nor integer output.
It preserves every input register, fractional/word width and reflection core.
The unchanged graph remains an admitted implementation baseline. Optimisation
is equally applicable to every financial method and coherent architecture.
"""
from fractions import Fraction as F
from pathlib import Path
import copy, hashlib, json
import integral_value_oracle as original


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def graph_description(graph):
    return {'q': graph.q, 'width': graph.width, 'inputs': graph.inputs, 'output': graph.output,
        'configuration': graph.configuration,
        'nodes': [{'name': name, 'kind': node.kind, 'args': list(node.args),
            'meta': node.meta.as_dict() if isinstance(node.meta, original.ScalarSpec) else str(node.meta),
            'bits': node.bits, 'output_bits': node.output_bits} for name, node in graph.nodes.items()]}


def _constant_value(graph, name, node, arguments, tables):
    """The exact original evaluate semantics, before final output-bit mask."""
    unit = 1 << graph.q
    mask = (1 << graph.width) - 1
    kind = node.kind
    if kind == 'constant':
        value = original.rounded(node.meta * unit)
    elif kind == 'add':
        value = arguments[0] + arguments[1]
    elif kind in ('sub', 'signed_sub'):
        value = arguments[0] - arguments[1]
    elif kind == 'positive_sub':
        value = max(0, arguments[0] - arguments[1])
    elif kind == 'mul':
        value = arguments[0] * arguments[1] >> graph.q
    elif kind == 'scale':
        integer, drop = original.scale_integer(node.meta, graph.q)
        value = arguments[0] * integer >> drop
    elif kind == 'uniform':
        value = (2 * arguments[0] + 1) * (unit >> (node.meta + 1))
    elif kind == 'minimum':
        value = min(arguments)
    elif kind == 'quantile':
        value = tables[original.table_key(name, node)][arguments[0]]
    elif kind == 'function':
        fractional = graph.q + original.span_log(node.meta.span) - node.bits
        if fractional < 0:
            raise ValueError('A negative interpolation remainder width is unsupported')
        index = arguments[0] >> fractional
        remainder = arguments[0] & ((1 << fractional) - 1)
        word = tables[original.table_key(name, node)][index]
        base = word & mask
        magnitude = (word >> graph.width) & mask
        decreasing = (word >> (2 * graph.width)) & 1
        adjustment = magnitude * remainder >> fractional
        value = base - adjustment if decreasing else base + adjustment
    elif kind in ('clamp', 'clamp_signed'):
        value = arguments[0]
        if kind == 'clamp_signed' and value >= (1 << (graph.width - 1)):
            value -= 1 << graph.width
        value = max(0, min(value, original.rounded(node.meta * unit)))
    else:
        raise ValueError('Unsupported original operation: ' + str(kind))
    return value & ((1 << node.output_bits) - 1)


def optimise(graph, tables):
    """Return an equivalent cloned graph, retained tables and rewrite record.

    No input is assumed constant, even when it becomes unused. Full-width
    known nodes are replaced by exactly dyadic constants, preserving names.
    Known narrow-output nodes propagate their value but are not replaced:
    the historical constant emitter allocates the full graph word width.
    Backwards liveness removes operations with no path to the final output.
    """
    if graph.width < graph.q + 1 or graph.output not in graph.nodes:
        raise ValueError('A finished financial graph with its output node is required')
    before = graph_description(graph)
    before_hash = digest(before)
    source_table_hash = digest({str(key): list(value) for key, value in sorted(tables.items(), key=lambda x: str(x[0]))})
    rewritten = copy.deepcopy(graph)
    known = {}
    rewrites = []
    skipped_narrow = []
    preceding = set(graph.inputs)
    for name, node in graph.nodes.items():
        if any(argument not in preceding for argument in node.args):
            raise ValueError('A dependency is missing or not topologically earlier')
        preceding.add(name)
        reason = None
        if node.kind == 'constant' or all(argument in known for argument in node.args):
            value = _constant_value(graph, name, node, [known[a] for a in node.args], tables)
            reason = 'all_arguments_are_exact_masked_constants'
        elif node.kind == 'scale' and original.scale_integer(node.meta, graph.q)[0] == 0:
            value = 0
            reason = 'exact_implemented_scale_integer_is_zero'
        elif node.kind == 'mul' and any(known.get(a) == 0 for a in node.args):
            value = 0
            reason = 'exact_masked_zero_multiplicand'
        elif node.kind in ('sub', 'signed_sub', 'positive_sub') and node.args[0] == node.args[1]:
            value = 0
            reason = 'identical_integer_arguments_subtract_to_zero'
        elif node.kind in ('clamp', 'clamp_signed') and original.rounded(node.meta * (1 << graph.q)) == 0:
            value = 0
            reason = 'clamp_cap_integer_is_zero'
        else:
            continue
        value &= (1 << node.output_bits) - 1
        known[name] = value
        if node.output_bits != graph.width:
            skipped_narrow.append(name)
            continue
        if node.kind == 'constant':
            continue
        replacement = original.IntegralNode('constant', (), meta=F(value, 1 << graph.q),
            bits=0, output_bits=node.output_bits)
        rewritten.nodes[name] = replacement
        rewritten.bounds[name] = (F(value, 1 << graph.q),) * 2
        rewrites.append({'node': name, 'original_kind': node.kind, 'original_args': list(node.args),
            'original_output_bits': node.output_bits, 'integer_value': value,
            'exact_dyadic_constant': str(F(value, 1 << graph.q)), 'reason': reason})
    live = set()
    stack = [graph.output]
    while stack:
        name = stack.pop()
        if name in graph.inputs or name in live:
            continue
        if name not in rewritten.nodes:
            raise ValueError('Missing rewritten dependency')
        live.add(name)
        stack.extend(rewritten.nodes[name].args)
    removed = [name for name in graph.nodes if name not in live]
    rewritten.nodes = {name: node for name, node in rewritten.nodes.items() if name in live}
    rewritten.bounds = {name: bounds for name, bounds in rewritten.bounds.items() if name in graph.inputs or name in live}
    rewritten.constant_names = {value: name for value, name in rewritten.constant_names.items() if name in live}
    if hasattr(rewritten, 'semantic_nodes'):
        rewritten.semantic_nodes = {name: value for name, value in rewritten.semantic_nodes.items() if value in live or value in graph.inputs}
    if hasattr(rewritten, 'table_certificates'):
        rewritten.table_certificates = {name: value for name, value in rewritten.table_certificates.items() if name in live}
    required_tables = {original.table_key(name, node) for name, node in rewritten.nodes.items() if node.kind in ('quantile', 'function')}
    retained_tables = {key: value for key, value in tables.items() if key in required_tables}
    report = {'schema': 'qbp-exact-financial-integer-optimisation-v1', 'status': 'exact integer equivalence by induction',
        'source_graph_sha256': before_hash, 'source_tables_sha256': source_table_hash,
        'rewritten_graph_sha256': digest(graph_description(rewritten)),
        'source_node_count': len(graph.nodes), 'retained_node_count': len(rewritten.nodes),
        'constant_rewrites': rewrites, 'dead_nodes_removed': removed,
        'known_narrow_nodes_left_unchanged': skipped_narrow,
        'original_input_registers_preserved': graph.inputs == rewritten.inputs,
        'word_width_preserved': graph.width == rewritten.width, 'fractional_width_preserved': graph.q == rewritten.q,
        'output_node_preserved': graph.output == rewritten.output,
        'reflection_core_qubits': sum(graph.inputs.values()) + graph.q + 1,
        'source_table_count': len(tables), 'retained_table_count': len(retained_tables),
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Exact implementation optimisation of the identical certified integer output. Original source graph/certificate and unchanged implementation baseline remain retained. No improved financial approximation or changed payoff/model is claimed.'}
    if graph_description(graph) != before:
        raise AssertionError('The optimiser modified its source graph')
    return rewritten, retained_tables, report


def compile_value(graph, tables, mode='exact7t', store_gates=True):
    optimised, retained, report = optimise(graph, tables)
    circuit, inputs, output = original.compile_value(optimised, retained, mode, store_gates)
    return circuit, inputs, output, report


def optimise_graph(graph, tables):
    """Root-driver adapter; original tables may be reused unchanged.

    The retained compiler only reads tables of live lookup nodes, so the
    caller may pass its unchanged original table dictionary when costing.
    """
    optimised, retained_tables, report = optimise(graph, tables)
    return optimised, report


def compile_flag(graph, tables, mode='exact7t', store_gates=True):
    optimised, retained, report = optimise(graph, tables)
    circuit, inputs, comparison, flag = original.compile_flag(optimised, retained, mode, store_gates)
    if len(circuit.reflection_core) != report['reflection_core_qubits']:
        raise AssertionError('The optimiser changed the original reflection core')
    return circuit, inputs, comparison, flag, report


def counted_composition(graph, tables, cache=None, mode='exact7t'):
    optimised, retained, report = optimise(graph, tables)
    profile = original.counted_composition(optimised, retained, cache, mode)
    if profile['reflection_core_qubits'] != report['reflection_core_qubits']:
        raise AssertionError('The counted optimiser changed the original reflection core')
    profile['integer_optimisation'] = report
    profile['implementation'] = 'exact_masked_constant_propagation_and_dead_operation_removal'
    profile['integer_optimisation_source_sha256'] = report['source_sha256']
    return profile
