"""Independent cross-check of root's block composition against a complete tape."""

import json
import hashlib
from pathlib import Path
import sys

ROOT_SCRIPT = Path(__file__).resolve().parents[1]/'encoding_resources.py'
sys.path.insert(0, str(ROOT_SCRIPT.parent))
import encoding_resources as e


def run_audit():
    records = []
    tables = {}
    for q, by, bz, bu, exp_bits, log_bits, cdf_bits in (
            (8, 2, 2, 2, 3, 3, 4), (10, 3, 3, 3, 3, 4, 4)):
        for kind in ('direct', 'interval', 'full', 'interval_kernel', 'full_kernel'):
            graph = e.build_graph(kind, q, by, bz, bu, exp_bits, log_bits, cdf_bits, 10)
            e.make_tables(graph, tables)
            complete, _, _, _ = e.compile_small(graph, tables)
            actual = complete.resources()
            composed = e.counted_composition(graph, tables, {})
            counts = composed['expanded_gate_counts']
            for label, key in (('X', 'x'), ('CNOT', 'cnot'), ('Toffoli', 'toffoli')):
                assert counts[label] == actual[key], (kind, q, label, counts[label], actual[key])
            assert composed['logical_T_count_exact_7T_Toffoli'] == actual['t_count_7t_decomposition']
            assert composed['logical_qubits_serial_upper_bound'] >= actual['wires'], (kind, q, 'qubits')
            assert composed['serial_all_gate_depth'] >= actual['all_gate_dependency_depth'], (kind, q, 'depth')
            assert composed['logical_T_depth_serial_3_layers_per_Toffoli'] >= actual['t_depth_3layer_schedule'], (kind, q, 'T depth')
            assert composed['includes_finite_uniform_preparation_H_gates'] == sum(graph.inputs.values()) + q
            # Separately materialise the body, then the clean comparator, and
            # compare the measured counts with the exact inverse multiplicity.
            c = e.Circuit()
            inputs = {name: c.alloc(bits) for name, bits in graph.inputs.items()}
            comparison, zero, flag = c.alloc(q), c.alloc(1), c.alloc(1)[0]
            regs = dict(inputs)
            measured_blocks = []
            for name, node in graph.nodes.items():
                before = c.resources()
                key = ('quantile', node.bits, q) if node.kind == 'quantile' else (node.meta, node.bits, q)
                regs[name] = e.emit_node(c, node, [regs[a] for a in node.args], q, graph.width, tables.get(key))
                after = c.resources()
                measured_blocks.append({label: after[key] - before[key] for label, key in
                                        (('X', 'x'), ('CNOT', 'cnot'), ('Toffoli', 'toffoli'))})
            before_comp = c.resources()
            c.compare_into(comparison + zero, regs[graph.output][:q + 1], flag)
            after_comp = c.resources()
            for actual_block, composed_block in zip(measured_blocks, composed['blocks']):
                assert actual_block == composed_block['counts'], (kind, q, 'block', actual_block, composed_block)
            for label, key in (('X', 'x'), ('CNOT', 'cnot'), ('Toffoli', 'toffoli')):
                independent_total = 2 * before_comp[key] + (after_comp[key] - before_comp[key])
                assert independent_total == actual[key], (kind, q, 'inverse multiplicity', label)
            records.append({
                'kind': kind, 'q': q, 'exact_gate_count_agreement': True,
                'per_node_gate_count_agreement': True, 'arithmetic_inverse_counted_twice': True,
                'clean_comparator_counted_once': True,
                'actual_counts': counts,
                'actual_complete_wires': actual['wires'],
                'composed_qubit_upper_bound': composed['logical_qubits_serial_upper_bound'],
                'actual_all_gate_depth': actual['all_gate_dependency_depth'],
                'composed_serial_depth_bound': composed['serial_all_gate_depth'],
                'actual_T_depth_schedule': actual['t_depth_3layer_schedule'],
                'composed_T_depth_bound': composed['logical_T_depth_serial_3_layers_per_Toffoli'],
                'finite_uniform_H_gates': sum(graph.inputs.values()) + q,
                'H_gates_excluded_from_X_CNOT_Toffoli_counts': True,
            })
    report = {
        'source': ROOT_SCRIPT.name,
        'encoding_source_sha256':hashlib.sha256(ROOT_SCRIPT.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'status': 'pass: ten complete circuits; gate counts exact and space/depth bounds conservative',
        'circuits': records,
        'limitation': 'The selected graphs have no aliased binary operand registers. Standalone block allocation presumes distinct operands; arbitrary Graph nodes with a repeated operand require a distinct cache key and emitted aliased operands for exact counts.',
    }
    dest = ROOT_SCRIPT.parent/'results'/'conditioning_composition_audit.json'
    dest.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run_audit()
