"""Targeted exact constant-table interpolation and lookup rule audit."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import integral_value_oracle as original
import boundary_integer_optimisation as z
from audit_boundary_integer_optimisation import interpret, check_clean


def fixture(kind):
    g = original.IntegralGraph(3, 1)
    g.width = 8
    g.inputs = {'original_unused_input': 3}
    g.configuration = {'kind': kind, 'fixture': 'constant packed-table interpolation'}
    g.bounds = {'original_unused_input': (F(0), F(7))}
    def node(name, op, args=(), meta=None, bits=0, output_bits=None):
        g.nodes[name] = original.IntegralNode(op, args, meta, bits, g.width if output_bits is None else output_bits)
        g.bounds[name] = (F(0), F(255, 8))
        return name
    tables = {}
    if kind.startswith('function'):
        argument = 8 if kind.endswith('endpoint') else 3
        a = node('known_argument', 'constant', meta=F(argument, 8))
        spec = original.ScalarSpec('fixture', F(0), F(1), (), (F(0), F(4)))
        fun = node('known_function', 'function', (a,), spec, 1)
        if kind == 'function_increasing':
            words = (2 | (4 << 8), 6 | (2 << 8), 8)
        else:
            words = (20 | (8 << 8) | (1 << 16), 12 | (4 << 8) | (1 << 16), 8)
        tables[original.table_key(fun, g.nodes[fun])] = words
        result = node('bounded_output', 'clamp', (fun,), F(1))
    elif kind == 'constant_quantile':
        a = node('known_label', 'constant', meta=F(1, 8))
        fun = node('known_quantile', 'quantile', (a,), F(2), 1)
        tables[original.table_key(fun, g.nodes[fun])] = (6, 7)
        result = node('bounded_output', 'clamp', (fun,), F(1))
    elif kind == 'constant_uniform':
        a = node('known_label', 'constant', meta=F(1, 8))
        label = node('narrow_known_label', 'minimum', (a, a), 1, output_bits=1)
        result = node('known_uniform', 'uniform', (label,), 1)
    elif kind == 'identical_positive_subtraction':
        word = node('full_word_argument', 'scale', ('original_unused_input',), F(1))
        result = node('known_zero', 'positive_sub', (word, word))
    else:
        raise ValueError(kind)
    g.output = result
    return g, tables


def main():
    records = []
    caches = {}
    for kind in ('function_increasing', 'function_decreasing', 'function_decreasing_endpoint',
                 'constant_quantile', 'constant_uniform', 'identical_positive_subtraction'):
        graph, tables = fixture(kind)
        optimised, retained, report = z.optimise(graph, tables)
        samples = [{'original_unused_input': k} for k in range(8)]
        outputs = [interpret(graph, sample, tables)[0] for sample in samples]
        if outputs != [interpret(optimised, sample, retained)[0] for sample in samples]:
            raise AssertionError('A constant-table fixture changed an integer output')
        for mode in ('exact7t', 'unitary8t'):
            baseline = check_clean(graph, tables, mode, samples, outputs, caches)
            rewritten = check_clean(optimised, retained, mode, samples, outputs, caches)
            if baseline['reflection_core_qubits'] != rewritten['reflection_core_qubits']:
                raise AssertionError('Constant-table rewrite changed the original reflection core')
            records.append({'kind': kind, 'architecture': mode, 'outputs': outputs,
                'baseline': baseline, 'optimised': rewritten, 'rewrite_record': report})
    result = {'schema': 'qbp-integer-constant-table-targeted-audit-v1', 'status': 'passed',
        'fixtures': 6, 'fixture_architecture_configurations': len(records), 'records': records,
        'audit_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'optimiser_sha256': hashlib.sha256((ROOT / 'boundary_integer_optimisation.py').read_bytes()).hexdigest(),
        'reused_audit_sha256': hashlib.sha256((ROOT / 'checks/audit_boundary_integer_optimisation.py').read_bytes()).hexdigest(),
        'scope': 'Constant interpolation in both directions and at the final endpoint, quantile and uniform lookup, and exact alias subtraction. Both original and optimised clean dirty-target circuits and literal inverse/accounting are checked. No financial model or monetary accuracy claim.'}
    (ROOT / 'results/boundary_integer_constant_table_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'fixtures', 'fixture_architecture_configurations')}))


if __name__ == '__main__':
    main()
