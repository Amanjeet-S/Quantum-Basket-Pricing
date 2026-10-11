"""Isolate bound changes at fixed prefix and verify exact emitted cost algebra."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import encoding_resources as e
import generalisation as study
from resource_study import propagate


def internal_nodes(table, bits):
    def visit(lo, hi, depth):
        if not depth or not any(table[lo:hi]):
            return 0
        middle = (lo + hi) // 2
        return 1 + visit(lo, middle, depth-1) + visit(middle, hi, depth-1)
    return visit(0, len(table), bits)


def main():
    resource_path = ROOT/'results/generalisation_resources.json'
    resources = json.loads(resource_path.read_text())
    for name, expected in resources['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected
    models = {m['model_id']: m for m in resources['models']}
    tables, cache, rows = {}, {}, []
    for base in resources['configurations']:
        if base['method'] != 'full_global':
            continue
        model = models[base['model_id']]
        selected = {}
        for mode in ('full_global', 'full_cap'):
            for bits in range(4, min(15, base['q']-1)+1):
                graph = study.build_graph(model, mode, base['q'], base['b_y'], 0, 0,
                    base['exp_address_bits'], base['log_address_bits'], bits, base['cap'])
                certificate = propagate(graph, base['cap'])
                if certificate['arithmetic_price_error_upper'] > .001:
                    continue
                study.make_tables(graph, tables)
                counts = e.counted_composition(graph, tables, cache)
                key = (f"general_cap:{model['lambda']}:{base['cap']}", bits, base['q'])
                selectors = internal_nodes(tables[key], bits)
                remainder = base['q']+2-bits
                block_formula = 2*selectors+2*remainder*(4*graph.width+3*remainder)+6*graph.width
                emitted = [b['counts']['Toffoli'] for b in counts['blocks']
                           if b['function'] == key[0]]
                assert emitted == [block_formula, block_formula]
                selected[mode] = {'bits': bits, 'remainder_bits': remainder,
                    'selector_internal_nodes': selectors, 'block_toffolis': block_formula,
                    'flag_toffolis': counts['expanded_gate_counts']['Toffoli'],
                    'T': counts['logical_T_count_exact_7T_Toffoli'],
                    'arithmetic_price_error_upper': certificate['arithmetic_price_error_upper']}
                break
            assert mode in selected
        a, b = selected['full_global'], selected['full_cap']
        saving = a['flag_toffolis']-b['flag_toffolis']
        assert saving == 4*(a['block_toffolis']-b['block_toffolis'])
        one_bit_formula = None
        if a['bits'] == b['bits']+1:
            f = a['remainder_bits']; width = base['q']+7
            one_bit_formula = 4*(2*(a['selector_internal_nodes']-b['selector_internal_nodes'])
                                -(8*width+12*f+6))
            assert one_bit_formula == saving
        rows.append({'case_id': base['case_id'], 'fixed_q': base['q'],
            'fixed_b_y': base['b_y'], 'fixed_exp_address_bits': base['exp_address_bits'],
            'fixed_log_address_bits': base['log_address_bits'], 'results': selected,
            'flag_toffoli_saving': saving, 'T_saving': 7*saving,
            'one_bit_coarsening_formula': one_bit_formula})
    assert len(rows) == 27
    output = {'scope': 'fixed-prefix causal comparison; minimum feasible probability width in 4..15; complete emitted counts; no global optimisation claim',
        'source_hashes': resources['source_hashes'],
        'resource_sha256': hashlib.sha256(resource_path.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'configurations': rows}
    path = ROOT/'results/generalisation_precision_effects.json'
    path.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'cases': len(rows), 'changed_cases': [r for r in rows if r['T_saving']]}, indent=2))


if __name__ == '__main__':
    main()
