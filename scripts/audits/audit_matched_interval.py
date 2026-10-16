"""Post-study matched-prefix diagnostic, kept distinct from the frozen search."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import encoding_resources as e
import generalisation as study
from resource_study import propagate
from audit_precision_effects import internal_nodes


def main():
    path = ROOT/'results/generalisation_resources.json'
    data = json.loads(path.read_text())
    for name, expected in data['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected
    models = {m['model_id']: m for m in data['models']}
    tables, cache, rows = {}, {}, []
    for base in data['configurations']:
        if base['method'] != 'interval':
            continue
        model = models[base['model_id']]
        counted = {}
        for mode in ('interval', 'full_global', 'full_cap'):
            graph = study.build_graph(model, mode, base['q'], base['b_y'], 0,
                base['b_u'] if mode == 'interval' else 0,
                base['exp_address_bits'], base['log_address_bits'],
                base['probability_address_bits'], base['cap'])
            certificate = propagate(graph, base['cap'])
            assert certificate['arithmetic_price_error_upper'] <= .001
            study.make_tables(graph, tables)
            counts = e.counted_composition(graph, tables, cache)
            if mode == 'interval':
                assert counts['expanded_gate_counts'] == base['resources']['expanded_gate_counts']
            counted[mode] = {'T': counts['logical_T_count_exact_7T_Toffoli'],
                'flag_toffolis': counts['expanded_gate_counts']['Toffoli'],
                'arithmetic_price_error_upper': certificate['arithmetic_price_error_upper']}
            if mode == 'interval':
                strike = next(b['counts']['Toffoli'] for b in counts['blocks']
                              if b['function'] == f"general_strike:{base['cap']}")
        q = base['q']; b = base['probability_address_bits']; w = q+7
        cdfkey = (f"general_cdf:{model['lambda']}", b, q)
        capkey = (f"general_cap:{model['lambda']}:{base['cap']}", b, q)
        cdf, cap = tables[cdfkey], tables[capkey]
        assert all(not row or upper for row, upper in zip(cap, cdf))
        i_cdf, i_cap = internal_nodes(cdf, b), internal_nodes(cap, b)
        assert i_cap <= i_cdf
        exact_gap = 2*strike+24*w+8*(i_cdf-i_cap)
        actual_gap = counted['interval']['flag_toffolis']-counted['full_global']['flag_toffolis']
        assert actual_gap == exact_gap
        assert counted['full_global']['T'] < counted['interval']['T']
        rows.append({'case_id': base['case_id'], 'q': q, 'b_y': base['b_y'],
            'exp_address_bits': base['exp_address_bits'],
            'log_address_bits': base['log_address_bits'], 'probability_address_bits': b,
            'results': counted, 'cdf_internal_nodes': i_cdf,
            'cap_internal_nodes': i_cap, 'strike_block_toffolis': strike,
            'exact_flag_toffoli_gap': exact_gap})
    assert len(rows) == 27
    output = {'status': 'post-study matched-prefix diagnostic; original prespecified outputs preserved',
        'question': 'Can the interval winner be replaced by a certified full kernel at identical prefix and table precision?',
        'scope': 'all 27 retained interval prefixes; not a search optimum or hindsight-prespecified comparison',
        'resource_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_hashes': data['source_hashes'],
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'configurations': rows}
    target = ROOT/'results/generalisation_matched_interval.json'
    target.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'cases': len(rows), 'exception': [r for r in rows
        if r['case_id'] == 'n8_lambda0.4_cap25']}, indent=2))


if __name__ == '__main__':
    main()
