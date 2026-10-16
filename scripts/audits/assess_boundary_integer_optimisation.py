"""Table-independent q28 structural assessment of the exact optimiser.

For the frozen production geometries there is no constant-dependent lookup.
Thus a no-rewrite result is independent of every scalar ordinate and every
function address allocation. This assessment does not compile or claim a
financial price certificate; the root comparison retains actual tables.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import boundary_value_oracle as o
import boundary_integer_optimisation as z


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    protocol = json.loads((ROOT / 'results/generalisation_protocol.json').read_text())
    preparations = json.loads((ROOT / 'results/higherdim_preparation.json').read_text())
    models = {m['model_id']: m for m in protocol['models']}
    records = []
    for preparation in preparations['records']:
        model = models[preparation['model_id']]
        domain = preparation['domain']
        for case in preparation['cases']:
            for method in ('boundary_r2', 'stable_gauss_r2', 'raw_direct_r2', 'positive_r1'):
                encodings = ('centred', 'positive_part', 'negative_part') if case['cap'] is not None and method != 'positive_r1' else ('centred',)
                for encoding in encodings:
                    graph = o.build_graph(model, case['cap'], q=28, b_y=10, beta_bits=8, table_bits=12,
                        outer_bound=6, U_star=F(domain['U_star']),
                        retained_bound=max(abs(F(v)) for v in domain['coordinate_bounds']),
                        kind=method, output_encoding=encoding)
                    # A lookup with known arguments would require a real table
                    # and raise here. No such lookup appears in these graphs.
                    optimised, retained, proof = z.optimise(graph, {})
                    unchanged = z.graph_description(graph) == z.graph_description(optimised)
                    if not unchanged or proof['constant_rewrites'] or proof['dead_nodes_removed']:
                        raise AssertionError('Production no-effect conclusion requires table/cost evaluation')
                    records.append({'model_id': model['model_id'], 'payoff': case['payoff'], 'cap': case['cap'],
                        'method': method, 'encoding': encoding, 'q': 28,
                        'source_node_count': len(graph.nodes), 'retained_node_count': len(optimised.nodes),
                        'identical_graph': unchanged, 'constant_rewrites': 0, 'dead_nodes_removed': 0,
                        'constant_dependent_lookup_reads': 0,
                        'source_graph_sha256': proof['source_graph_sha256'],
                        'scope': 'Structural zero effect independent of scalar table ordinates and address allocations. No original tables, compilation, complete financial budget or quantum execution are claimed in this assessment.'})
    result = {'schema': 'qbp-boundary-integer-optimisation-production-structure-v1', 'status': 'passed',
        'records': records, 'production_graphs': len(records), 'zero_effect_graphs': sum(r['identical_graph'] for r in records),
        'fractional_width': 28, 'outer_input_bits': 10,
        'scalar_precision_scope': 'All source address allocations preserving the same q28 arithmetic graph. The pass does not test globally flat tables; absent constant-dependent lookups, address width and table contents cannot create a rewrite.',
        'source_sha256': sha(__file__), 'optimiser_source_sha256': sha(ROOT / 'boundary_integer_optimisation.py'),
        'input_sha256': {p: sha(ROOT / p) for p in ['results/generalisation_protocol.json', 'results/higherdim_preparation.json', 'boundary_value_oracle.py']},
        'interpretation': 'The exact pass has zero structural effect on the frozen q28 p3/p7 graphs. Coarse-fixture savings must not be reported as production savings. The unchanged implementation remains fully admitted.'}
    (ROOT / 'results/boundary_integer_optimisation_production_structure.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k in ('status', 'production_graphs', 'zero_effect_graphs')}))


if __name__ == '__main__':
    main()
