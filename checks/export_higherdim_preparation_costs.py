"""Retain direct source-level operation counts of shared preparation."""
from pathlib import Path
import hashlib, json
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    path = ROOT / 'results/higherdim_preparation.json'
    data = json.loads(path.read_text())
    rows = []
    for record in data['records']:
        n = record['moments']['assets']
        for case in record['cases']:
            capped = case['payoff'] == 'capped'
            rows.append({'model_id': record['model_id'], 'payoff': case['payoff'], 'cap': case['cap'],
                'coefficient_kernel_operation_counts': {'exponential_calls': 2 * n + 2,
                    'Gaussian_CDF_calls': 2 * n + 6 if capped else n + 3,
                    'logarithm_calls': n + 2 if capped else 0, 'square_root_calls': 1},
                'count_basis': 'Direct source-level tally of boundary_grid_certificate.coefficient_certificate, including repeated arguments. No scalar caching claimed for coefficients.'})
    result = {'schema': 'qbp-higherdim-preparation-cost-ledger-v1', 'status': 'retained explicit shared work',
        'preparation_sha256': sha(path), 'source_sha256': sha(__file__),
        'coefficient_source_sha256': sha(ROOT / 'boundary_grid_certificate.py'),
        'coefficient_cases': rows, 'totals': {
            'coefficient_exponential_calls': sum(row['coefficient_kernel_operation_counts']['exponential_calls'] for row in rows),
            'coefficient_Gaussian_CDF_calls': sum(row['coefficient_kernel_operation_counts']['Gaussian_CDF_calls'] for row in rows),
            'coefficient_logarithm_calls': sum(row['coefficient_kernel_operation_counts']['logarithm_calls'] for row in rows),
            'coefficient_square_root_calls': sum(row['coefficient_kernel_operation_counts']['square_root_calls'] for row in rows),
            'multinomial_and_correlated_gap_exponential_evaluations': sum(row['moments']['additional_multinomial_exponential_evaluations'] for row in data['records']),
            'vertex_exponential_evaluations': sum(row['domain']['additional_vertex_exponential_evaluations'] for row in data['records'])},
        'shared_MGF_profile': data['shared_MGF_profile'],
        'finite_MGF_evaluation_policy': {
            'nodal_exponential_evaluations': 0,
            'remainder_exponential_evaluations': data['shared_MGF_profile']['remainder_exponential_evaluations'],
            'empirical_even_power_multiplications': data['shared_MGF_profile']['empirical_even_power_multiplications'],
            'empirical_even_power_recurrence_divisions': data['shared_MGF_profile']['empirical_even_power_multiplications'],
            'polynomial_coefficient_evaluations': data['shared_MGF_profile']['polynomial_coefficient_evaluations'],
            'inverse_error_function_evaluations': data['shared_MGF_profile']['node_certificate']['inverse_error_function_evaluations'],
            'node_generation_Gaussian_erf_evaluations': 1,
            'node_generation_square_root_evaluations': 1,
            'recurrence_reference': 'higherdim_preparation_proof.tex, finite empirical exponential moments proposition. Actual dyadic powers are accumulated once by the factorial-divided recurrence. Each distinct rate is evaluated by Horner recurrence plus a certified positive remainder; no rate*N nodal exponentials are omitted from the count.'},
        'comparison_policy': 'Exactly the same paid moments, coefficients, correlated gap, norm and finite-domain information enters all competitors. Existing shared derivative-norm generation has its own retained certificates and is reused here.'}
    (ROOT / 'results/higherdim_preparation_costs.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
