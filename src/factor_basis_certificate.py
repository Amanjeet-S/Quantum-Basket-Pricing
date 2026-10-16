"""Exact finite-pool factor-basis certificates, with no financial price input.

This module precompiles a rational rotation into classical GBM coefficients.
It does not rotate the existing finite tensor grid or add a quantum unitary.
Register before producing the full deterministic screen. Earlier loading-only
proxy diagnostics motivated the direction, so registration is not blind.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import argparse, copy, hashlib, json, sys, time

ROOT = Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from flint import ctx
import tensor_transfer as tt
from boundary_grid_certificate import parameters, serial, upper
from boundary_derivative_certificate import payoff_norms

DEFAULT_RESULTS = Path(__file__).resolve().parent / 'results'
DEPENDENCIES = ('tensor_transfer.py', 'higherdim_preparation.py',
    'boundary_grid_certificate.py', 'boundary_derivative_certificate.py',
    'integral_value_oracle.py', 'continuous_quadrature.py', 'taylor_pricing.py')
INPUTS = ('results/four_survivor_quantum_protocol.json',
    'results/boundary_derivative_norms.json', 'results/gaussian_derivative_norms.json')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
        separators=(',', ':'), default=str).encode()).hexdigest()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=str) + '\n')

def candidates(model):
    p = model['outer_dimension']
    if type(p) is not int or p < 1:
        raise ValueError('Positive factor dimension required')
    out = [{'id': 'identity', 'dimension': p, 'plane': None,
            'cosine': '1', 'sine': '0'}]
    for a in range(p):
        for b in range(a + 1, p):
            for sign, word in ((-1, 'minus'), (1, 'plus')):
                out.append({'id': f'plane_{a}_{b}_{word}', 'dimension': p,
                    'plane': [a, b], 'cosine': '3/5', 'sine': str(sign * F(4, 5))})
    return out

def rotation_matrix(candidate):
    p = candidate['dimension']
    out = [[F(int(a == b)) for b in range(p)] for a in range(p)]
    if candidate['plane'] is not None:
        a, b = candidate['plane']
        c, s = F(candidate['cosine']), F(candidate['sine'])
        if not 0 <= a < b < p or c * c + s * s != 1:
            raise ValueError('Invalid exact rational plane rotation')
        out[a][a] = out[b][b] = c
        out[a][b] = -s
        out[b][a] = s
    elif candidate['id'] != 'identity' or F(candidate['cosine']) != 1 or F(candidate['sine']):
        raise ValueError('Invalid identity descriptor')
    return tuple(map(tuple, out))

def multiply_rows(rows, matrix):
    p = len(matrix)
    if any(len(row) != p for row in rows):
        raise ValueError('Loading dimension differs from rotation')
    return tuple(tuple(sum(F(row[k]) * matrix[k][j] for k in range(p))
        for j in range(p)) for row in rows)

def rotate_model(model, candidate):
    if candidate not in candidates(model):
        raise ValueError('Candidate is outside the declared pool')
    result = copy.deepcopy(model)
    rows = multiply_rows(model['loadings'], rotation_matrix(candidate))
    result['loadings'] = [list(map(str, row)) for row in rows]
    if candidate['id'] != 'identity':
        result['model_id'] = model['model_id'] + '__' + candidate['id']
    # The financial covariance, drifts and all other fields are deliberately
    # identical. This exact validator rejects any drift or dimension mismatch.
    parameters(result, False)
    return result

def gram(rows):
    return tuple(tuple(sum(a * b for a, b in zip(x, y)) for y in rows)
        for x in rows)

def string_matrix(matrix):
    return [list(map(str, row)) for row in matrix]

def invariant_certificate(original, rotated, candidate):
    w, d, s, mg, rg, lam = parameters(original, False)
    ww, dd, ss, mmg, rrg, llam = parameters(rotated, False)
    O = rotation_matrix(candidate)
    rows = tuple(tuple(map(F, row)) for row in original['loadings'])
    rrows = tuple(tuple(map(F, row)) for row in rotated['loadings'])
    p = len(rg)
    columns = tuple(tuple(O[k][j] for k in range(p)) for j in range(p))
    identity = tuple(tuple(F(int(a == b)) for b in range(p)) for a in range(p))
    if gram(O) != identity or gram(columns) != identity:
        raise ArithmeticError('Both exact orthogonality identities failed')
    outer = gram(rows)
    if outer != gram(rrows):
        raise ArithmeticError('Asset covariance changed')
    full = tuple(tuple(v + lam * lam for v in row) for row in outer)
    if full != tuple(tuple(map(F, row)) for row in original['covariance']):
        raise ArithmeticError('Original stated covariance is inconsistent')
    if full != tuple(tuple(map(F, row)) for row in rotated['covariance']):
        raise ArithmeticError('Rotated stated covariance changed')
    if (w, d, mg, lam) != (ww, dd, mmg, llam):
        raise ArithmeticError('Weights, centred drifts or common factor changed')
    centred_gram = gram(s)
    geometric_variance = sum(v * v for v in rg)
    cross = tuple(sum(a * b for a, b in zip(rg, si)) for si in s)
    if centred_gram != gram(ss) or geometric_variance != sum(v*v for v in rrg):
        raise ArithmeticError('Centred or geometric covariance changed')
    if cross != tuple(sum(a*b for a,b in zip(rrg, si)) for si in ss):
        raise ArithmeticError('Geometric-centred covariance changed')
    if tuple(rrg) != multiply_rows((rg,), O)[0] or ss != multiply_rows(s, O):
        raise ArithmeticError('Weighted geometric/centred loadings did not rotate')
    dispersion = sum(wi * sum(v*v for v in si) for wi, si in zip(w, s))
    if dispersion != F(original['weighted_residual_dispersion']):
        raise ArithmeticError('Stated residual dispersion is inconsistent')
    if geometric_variance != F(original['geometric_outer_variance']):
        raise ArithmeticError('Stated geometric variance is inconsistent')
    return {'status': 'passed exact rational identities',
        'rotation_matrix': string_matrix(O), 'both_orthogonality_identities': True,
        'asset_outer_gram_sha256': digest(string_matrix(outer)),
        'full_covariance_sha256': digest(string_matrix(full)),
        'centred_loading_gram_sha256': digest(string_matrix(centred_gram)),
        'geometric_centred_cross_products': list(map(str, cross)),
        'geometric_outer_variance': str(geometric_variance),
        'weighted_residual_dispersion': str(dispersion),
        'effective_geometric_variance': str(geometric_variance + lam * lam),
        'unchanged_log_drifts': original['log_drifts'],
        'unchanged_geometric_mean_log_drift': str(mg),
        'all_continuous_U_moments_invariant': True,
        'continuous_geometric_coefficients_invariant': True,
        'proof_basis': 'Exact Gram and cross products determine every Gaussian exponential monomial; the supporting proof covers every finite nonnegative integer moment order.',
        'finite_grid_invariance_claimed': False}

def direct_first_transfer(model, transport, norms):
    loadings = tuple(tuple(map(F, row)) for row in model['loadings'])
    _, _, _, _, rg, _ = parameters(model, False)
    p = len(rg)
    M1, M2 = F(norms['M1_upper']), F(norms['M2_upper'])
    b1, b2 = upper(transport.rho1(0)), upper(transport.rho2mass)
    other = upper(transport.mass ** (p - 1))
    linear = quadratic = F(0)
    geometry1 = geometry2 = F(0)
    rows = []
    for j in range(p):
        A, r = max(abs(row[j]) for row in loadings), rg[j]
        B1, B2 = A + abs(r), A * A + r * r
        D1, D2 = M1 * B1, (M1 + M2) * B2
        ell, quad = other * b1 * D1, other * b2 * D2
        linear += ell
        quadratic += quad
        geometry1 += B1
        geometry2 += B2
        rows.append({'coordinate': j, 'full_asset_loading_max': str(A),
            'geometric_loading': str(r), 'D1_upper': str(D1), 'D2_upper': str(D2),
            'linear_error_upper': str(ell), 'quadratic_error_upper': str(quad)})
    tail = upper(1 - transport.mass ** p)
    return {'model_id': model['model_id'], 'payoff': 'capped', 'cap': '1', 'order': 1,
        'linear_transport_upper': str(linear), 'quadratic_transport_upper': str(quadratic),
        'tail_upper': str(tail), 'financial_transfer_upper': str(linear + quadratic + tail),
        'mass_interval': serial(transport.mass ** p), 'coordinates': rows,
        'rho1_mass_upper': str(b1), 'rho2_mass_upper': str(b2),
        'other_coordinate_mass_upper': str(other),
        'geometric_B1': str(geometry1), 'geometric_B2': str(geometry2),
        'scope': 'Full bounded capped first-residual derivative transfer, with exact rational loading geometry and outward raw unnormalised Gaussian cell masses. No financial price input.'}

def source_hashes(repo):
    return {name: sha(repo/name) for name in DEPENDENCIES}

def input_hashes(repo):
    return {name: sha(repo/name) for name in INPUTS}

def register(repo, results):
    path = results/'factor_basis_protocol.json'
    if path.exists():
        raise FileExistsError('Preserve existing registration')
    models = json.loads((repo/INPUTS[0]).read_text())['models']
    protocol = {'schema': 'qbp-exact-factor-basis-protocol-v1',
        'registered_before_fresh_full_certificate_execution': True,
        'registration_time_utc': datetime.now(timezone.utc).isoformat(),
        'research_record_date': '2026-10-10',
        'earlier_proxy_diagnostics_known': True,
        'registration_scope': 'Freeze the full exact covariance-preserving rational candidate pool and deterministic selection before fresh full transfer certificates. Earlier loading-only proxy diagnostics informed this direction; this is not blind registration.',
        'models': models, 'candidate_pools': {m['model_id']: candidates(m) for m in models},
        'candidate_count': sum(len(candidates(m)) for m in models),
        'rotation_convention': 'Rprime=R O, with O[a,a]=O[b,b]=3/5, O[a,b]=-s, O[b,a]=s and s in {-4/5,4/5}',
        'cap': '1', 'financial_law': {'L': '6', 'b_y': 10, 'q': 28, 'precision_bits': 512},
        'selection': 'Exact minimum direct R1 financial_transfer_upper, then lexical candidate ID. Identity retained. No quantum resource or price answer enters this screen.',
        'resource_pool': 'For each case, identity and the selected direct-transfer rotation; duplicate identity omitted. Full quantum resources determine deployment only within this declared pair.',
        'quantum_comparison_scope': 'No cost optimisation over the complete analytic rotation pool and no global basis optimality claim.',
        'classical_information_policy': 'The same exact rotated financial representation, deterministic certificates and finite-grid preparation are available to every classical comparator. Continuous price identity is proved; earlier identity classical price certificates remain valid for the common continuous target.',
        'data_origin': 'Exact rational synthetic GBM parameters. No LSEG observations or derived market data.',
        'financial_price_answer_excluded': True,
        'source_sha256': sha(__file__), 'dependency_source_sha256': source_hashes(repo),
        'input_sha256': input_hashes(repo)}
    write_json(path, protocol)
    print('Registered', protocol['candidate_count'], 'full analytic candidates', flush=True)

def execute(repo, results):
    protocol_path = results/'factor_basis_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    expected = (protocol['source_sha256'], protocol['dependency_source_sha256'], protocol['input_sha256'])
    if expected != (sha(__file__), source_hashes(repo), input_hashes(repo)):
        raise RuntimeError('Frozen source or input changed')
    for name in ('factor_basis_certificates.json', 'factor_basis_summary.json'):
        if (results/name).exists():
            raise FileExistsError('Preserve existing full certificate outcomes')
    start = time.perf_counter()
    ctx.prec = 512
    transport = tt.TransportMGF(6, 10, 28)
    records, selected, summary = [], [], []
    for model in protocol['models']:
        norms = payoff_norms(model['lambda'], F(1))
        family = []
        for candidate in protocol['candidate_pools'][model['model_id']]:
            rotated = rotate_model(model, candidate)
            invariants = invariant_certificate(model, rotated, candidate)
            transfer = direct_first_transfer(rotated, transport, norms)
            row = {'original_model_id': model['model_id'], 'candidate_id': candidate['id'],
                'candidate': candidate, 'rotated_model': rotated, 'invariants': invariants,
                'direct_global_transfer': transfer,
                'screen_status': 'passed geometry and full direct transfer',
                'resource_selection_used_price_answer': False}
            family.append(row)
            records.append(row)
        winner = min(family, key=lambda x: (F(x['direct_global_transfer']['financial_transfer_upper']), x['candidate_id']))
        identity = next(row for row in family if row['candidate_id'] == 'identity')
        selected.append(copy.deepcopy(winner))
        b = identity['direct_global_transfer']
        q = winner['direct_global_transfer']
        summary.append({'model_id': model['model_id'], 'analytic_candidates': len(family),
            'selected_candidate': winner['candidate_id'],
            'identity_transfer': b['financial_transfer_upper'], 'selected_transfer': q['financial_transfer_upper'],
            'transfer_reduction_fraction': str(1-F(q['financial_transfer_upper'])/F(b['financial_transfer_upper'])),
            'B1_reduction_fraction': str(1-F(q['geometric_B1'])/F(b['geometric_B1'])),
            'B2_reduction_fraction': str(1-F(q['geometric_B2'])/F(b['geometric_B2'])),
            'identity_retained_for_full_resource_comparison': True,
            'screen_choice_is_global_resource_optimum': False})
        print(model['model_id'], winner['candidate_id'], 'transfer', float(F(q['financial_transfer_upper'])),
            'reduction', float(F(summary[-1]['transfer_reduction_fraction'])), flush=True)
    if expected != (sha(__file__), source_hashes(repo), input_hashes(repo)):
        raise RuntimeError('Frozen source or input changed during execution')
    result = {'schema': 'qbp-exact-factor-basis-certificates-v1', 'status': 'completed with every analytic candidate retained',
        'protocol_sha256': sha(protocol_path), 'records': records, 'selected_models': selected,
        'source_sha256': sha(__file__), 'dependency_source_sha256': protocol['dependency_source_sha256'],
        'input_sha256': protocol['input_sha256'], 'transport_profile': transport.profile(),
        'actual_financial_node_integers': [c['actual_node_integer'] for c in transport.cells],
        'elapsed_seconds_diagnostic': time.perf_counter()-start,
        'financial_price_evaluations': 0,
        'scope': 'Continuous financial law preserved exactly, new finite product representation certified separately. Full quantum costs and classical prices are not produced by this analytic screen.'}
    write_json(results/'factor_basis_certificates.json', result)
    write_json(results/'factor_basis_summary.json', {'schema': 'qbp-exact-factor-basis-summary-v1',
        'status': 'passed', 'certificate_sha256': sha(results/'factor_basis_certificates.json'),
        'records': summary, 'analytic_candidates_retained': len(records)})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=ROOT)
    parser.add_argument('--results', type=Path, default=DEFAULT_RESULTS)
    parser.add_argument('action', choices=('register', 'run'))
    args = parser.parse_args()
    (register if args.action == 'register' else execute)(args.repo, args.results)

if __name__ == '__main__':
    main()
