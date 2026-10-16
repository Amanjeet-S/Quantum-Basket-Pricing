"""Independent one-sided high-address exclusion for the fixed literal graph.

This is not an algorithmic quantum lower bound. No producer, graph builder,
controller or earlier exclusion auditor is imported.
"""
import argparse, gzip, hashlib, json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True

def read(p):
    p = Path(p)
    return json.load(gzip.open(p, 'rt')) if p.suffix == '.gz' else json.loads(p.read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def digest(v):
    return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()

def root_upper(v):
    v = F(v); scale = 1 << 100
    integer = -(-(v*scale*scale).numerator//(v*scale*scale).denominator)
    root = math.isqrt(integer)
    if root*root < integer:
        root += 1
    return F(root, scale)

def failure(n, success):
    return sum(F(math.comb(n, k))*(1-success)**k*success**(n-k)
               for k in range((n+1)//2, n+1))

def fixed_record(row, mode):
    sem = row['graph_semantic_nodes']; omit = {sem['PG_exp'], sem['PA_exp']}
    profile = row['profiles'][mode]['centred']
    blocks = [b for b in profile['blocks'] if b['node'] not in omit]
    exp = [b for b in profile['blocks'] if b['node'] in omit]
    assert len(exp) == 2 and all(b['lookup_invocations'] == 2 for b in exp)
    assert all(b['kind'] == 'function' and b['function']['name'] == 'exp' and
               F(b['function']['lo']) == 0 and F(b['function']['span']) == F(1, 64) for b in exp)
    counts = profile['comparator_primitive_counts']
    t = 2*sum(b['logical_T_count'] for b in blocks)+counts.get('T', 0)+counts.get('Tdg', 0)
    assert t == profile['logical_T_count']-2*sum(b['logical_T_count'] for b in exp)
    graph = [n for n in row['graph_nodes'] if n['name'] not in omit]
    assert profile['integer_word_width'] == row['graph_configuration']['q']+row['graph_configuration']['headroom_bits']
    assert all(n['output_bits'] == profile['integer_word_width'] for n in row['graph_nodes'])
    return {'fixed_T_count': t, 'fixed_blocks_sha256': digest(blocks),
        'fixed_graph_sha256': digest(graph), 'comparator_sha256': digest(counts),
        'integer_word_width': profile['integer_word_width'],
        'fixed_block_count': len(blocks), 'private_lookups_per_exp_node': 2,
        'complete_lookups_per_exp_side': 4}

def main(args):
    repo, directory = Path(args.repo), Path(args.results)
    resource_path = directory/'asymmetric_raw_payoff_resources.json.gz'
    protocol_path = directory/'asymmetric_raw_payoff_protocol.json'
    data, protocol = read(resource_path), read(protocol_path)
    assert data['protocol_sha256'] == sha(protocol_path)
    assert protocol['q'] == 28 and protocol['payoff_exp_span'] == '1/64'
    assert protocol['inner_normal_bits'] == [8, 10, 12]
    # Exact source slices in the proof identify the fixed compiler whose
    # private clean lookup structure is being bounded.
    source_bindings = {n: sha(repo/n) for n in ('integral_value_oracle.py', 'architecture.py',
        'reversible.py', 'qpe_estimator.py', 'certified_scalar.py',
        'correction_information_priors.py', 'raw_payoff_oracle_tight.py')}
    assert source_bindings['raw_payoff_oracle_tight.py'] == protocol['retained_raw_builder_source_sha256']
    for n in ('integral_value_oracle.py', 'architecture.py', 'reversible.py',
              'qpe_estimator.py', 'certified_scalar.py', 'correction_information_priors.py'):
        assert source_bindings[n] == protocol['dependencies_sha256'][n]
    # Check that the actual 252 fresh emitted profiles share the entire
    # non-exp graph and every retained non-exp block at fixed model/normal.
    contexts = {}; fresh_profile_checks = 0
    for row in data['records']:
        for mode in ('exact7t', 'unitary8t'):
            r = fixed_record(row, mode)
            key = row['original_model_id'], row['normal_bits'], mode
            if key in contexts:
                assert contexts[key] == r
            else:
                contexts[key] = r
            fresh_profile_checks += 1
    assert len(contexts) == 24
    # Referenced diagonal graphs remain exact evidence, with no re-emission.
    old = read(repo/'results/raw_payoff_tight/raw_payoff_resources.json.gz')
    expanded = read(repo/'results/expanded_raw_payoff/expanded_raw_payoff_resources.json.gz')
    previous = {(r['original_model_id'], r['normal_bits'], r['payoff_exp_bits']): r
                for r in old['records']+expanded['records']}
    reused_profile_checks = 0
    for binding in protocol['reused_diagonal_bindings']:
        row = previous[binding['model_id'], binding['normal_bits'], binding['PG_exp_bits']]
        assert digest(row) == binding['pointwise_record_sha256']
        for mode in ('exact7t', 'unitary8t'):
            assert fixed_record(row, mode) == contexts[binding['model_id'], binding['normal_bits'], mode]
            reused_profile_checks += 1
    pi = F(355, 113); alpha = F(1, 3600); s0 = 8/pi**2
    failures = {n: failure(n, s0) for n in range(1, 26, 2)}
    assert all(v > alpha for n, v in failures.items() if n < 25)
    assert failures[25] <= alpha
    records = []
    for comparison in data['selected_comparisons']:
        mid = comparison['model_id']
        rows = [r for r in data['records'] if r['original_model_id'] == mid]
        reference = rows[0]; bud = reference['complete_deterministic_budget']
        prior = reference['coefficient_informed_prior']
        fields = ('coefficient_radius_upper', 'financial_transfer_upper', 'external_factor_upper', 'external_factor_radius')
        assert all(all(r['complete_deterministic_budget'][k] == bud[k] for k in fields) for r in rows)
        assert all(r['coefficient_informed_prior']['finite_ideal_mean_lower'] == prior['finite_ideal_mean_lower'] and
                   r['coefficient_informed_prior']['finite_ideal_mean_upper'] == prior['finite_ideal_mean_upper'] for r in rows)
        A = F(bud['external_factor_upper'])
        l, u = F(prior['finite_ideal_mean_lower']), F(prior['finite_ideal_mean_upper'])
        d0 = F(bud['coefficient_radius_upper'])+F(bud['financial_transfer_upper'])+F(bud['external_factor_radius'])
        remaining = F(1, 100)-d0; eps = remaining/A
        assert 0 <= l < u <= 1 and A > 0 and remaining > 0
        radius = A*(u-l)/2
        assert radius > remaining and u > 2*eps
        variance = u*(1-u) if u <= F(1, 2) else F(1, 4)
        root = root_upper(variance); M = 2
        tested = []
        while True:
            error = 2*pi*root/M+pi*pi/(M*M)+F(1, 10**60)
            tested.append({'M': M, 'amplitude_plus_postprocess_upper': str(error), 'admitted': error <= eps})
            if error <= eps:
                break
            M *= 2
        assert M in (256, 512)
        fixed = {mode: {str(b): contexts[mid, b, mode] for b in (8, 10, 12)}
                 for mode in ('exact7t', 'unitary8t')}
        minima = {mode: min(x['fixed_T_count'] for x in bynormal.values()) for mode, bynormal in fixed.items()}
        assert all(fixed['exact7t'][str(b)]['fixed_T_count'] >= fixed['unitary8t'][str(b)]['fixed_T_count'] for b in (8, 10, 12))
        incumbent = comparison['raw_T']
        assert comparison['selected_PG_exp_bits'] == comparison['selected_PA_exp_bits'] == 9
        floors = {}
        for mode, coeff in (('exact7t', 56), ('unitary8t', 32)):
            floors[mode] = {str(b): 25*(2*M-1)*(minima[mode]+coeff*((1 << b)+b)) for b in range(15, 23)}
            assert all(t > incumbent for t in floors[mode].values())
        records.append({'model_id': mid, 'ideal_probability_lower': str(l), 'ideal_probability_upper': str(u),
            'zero_eta_deterministic_error': str(d0), 'optimistic_probability_tolerance': str(eps),
            'two_sided_prior_radius': str(radius), 'two_sided_zero_shot_excluded': True,
            'one_sided_zero_shot_excluded': True, 'optimistic_variance_upper': str(variance),
            'dyadic_square_root_upper': str(root), 'supported_M_checks': tested, 'M0': M, 'r0': 25,
            'fixed_profile_bindings': fixed, 'fixed_T_minima': minima,
            'one_side_complete_T_floors': floors, 'PG9_PA9_incumbent_T': incumbent,
            'unitary_width15_T_floor': floors['unitary8t']['15'],
            'unitary_width15_gap_above_incumbent': floors['unitary8t']['15']-incumbent,
            'either_side_width15_through22_excluded': True})
    result = {'schema': 'qbp-asymmetric-one-side-high-width-exclusion-v1', 'status': 'passed',
        'auditor_sha256': sha(__file__), 'resources_sha256': sha(resource_path),
        'protocol_sha256': sha(protocol_path), 'literal_sources_sha256': source_bindings,
        'fresh_fixed_profile_invariance_checks': fresh_profile_checks,
        'referenced_diagonal_fixed_profile_invariance_checks': reused_profile_checks,
        'contexts': len(contexts), 'per_side_selector_coefficients': {'exact7t': 56, 'unitary8t': 32},
        'maximum_single_shot_certified_success': str(s0),
        'exact_majority_failures': {str(n): str(p) for n, p in failures.items()},
        'alpha': str(alpha), 'minimum_odd_repetitions': 25, 'records': records,
        'address_domain': {'each_side_feasible_widths': list(range(1, 23)),
            'excluded_if_either_side_in': list(range(15, 23)),
            'other_side_scope': 'Any feasible independently addressed width1 through22 in the unchanged literal graph; the other side and all exponential arithmetic are discarded in the lower bound.',
            'unresolved': 'Asymmetric widths13 and14, omitted widths below8, other normal widths, fixed-point precisions, optimised arithmetic sharing, other lookup architectures and estimators.'},
        'scope': 'Prospective exclusion within the retained unary-look-up compiler and canonical-QPE certification controller. No registered numerical outcome was changed or new quantum execution performed. No global precision optimum, generic query lower bound, physical advantage or theoretical priority claim.'}
    output = Path(args.output)
    if output.exists():
        raise FileExistsError('Preserve prior retained exclusion')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print('One-side high-width exclusion passed:', fresh_profile_checks+reused_profile_checks, 'fixed profiles;',
          [r['unitary_width15_T_floor'] for r in records])

if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--repo', required=True)
    p.add_argument('--results', required=True); p.add_argument('--output', required=True)
    main(p.parse_args())
