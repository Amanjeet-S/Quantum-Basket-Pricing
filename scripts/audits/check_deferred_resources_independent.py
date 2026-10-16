"""Own formula audit of fair deferred/cancellation-optimised costs.

The new producer and all financial confidence/resource frontends are excluded.
Existing independently derived QPE-plan and count audit helpers are reused.
No fresh full financial tape emission is claimed by this count transformation.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import ast
import gzip
import hashlib
import importlib.util
import json
import sys
sys.dont_write_bytecode = True


def read(path):
    path = Path(path)
    return json.load(gzip.open(path, 'rt')) if path.suffix == '.gz' else json.loads(path.read_text())


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def inverse(counts):
    names = {'T': 'Tdg', 'Tdg': 'T', 'S': 'Sdg', 'Sdg': 'S'}
    out = Counter()
    for op, value in counts.items():
        out[names.get(op, op)] += value
    return out


def controlled_reflection_source_guard(repo):
    """Check the frozen control-last conjunction source and its clean domain.

    This constructs/replays Boolean phase ladders independently and imports
    no gate emitter. Literal Clifford/T phases are already certified by the
    retained primitive audits; the check here concerns actual scratch use.
    """
    source = (repo / 'quantum_gates.py').read_text()
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == 'zero_phase')
    body = ast.get_source_segment(source, node)
    for text in ("controls=core+(() if control is None else (control,))",
                 "('CCX',controls[0],controls[1],work[0])",
                 "('CCX',work[j-2],controls[j],work[j-1])",
                 "('CZ',work[len(controls)-3],controls[-1])",
                 'gates+=body[::-1]'):
        assert text in body
    cases = 0
    for n in range(2, 9):
        core = tuple(range(n))
        control = n
        work = tuple(range(n + 1, 2 * n))
        controls = core + (control,)
        gates = [('X', w) for w in core]
        conjunction = [('CCX', controls[0], controls[1], work[0])]
        for j in range(2, len(controls) - 1):
            conjunction.append(('CCX', work[j - 2], controls[j], work[j - 1]))
        gates += conjunction + [('CZ', work[len(controls) - 3], controls[-1])] + list(reversed(conjunction))
        gates += [('X', w) for w in core]
        assert sum(g[0] == 'CCX' for g in gates) == 2 * (n - 1)
        for core_value in range(1 << n):
            for control_value in (0, 1):
                initial = core_value | (control_value << control)
                value, phase = initial, 1
                for op, *w in gates:
                    if op == 'X':
                        value ^= 1 << w[0]
                    elif op == 'CCX' and (value & (1 << w[0])) and (value & (1 << w[1])):
                        value ^= 1 << w[2]
                    elif op == 'CZ' and (value & (1 << w[0])) and (value & (1 << w[1])):
                        phase *= -1
                assert value == initial
                assert phase == (-1 if control_value and core_value == 0 else 1)
                cases += 1
    return {'status': 'passed', 'clean_domain_basis_cases': cases,
            'checked_core_widths': list(range(2, 9)),
            'zero_phase_source_sha256': hashlib.sha256(body.encode()).hexdigest(),
            'actual_controlled_tape_contract_is_restricted_to_zero_work': True,
            'off_domain_controlled_uncontrolled_equality_claim': False}


def own_clean_extension(dirty, clean, profile):
    reps = dirty['repetitions']
    suffix = Counter(profile['terminal_inverse_body_counts'])
    per = Counter(dirty['per_repetition_metrics'])
    per['T'] += suffix['T'] + suffix['Tdg']
    per['expanded_primitive_operations'] += sum(suffix.values())
    per['explicit_T_depth_serial_upper'] += profile['terminal_inverse_body_explicit_T_depth_serial_upper']
    per['conservative_5layer_T_depth_serial_upper'] += profile['terminal_inverse_body_conservative_5layer_T_depth_serial_upper']
    for op, value in suffix.items():
        per['gate_' + op] += value
    assert Counter(clean['per_repetition_metrics']) == per
    assert Counter(clean['metrics']) == Counter({k: reps * v for k, v in per.items()})
    assert clean['peak_logical_qubits'] == dirty['peak_logical_qubits']
    assert not clean['work_dirty_before_terminal_reset']
    assert clean['terminal_cleanup_T_per_repetition'] == suffix['T'] + suffix['Tdg']


def run(repo, directory, output, proof):
    repo, directory = Path(repo), Path(directory)
    resources_path = directory / 'deferred_cleanup_resources.json.gz'
    protocol_path = directory / 'deferred_cleanup_protocol.json'
    data = read(resources_path)
    protocol = read(protocol_path)
    old = read(repo / 'results/common_shift/common_shift_resources.json.gz')
    old_half = read(repo / 'results/common_shift/common_shift_poolwise_fallback.json.gz')
    assert data['registered_configurations'] == 416
    assert len(data['records']) == 808 and len(data['rejections']) == 12
    expected_admitted = {(i, architecture) for i, r in enumerate(old['records'])
                         if r['strategies'] for architecture in ('exact7t', 'unitary8t')}
    observed_admitted = {(r['record_index'], r['architecture']) for r in data['records']}
    assert expected_admitted == observed_admitted and len(observed_admitted) == 808
    expected_rejected = {i for i, r in enumerate(old['records']) if not r['strategies']}
    assert {r['record_index'] for r in data['rejections']} == expected_rejected
    for rejected in data['rejections']:
        prior = old['records'][rejected['record_index']]
        assert rejected['record_sha256'] == digest(prior)
        assert rejected['retained_rejection'] == prior['rejection']
    assert data['checked_proof_sha256'] == protocol['checked_proof_sha256']
    assert hashlib.sha256(Path(proof).read_bytes()).hexdigest() == protocol['checked_proof_sha256']
    assert hashlib.sha256((repo / 'deferred_cleanup_comparison.py').read_bytes()).hexdigest() == protocol['producer_source_sha256']
    for name, expected in protocol['inputs_sha256'].items():
        assert hashlib.sha256((repo / name).read_bytes()).hexdigest() == expected
    helper_path = repo / 'check_expanded_precision_resources.py'
    spec = importlib.util.spec_from_file_location('prior_independent_qpe_audit', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    reflection_guard = controlled_reflection_source_guard(repo)
    rotations = read(repo / 'results/qft_rotation_tapes.json')['rotations']
    bounds_raw = read(repo / 'results/qft_rotation_certificates.json')
    bounds = {c['key']: F(c['operator_norm_upper_rational']) for c in bounds_raw['certificates']}
    source_cache, checked_profiles = {}, set()
    half_index = {(r['record_index'], r['architecture']): r for r in old_half['records']}
    results, plan_checks, charge_checks = [], 0, 0
    for record in data['records']:
        index, architecture = record['record_index'], record['architecture']
        prior = old['records'][index]
        assert record['record_sha256'] == digest(prior)
        for field in ('model_id', 'method', 'source_path', 'source_row_index',
                      'source_row_sha256', 'profiles_sha256'):
            assert record[field] == prior[field]
        expected_configuration = {k: prior[k] for k in
            ('normal_bits', 'PG_exp_bits', 'PA_exp_bits', 'PG_bits', 'PA_bits') if k in prior}
        assert record['configuration'] == expected_configuration
        assert record['probability_certificate_sha256'] == digest(prior['probability_certificate'])
        assert record['complete_deterministic_budget_sha256'] == digest(prior['complete_deterministic_budget'])
        assert record['coefficient_informed_prior_sha256'] == digest(prior['coefficient_informed_prior'])
        source = prior['source_path']
        if source not in source_cache:
            source_cache[source] = read(repo / 'results' / source)
        material = source_cache[source]
        base = material['rows' if prior['method'] == 'shared_threshold' else 'records'][prior['source_row_index']]
        assert digest(base) == prior['source_row_sha256']
        assert digest(base['profiles']) == prior['profiles_sha256']
        clean = base['profiles'][architecture]['centred']
        profile = data['profiles'][record['deferred_profile_sha256']]
        assert digest(profile) == record['deferred_profile_sha256']
        assert profile['proof_sha256'] == protocol['checked_proof_sha256']
        assert profile['source_clean_profile_sha256'] == digest(clean)
        assert profile['source_blocks_sha256'] == digest(clean['blocks'])
        assert profile['source_emitter_sha256'] == clean['source_sha256']
        fw = Counter(clean['primitive_counts_forward_body'])
        predicate = Counter(clean.get('predicate_primitive_counts', clean.get('comparator_primitive_counts')))
        H = Counter({'H': clean['uniform_H']})
        assert Counter(clean['primitive_counts']) == fw + fw + predicate + H
        assert Counter(profile['primitive_counts']) == fw + predicate + H
        assert Counter(profile['terminal_inverse_body_counts']) == inverse(fw)
        assert profile['logical_T_count'] == fw['T'] + fw['Tdg'] + predicate['T'] + predicate['Tdg']
        assert profile['terminal_inverse_body_T'] == fw['T'] + fw['Tdg']
        for field in ('reflection_core_qubits', 'uniform_H', 'logical_qubits_upper',
                      'explicit_T_depth_serial_upper', 'conservative_5layer_T_depth_serial_upper'):
            assert profile[field] == clean[field]
        assert profile['terminal_inverse_body_explicit_T_depth_serial_upper'] == sum(
            b['explicit_T_depth_serial_upper'] for b in clean['blocks'])
        assert profile['terminal_inverse_body_conservative_5layer_T_depth_serial_upper'] == sum(
            b['conservative_5layer_T_depth_serial_upper'] for b in clean['blocks'])
        assert profile['work_dirty_after_preparation']
        assert not profile['private_work_in_zero_reflection_core']
        assert profile['actual_reflection_tape_borrows_clean_work']
        assert not profile['fresh_whole_preparation_emission_claim']
        checked_profiles.add(record['deferred_profile_sha256'])
        new_result = record['deferred_B_full_result']
        flag = new_result['flags'][0]
        old_result = next(s['result'] for s in prior['strategies'] if s['architecture'] == architecture)
        old_flag = old_result['flags'][0]
        for field in ('plan', 'fourier_profile', 'probability_tolerance', 'prior_upper'):
            assert digest(flag[field]) == digest(old_flag[field])
        for field in ('deterministic_error_upper', 'statistical_price_error_upper',
                      'total_price_error_upper', 'failure_probability_upper',
                      'coefficient_informed_probability_interval'):
            assert digest(new_result[field]) == digest(old_result[field])
        plan, qft = flag['plan'], flag['fourier_profile']
        helper.ALPHA = F(1, 3600)
        helper.own_plan(plan, F(flag['probability_tolerance']), F(flag['prior_upper']), qft, rotations, bounds)
        helper.own_charge(flag['complete_charge'], profile, qft, plan)
        own_clean_extension(flag['complete_charge'], record['cancellation_optimised_clean_full_charge'], profile)
        hp = record['half_plan']
        helper.ALPHA = F(1, 7200)
        helper.own_plan(hp, F(flag['probability_tolerance']), F(flag['prior_upper']), qft, rotations, bounds)
        helper.own_charge(record['deferred_B_half_charge'], profile, qft, hp)
        own_clean_extension(record['deferred_B_half_charge'], record['cancellation_optimised_clean_half_charge'], profile)
        assert (plan['repetitions'], hp['repetitions']) == (25, 27)
        assert hp['M'] == plan['M']
        f = fw['T'] + fw['Tdg']
        naive = old_result['T_count']
        dirty_T = flag['complete_charge']['metrics']['T']
        optimised_T = record['cancellation_optimised_clean_full_charge']['metrics']['T']
        assert naive - dirty_T == plan['repetitions'] * (2 * plan['M'] - 1) * f
        assert optimised_T - dirty_T == plan['repetitions'] * f
        assert record['previous_naive_clean_T_half_alpha_fallback'] == half_index[index, architecture]['T_half_alpha_fallback']
        for name, full, half in (
            ('deferred_B', flag['complete_charge'], record['deferred_B_half_charge']),
            ('cancellation_optimised_clean', record['cancellation_optimised_clean_full_charge'], record['cancellation_optimised_clean_half_charge'])
        ):
            canonical, fallback = full['metrics']['T'], half['metrics']['T']
            assert F(fallback, canonical) == F(27, 25)
            assert record[name] == {'T_full_alpha': canonical, 'T_half_alpha_fallback': fallback,
                                    'pilot_T_budget': fallback // 10, 'T_worst_history_upper': fallback + fallback // 10}
        plan_checks += 2
        charge_checks += 4
        results.append(record)
    # Poolwise selection is independently performed, with its declared tie rule.
    def tie(r):
        c = r['configuration']
        return (c.get('normal_bits', 0), c.get('PG_exp_bits', c.get('PG_bits', 0)),
                c.get('PA_exp_bits', c.get('PA_bits', 0)), r['architecture'], r['record_index'])
    comparisons = []
    for case in data['selected_comparisons']:
        model = case['model_id']
        item = {'model_id': model, 'variants': {}}
        for variant in ('deferred_B', 'cancellation_optimised_clean'):
            minima = {}
            for method, value in case['variants'][variant]['methods'].items():
                group = [r for r in results if r['model_id'] == model and r['method'] == method]
                best = min(group, key=lambda r: (r[variant]['T_full_alpha'], tie(r)))
                worst = min(group, key=lambda r: (r[variant]['T_worst_history_upper'], tie(r)))
                assert best['record_index'] == value['canonical_record_index']
                assert best['architecture'] == value['canonical_architecture']
                assert best['configuration'] == value['canonical_configuration']
                assert best[variant]['T_full_alpha'] == value['T_full_alpha']
                assert worst['record_index'] == value['worst_history_record_index']
                assert worst[variant]['T_worst_history_upper'] == value['T_worst_history_upper']
                minima[method] = best[variant]['T_full_alpha']
            raw, numerical = ('raw_payoff', 'numerical_conditional_price')
            rival = min((raw, numerical), key=lambda m: (minima[m], m))
            assert rival == case['variants'][variant]['strongest_rival']
            saving = F(minima[rival] - minima['shared_threshold'], minima[rival])
            assert str(saving) == case['variants'][variant]['saving_against_strongest_canonical_exact']
            item['variants'][variant] = {'canonical_minima': minima, 'strongest_rival': rival,
                                         'shared_threshold_saving_exact': str(saving)}
        comparisons.append(item)
    out = {'schema': 'qbp-deferred-resources-independent-review-v1', 'status': 'passed',
           'checked_architecture_profiles': len(checked_profiles), 'checked_qpe_plans': plan_checks,
           'checked_complete_variant_charges': charge_checks, 'retained_rejections': len(data['rejections']),
           'registered_settings': 416, 'selected_comparisons': comparisons,
           'actual_controlled_reflection_source_guard': reflection_guard,
           'resources_sha256': hashlib.sha256(resources_path.read_bytes()).hexdigest(),
           'protocol_sha256': hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
           'independent_checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'proof_sha256': protocol['checked_proof_sha256'],
           'producer_sha256_read_without_import': protocol['producer_source_sha256'],
           'producer_imported': False, 'confidence_or_charging_frontend_imported': False,
           'reused_independent_helper': str(helper_path),
           'reused_independent_helper_sha256': hashlib.sha256(helper_path.read_bytes()).hexdigest(),
           'fresh_whole_financial_tape_emission_claim': False,
           'scope': 'Own integer primitive/reflection/QFT count formulas, exact profile decomposition, unchanged financial priors and canonical/half-alpha plans, all pool minima. Actual financial primitive correctness reuses retained earlier independent audits.'}
    Path(output).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'selected_comparisons'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', required=True)
    parser.add_argument('--results', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--proof', required=True)
    args = parser.parse_args()
    run(args.repo, args.results, args.output, args.proof)
