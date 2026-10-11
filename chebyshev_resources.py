"""Bounded Chebyshev sampling on retained, complete financial preparations.

This is a source-matched finite prototype, not an implementation of every
heuristic in published ChebAE and not a new financial model.  No derivative
price or market observation enters the diagnostics.  Register before running.
"""
from pathlib import Path
from collections import Counter
from dataclasses import asdict
from fractions import Fraction as F
from datetime import datetime, timezone
import argparse, hashlib, importlib, json, random, statistics, sys

HERE = Path(__file__).resolve().parent
DEGREE_CAP = 2047
LOOK_CAP = 12
BATCH_SHOTS = 100
SEEDS = tuple(range(20))
ALPHA = F(1, 3600)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str) + '\n')


def clean(counter):
    return {k: int(v) for k, v in counter.items() if v}


def retained_cases(repo):
    data = json.loads((repo / 'results/four_survivor_quantum_resources.json').read_text())
    cases = []
    for selected in data['selected_prices']:
        choice = selected['methods']['positive_r1']
        row = data['records'][choice['candidate_index']]
        selected_strategy = choice['selected']
        assert selected_strategy['architecture'] == 'unitary8t'
        assert selected_strategy['prior_information'] == 'coefficient_informed_direct_minimum'
        assert selected_strategy['strategy'] == 'positive'
        result = selected_strategy['result']
        assert len(result['flags']) == 1
        flag = result['flags'][0]
        profile = row['profiles']['unitary8t']['centred']
        assert digest(profile) == flag['financial_profile_sha256']
        prior = result['correction_prior']
        cases.append({'model_id': selected['model_id'], 'candidate_index': row['candidate_index'],
            'financial_preparation_profile': profile, 'financial_preparation_profile_sha256': digest(profile),
            'probability_lower': prior['probability_lower'], 'probability_upper': prior['probability_upper'],
            'probability_tolerance': flag['probability_tolerance'],
            'monetary_signal_factor_upper': result['monetary_signal_factor_upper'],
            'deterministic_price_error_upper': result['deterministic_error_upper'],
            'previous_complete_QPE_T': result['T_count'],
            'prior_scope': 'Checked coefficient-informed direct-minimum bounds, without a classical price answer.'})
    assert len(cases) == 4
    return cases


def load_dependencies(repo):
    # Financial files are read, never modified or imported as producer APIs.
    sys.path.insert(0, str(repo))
    from reflections import zero_reflection
    from quantum_gates import expand_toffoli
    from integral_pricing_budget import compiled_extended_plan, charge_full_oracle
    return zero_reflection, expand_toffoli, compiled_extended_plan, charge_full_oracle


def emitted_reflection(core, zero_reflection, expand_toffoli):
    ideal = zero_reflection(core)
    tape = expand_toffoli(ideal)
    counts = Counter(g[0] for g in tape)
    ccx = sum(g[0] == 'CCX' for g in ideal)
    return {'core_qubits': core, 'logical_qubits': 2 * core - 2,
        'ideal_gate_list': ideal, 'chronological_primitive_gate_list': tape,
        'chronological_tape_sha256': digest(tape), 'primitive_counts': clean(counts),
        'T_count': counts['T'] + counts['Tdg'], 'Toffoli_count': ccx,
        'explicit_T_depth_serial_upper': 4 * ccx,
        'conservative_5layer_T_depth_serial_upper': 5 * ccx,
        'scratch_clean_on_invariant_inputs': True}


def build_profile(preparation, reflection, degree):
    if isinstance(degree, bool) or not isinstance(degree, int) or not 1 <= degree <= DEGREE_CAP:
        raise ValueError('Sampler degree must be an integer from 1 through 2047')
    core = preparation['reflection_core_qubits']
    assert core == reflection['core_qubits']
    forward = Counter(preparation['primitive_counts'])
    adjoint = Counter()
    swap = {'T': 'Tdg', 'Tdg': 'T', 'S': 'Sdg', 'Sdg': 'S'}
    for opcode, count in forward.items():
        adjoint[swap.get(opcode, opcode)] += count
    k = degree // 2
    if degree % 2:
        fwd_calls, adj_calls, zero_calls, good_calls = k + 1, k, k, k
        observable = 'Good flag equals one'
        symbolic = ['A', {'repeat': k, 'chronological_body': ['S_good', 'A_adjoint', 'S_zero', 'A']}]
    else:
        fwd_calls, adj_calls, zero_calls, good_calls = k, k, k - 1, k
        observable = 'All invariant core input wires equal zero'
        symbolic = ['A', {'repeat': k - 1, 'chronological_body': ['S_good', 'A_adjoint', 'S_zero', 'A']}, 'S_good', 'A_adjoint']
    counts = Counter({g: fwd_calls * n for g, n in forward.items()})
    counts.update({g: adj_calls * n for g, n in adjoint.items()})
    counts.update({g: zero_calls * n for g, n in reflection['primitive_counts'].items()})
    counts['Z'] += good_calls
    peak = max(preparation['logical_qubits_upper'], reflection['logical_qubits'])
    T = counts['T'] + counts['Tdg']
    explicit = degree * preparation['explicit_T_depth_serial_upper'] + zero_calls * reflection['explicit_T_depth_serial_upper']
    five = degree * preparation['conservative_5layer_T_depth_serial_upper'] + zero_calls * reflection['conservative_5layer_T_depth_serial_upper']
    metrics = {'T': T, 'measurements': peak, 'reset_X': peak, 'feedback_CZ': 0,
        'expanded_primitive_operations': sum(counts.values()) + 2 * peak,
        'uniform_H': degree * preparation['uniform_H'],
        'explicit_T_depth_serial_upper': explicit,
        'conservative_5layer_T_depth_serial_upper': five}
    metrics.update({f'gate_{op}': n for op, n in counts.items() if n})
    return {'degree': degree, 'parity': 'odd' if degree % 2 else 'even',
        'chronological_macro_template': symbolic, 'observable': observable,
        'signal_identity': 'T_d(sqrt(p))^2 = cos(d arccos(sqrt(p)))^2',
        'theta_signal_identity': 'sin(d theta)^2' if degree % 2 else 'cos(d theta)^2',
        'A_forward_calls': fwd_calls, 'A_adjoint_calls': adj_calls,
        'zero_reflection_calls': zero_calls, 'good_Z_calls': good_calls,
        'primitive_counts': clean(counts), 'per_shot_metrics': metrics,
        'peak_logical_qubits': peak,
        'full_observation_reset_scope': 'Every allocated wire is measured and conditionally X-reset after each fresh shot; the zero predicate is classical postprocessing of those observations.',
        'global_phase_convention': 'Q=-A S_zero A_adjoint S_good. Emitted uncontrolled Q omits its physically immaterial global minus; all probabilities and counts are invariant under that omission.',
        'financial_preparation_sha256': digest(preparation),
        'full_adjoint_opcode_rule': swap,
        'financial_cleanup_scope': 'The retained preparation already includes freshly emitted Boolean-body cleanup; the adjoint here reverses the full preparation and swaps phase opcodes.'}


def fallback(case, compiled_extended_plan, charge_full_oracle):
    tolerance = F(case['probability_tolerance'])
    prior = F(case['probability_upper'])
    setup, fourier = compiled_extended_plan(tolerance, ALPHA / 2, prior)
    charge = charge_full_oracle(setup, case['financial_preparation_profile'], fourier)
    assert not setup.zero_shots
    return {'plan': asdict(setup), 'fourier_profile': fourier, 'complete_charge': charge,
        'T_count': charge['metrics']['T'], 'alpha': str(ALPHA / 2),
        'probability_tolerance': str(tolerance), 'prior_upper': str(prior),
        'financial_preparation_sha256': case['financial_preparation_profile_sha256'],
        'scope': 'Fresh canonical QPE fallback with its own half of the failure allowance, compiled Fourier norm and full initial preparation, controlled iterations, measurements and reinitialisations.'}


def source_manifest(repo):
    # Capture every imported local implementation involved in fallback synthesis.
    sources = {'chebyshev_resources.py': sha(__file__),
        'chebyshev_controller.py': sha(HERE / 'chebyshev_controller.py')}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            p = Path(path)
            if p.parent == repo and p.suffix == '.py':
                sources[p.name] = sha(p)
    return dict(sorted(sources.items()))


def register(repo, output):
    output.mkdir(parents=True, exist_ok=True)
    target = output / 'chebyshev_protocol.json'
    if target.exists():
        raise FileExistsError('Protocol already exists; preserve the executed record')
    dependencies = load_dependencies(repo)
    load_controller()
    protocol = {'schema': 'qbp-bounded-chebyshev-protocol-v1',
        'registered_time_utc': datetime.now(timezone.utc).isoformat(),
        'registered_before_new_diagnostics': True,
        'research_record_date': '2026-10-09',
        'scope': 'Fixed finite adaptation of published ChebAE sampling geometry with exact-rational sequential confidence, complete retained pricing preparations and a canonical QPE fallback. No theoretical novelty or physical advantage asserted.',
        'financial_input': 'results/four_survivor_quantum_resources.json',
        'input_sha256': {'results/four_survivor_quantum_resources.json': sha(repo / 'results/four_survivor_quantum_resources.json'),
                         'results/qft_rotation_tapes.json': sha(repo / 'results/qft_rotation_tapes.json'),
                         'results/qft_rotation_certificates.json': sha(repo / 'results/qft_rotation_certificates.json')},
        'source_sha256': source_manifest(repo), 'cases': retained_cases(repo),
        'alpha_per_complete_price': str(ALPHA), 'alpha_pilot': str(ALPHA / 2),
        'alpha_QPE_fallback': str(ALPHA / 2), 'alpha_per_look': str(ALPHA / (2 * LOOK_CAP)),
        'simultaneous_four_price_union_failure_upper': str(4 * ALPHA),
        'max_looks': LOOK_CAP, 'degree_cap': DEGREE_CAP, 'fresh_batch_shots': BATCH_SHOTS,
        'pilot_T_cap_rule': 'floor(fresh half-alpha canonical QPE fallback T / 10)',
        'whole_history_T_cap_rule': 'pilot cap plus the complete fallback T, including discarded pilot shots',
        'interval_method': 'Fixed two-sided Clopper-Pearson binomial inversion with outward 64-bit dyadic endpoints, separately certified by exact rational tail evaluations; bounded Hoeffding retained as a reference.',
        'degree_rule': 'Largest currently monotone degree in 1 through 2047; subsequent increases require at least doubling. Full current interval is retained; no pooling of distinct fresh batches.',
        'diagnostics': {'prior_mesh_fractions': ['0', '1/4', '1/2', '3/4', '1'],
            'fixed_seed_indices': list(SEEDS), 'seed_rule': 'base 1729 + 1000000*case_index + 1000*mesh_index + seed_index',
            'signal_evaluation': '256-bit Arb enclosure; q64 floor used only by the reproducible pseudorandom Bernoulli simulator',
            'bernoulli_draw_rule': 'random.Random(seed).getrandbits(64) < floor(2^64 signal_probability)',
            'interpretation': 'Synthetic probability-prior diagnostics, not a sampled financial price, physical quantum execution or a confidence validation of randomness.'},
        'data_origin': 'Only previously retained synthetic model certificates, quantum resource profiles and compiled Fourier tapes; no LSEG observations, derived LSEG data or classical price answers.'}
    write(target, protocol)
    return protocol


def load_controller():
    sys.path.insert(0, str(HERE))
    return importlib.import_module('chebyshev_controller')


def verify_protocol(repo, output, protocol):
    assert protocol['source_sha256'] == source_manifest(repo)
    for name, checksum in protocol['input_sha256'].items():
        assert sha(repo / name) == checksum
    assert protocol['cases'] == retained_cases(repo)


def run(repo, output):
    zero_reflection, expand_toffoli, compiled_extended_plan, charge_full_oracle = load_dependencies(repo)
    controller = load_controller()
    protocol = json.loads((output / 'chebyshev_protocol.json').read_text())
    verify_protocol(repo, output, protocol)
    from flint import arb, ctx
    ctx.prec = 256
    cases = []
    for ci, case in enumerate(protocol['cases']):
        financial = case['financial_preparation_profile']
        reflection = emitted_reflection(financial['reflection_core_qubits'], zero_reflection, expand_toffoli)
        fb = fallback(case, compiled_extended_plan, charge_full_oracle)
        cap = fb['T_count'] // 10
        costs = {d: BATCH_SHOTS * build_profile(financial, reflection, d)['per_shot_metrics']['T'] for d in range(1, DEGREE_CAP + 1)}
        setup = controller.plan(F(case['probability_lower']), F(case['probability_upper']), F(case['probability_tolerance']), ALPHA, cap, costs)
        diagnostics = []
        used = {1, 2, 3, 4, DEGREE_CAP}
        plo, phi = F(case['probability_lower']), F(case['probability_upper'])
        for mi, fraction in enumerate(map(F, protocol['diagnostics']['prior_mesh_fractions'])):
            truth = plo + fraction * (phi - plo)
            truth_ball = arb(truth.numerator) / truth.denominator
            theta = truth_ball.sqrt().asin()
            signals = {}
            for seed_index in SEEDS:
                seed = 1729 + 1000000 * ci + 1000 * mi + seed_index
                rng = random.Random(seed)
                raw_batches = []
                def count_ones(degree, n):
                    if degree not in signals:
                        value = (degree * theta).sin() ** 2 if degree % 2 else (degree * theta).cos() ** 2
                        scaled = value * (1 << 64)
                        lower, upper = scaled.lower(), scaled.upper()
                        if int(lower.floor().unique_fmpz()) != int(upper.floor().unique_fmpz()):
                            raise ArithmeticError('Signal q64 floor straddles a boundary')
                        threshold = int(lower.floor().unique_fmpz())
                        signals[degree] = {'probability_q64_numerator': threshold,
                            'probability_q64_denominator': 1 << 64,
                            'Arb_enclosure': str(value),
                            'simulator_rounding_absolute_upper': '1/18446744073709551616'}
                    threshold = signals[degree]['probability_q64_numerator']
                    words = [rng.getrandbits(64) for _ in range(n)]
                    successes = sum(word < threshold for word in words)
                    raw_batches.append({'degree': degree, 'shots': n, 'successes': successes,
                        'seeded_u64_words': words, 'batch_T_count': costs[degree],
                        'signal': signals[degree]})
                    used.add(degree)
                    return successes
                estimate = controller.estimate(count_ones, setup)
                # A fallback supplies the theorem's canonical guaranteed radius;
                # this diagnostic records its complete cost, without generating
                # a phase observation or any financial answer.
                needs = estimate['status'] != 'certified'
                total = estimate['pilot_T_count'] + (fb['T_count'] if needs else 0)
                interval = estimate['probability_interval']
                covers = interval is not None and F(interval[0]) <= truth <= F(interval[1])
                diagnostics.append({'mesh_index': mi, 'mesh_fraction': str(fraction),
                    'truth_probability': str(truth), 'seed_index': seed_index, 'seed': seed,
                    'controller_result': estimate, 'raw_batches': raw_batches,
                    'pilot_interval_covers_truth': covers,
                    'fallback_invoked_in_resource_history': needs,
                    'actual_history_T_count': total,
                    'fallback_execution_scope': 'Cost charged when required; no QPE outcome simulated or price returned.' if needs else 'Not required in this synthetic history.'})
        totals = [x['actual_history_T_count'] for x in diagnostics]
        stops = sum(x['controller_result']['status'] == 'certified' for x in diagnostics)
        cases.append({'model_id': case['model_id'], 'retained_financial_information': case,
            'uncontrolled_zero_reflection': reflection,
            'half_alpha_QPE_fallback': fb, 'pilot_T_cap': cap,
            'whole_history_T_cap': cap + fb['T_count'],
            'controller_setup': controller.rational_record(setup),
            'batch_cost_by_degree_sha256': digest(costs),
            'sampler_profiles': {str(d): build_profile(financial, reflection, d) for d in sorted(used)},
            'diagnostics': diagnostics,
            'diagnostic_summary': {'histories': len(diagnostics), 'pilot_certified': stops,
                'fallback_charged': len(diagnostics) - stops,
                'observed_history_T_min': min(totals), 'observed_history_T_median': statistics.median(totals),
                'observed_history_T_max': max(totals),
                'coverage_diagnostic_failures': sum(not x['pilot_interval_covers_truth'] for x in diagnostics if x['controller_result']['status'] == 'certified'),
                'scope': 'Descriptive fixed synthetic mesh and seeds only; these observations do not certify best-case or expected complexity.'}})
        print(case['model_id'], 'pilot certified', stops, 'of', len(diagnostics), 'max T', max(totals), flush=True)
    verify_protocol(repo, output, protocol)
    result = {'schema': 'qbp-bounded-chebyshev-resource-diagnostics-v1',
        'status': 'completed', 'protocol_sha256': sha(output / 'chebyshev_protocol.json'),
        'source_sha256': sha(__file__), 'controller_source_sha256': sha(HERE / 'chebyshev_controller.py'),
        'input_sha256': protocol['input_sha256'], 'cases': cases,
        'alpha_per_complete_price': str(ALPHA), 'simultaneous_four_price_union_failure_upper': str(4 * ALPHA),
        'scope': protocol['scope'], 'data_origin': protocol['data_origin']}
    write(output / 'chebyshev_resources.json', result)
    summary = {'schema': 'qbp-bounded-chebyshev-resource-summary-v1',
        'result_sha256': sha(output / 'chebyshev_resources.json'), 'cases': [{
            'model_id': c['model_id'], 'previous_complete_QPE_T': c['retained_financial_information']['previous_complete_QPE_T'],
            'half_alpha_QPE_fallback_T': c['half_alpha_QPE_fallback']['T_count'],
            'pilot_T_cap': c['pilot_T_cap'], 'whole_history_T_cap': c['whole_history_T_cap'],
            'diagnostic_summary': c['diagnostic_summary']} for c in cases],
        'scope': protocol['scope'], 'data_origin': protocol['data_origin']}
    write(output / 'chebyshev_resource_summary.json', summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=HERE)
    parser.add_argument('--output', type=Path, default=HERE / 'results')
    parser.add_argument('action', choices=['register', 'run'])
    args = parser.parse_args()
    repo, output = args.repo.resolve(), args.output.resolve()
    register(repo, output) if args.action == 'register' else run(repo, output)
