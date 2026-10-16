"""Matched finite Chebyshev challenge with paired common-signal diagnostics.

The common signal is stipulated within independently proved priors; it is
not a price answer or an assertion that distinct implemented flags coincide.
All canonical inputs are read only. Register before running new outcomes.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from dataclasses import asdict
from datetime import datetime, timezone
import argparse, gzip, hashlib, importlib, json, random, statistics, sys

HERE = Path(__file__).resolve().parent
ALPHA = F(1, 3600)
METHODS = ('shared_threshold', 'raw_payoff', 'numerical_conditional')
MESH = tuple(map(F, ('0', '1/4', '1/2', '3/4', '1')))
SEEDS = tuple(range(20))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def digest(v):
    return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'),
                                    default=str).encode()).hexdigest()


def read(p):
    raw = Path(p).read_bytes()
    return json.loads(gzip.decompress(raw) if str(p).endswith('.gz') else raw)


def write(p, value):
    p = Path(p)
    if p.exists():
        raise FileExistsError('Preserve executed evidence: '+str(p))
    p.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(value, indent=2, default=str)+'\n').encode()
    p.write_bytes(gzip.compress(raw, mtime=0) if str(p).endswith('.gz') else raw)


def load(repo):
    sys.path.insert(0, str(repo))
    return importlib.import_module('chebyshev_resources'), importlib.import_module('chebyshev_controller')


def source_manifest(repo):
    result = {Path(__file__).name: sha(__file__)}
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            p = Path(path)
            if p.parent == repo and p.suffix == '.py':
                result[p.name] = sha(p)
    return dict(sorted(result.items()))


def current_cases(repo):
    paths = ('results/asymmetric_raw_payoff/asymmetric_raw_payoff_resources.json.gz',
             'results/shared_threshold_resources.json', 'results/strengthened_r1_resources.json')
    data = [read(repo/p) for p in paths]
    profiles = {}
    for obj in data:
        for row in obj.get('records', obj.get('rows', [])):
            for architecture, flags in row['profiles'].items():
                for flag, profile in flags.items():
                    profiles[digest(profile)] = profile
    cases = []
    keys = ('threshold_selected_strategy', 'selected_raw_strategy', 'numerical_selected_strategy')
    for comparison in data[0]['selected_comparisons']:
        methods = {}
        for method, key in zip(METHODS, keys):
            strategy = comparison[key]
            assert strategy['architecture'] == 'unitary8t'
            result = strategy['result']
            assert len(result['flags']) == 1
            flag = result['flags'][0]
            profile = profiles[flag['financial_profile_sha256']]
            assert digest(profile) == flag['financial_profile_sha256']
            prior = result['correction_prior']
            methods[method] = {'architecture': 'unitary8t',
                'financial_preparation_profile': profile,
                'financial_preparation_profile_sha256': digest(profile),
                'probability_lower': prior['probability_lower'],
                'probability_upper': prior['probability_upper'],
                'finite_ideal_mean_lower': prior['finite_ideal_mean_lower'],
                'finite_ideal_mean_upper': prior['finite_ideal_mean_upper'],
                'probability_tolerance': flag['probability_tolerance'],
                'monetary_signal_factor_upper': result['monetary_signal_factor_upper'],
                'deterministic_price_error_upper': result['deterministic_error_upper'],
                'full_alpha_canonical_T': result['T_count'],
                'selected_strategy_sha256': digest(strategy),
                'target_monetary_radius': '1/100', 'alpha': str(ALPHA)}
        cases.append({'model_id': comparison['model_id'], 'methods': methods})
    return cases, {p: sha(repo/p) for p in paths}


def supplied_cases(path):
    data = read(path)
    assert data['schema'] == 'qbp-matched-chebyshev-inputs-v1'
    cases = data['cases']
    assert len(cases) == 4
    for case in cases:
        assert tuple(case['methods']) == METHODS
        for method in METHODS:
            m = case['methods'][method]
            assert m['architecture'] == 'unitary8t'
            assert digest(m['financial_preparation_profile']) == m['financial_preparation_profile_sha256']
            assert m['target_monetary_radius'] == '1/100' and F(m['alpha']) == ALPHA
            assert F(m['monetary_signal_factor_upper'])*F(m['probability_tolerance'])+F(m['deterministic_price_error_upper']) <= F(1, 100)
    return cases, {str(Path(path).resolve()): sha(path)}


def common_signal_interval(case):
    lower = max(max(F(m['probability_lower']), F(m['finite_ideal_mean_lower']))
                for m in case['methods'].values())
    upper = min(min(F(m['probability_upper']), F(m['finite_ideal_mean_upper']))
                for m in case['methods'].values())
    assert 0 <= lower <= upper <= 1
    return lower, upper


def annotate_deferred_sampler_metadata(resources):
    """Correct historical clean-work wording; preserve all numerical logic."""
    original = resources.build_profile
    def annotated(preparation, reflection, degree):
        assert preparation['work_dirty_after_preparation']
        profile = original(preparation, reflection, degree)
        profile['financial_cleanup_scope'] = (
            'The supplied preparation is B=P U H, with correlated arithmetic work. '
            'Its exact adjoint restores zero work before the same actual reflection tape. '
            'Odd terminal work remains dirty until full-allocation observation/reset; '
            'even core-zero observation occurs in the restored clean-work invariant. '
            'The full original allocated-wire measurement/reset charge is retained.')
        profile['preparation_symbol'] = 'B'
        profile['preparation_conjugacy_proof_sha256'] = preparation['proof_sha256']
        profile['financial_preparation_profiles_reused'] = True
        return profile
    resources.build_profile = annotated


def register(repo, output, inputs):
    resources, controller = load(repo)
    resources.load_dependencies(repo)
    cases, evidence = supplied_cases(inputs) if inputs else current_cases(repo)
    for case in cases:
        case['common_signal_interval'] = list(map(str, common_signal_interval(case)))
    protocol = {'schema': 'qbp-matched-chebyshev-protocol-v1',
        'registered_time_utc': datetime.now(timezone.utc).isoformat(),
        'registered_before_new_outcomes': True,
        'prior_outcome_disclosure': 'All previous clean canonical and matched Chebyshev diagnostics, the cleanup-deferral identities, and every B/optimised-clean canonical and half-alpha cost are already known. No fresh B diagnostic history has yet run. This is a disclosed sequential matched resource challenge using the same prior mesh, seeds, controller, confidence allocation, limits and pilot-cap ratio.',
        'source_sha256': source_manifest(repo), 'evidence_sha256': evidence,
        'input_mode': 'supplied' if inputs else 'current',
        'input_path': str(Path(inputs).resolve()) if inputs else None,
        'cases': cases, 'methods': METHODS, 'alpha_per_alternative_price': str(ALPHA),
        'alpha_pilot': str(ALPHA/2), 'alpha_fallback': str(ALPHA/2),
        'alpha_per_fresh_look': str(ALPHA/24),
        'all_twelve_alternative_union_failure_upper': str(12*ALPHA),
        'one_selected_method_on_four_cases_union_failure_upper': str(4*ALPHA),
        'max_looks': 12, 'batch_shots': 100, 'degree_cap': 2047,
        'pilot_cap': 'floor(fresh half-alpha canonical complete T/10)',
        'whole_history_cap': 'pilot cap plus fresh half-alpha canonical complete T',
        'mesh_fractions': list(map(str, MESH)), 'seed_indices': SEEDS,
        'seed_rule': '1729+1000000*case_index+1000*mesh_index+seed_index',
        'paired_words_rule': 'Precompute12 fresh100-word batches for every seed. Method look j uses the identical batch j; method degrees and confidence decisions retain their own prior/tolerance.',
        'word_record': 'Seed, batch-index and SHA-256 retain exact deterministic word replay; full words are reconstructed rather than duplicated.',
        'signal_rule': '256-bit Arb parity-correct signal; q64 floor only for descriptive Bernoulli simulation.',
        'primary_comparison': 'Complete worst-history cap at the same monetary guarantee and confidence. Mesh means/medians are descriptive only.',
        'scope': 'Matched finite bounded adaptation, not unchanged published ChebAE or every strongest quantum estimator. Common ideal-signal diagnostic probabilities are stipulated, not actual financial flags or classical price answers.',
        'data_origin': 'Previously retained synthetic certificates and complete logical profiles only. No licensed data or physical quantum shots.'}
    write(output/'matched_chebyshev_protocol.json', protocol)
    return protocol


def verify(repo, protocol):
    assert source_manifest(repo) == protocol['source_sha256']
    for name, checksum in protocol['evidence_sha256'].items():
        path = Path(name) if Path(name).is_absolute() else repo/name
        assert sha(path) == checksum


def run(repo, output):
    resources, controller = load(repo)
    zero, expand, compiled, charge = resources.load_dependencies(repo)
    protocol = read(output/'matched_chebyshev_protocol.json')
    annotate_deferred_sampler_metadata(resources)
    verify(repo, protocol)
    from flint import arb, ctx
    ctx.prec = 256
    all_cases = []
    for ci, case in enumerate(protocol['cases']):
        lower, upper = map(F, case['common_signal_interval'])
        truths = [lower+f*(upper-lower) for f in MESH]
        signal_cache = {}
        setups, method_results = {}, {}
        for method in METHODS:
            info = case['methods'][method]
            financial = info['financial_preparation_profile']
            reflection = resources.emitted_reflection(financial['reflection_core_qubits'], zero, expand)
            fb = resources.fallback(info, compiled, charge)
            cap = fb['T_count']//10
            aT = financial['primitive_counts'].get('T', 0)+financial['primitive_counts'].get('Tdg', 0)
            rT = reflection['T_count']
            costs = {d: 100*(d*aT+(d//2 if d%2 else d//2-1)*rT) for d in range(1, 2048)}
            setup = controller.plan(F(info['probability_lower']), F(info['probability_upper']),
                F(info['probability_tolerance']), ALPHA, cap, costs)
            setups[method] = setup
            method_results[method] = {'financial_preparation_profile_sha256': info['financial_preparation_profile_sha256'],
                'uncontrolled_zero_reflection': reflection, 'half_alpha_QPE_fallback': fb,
                'full_alpha_canonical_T': info['full_alpha_canonical_T'],
                'pilot_T_cap': cap, 'whole_history_T_cap': cap+fb['T_count'],
                'controller_setup': controller.rational_record(setup),
                'batch_cost_by_degree_sha256': digest(costs),
                'used_degrees': {1, 2, 3, 4, 2047}, 'diagnostics': [],
                'financial_monetary_error_check': str(F(info['deterministic_price_error_upper'])+F(info['monetary_signal_factor_upper'])*F(info['probability_tolerance']))}
        for mi, truth in enumerate(truths):
            theta = (arb(truth.numerator)/truth.denominator).sqrt().asin()
            for si in SEEDS:
                seed = 1729+1000000*ci+1000*mi+si
                rng = random.Random(seed)
                words = [[rng.getrandbits(64) for _ in range(100)] for _ in range(12)]
                for method in METHODS:
                    record = method_results[method]
                    setup = setups[method]
                    batches = []
                    def count_ones(degree, n):
                        assert n == 100
                        look = len(batches)
                        key = (mi, degree)
                        if key not in signal_cache:
                            value = (degree*theta).sin()**2 if degree%2 else (degree*theta).cos()**2
                            scaled = value*(1<<64)
                            lo = int(scaled.lower().floor().unique_fmpz())
                            hi = int(scaled.upper().floor().unique_fmpz())
                            if lo != hi:
                                raise ArithmeticError('A simulator q64 floor crosses an Arb enclosure')
                            signal_cache[key] = {'probability_q64_numerator': lo,
                                'probability_q64_denominator': 1<<64, 'Arb_enclosure': str(value),
                                'simulator_rounding_absolute_upper': '1/18446744073709551616'}
                        signal = signal_cache[key]
                        successes = sum(w < signal['probability_q64_numerator'] for w in words[look])
                        batches.append({'look': look+1, 'degree': degree, 'shots': n,
                            'successes': successes, 'paired_seeded_words_sha256': digest(words[look]),
                            'signal': signal, 'batch_T_count': setup.batch_cost(degree)})
                        record['used_degrees'].add(degree)
                        return successes
                    estimate = controller.estimate(count_ones, setup)
                    requires = estimate['status'] != 'certified'
                    total = estimate['pilot_T_count']+(record['half_alpha_QPE_fallback']['T_count'] if requires else 0)
                    assert total <= record['whole_history_T_cap']
                    interval = estimate['probability_interval']
                    record['diagnostics'].append({'mesh_index': mi, 'mesh_fraction': str(MESH[mi]),
                        'truth_probability': str(truth), 'seed_index': si, 'seed': seed,
                        'controller_result': controller.rational_record(estimate), 'raw_batches': batches,
                        'pilot_interval_covers_stipulated_truth': F(interval[0]) <= truth <= F(interval[1]),
                        'fallback_invoked_in_resource_history': requires, 'actual_history_T_count': total,
                        'fallback_execution_scope': 'Full cost charged when required; no phase outcome or price simulated.'})
        for method in METHODS:
            record = method_results[method]
            info = case['methods'][method]
            financial = info['financial_preparation_profile']
            reflection = record['uncontrolled_zero_reflection']
            record['sampler_profiles'] = {str(d): resources.build_profile(financial, reflection, d)
                for d in sorted(record.pop('used_degrees'))}
            totals = [h['actual_history_T_count'] for h in record['diagnostics']]
            stops = sum(h['controller_result']['status'] == 'certified' for h in record['diagnostics'])
            mean = F(sum(totals), len(totals))
            pilot_mean = F(sum(h['controller_result']['pilot_T_count'] for h in record['diagnostics']), len(totals))
            tf, t0 = record['half_alpha_QPE_fallback']['T_count'], record['full_alpha_canonical_T']
            rate = F(stops, len(totals))
            assert mean == pilot_mean+(1-rate)*tf
            threshold = 1-F(t0, tf)+pilot_mean/tf
            assert (mean < t0) == (rate > threshold)
            record['diagnostic_summary'] = {'histories': len(totals), 'pilot_certified': stops,
                'fallback_charged': len(totals)-stops, 'observed_history_T_min': min(totals),
                'observed_history_T_median': str(F(str(statistics.median(totals)))),
                'observed_history_T_max': max(totals), 'observed_history_T_mean': str(mean),
                'descriptive_mean_pilot_T': str(pilot_mean), 'descriptive_exact_break_even_threshold': str(threshold),
                'descriptive_success_fraction': str(rate),
                'coverage_diagnostic_failures': sum(not h['pilot_interval_covers_stipulated_truth'] for h in record['diagnostics'] if h['controller_result']['status'] == 'certified'),
                'scope': 'Common stipulated prior mesh and paired seeds only; not expected cost at an actual financial flag.'}
        all_cases.append({'model_id': case['model_id'], 'common_signal_interval': case['common_signal_interval'],
            'methods': method_results})
        print(case['model_id'], [(m, method_results[m]['diagnostic_summary']['pilot_certified']) for m in METHODS], flush=True)
    verify(repo, protocol)
    result = {'schema': 'qbp-matched-chebyshev-results-v1', 'status': 'completed',
        'protocol_sha256': sha(output/'matched_chebyshev_protocol.json'),
        'source_sha256': source_manifest(repo), 'cases': all_cases,
        'scope': protocol['scope'], 'data_origin': protocol['data_origin']}
    write(output/'matched_chebyshev_resources.json.gz', result)
    summary = {'schema': 'qbp-matched-chebyshev-summary-v1', 'status': 'completed',
        'resources_sha256': sha(output/'matched_chebyshev_resources.json.gz'),
        'cases': [{'model_id': c['model_id'], 'methods': {m: {k: c['methods'][m][k] for k in
            ('full_alpha_canonical_T', 'pilot_T_cap', 'whole_history_T_cap', 'diagnostic_summary')}
            for m in METHODS}} for c in all_cases],
        'paired_histories': 400, 'alternative_histories': 1200,
        'scope': protocol['scope'], 'data_origin': protocol['data_origin']}
    write(output/'matched_chebyshev_summary.json', summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('register', 'run'))
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inputs', type=Path)
    args = parser.parse_args()
    repo, output = args.repo.resolve(), args.output.resolve()
    register(repo, output, args.inputs) if args.action == 'register' else run(repo, output)
