"""Independent moment, coefficient, interval, estimator and selection audit.

No producer prior/plan function is called. The original scalar/circuit pools
are read unchanged; their separate audits establish their base certificates.
This audit recomputes the new information and every retained/refined strategy.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse, gzip, hashlib, json, sys, time
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))
from artifact_archive import archive_sha256, load_json, iter_top_level_array
from audit_higherdim_precision import check_plan, charge_check
EPS, ALPHA = F(1, 100), F(1, 3600)
KINDS = ('boundary_r2', 'stable_gauss_r2', 'raw_direct_r2', 'positive_r1')


def sha(path):
    path = Path(path)
    return archive_sha256(path) if path.suffix == '.gz' else hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def read(path):
    path = Path(path)
    return load_json(path) if path.suffix == '.gz' else json.loads(path.read_text())


def lower(record):
    return F(record['lower_exact_fraction'])


def upper(record):
    return F(record['upper_exact_fraction'])


def ball(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def endpoints(value):
    result = []
    for endpoint in (value.lower(), value.upper()):
        mantissa, exponent = map(int, endpoint.man_exp())
        result.append(F(mantissa * (1 << exponent)) if exponent >= 0 else F(mantissa, 1 << -exponent))
    return tuple(result)


def serial(value):
    lo, hi = endpoints(value)
    return {'lower_exact_fraction': str(lo), 'upper_exact_fraction': str(hi)}


def geometry(model):
    w = tuple(map(F, model['weights']))
    rows = tuple(tuple(map(F, row)) for row in model['loadings'])
    mu = tuple(map(F, model['log_drifts']))
    lam = F(model['lambda'])
    p = model['outer_dimension']
    assert len(w) == model['assets'] and all(len(row) == p for row in rows)
    assert F(model['spot']) == F(model['strike']) == 100 and F(model['maturity']) == 1
    covariance = tuple(tuple(map(F, row)) for row in model['covariance'])
    assert all(covariance[i][j] == lam * lam + sum(a * b for a, b in zip(rows[i], rows[j])) for i in range(len(w)) for j in range(len(w)))
    assert min(w) > 0 and sum(w) == 1 and lam > 0
    rate = F(model['risk_free_rate']) - F(model['dividend_yield'])
    for drift, row in zip(mu, rows):
        assert drift + (lam * lam + sum(v * v for v in row)) / 2 == rate
    mg = sum(a * b for a, b in zip(w, mu))
    rg = tuple(sum(w[i] * rows[i][j] for i in range(len(w))) for j in range(p))
    d = tuple(x - mg for x in mu)
    s = tuple(tuple(x - y for x, y in zip(row, rg)) for row in rows)
    assert sum(a * b for a, b in zip(w, d)) == 0
    assert all(sum(w[i] * s[i][j] for i in range(len(w))) == 0 for j in range(p))
    return w, d, s, mg, rg, lam, rate


def independent_U2(model):
    """Full ordered-pair Gaussian expansion, not producer compositions."""
    ctx.prec = 832
    w, d, s, mg, rg, lam, rate = geometry(model)
    first = sum((ball(wi) * ball(di + sum(v * v for v in si) / 2).exp() for wi, di, si in zip(w, d, s)), arb(0))
    second = arb(0)
    for i in range(len(w)):
        for j in range(len(w)):
            combined = tuple(a + b for a, b in zip(s[i], s[j]))
            exponent = d[i] + d[j] + sum(v * v for v in combined) / 2
            second += ball(w[i] * w[j]) * ball(exponent).exp()
    return second - 2 * first + 1


def independent_c0(model, cap):
    """Direct single-lognormal call spread, independent of coefficient builder."""
    ctx.prec = 832
    w, d, s, mg, rg, lam, rate = geometry(model)
    variance = lam * lam + sum(v * v for v in rg)
    sigma = ball(variance).sqrt()
    sqrt2 = arb(2).sqrt()
    def Phi(z):
        return (1 + (z / sqrt2).erf()) / 2
    forward = 100 * ball(mg + variance / 2).exp()
    if cap is None:
        return forward * Phi(ball(mg + variance) / sigma) - 100 * Phi(ball(mg) / sigma)
    width = ball(1 + F(cap, 100)).log()
    return forward * (Phi(ball(mg + variance) / sigma) - Phi((ball(mg + variance) - width) / sigma)) - 100 * Phi(ball(mg) / sigma) + (100 + cap) * Phi((ball(mg) - width) / sigma)


def independent_c1(model, cap):
    """Forward-price Gaussian tilt, avoiding producer prefactor products."""
    ctx.prec = 832
    w, d, s, mg, rg, lam, rate = geometry(model)
    variance = lam * lam + sum(v * v for v in rg)
    sigma = ball(variance).sqrt()
    sqrt2 = arb(2).sqrt()
    def Phi(z):
        return (1 + (z / sqrt2).erf()) / 2
    width = None if cap is None else ball(1 + F(cap, 100)).log()
    def probability(mean):
        value = Phi((ball(mean) + ball(variance)) / sigma)
        if width is not None:
            value -= Phi((ball(mean) - width + ball(variance)) / sigma)
        return value
    weighted = arb(0)
    for wi, di, si in zip(w, d, s):
        shifted = mg + sum(a * b for a, b in zip(rg, si))
        assert di + sum(v * v for v in si) / 2 + shifted + variance / 2 == rate
        weighted += ball(wi) * probability(shifted)
    return 100 * ball(rate).exp() * weighted - 100 * ball(mg + variance / 2).exp() * probability(mg)


def key(record):
    return record['model_id'], record['payoff'], None if record['cap'] is None else int(F(record['cap']))


def contained(fresh, retained):
    lo, hi = endpoints(fresh)
    assert lower(retained) <= lo <= hi <= upper(retained), 'Fresh Gaussian value is not contained by the retained enclosure'


def profile_family(record, strategy):
    family = record['profiles']
    implementation = strategy.get('implementation', 'original')
    if strategy['architecture'] in family:
        assert implementation == 'original'
        return family[strategy['architecture']]
    return family[implementation][strategy['architecture']]


def first_interval(base, coefficient, beta, mass, norm, moment_record):
    """Reconstruct each rational interval endpoint without producer APIs."""
    cap = key(base)[2]
    cfg = base['base_graph_certificate']['configuration']
    B = F(cfg['monetary_signal_normaliser'])
    eta = F(base['base_graph_certificate']['probability_error_upper'])
    previous = base['correction_prior']
    remainder = F(norm['M2_upper']) * upper(moment_record['continuous_U2_interval']) / 2
    c1 = coefficient['c1']
    lq = lower(c1) - (remainder if cap is not None else F(0))
    uq = upper(c1) + remainder
    ln = max(F(0), lq - beta) / upper(mass)
    un = min(B, (uq + beta) / lower(mass), F(previous['ideal_mean_upper']))
    assert 0 <= ln <= un <= B
    pl = max(F(0), ln / B - eta)
    pu = min(F(1), F(previous['probability_upper']), un / B + eta)
    assert 0 <= pl <= pu <= 1
    return {'probability_lower': pl, 'probability_upper': pu, 'finite_ideal_mean_lower': ln,
        'finite_ideal_mean_upper': un, 'ideal_mean_upper': un,
        'continuous_R2_mean_absolute_upper': remainder,
        'continuous_R1_mean_lower': lq, 'continuous_R1_mean_upper': uq,
        'financial_transfer_upper': beta, 'arithmetic_probability_enlargement': eta}


def audit_zero(base, result, lo, hi, det):
    """Exact midpoint and complete price-envelope calculation."""
    budget = base['complete_deterministic_budget']
    cfg = base['base_graph_certificate']['configuration']
    B = F(cfg['monetary_signal_normaliser'])
    ml, mh = lower(budget['mass_interval']), upper(budget['mass_interval'])
    coefficient = budget['coefficient_interval']
    cl, ch = lower(coefficient), upper(coefficient)
    cbias = F(budget['financial_transfer_upper']) + mh * (F(budget['auxiliary_components']['bias_upper']) + F(budget['arithmetic_before_mass_upper']))
    if lo is not None:
        point = (cl + ch) / 2 + (ml + mh) * B * (lo + hi) / 4
        lower_price = cl + ml * B * lo - cbias
        upper_price = ch + mh * B * hi + cbias
        radius = mh * B * (hi - lo) / 2
    else:
        point = (cl + ch) / 2
        lower_price = cl - mh * B * hi - cbias
        upper_price = ch + mh * B * hi + cbias
        radius = mh * B * hi
    exact_midpoint_error = max(abs(point - lower_price), abs(point - upper_price))
    assert exact_midpoint_error <= det + radius <= EPS
    assert F(result['certified_correction_radius']) == radius
    assert F(result['total_price_error_upper']) == det + radius
    assert result['T_count'] == result['quantum_shots'] == 0
    return {'undiscounted_price_midpoint': str(point), 'exact_complete_midpoint_error_upper': str(exact_midpoint_error),
        'claimed_error_upper': str(det + radius)}


def audit(root, out):
    root = Path(root)
    start = time.perf_counter()
    auditor_start_sha256 = sha(__file__)
    archive_adapter_start_sha256 = sha(root / 'artifact_archive.py')
    refined_path = root / 'results/complete_information_comparison.json.gz'
    refined_start_sha256 = sha(refined_path)
    refined = read(refined_path)
    assert refined['source_sha256'] == sha(root / 'complete_information_comparison.py')
    for name, expected in refined['dependency_source_sha256'].items():
        assert sha(root / name) == expected
    for name, expected in refined['input_sha256'].items():
        assert sha(root / 'results' / name) == expected
    allowed = {'boundary_precision_comparison.json.gz', 'higherdim_precision_comparison.json.gz',
        'generalisation_protocol.json', 'higherdim_preparation.json', 'boundary_financial_grid_transfer_summary.json',
        'tensor_first_transfer.json', 'boundary_derivative_norms.json', 'gaussian_derivative_norms.json',
        'qft_rotation_tapes.json', 'qft_rotation_certificates.json'}
    assert set(refined['input_sha256']) == allowed
    assert F(refined['monetary_radius']) == EPS and F(refined['alpha_per_selected_price']) == ALPHA
    assert F(refined['union_over_selected36_failure_upper']) == 36 * ALPHA == F(1, 100)
    models = {m['model_id']: m for m in read(root / 'results/generalisation_protocol.json')['models']}
    moments = {r['model_id']: r for r in refined['continuous_ratio_moment_certificates']}
    assert set(moments) == set(models) and len(models) == 9
    moment_checks = []
    for mid, model in models.items():
        fresh = independent_U2(model)
        contained(fresh, moments[mid]['continuous_U2_interval'])
        moment_checks.append({'model_id': mid, 'ordered_pair_terms': model['assets'] ** 2,
            'fresh_U2_interval': serial(fresh), 'retained_moment_record_sha256': digest(moments[mid])})
    lower_family = read(root / 'results/boundary_financial_grid_transfer_summary.json')
    prep = read(root / 'results/higherdim_preparation.json')
    coefficients = {key(c): c for c in lower_family['coefficients']}
    coefficients.update({key(c): c['coefficients'] for r in prep['records'] for c in r['cases']})
    assert len(coefficients) == 36
    c0_checks, c1_checks = [], []
    for case, coefficient in coefficients.items():
        fresh_c0 = independent_c0(models[case[0]], case[2])
        contained(fresh_c0, coefficient['c0'])
        c0_checks.append({'model_id': case[0], 'payoff': case[1], 'cap': case[2], 'fresh_direct_lognormal_c0_interval': serial(fresh_c0)})
        fresh = independent_c1(models[case[0]], case[2])
        contained(fresh_c0 + fresh, coefficient['coefficient'])
        contained(fresh, coefficient['c1'])
        c1_checks.append({'model_id': case[0], 'payoff': case[1], 'cap': case[2], 'fresh_forward_tilt_c1_interval': serial(fresh)})
    financial = {key(r): r for r in lower_family['records'] if r['order'] == 1}
    financial.update({key(r): r for r in read(root / 'results/tensor_first_transfer.json')['records'] if models[r['model_id']]['outer_dimension'] in (3, 7)})
    one = prep['shared_MGF_profile']['node_certificate']['mass_interval']
    rotationdata = read(root / 'results/qft_rotation_tapes.json')['rotations']
    rotationbounds = {r['key']: F(r['operator_norm_upper_rational']) for r in read(root / 'results/qft_rotation_certificates.json')['certificates']}
    norms = {(F(r['lambda']), None if r['cap'] is None else int(F(r['cap']))): r for r in read(root / 'results/boundary_derivative_norms.json')['records']}
    rows_by_family = {}
    for index, row in enumerate(refined['records']):
        rows_by_family.setdefault(row['family_archive'], []).append((index, row))
    assert set(rows_by_family) == {'boundary_precision_comparison.json.gz', 'higherdim_precision_comparison.json.gz'}
    interval_checks = retained_strategies = quantum_flags = quantum_strategies = zero_strategies = rejected_strategies = 0
    groups = {}
    selected_per_case = {}
    zero_midpoints = []
    profile_hashes = {}
    family_hashes = {}
    for family_name in ('boundary_precision_comparison.json.gz', 'higherdim_precision_comparison.json.gz'):
        family_path = root / 'results' / family_name
        family_hashes[family_name] = sha(family_path)
        rows = rows_by_family[family_name]
        base_records = iter_top_level_array(family_path, 'records')
        for (refined_index, row), (candidate_index, base) in zip(rows, enumerate(base_records), strict=True):
            # Cache only within a retained record; object ids can be reused after streaming.
            profile_hashes = {}
            assert row['candidate_index'] == candidate_index and row['family_archive_sha256'] == family_hashes[family_name]
            assert key(row) == key(base) and row['kind'] == base['kind']
            assert row['node_address_bits'] == base['node_address_bits'] and row['allocation_limits'] == base['allocation_limits']
            assert row['base_graph_certificate_sha256'] == base['base_graph_certificate']['precision_search_certificate_sha256']
            original_rows = []
            for strategy in base['strategies']:
                retained = dict(strategy)
                retained.setdefault('implementation', 'original')
                retained['prior_information'] = 'previous_proved_prior'
                original_rows.append(retained)
            assert row['strategies'][:len(original_rows)] == original_rows
            retained_strategies += len(original_rows)
            expected_prior = None
            has_refinement = base['kind'] == 'positive_r1' and base['complete_deterministic_budget'] is not None
            if has_refinement:
                case = key(base)
                supplied = row['coefficient_informed_prior']
                mass = supplied['external_mass_interval']
                p = models[case[0]]['outer_dimension']
                if p == 1:
                    assert mass == financial[case]['mass_interval']
                else:
                    assert 0 < lower(mass) <= lower(one) ** p <= upper(one) ** p <= upper(mass)
                    assert upper(mass) - lower(mass) < F(1, 10 ** 140)
                beta = F(financial[case]['financial_transfer_upper'])
                expected_prior = first_interval(base, coefficients[case], beta, mass, norms[F(models[case[0]]['lambda']), case[2]], moments[case[0]])
                for field, value in expected_prior.items():
                    assert F(supplied[field]) == value, (field, case)
                assert supplied['continuous_c1_interval'] == coefficients[case]['c1']
                assert supplied['continuous_R2_nonnegative'] == (case[2] is None)
                assert supplied['continuous_moment_record_sha256'] == digest(moments[case[0]])
                assert supplied['previous_probability_upper'] == base['correction_prior']['probability_upper']
                assert supplied['previous_ideal_mean_upper'] == base['correction_prior']['ideal_mean_upper']
                assert supplied['native_positive'] is True
                assert len(row['strategies']) == 2 * len(original_rows)
                interval_checks += 1
            else:
                assert len(row['strategies']) == len(original_rows) and 'coefficient_informed_prior' not in row
            for strategy_index, strategy in enumerate(row['strategies']):
                is_new = strategy['prior_information'] == 'coefficient_informed_interval'
                if is_new:
                    assert has_refinement and expected_prior is not None
                    old = base['strategies'][strategy_index - len(original_rows)]
                    old_headers = {k: v for k, v in old.items() if k != 'result'}
                    old_headers.setdefault('implementation', 'original')
                    old_headers['prior_information'] = 'coefficient_informed_interval'
                    assert {k: v for k, v in strategy.items() if k != 'result'} == old_headers
                    prior = row['coefficient_informed_prior']
                else:
                    assert strategy['prior_information'] == 'previous_proved_prior'
                    prior = base['correction_prior']
                result = strategy['result']
                budget = base['complete_deterministic_budget']
                det = F(budget['total_deterministic_undiscounted_price_error_upper'])
                cfg = base['base_graph_certificate']['configuration']
                B = F(cfg['monetary_signal_normaliser'])
                C = F(budget['external_factor_upper']) * B
                remaining = EPS - det
                positive = cfg['output_encoding'] == 'positive'
                pu = F(prior['probability_upper'])
                pl = F(prior['probability_lower']) if is_new else (F(0) if positive else None)
                radius = C * (pu - pl) / 2 if pl is not None else C * pu
                if remaining <= 0:
                    assert result['kind'] == 'rejected'
                    rejected_strategies += 1
                    continue
                if radius <= remaining:
                    assert result['kind'] == 'deterministic_correction_prior'
                    checked = audit_zero(base, result, pl, pu, det)
                    if is_new:
                        assert F(result['implemented_probability_midpoint']) == (pl + pu) / 2
                        assert result['implemented_probability_interval'] == {'lower': str(pl), 'upper': str(pu)}
                    zero_midpoints.append({'refined_candidate_index': refined_index,
                        'strategy_index': strategy_index, 'prior_information': strategy['prior_information'], **checked})
                    zero_strategies += 1
                    continue
                assert result['kind'] == 'complete_quantum_resource_plan'
                encoding = strategy['strategy']
                labels = ['centred'] if positive or encoding == 'centred' else ['positive_part', 'negative_part']
                target = remaining / C if positive else remaining / (2 * C)
                priors = [pu] * len(labels) if positive or encoding == 'split' else [F(1)]
                alphas = [ALPHA / len(labels)] * len(labels)
                decode = 1 if positive or encoding == 'split' else 2
                assert [f['flag'] for f in result['flags']] == labels
                stat = failure = F(0)
                T = shots = peak = 0
                selected_profile_family = profile_family(base, strategy)
                for label, prior_upper, alpha, flag in zip(labels, priors, alphas, result['flags']):
                    profile = selected_profile_family[label]
                    profile_key = id(profile)
                    if profile_key not in profile_hashes:
                        profile_hashes[profile_key] = digest(profile)
                    assert flag['financial_profile_sha256'] == profile_hashes[profile_key]
                    assert F(flag['probability_tolerance']) == target and F(flag['prior_upper']) == prior_upper
                    check_plan(flag['plan'], target, alpha, prior_upper, flag['fourier_profile'], rotationdata, rotationbounds)
                    charge_check(flag['complete_charge'], profile, flag['fourier_profile'], flag['plan'])
                    stat += F(flag['plan']['amplitude_error_upper']) + F(flag['plan']['postprocess_allowance'])
                    failure += F(flag['plan']['binomial_failure_upper'])
                    T += flag['complete_charge']['metrics']['T']
                    shots += flag['plan']['repetitions']
                    peak = max(peak, flag['complete_charge']['peak_logical_qubits'])
                    quantum_flags += 1
                stat *= C * decode
                assert F(result['statistical_price_error_upper']) == stat
                assert F(result['total_price_error_upper']) == det + stat <= EPS
                assert F(result['failure_probability_upper']) == failure <= ALPHA
                assert result['T_count'] == T and result['quantum_shots'] == shots and result['logical_qubits_upper'] == peak
                if is_new:
                    assert result['coefficient_informed_probability_interval'] == {'lower': str(pl), 'upper': str(pu)}
                quantum_strategies += 1
            eligible = [s for s in row['strategies'] if s['result']['kind'] != 'rejected']
            winner = min(eligible, key=lambda s: (s['result']['T_count'], s['implementation'], s['architecture'], s['strategy'], s['prior_information'])) if eligible else None
            assert row['selected_strategy'] == winner
            groups.setdefault((*key(row), row['kind']), []).append((refined_index, row))
        del base_records
    assert len(groups) == 36 * len(KINDS)
    assert len(refined['selected_prices']) == 36
    selected_cases = set()
    selected_quantum = selected_zero = shared_stops = 0
    shared_screen_checks = []
    for row in refined['selected_prices']:
        case = key(row)
        assert case not in selected_cases and case in coefficients
        selected_cases.add(case)
        base_interval = coefficients[case]['coefficient']
        M2 = F(norms[F(models[case[0]]['lambda']), case[2]]['M2_upper'])
        remainder = M2 * upper(moments[case[0]]['continuous_U2_interval']) / 2
        sl = max(F(0), lower(base_interval) - (remainder if case[2] is not None else F(0)))
        sh = upper(base_interval) + remainder
        if case[2] is not None:
            sh = min(F(case[2]), sh)
        assert 0 <= sl <= sh
        screen = row['shared_continuous_screen']
        assert key(screen) == case
        assert screen['price_interval'] == {'lower_exact_fraction': str(sl), 'upper_exact_fraction': str(sh)}
        assert F(screen['price_midpoint']) == (sl + sh) / 2
        assert F(screen['price_radius_upper']) == (sh - sl) / 2
        assert F(screen['monetary_radius']) == EPS
        stops = (sh - sl) / 2 <= EPS
        assert screen['stops_without_quantum'] is stops
        if stops:
            shared_stops += 1
        shared_screen_checks.append({'model_id': case[0], 'payoff': case[1], 'cap': case[2],
            'price_midpoint': str((sl + sh) / 2), 'price_radius_upper': str((sh - sl) / 2),
            'stops_without_quantum': stops})
        for kind in KINDS:
            eligible = [(index, r) for index, r in groups[(*case, kind)] if r['selected_strategy']]
            winner = min(eligible, key=lambda pair: (pair[1]['selected_strategy']['result']['T_count'], pair[0])) if eligible else None
            selected = row['methods'][kind]
            assert (selected is None) == (winner is None)
            if winner:
                index, candidate = winner
                assert selected['refined_candidate_index'] == index
                assert selected['family_archive'] == candidate['family_archive'] and selected['candidate_index'] == candidate['candidate_index']
                if stops:
                    assert selected['conditional_estimator_selection'] == candidate['selected_strategy']
                    stopped = selected['selected']
                    assert stopped['implementation'] == 'shared_continuous_screen'
                    assert stopped['architecture'] is None and stopped['strategy'] == 'analytical_price_interval'
                    assert stopped['prior_information'] == 'shared_low_order_control'
                    result = stopped['result']
                    assert result['kind'] == 'deterministic_correction_prior'
                    assert result['T_count'] == result['quantum_shots'] == 0
                    assert F(result['total_price_error_upper']) == (sh - sl) / 2 <= EPS
                    assert result['certified_price_interval'] == screen['price_interval']
                    assert F(result['price_midpoint']) == (sl + sh) / 2
                else:
                    assert 'conditional_estimator_selection' not in selected
                    assert selected['selected'] == candidate['selected_strategy']
                if selected['selected']['result']['T_count'] == 0:
                    selected_zero += 1
                else:
                    selected_quantum += 1
        available = [k for k in KINDS if row['methods'][k]]
        conditional_winner = min(available, key=lambda kind: row['methods'][kind].get('conditional_estimator_selection', row['methods'][kind]['selected'])['result']['T_count']) if available else None
        assert row['conditional_quantum_resource_winner'] == conditional_winner
        winner = min(available, key=lambda kind: row['methods'][kind]['selected']['result']['T_count']) if available else None
        assert row['quantum_resource_winner'] == winner
    assert selected_cases == set(coefficients)
    for family_name, expected in family_hashes.items():
        assert sha(root / 'results' / family_name) == expected
    assert sha(refined_path) == refined_start_sha256
    assert sha(__file__) == auditor_start_sha256
    assert sha(root / 'artifact_archive.py') == archive_adapter_start_sha256
    assert sha(root / 'complete_information_comparison.py') == refined['source_sha256']
    for name, expected in refined['dependency_source_sha256'].items():
        assert sha(root / name) == expected
    for name, expected in refined['input_sha256'].items():
        assert sha(root / 'results' / name) == expected
    result = {'schema': 'qbp-complete-information-independent-audit-v1', 'status': 'passed',
        'continuous_U2_ordered_pair_models': len(moment_checks), 'independent_direct_lognormal_c0_contracts': len(c0_checks), 'independent_forward_tilt_c1_contracts': len(c1_checks),
        'shared_continuous_screen_cases_checked': len(shared_screen_checks), 'shared_continuous_analytical_stops': shared_stops,
        'coefficient_informed_candidate_intersections': interval_checks,
        'unchanged_previous_strategy_records': retained_strategies,
        'complete_quantum_strategies_checked': quantum_strategies, 'complete_QPE_flag_plans_checked': quantum_flags,
        'exact_analytical_midpoint_strategy_records_checked': zero_strategies,
        'rejected_strategy_records_checked': rejected_strategies,
        'candidate_selection_ties_checked': len(refined['records']), 'method_selection_ties_checked': 36 * len(KINDS),
        'selected_price_cases_checked': 36, 'selected_method_quantum_plans': selected_quantum,
        'selected_method_analytical_stops': selected_zero,
        'moment_checks': moment_checks, 'c0_checks': c0_checks, 'c1_checks': c1_checks, 'shared_screen_checks': shared_screen_checks, 'exact_analytical_midpoints': zero_midpoints,
        'comparison_sha256': refined_start_sha256, 'source_sha256': auditor_start_sha256,
        'portable_archive_adapter_sha256': archive_adapter_start_sha256,
        'reused_independent_QPE_audit_sha256': sha(root / 'checks/audit_higherdim_precision.py'),
        'producer_source_sha256': {p: sha(root / p) for p in ('correction_information_priors.py', 'complete_information_comparison.py')},
        'input_sha256': refined['input_sha256'],
        'precision_bits': 832, 'elapsed_seconds_diagnostic': time.perf_counter() - start,
        'scope': 'Independent Gaussian ordered-pair moments, direct lognormal c0 and alternate forward-tilt c1 coefficients, all shared continuous screens and final analytical overrides, every strengthened finite/implemented interval, preserved old strategies, exact rational analytical midpoint errors, canonical binary QPE/finite-binomial confidence, complete physical logical charges and all declared selection ties. Scalar/circuit pool correctness is covered by separately retained base audits; no financial price-answer input or theoretical-priority claim.'}
    Path(out).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if isinstance(v, (str, int, float)) and k not in ('scope', 'comparison_sha256', 'source_sha256', 'reused_independent_QPE_audit_sha256')}))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', default=str(ROOT))
    parser.add_argument('--out', default=str(ROOT / 'results/complete_information_independent_audit.json'))
    args = parser.parse_args()
    audit(args.root, args.out)
