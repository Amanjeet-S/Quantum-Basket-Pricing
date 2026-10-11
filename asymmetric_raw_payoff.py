"""Disclosed, separately registered PG/PA precision challenge.

Old diagonal evidence is referenced exactly, rather than regenerated. Novel
settings are freshly certified and emitted. Positive projection is retained
because distinct exponential tables need not preserve endpoint ordering.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime, timezone
import argparse, copy, gzip, hashlib, importlib.util, json, sys, time
sys.dont_write_bytecode = True

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def digest(v):
    return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()

def read(p):
    p = Path(p)
    return json.load(gzip.open(p, 'rt')) if p.suffix == '.gz' else json.loads(p.read_text())

def write(p, v):
    p = Path(p)
    if p.exists():
        raise FileExistsError('Preserve prior evidence: '+str(p))
    if p.suffix == '.gz':
        with p.open('wb') as file:
            with gzip.GzipFile(fileobj=file, mode='wb', mtime=0) as stream:
                encoder = json.JSONEncoder(sort_keys=True, separators=(',', ':'), default=str)
                for block in encoder.iterencode(v):
                    stream.write(block.encode())
                stream.write(b'\n')
    else:
        p.write_text(json.dumps(v, indent=2, default=str)+'\n')

def module(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def configure(args):
    repo, out = Path(args.repo).resolve(), Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    raw = module('retained_raw_builder', repo/'raw_payoff_oracle_tight.py')
    raw.configure(repo, repo/'results/factor_basis_quantum_resources.json',
                  repo/'results/shared_threshold_resources.json',
                  repo/'results/strengthened_r1_resources.json', out)
    smooth = module('retained_smoothing', repo/'raw_payoff_smoothing_amendment.py')
    expanded = module('retained_expanded_helpers', repo/'expanded_raw_payoff.py')
    return repo, out, raw, smooth, expanded

def old_paths(repo):
    names = (
        'raw_payoff_tight/raw_payoff_resources.json.gz',
        'raw_payoff_tight/raw_payoff_scalar_ordinates.json.gz',
        'raw_payoff_tight/raw_payoff_protocol.json',
        'raw_payoff_smoothing/raw_payoff_smoothing_resources.json',
        'raw_payoff_smoothing/raw_payoff_smoothing_protocol.json',
        'raw_payoff_smoothing/raw_payoff_smoothing_audit.json',
        'expanded_raw_payoff/expanded_raw_payoff_resources.json.gz',
        'expanded_raw_payoff/expanded_raw_payoff_scalar_ordinates.json.gz',
        'expanded_raw_payoff/expanded_raw_payoff_protocol.json')
    return {n: repo/'results'/n for n in names}

def signatures(args, raw, repo):
    return {'producer_source_sha256': sha(__file__),
            'checked_asymmetric_proof_sha256': sha(args.proof),
            'retained_raw_builder_source_sha256': sha(repo/'raw_payoff_oracle_tight.py'),
            'retained_smoothing_source_sha256': sha(repo/'raw_payoff_smoothing_amendment.py'),
            'retained_expanded_helpers_sha256': sha(repo/'expanded_raw_payoff.py'),
            'old_evidence_sha256': {k: sha(p) for k, p in old_paths(repo).items()},
            **{k: v for k, v in raw.signatures().items() if k != 'source_sha256'}}

def known_records(repo):
    raw_old = read(repo/'results/raw_payoff_tight/raw_payoff_resources.json.gz')
    smooth_old = read(repo/'results/raw_payoff_smoothing/raw_payoff_smoothing_resources.json')
    expanded = read(repo/'results/expanded_raw_payoff/expanded_raw_payoff_resources.json.gz')
    old = {}
    for point, amended in zip(raw_old['records'], smooth_old['records']):
        assert point['original_model_id'] == amended['model_id']
        key = point['original_model_id'], point['normal_bits'], point['payoff_exp_bits']
        old[key] = (point, amended, 'raw_payoff_tight/raw_payoff_resources.json.gz',
                    'raw_payoff_smoothing/raw_payoff_smoothing_resources.json')
    for point, amended in zip(expanded['records'], expanded['amended_records']):
        assert point['original_model_id'] == amended['model_id']
        key = point['original_model_id'], point['normal_bits'], point['payoff_exp_bits']
        assert key not in old
        old[key] = (point, amended, 'expanded_raw_payoff/expanded_raw_payoff_resources.json.gz',
                    'expanded_raw_payoff/expanded_raw_payoff_resources.json.gz')
    assert len(old) == 96
    return old, expanded['selected_comparisons']

def register(args):
    repo, out, raw, smooth, expanded = configure(args)
    old, previous = known_records(repo)
    threshold = read(repo/'results/shared_threshold_resources.json')
    audit = read(repo/'results/raw_payoff_smoothing/raw_payoff_smoothing_audit.json')
    assert audit['status'] == 'passed'
    bindings = []
    for mid, bn, be in sorted(old):
        if be not in (8, 10, 11, 12):
            continue
        point, amended, pp, ap = old[mid, bn, be]
        bindings.append({'model_id': mid, 'normal_bits': bn, 'PG_exp_bits': be,
            'PA_exp_bits': be, 'pointwise_source': pp, 'amended_source': ap,
            'pointwise_record_sha256': digest(point), 'amended_record_sha256': digest(amended),
            'graph_sha256': digest(point['graph_nodes']),
            'configuration_sha256': digest(point['graph_configuration']),
            'profiles_sha256': digest(point['profiles'])})
    assert len(bindings) == 48
    protocol = {'schema': 'qbp-asymmetric-raw-payoff-protocol-v1',
        'registration_utc': datetime.now(timezone.utc).isoformat(),
        'research_record_date': '2026-10-10', **signatures(args, raw, repo),
        'registered_before_novel_scalar_and_profile_outcomes': True,
        'known_previous_outcomes': {'diagonal_pool_configurations': 96,
            'diagonal_exp_widths': [4, 6, 8, 10, 11, 12, 13, 14],
            'selected_comparisons': previous,
            'disclosure': 'All prior 96 raw-payoff settings, the smoothing amendment, numerical conditional-price and shared-threshold outcomes are known. This is a disclosed sequential asymmetric precision challenge, not a blind initial registration.'},
        'pairs': [[r['original_model_id'], r['candidate_id']] for r in threshold['rows']],
        'inner_normal_bits': [8, 10, 12], 'PG_exp_address_bits': [8, 9, 10, 11, 12],
        'PA_exp_address_bits': [8, 9, 10, 11, 12],
        'architectures': ['exact7t', 'unitary8t'],
        'configurations': 300, 'reused_diagonal_configurations': 48,
        'novel_configurations': 252, 'reused_diagonal_bindings': bindings,
        'q': 28, 'outer_quantile_bits': 10, 'outer_scalar_table_bits': 12,
        'outer_L': '6', 'inner_normal_L': '6', 'K': '100', 'H': '1',
        'payoff_exp_span': '1/64', 'payoff_exp_output_range': ['1', '64/63'],
        'representation': 'Retain the exact common-normal graph and final clamp(positive_sub(PA,PG),1). Independent PG/PA exp tables may cross even at ordered represented endpoints, so the unprojected represented difference is not asserted nonnegative. The exact common-normal payoff difference remains in[0,1] and the projection is1-Lipschitz.',
        'error_rule': 'eta=min(pointwise,min(1,100(e_exp_G+e_exp_A)+202E_s+2D_N+(nu0/lambda)(2(e_m+e_C+e_O)+e_Delta))). Use the separate intrinsic error of each exp map. Quantile rounding is charged through D_N only in the smoothing term.',
        'graph_policy': 'For novel settings build the retained raw graph, set PG_exp and PA_exp address widths independently and retain all other nodes, guards, unsigned widths and projection. Same-width reused graphs/configurations/profiles are referenced byte-identically.',
        'scalar_policy': 'Every novel setting receives freshly outward certified scalar ordinates and packed table data using the retained exact scalar constructor. Identical certificates and word widths are deduplicated by exact hash, never by approximate cost.',
        'resource_policy': 'Every arithmetic-admitted novel graph receives fresh counted primitive composition in both architectures, including cleanup, comparator, state preparation, reflection, Fourier circuit and reset costs. Build both pointwise and minimum pointwise/smoothing price priors and full canonical QPE plans.',
        'financial_policy': 'Same established synthetic no-arbitrage model, factor basis, finite laws, smooth outer transfer, coefficient and price-independent ideal mean information. Monetary radius1/100 and alpha1/3600. No finance price answers supplied.',
        'selection': 'For each model minimise complete canonical QPE T cap over exactly the300 asymmetric settings and both architectures. Tie-break by normal bits, PG bits, PA bits and architecture. Then compare shared threshold to the minimum of this winner, the retained96 diagonal winner and the strengthened numerical rival. Retain negative or zero savings.',
        'failure_policy': 'Retain every production failure and exhausted deterministic budget. No adaptive repair, silent exclusion or pool widening after new outcomes.',
        'boundary_policy': 'Report winning coordinates at8 or12 as boundaries of this separately registered finite asymmetric pool. No precision global optimum or exclusion of other asymmetric widths is claimed.',
        'data_origin': 'Synthetic exact rational models and certified scalar data only. No LSEG observations, paid data, physical quantum executions or quantum price samples.'}
    write(out/'asymmetric_raw_payoff_protocol.json', protocol)
    print('Registered300 settings before all252 novel outcomes', flush=True)

def build_graph(raw, model, domain, bn, bg, ba):
    g = raw.build_graph(model, domain, bn, bg)
    g.nodes[g.semantic_nodes['PG_exp']].bits = bg
    g.nodes[g.semantic_nodes['PA_exp']].bits = ba
    g.configuration.pop('payoff_exp_bits')
    g.configuration.update({'PG_exp_bits': bg, 'PA_exp_bits': ba,
        'positive_projection_retained': True,
        'pricing_representation': 'common_normal_raw_capped_payoff_difference_independent_exp_addresses'})
    # Address widths change table selectors only. Preserve the already checked
    # word width and all node output widths rather than refinalising the graph.
    return g

def make_amended(base, smooth, archive, nu0, frontend, basis):
    binding = {'model_id': base['original_model_id'], 'candidate_id': base['candidate_id'],
        'normal_bits': base['normal_bits'], 'PG_exp_bits': base['PG_exp_bits'],
        'PA_exp_bits': base['PA_exp_bits'], 'graph_sha256': digest(base['graph_nodes']),
        'configuration_sha256': digest(base['graph_configuration']),
        'profile_sha256': digest(base['profiles'])}
    if base.get('rejection') and not base.get('complete_deterministic_budget'):
        return {**binding, 'pointwise_row_sha256': digest(base),
                'rejection': copy.deepcopy(base['rejection']), 'strategies': []}
    term = smooth.mapping(base, archive, nu0)
    eta = F(term['selected_probability_error_upper'])
    budget = copy.deepcopy(base['complete_deterministic_budget'])
    A, D = F(budget['external_factor_upper']), F(term['inner_normal_CDF_error_upper'])
    det = F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*eta+F(budget['external_factor_radius'])
    budget.update(arithmetic_before_mass_upper=str(min(F(base['probability_certificate']['arithmetic_probability_error_upper']), F(term['smoothing_arithmetic_error_upper']))),
        total_deterministic_undiscounted_price_error_upper=str(det),
        auxiliary_components={'bias_upper': str(2*D), 'inner_normal_payoff_difference_bias_upper': str(2*D)},
        certificate_amendment='smoothing_first_uniform_conditional_independent_exp_tables')
    previous = basis['coefficient_informed_prior']
    prior = {'probability_lower': str(max(F(0), F(previous['finite_ideal_mean_lower'])-eta)),
        'probability_upper': str(min(F(1), F(previous['finite_ideal_mean_upper'])+eta)),
        'native_positive': True, 'finite_ideal_mean_lower': previous['finite_ideal_mean_lower'],
        'finite_ideal_mean_upper': previous['finite_ideal_mean_upper'],
        'previous_mean_information_sha256': digest(previous),
        'information_scope': 'Unchanged price-independent finite ideal R1 mean interval, enlarged by minimum pointwise/smoothing conditional projected-payoff error.'}
    cert = {'status': 'passed', 'configuration': base['graph_configuration'],
        'probability_error_upper': str(eta), 'monetary_value_error_upper': str(eta),
        'probability_actual_bounds': base['probability_certificate']['probability_actual_bounds'],
        'retained_pointwise_certificate_sha256': digest(base['probability_certificate']),
        'smoothing_derivation': term,
        'scope': 'Uniform conditional inner-law expectation error; two independently addressed exp errors survive the retained positive projection.'}
    row = {**binding, 'pointwise_row_sha256': digest(base), 'probability_certificate': cert,
        'complete_deterministic_budget': budget, 'coefficient_informed_prior': prior, 'strategies': []}
    if det >= F(1, 100):
        row['rejection'] = {'reason': 'Amended deterministic monetary tolerance exhausted', 'deterministic_error_upper': str(det)}
    else:
        for mode in ('exact7t', 'unitary8t'):
            try:
                result = frontend.informed_complete_plan('positive', budget, cert, prior, base['profiles'][mode])
            except (ArithmeticError, ValueError, AssertionError) as e:
                result = {'kind': 'production_failure', 'reason': type(e).__name__, 'message': str(e)}
            row['strategies'].append({'architecture': mode, 'result': result})
        valid = [s for s in row['strategies'] if s['result']['kind'] in ('complete_quantum_resource_plan', 'deterministic_correction_prior')]
        if valid:
            row['selected_strategy'] = min(valid, key=lambda s: (s['result']['T_count'], s['architecture']))
        else:
            row['rejection'] = {'reason': 'No valid amended architecture plan'}
    return row

def selected_summary(pool, previous, threshold, numerical):
    comparisons = []
    for mid in dict.fromkeys(r['model_id'] for r in pool):
        group = [r for r in pool if r['model_id'] == mid]
        admitted = [r for r in group if r.get('selected_strategy')]
        old = next(r for r in previous if r['model_id'] == mid)
        num = next(r['selected_numerical_strategy'] for r in numerical['selected_comparisons'] if r['original_model_id'] == mid)
        event = next(r['selected_strategy'] for r in threshold['rows'] if r['original_model_id'] == mid)
        row = {'model_id': mid, 'admitted_configurations': len(admitted),
            'rejected_configurations': len(group)-len(admitted),
            'previous96_raw_T': old['raw_T'], 'previous96_selected_normal_bits': old['selected_normal_bits'],
            'previous96_selected_exp_bits': old['selected_exp_bits'],
            'previous96_comparison_sha256': digest(old),
            'numerical_T': num['result']['T_count'], 'numerical_selected_strategy': num,
            'shared_threshold_T': event['result']['T_count'], 'threshold_selected_strategy': event}
        if admitted:
            win = min(admitted, key=lambda r: (r['selected_strategy']['result']['T_count'],
                r['normal_bits'], r['PG_exp_bits'], r['PA_exp_bits'], r['selected_strategy']['architecture']))
            T = win['selected_strategy']['result']['T_count']
            rivals = {'asymmetric300_raw_payoff': T, 'retained96_raw_payoff': old['raw_T'],
                      'strengthened_numerical_conditional_price': num['result']['T_count']}
            strongest = min(rivals, key=lambda k: (rivals[k], k))
            S, E = rivals[strongest], event['result']['T_count']
            row.update(selected_normal_bits=win['normal_bits'], selected_PG_exp_bits=win['PG_exp_bits'],
                selected_PA_exp_bits=win['PA_exp_bits'], selected_evidence_generation=win['evidence_generation'],
                selected_record_sha256=win['record_sha256'], selected_raw_strategy=win['selected_strategy'],
                raw_T=T, rival_T_counts=rivals, strongest_registered_rival=strongest,
                strongest_registered_rival_T=S, event_saving_against_asymmetric_exact=str(F(T-E, T)),
                event_saving_against_strongest_exact=str(F(S-E, S)),
                event_saving_against_strongest_per_cent_diagnostic=float(100*F(S-E, S)),
                selected_PG_boundary=win['PG_exp_bits'] in (8, 12),
                selected_PA_boundary=win['PA_exp_bits'] in (8, 12))
        comparisons.append(row)
    return comparisons

def run(args):
    repo, out, raw, smooth, expanded = configure(args)
    protocol_path = Path(args.protocol).resolve() if args.protocol else out/'asymmetric_raw_payoff_protocol.json'
    protocol = read(protocol_path)
    frozen = signatures(args, raw, repo)
    for k, v in frozen.items():
        assert protocol[k] == v, (k, 'Frozen input changed')
    assert protocol['inner_normal_bits'] == [8, 10, 12]
    assert protocol['PG_exp_address_bits'] == protocol['PA_exp_address_bits'] == [8, 9, 10, 11, 12]
    old, previous = known_records(repo)
    factor = read(repo/'results/factor_basis_quantum_resources.json')
    threshold = read(repo/'results/shared_threshold_resources.json')
    numerical = read(repo/'results/strengthened_r1_resources.json')
    norms = read(repo/'results/gaussian_derivative_norms.json')
    nu0 = F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order'] == 0 and r['passed']))
    rows, amended_rows, bindings, archive = [], [], [], []
    seen, cache, block_cache = set(), {}, {}
    start = time.monotonic()
    for mid, candidate in protocol['pairs']:
        prep = next(r for r in factor['preparations'] if r['original_model_id'] == mid and r['candidate_id'] == candidate)
        basis = next(r for r in factor['records'] if r['original_model_id'] == mid and r['factor_basis_candidate_id'] == candidate and r['kind'] == 'positive_r1')
        for bn, bg, ba in product(protocol['inner_normal_bits'], protocol['PG_exp_address_bits'], protocol['PA_exp_address_bits']):
            if bg == ba and (mid, bn, bg) in old:
                point, amended, pp, ap = old[mid, bn, bg]
                registered = next(r for r in protocol['reused_diagonal_bindings'] if (r['model_id'], r['normal_bits'], r['PG_exp_bits']) == (mid, bn, bg))
                assert registered['pointwise_record_sha256'] == digest(point)
                assert registered['amended_record_sha256'] == digest(amended)
                bind = {**registered, 'candidate_id': candidate, 'evidence_generation': 'retained_diagonal',
                        'record_sha256': digest(amended), 'admitted': bool(amended.get('selected_strategy'))}
                if amended.get('selected_strategy'):
                    bind['selected_strategy'] = amended['selected_strategy']
                else:
                    bind['rejection'] = copy.deepcopy(amended['rejection'])
                bindings.append(bind)
                print(mid, bn, bg, ba, 'retained_diagonal_T', bind.get('selected_strategy', {}).get('result', {}).get('T_count'), flush=True)
                continue
            point = {'original_model_id': mid, 'candidate_id': candidate, 'normal_bits': bn,
                'PG_exp_bits': bg, 'PA_exp_bits': ba, 'profiles': {}, 'strategies': []}
            index = {}
            try:
                g = build_graph(raw, prep['model'], prep['domain'], bn, bg, ba)
                tables, ordinates = raw.o.make_tables(g, cache, True)
                body = raw.c.certify_graph(g, tables, ordinates)
                cert = raw.probability_certificate(g, body, ordinates, nu0)
                cert.update(positivity='Final positive subtraction and clamp project onto[0,1]. Distinct exp maps are not asserted to preserve represented payoff ordering.',
                    decoded_observable='Projection of the represented difference of two capped payoffs on the same normal draw; the exact target difference is nonnegative.')
                point.update(model=prep['model'], graph_configuration=g.configuration,
                    graph_semantic_nodes=g.semantic_nodes,
                    graph_nodes=[{'name': n, 'kind': v.kind, 'args': list(v.args),
                        'meta': v.meta.as_dict() if isinstance(v.meta, raw.o.ScalarSpec) else str(v.meta),
                        'bits': v.bits, 'output_bits': v.output_bits,
                        'ideal_bounds': list(map(str, g.bounds[n]))} for n, v in g.nodes.items()],
                    probability_certificate=cert)
                for n, cc in ordinates.items():
                    key = digest(cc), g.width
                    if key not in seen:
                        seen.add(key)
                        table = tables[raw.o.table_key(n, g.nodes[n])]
                        archive.append({'certificate_sha256': key[0], 'certificate': cc,
                            'packed_table': table, 'packed_table_sha256': digest(table), 'graph_word_width': g.width})
                    if cc.get('name') == 'exp' and F(cc['lo']) == 0 and F(cc['span']) == F(1, 64):
                        index['exp', F(cc['lo']), F(cc['span']), cc['table_bits'], cc['q']] = cc
                if cert['status'] != 'passed':
                    point['rejection'] = {'reason': cert['reason']}
                else:
                    eta = F(cert['probability_error_upper'])
                    bud = copy.deepcopy(basis['complete_deterministic_budget'])
                    A, D = F(bud['external_factor_upper']), F(cert['inner_normal_CDF_error_upper'])
                    det = F(bud['coefficient_radius_upper'])+F(bud['financial_transfer_upper'])+A*eta+F(bud['external_factor_radius'])
                    bud.update(arithmetic_before_mass_upper=cert['arithmetic_probability_error_upper'],
                        auxiliary_components={'inner_normal_payoff_difference_bias_upper': str(2*D), 'bias_upper': str(2*D)},
                        total_deterministic_undiscounted_price_error_upper=str(det))
                    previous_prior = basis['coefficient_informed_prior']
                    prior = {'probability_lower': str(max(F(0), F(previous_prior['finite_ideal_mean_lower'])-eta)),
                        'probability_upper': str(min(F(1), F(previous_prior['finite_ideal_mean_upper'])+eta)),
                        'native_positive': True, 'finite_ideal_mean_lower': previous_prior['finite_ideal_mean_lower'],
                        'finite_ideal_mean_upper': previous_prior['finite_ideal_mean_upper'],
                        'previous_mean_information_sha256': digest(previous_prior),
                        'information_scope': 'Identical price-independent finite ideal R1 mean interval, enlarged by full projected raw-payoff pointwise conditional error.'}
                    point.update(complete_deterministic_budget=bud, coefficient_informed_prior=prior)
                    for mode in protocol['architectures']:
                        try:
                            profile = raw.o.counted_composition(g, tables, block_cache, mode)
                            point['profiles'][mode] = {'centred': profile}
                            result = ({'kind': 'deterministic_budget_exhausted', 'deterministic_error_upper': str(det)}
                                if det >= F(1, 100) else raw.informed.informed_complete_plan('positive', bud, cert, prior, point['profiles'][mode]))
                        except (ArithmeticError, ValueError, AssertionError) as e:
                            result = {'kind': 'production_failure', 'reason': type(e).__name__, 'message': str(e)}
                        point['strategies'].append({'architecture': mode, 'result': result})
                    valid = [s for s in point['strategies'] if s['result']['kind'] in ('complete_quantum_resource_plan', 'deterministic_correction_prior')]
                    if valid:
                        point['selected_strategy'] = min(valid, key=lambda s: (s['result']['T_count'], s['architecture']))
                    else:
                        point['rejection'] = {'reason': 'Pointwise budget exhausted or no valid emitted architecture plan'}
            except (ArithmeticError, ValueError, AssertionError) as e:
                point['rejection'] = {'reason': 'Explicit production failure', 'type': type(e).__name__, 'message': str(e)}
            rows.append(point)
            if point.get('graph_nodes'):
                try:
                    amended = make_amended(point, smooth, index, nu0, raw.informed, basis)
                except (ArithmeticError, ValueError, AssertionError) as e:
                    amended = {'model_id': mid, 'candidate_id': candidate, 'normal_bits': bn, 'strategies': [],
                        'pointwise_row_sha256': digest(point), 'rejection': {'reason': 'Explicit smoothing production failure', 'type': type(e).__name__, 'message': str(e)}}
            else:
                amended = {'model_id': mid, 'candidate_id': candidate, 'normal_bits': bn, 'strategies': [],
                    'pointwise_row_sha256': digest(point), 'rejection': copy.deepcopy(point['rejection'])}
            amended.update(PG_exp_bits=bg, PA_exp_bits=ba)
            amended_rows.append(amended)
            bind = {'model_id': mid, 'candidate_id': candidate, 'normal_bits': bn,
                'PG_exp_bits': bg, 'PA_exp_bits': ba, 'evidence_generation': 'novel_asymmetric',
                'pointwise_record_sha256': digest(point), 'record_sha256': digest(amended),
                'admitted': bool(amended.get('selected_strategy'))}
            if amended.get('selected_strategy'):
                bind['selected_strategy'] = amended['selected_strategy']
            else:
                bind['rejection'] = copy.deepcopy(amended['rejection'])
            bindings.append(bind)
            print(mid, bn, bg, ba, 'pointwise_T', point.get('selected_strategy', {}).get('result', {}).get('T_count'),
                'amended_T', amended.get('selected_strategy', {}).get('result', {}).get('T_count'), 'rejection', amended.get('rejection'), flush=True)
    assert len(rows) == len(amended_rows) == 252 and len(bindings) == 300
    assert len({(r['model_id'], r['normal_bits'], r['PG_exp_bits'], r['PA_exp_bits']) for r in bindings}) == 300
    archive_path = out/'asymmetric_raw_payoff_scalar_ordinates.json.gz'
    write(archive_path, archive)
    comparisons = selected_summary(bindings, previous, threshold, numerical)
    for k, v in signatures(args, raw, repo).items():
        assert protocol[k] == v, (k, 'Frozen input changed during execution')
    write(out/'asymmetric_raw_payoff_resources.json.gz', {'schema': 'qbp-asymmetric-raw-payoff-resources-v1',
        'status': 'completed retaining every pointwise and smoothing outcome', **frozen,
        'protocol_sha256': sha(protocol_path), 'configurations': 300,
        'reused_diagonal_configurations': 48, 'novel_configurations': 252,
        'records': rows, 'amended_records': amended_rows, 'combined_pool_bindings': bindings,
        'selected_comparisons': comparisons, 'scalar_archive_file': archive_path.name,
        'scalar_archive_sha256': sha(archive_path), 'unique_scalar_certificates_and_wordwidths': len(archive),
        'emitted_block_profiles': len(block_cache), 'generation_seconds_diagnostic': time.monotonic()-start,
        'scope': 'Finite sequential asymmetric precision challenge, checked sufficient monetary error and complete canonical QPE logical caps only. No global optimum, hardware quantum prices, practical advantage or theoretical priority.'})
    for r in comparisons:
        print('Pool300', r['model_id'], 'winner', r.get('selected_normal_bits'), r.get('selected_PG_exp_bits'),
            r.get('selected_PA_exp_bits'), 'raw_T', r.get('raw_T'), 'strongest_T', r.get('strongest_registered_rival_T'),
            'event_saving_per_cent', r.get('event_saving_against_strongest_per_cent_diagnostic'), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=('register', 'run', 'replay'))
    parser.add_argument('--repo', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--proof', required=True)
    parser.add_argument('--protocol')
    args = parser.parse_args()
    if args.mode == 'register':
        if args.protocol:
            raise ValueError('Registration writes its protocol into output')
        register(args)
    else:
        if args.mode == 'replay' and not args.protocol:
            raise ValueError('Replay requires the exact frozen protocol')
        run(args)
