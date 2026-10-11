"""Independent portable resource audit of the expanded raw-payoff pool.

The producer and all confidence/charging frontends are excluded. Established
low-level financial node emitters are reused, explicitly, to emit fresh
chronological blocks and winning preparations. Exact retained data are used.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from math import comb, isqrt
import argparse, gzip, hashlib, json, sys, time
sys.dont_write_bytecode = True
PI = F(355, 113)
POST = F(1, 10**60)
EPS = F(1, 100)
ALPHA = F(1, 3600)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def canonical(v):
    return json.dumps(v, sort_keys=True, separators=(',', ':'), default=str).encode()

def digest(v):
    return hashlib.sha256(canonical(v)).hexdigest()

def read(p):
    p = Path(p)
    return json.load(gzip.open(p, 'rt')) if p.suffix == '.gz' else json.loads(p.read_text())

def sqrt_cap(x):
    y = x*(1 << 200)
    k = (y.numerator+y.denominator-1)//y.denominator
    r = isqrt(k)
    return F(r+(r*r < k), 1 << 100)

def binomial_failure(n, success):
    return sum((F(comb(n, j))*(1-success)**j*success**(n-j)
                for j in range((n+1)//2, n+1)), F(0))

def own_qft(bits, rotations, bounds):
    tape = []
    phase = norm = F(0)
    for j in range(bits//2):
        k = bits-1-j
        tape.extend([('CX', j, k), ('CX', k, j), ('CX', j, k)])
    for j in range(bits-1, -1, -1):
        for k in range(bits-1, j, -1):
            d = k-j
            if d >= 32:
                norm += PI/(1 << d)
                continue
            key = d+1
            minus, plus = rotations[f'-1:{key}'], rotations[f'1:{key}']
            r, s = minus['chronological_gates'], plus['chronological_gates']
            local = [(op, k) for op in r]+[(op, j) for op in r]+[('CX', k, j)]
            local += [(op, j) for op in s]+[('CX', k, j)]
            if not r and not s:
                local = []
            tape.extend(local)
            phase += -F(1, 1 << (d+2))+2*F(minus['phase_pi'])+F(plus['phase_pi'])
            norm += 2*bounds[f'-1:{key}']+bounds[f'1:{key}']
        tape.append(('H', j))
    return tape, phase, norm

def own_plan(plan, target, prior, qft, rotations, bounds):
    assert F(plan['epsilon']) == target and F(plan['alpha']) == ALPHA
    assert F(plan['prior_upper']) == prior
    var = prior*(1-prior) if prior <= F(1, 2) else F(1, 4)
    assert F(plan['variance_upper']) == var
    M, bits, reps = plan['M'], plan['phase_register_bits'], plan['repetitions']
    assert M == 1 << bits and reps > 0 and reps % 2 == 1
    error = 2*PI*sqrt_cap(var)/M+PI**2/M**2
    assert F(plan['amplitude_error_upper']) == error
    assert F(plan['postprocess_allowance']) == POST and error+POST <= target
    if M > 2:
        assert 2*PI*sqrt_cap(var)/(M//2)+PI**2/F((M//2)**2)+POST > target
    tape, phase, norm = own_qft(bits, rotations, bounds)
    counts = Counter(op[0] for op in tape)
    assert counts == Counter(qft['gate_counts'])
    assert qft['T_count'] == counts['T']+counts['Tdg']
    assert F(qft['norm_error_upper']) == norm
    assert F(qft['global_phase_correction_pi']) == phase
    assert qft['chronological_tape_sha256'] == hashlib.sha256(json.dumps(tuple(tape), separators=(',', ':')).encode()).hexdigest()
    scale = 1 << 80
    z = norm*scale
    rounded = F((z.numerator+z.denominator-1)//z.denominator, scale)
    assert F(plan['qft_norm_error_allowance']) == rounded
    success = 8/PI**2-rounded
    assert success > F(1, 2) and F(plan['good_probability_lower']) == success
    failure = binomial_failure(reps, success)
    assert F(plan['binomial_failure_upper']) == failure <= ALPHA
    if reps > 1:
        assert binomial_failure(reps-2, success) > ALPHA
    return error+POST, failure

def own_charge(charge, p, qft, plan):
    ap, n = Counter(p['primitive_counts']), p['reflection_core_qubits']
    ccx = 2*(n-1)
    refl = Counter({'X': 2*n, 'CZ': 2, 'Z': 1, 'T': 4*ccx,
                    'Tdg': 3*ccx, 'H': 2*ccx, 'CX': 6*ccx})
    r = charge['controlled_reflection']
    assert Counter(r['gate_counts']) == refl and r['T'] == 7*ccx
    assert r['logical_qubits'] == 2*n
    assert r['explicit_T_depth_serial_upper'] == 4*ccx
    assert r['conservative_5layer_T_depth_serial_upper'] == 5*ccx
    assert r['expanded_primitive_operations'] == sum(refl.values())
    assert r['global_minus_Z_charged'] and r['controlled_good_CZ_charged']
    adj = Counter()
    inverse = {'T': 'Tdg', 'Tdg': 'T', 'S': 'Sdg', 'Sdg': 'S'}
    for op, value in ap.items():
        adj[inverse.get(op, op)] += value
    Q = ap+adj+refl
    bits, M, reps = plan['phase_register_bits'], plan['M'], plan['repetitions']
    peak = max(p['logical_qubits_upper'], 2*n)+bits+qft.get('extra_ancillas', 0)
    per_counts = ap+Counter({op: (M-1)*v for op, v in Q.items()})+Counter(qft['gate_counts'])
    per_counts['H'] += bits
    at, qt = ap['T']+ap['Tdg'], Q['T']+Q['Tdg']
    per = {'T': at+(M-1)*qt+qft['T_count'], 'measurements': peak,
           'reset_X': peak, 'feedback_CZ': 0,
           'expanded_primitive_operations': sum(per_counts.values())+2*peak,
           'uniform_H': p['uniform_H']*(2*M-1)+bits,
           'explicit_T_depth_serial_upper': p['explicit_T_depth_serial_upper']*(2*M-1)+(M-1)*4*ccx+qft['T_count'],
           'conservative_5layer_T_depth_serial_upper': p['conservative_5layer_T_depth_serial_upper']*(2*M-1)+(M-1)*5*ccx+qft['T_count']}
    per.update({f'gate_{op}': value for op, value in per_counts.items()})
    assert Counter(charge['per_repetition_metrics']) == Counter(per)
    assert Counter(charge['metrics']) == Counter({k: reps*v for k, v in per.items()})
    assert charge['peak_logical_qubits'] == peak
    assert charge['complete_oracle_T'] == at and charge['controlled_Q_T'] == qt
    assert charge['controlled_Q_calls'] == reps*(M-1)
    assert charge['phase_register_measurements'] == reps*bits
    assert charge['complete_reinitialisations'] == reps

def archive_index(path):
    index = {}
    certificates = 0
    for a in read(path):
        c, w = a['certificate'], a['graph_word_width']
        assert digest(c) == a['certificate_sha256']
        if c['schema'] == 'qbp-integral-quantile-v1':
            key = ('quantile', c['b_y'], c['L'])
            own = tuple(n['offset_integer'] for n in c['nodes'])
        else:
            key = ('scalar', c['name'], c['lo'], c['span'], tuple(sorted(c['params'].items())), c['table_bits'])
            vals = [n['integer'] for n in c['nodes']]
            own = []
            for j, v in enumerate(vals):
                d = vals[j+1]-v if j+1 < len(vals) else 0
                assert 0 <= v < 1 << w and abs(d) < 1 << w
                own.append(v | (abs(d) << w) | ((d < 0) << (2*w)))
            own = tuple(own)
        assert own == tuple(a['packed_table']) and digest(own) == a['packed_table_sha256']
        if (key, w) in index:
            assert index[key, w] == (own, c)
        index[key, w] = own, c
        certificates += 1
    return index, certificates

def own_graph(row, e, index):
    g = e.IntegralGraph(28, 12)
    for j in range(row['model']['outer_dimension']):
        g.inp('iy'+str(j), 10)
    g.inp('inner_normal_label', row['normal_bits'])
    for n in row['graph_nodes']:
        kind, meta = n['kind'], n['meta']
        if kind == 'function':
            meta = e.ScalarSpec(meta['name'], F(meta['lo']), F(meta['span']),
                 tuple(sorted((k, F(v)) for k, v in meta['parameters'].items())),
                 tuple(map(F, meta['output_bounds'])))
        elif kind in ('constant', 'scale', 'clamp', 'clamp_signed', 'quantile'):
            meta = F(meta)
        elif meta == 'None':
            meta = None
        assert g.op(kind, *n['args'], meta=meta, bits=n['bits'], bounds=n['ideal_bounds']) == n['name']
    g.semantic_nodes = row['graph_semantic_nodes']
    g.configuration = dict(row['graph_configuration'])
    g.finish(g.semantic_nodes['probability'])
    assert g.configuration == row['graph_configuration']
    assert all(g.nodes[n['name']].output_bits == n['output_bits'] for n in row['graph_nodes'])
    tables = {}
    for name, n in g.nodes.items():
        if n.kind == 'quantile':
            key = ('quantile', n.bits, str(n.meta))
        elif n.kind == 'function':
            key = ('scalar', n.meta.name, str(n.meta.lo), str(n.meta.span),
                   tuple(sorted((k, str(v)) for k, v in n.meta.parameters)), n.bits)
        else:
            continue
        tables[e.table_key(name, n)] = index[key, g.width][0]
    return g, tables

def profile_check(g, tables, p, e, cache):
    assert p['configuration'] == g.configuration
    assert p['integer_word_width'] == g.width and p['integer_fractional_bits'] == 28
    assert p['comparison_bits'] == 28 and p['maximum_probability_integer'] == 1 << 28
    assert p['reflection_core_qubits'] == sum(g.inputs.values())+29
    assert p['uniform_H'] == sum(g.inputs.values())+28 and p['endpoint_high_zero_bits'] == 1
    assert [b['node'] for b in p['blocks']] == list(g.nodes)
    fw, raw = Counter(), Counter()
    depth = five = all_depth = lookups = 0
    retained = sum(g.inputs.values())+30
    peak = retained
    for b, (name, node) in zip(p['blocks'], g.nodes.items()):
        widths = tuple(g.inputs[a] if a in g.inputs else g.nodes[a].output_bits for a in node.args)
        distinct = tuple(dict.fromkeys(node.args))
        aliases = tuple(distinct.index(a) for a in node.args)
        table = tables.get(e.table_key(name, node)) if node.kind in ('function', 'quantile') else None
        table_hash = hashlib.sha256(json.dumps(table, separators=(',', ':')).encode()).hexdigest() if table is not None else None
        assert b['table_sha256'] == table_hash
        assert b['argument_widths'] == list(widths) and b['argument_alias_pattern'] == list(aliases)
        assert b['kind'] == node.kind and b['output_bits'] == node.output_bits
        assert b['function'] == (node.meta.as_dict() if node.kind == 'function' else None)
        meta = node.meta.as_dict() if node.kind == 'function' else str(node.meta)
        key = (p['architecture'] if table is not None else 'shared', node.kind,
               json.dumps(meta, sort_keys=True), node.bits, node.output_bits,
               g.q, g.width, widths, aliases, table_hash)
        if key not in cache:
            circuit = e.old.vo.ValueCircuit(store_gates=False, mode=p['architecture'])
            regs = [circuit.alloc(widths[node.args.index(a)]) for a in distinct]
            incoming = sum(map(len, regs))
            args = [regs[i] for i in aliases]
            out = e.old.emit_node(circuit, node, args, g.q, g.width, table)
            m = circuit.resources()
            assert len(out) == node.output_bits and m['active_wires_at_end']-incoming == len(out)
            cache[key] = m, m['wires']-incoming
        m, extra = cache[key]
        counts = Counter(m['clifford_T_expanded_gate_counts'])
        assert counts == Counter(b['primitive_counts'])
        assert b['logical_T_count'] == m['logical_T_count']
        assert b['explicit_T_depth_serial_upper'] == m['explicit_7T_sequence_T_depth']
        assert b['conservative_5layer_T_depth_serial_upper'] == m['conservative_5layer_CCX_T_depth_upper']
        assert b['lookup_invocations'] == m['lookup_macro_invocations']
        assert b['extra_qubits_upper'] == extra
        fw.update(counts)
        raw.update(m['primitive_gate_counts'])
        depth += m['explicit_7T_sequence_T_depth']
        five += m['conservative_5layer_CCX_T_depth_upper']
        all_depth += m['explicit_Clifford_T_all_gate_dependency_depth']
        lookups += m['lookup_macro_invocations']
        peak = max(peak, retained+extra)
        retained += node.output_bits
    comp = e.old.vo.ValueCircuit(store_gates=False, mode=p['architecture'])
    a, b = comp.alloc(29), comp.alloc(29)
    flag = comp.alloc(1)[0]
    comp.compare_into(a, b, flag)
    cm = comp.resources()
    assert Counter(cm['clifford_T_expanded_gate_counts']) == Counter(p['comparator_primitive_counts'])
    full = Counter({op: 2*v for op, v in fw.items()})+Counter(cm['clifford_T_expanded_gate_counts'])
    full['H'] += p['uniform_H']
    assert Counter(p['primitive_counts_forward_body']) == fw and Counter(p['pre_expansion_forward_counts']) == raw
    assert Counter(p['primitive_counts']) == full
    assert p['logical_T_count'] == full['T']+full['Tdg']
    assert p['explicit_T_depth_serial_upper'] == 2*depth+cm['explicit_7T_sequence_T_depth']
    assert p['conservative_5layer_T_depth_serial_upper'] == 2*five+cm['conservative_5layer_CCX_T_depth_upper']
    assert p['explicit_all_gate_depth_serial_upper'] == 1+2*all_depth+cm['explicit_Clifford_T_all_gate_dependency_depth']
    assert p['lookup_invocations_forward_and_cleanup'] == 2*lookups
    peak = max(peak, retained+cm['wires']-59)
    assert p['logical_qubits_upper'] == peak
    return len(cache)

def whole_preparation(row, g, tables, e, qg):
    bindings = []
    for mode in ('exact7t', 'unitary8t'):
        circuit, ins, compare, flag = e.compile_flag(g, tables, mode, store_gates=False)
        metrics = circuit.resources()
        p = row['profiles'][mode]['centred']
        counts = Counter(metrics['clifford_T_expanded_gate_counts'])
        counts['H'] += len(circuit.uniform_preparation_gates)
        assert counts == Counter(p['primitive_counts'])
        core = tuple(w for reg in ins.values() for w in reg)+tuple(compare)+(flag,)
        assert core == circuit.reflection_core and len(core) == p['reflection_core_qubits']
        excluded = flag-1
        assert excluded not in core and len(compare) == 28
        assert metrics['wires'] <= p['logical_qubits_upper']
        work = tuple(w for w in range(metrics['wires']) if w not in core and w != excluded)[:len(core)-1]
        assert len(work) == len(core)-1
        reflection = qg.expand_toffoli(qg.grover([], core, flag, work, control=metrics['wires']))
        ccx = 2*(len(core)-1)
        assert Counter(op[0] for op in reflection) == Counter({'X': 2*len(core), 'CZ': 2, 'Z': 1, 'T': 4*ccx, 'Tdg': 3*ccx, 'H': 2*ccx, 'CX': 6*ccx})
        bindings.append({'model_id': row['original_model_id'], 'normal_bits': row['normal_bits'],
           'exp_bits': row['payoff_exp_bits'], 'architecture': mode, 'flag': flag,
           'excluded_high_zero': excluded, 'reflection_core_qubits': len(core),
           'chronological_emitted_peak': metrics['wires'], 'serial_profile_peak_upper': p['logical_qubits_upper'],
           'complete_primitive_counts_verified': True, 'actual_core_and_reflection_verified': True})
        print('Fresh whole preparation verified', row['original_model_id'], row['normal_bits'], row['payoff_exp_bits'], mode, flush=True)
    return bindings

def budget_and_plan(row, profiles, baseline, rotations, bounds):
    """Financial decode and probability prior from the exact retained fields."""
    if 'complete_deterministic_budget' not in row:
        assert row.get('rejection') and not row.get('selected_strategy')
        return 0, 0, 1
    budget, cert, prior = row['complete_deterministic_budget'], row['probability_certificate'], row['coefficient_informed_prior']
    old_budget, old_prior = baseline['complete_deterministic_budget'], baseline['coefficient_informed_prior']
    for k in ('coefficient_interval', 'coefficient_radius_upper', 'mass_interval',
              'external_factor_midpoint', 'external_factor_upper', 'external_factor_radius',
              'financial_transfer_upper', 'monetary_probability_decode_factor'):
        assert budget[k] == old_budget[k]
    assert F(budget['monetary_probability_decode_factor']) == 1
    cfg = cert['configuration']
    assert F(cfg['monetary_signal_normaliser']) == F(cfg['positive_R1_normaliser']) == 1
    assert cfg['output_encoding'] == 'positive' and F(cfg['cap']) == 1
    assert F(cfg['K']) == 100 and F(cfg['inner_normal_L']) == 6
    eta, A = F(cert['probability_error_upper']), F(budget['external_factor_upper'])
    assert 0 <= eta <= 1 and A == F(budget['mass_interval']['upper_exact_fraction'])
    alo = F(budget['mass_interval']['lower_exact_fraction'])
    assert F(budget['external_factor_midpoint']) == (A+alo)/2
    assert F(budget['external_factor_radius']) == (A-alo)/2
    cl, ch = (F(budget['coefficient_interval'][k]) for k in ('lower_exact_fraction', 'upper_exact_fraction'))
    assert F(budget['coefficient_radius_upper']) == (ch-cl)/2
    det = F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*eta+F(budget['external_factor_radius'])
    assert det == F(budget['total_deterministic_undiscounted_price_error_upper'])
    ideal_lo, ideal_hi = F(old_prior['finite_ideal_mean_lower']), F(old_prior['finite_ideal_mean_upper'])
    assert F(prior['finite_ideal_mean_lower']) == ideal_lo and F(prior['finite_ideal_mean_upper']) == ideal_hi
    assert prior['native_positive']
    assert F(prior['probability_lower']) == max(F(0), ideal_lo-eta)
    assert F(prior['probability_upper']) == min(F(1), ideal_hi+eta)
    # Original tight rows reference the basis prior, while retained smoothing
    # rows reference that same earlier information through a separate hash.
    assert prior['previous_mean_information_sha256'] == digest(old_prior)
    plans = zeros = 0
    valid = []
    for strategy in row['strategies']:
        result = strategy['result']
        if result['kind'] == 'deterministic_budget_exhausted':
            assert det >= EPS and F(result['deterministic_error_upper']) == det
            continue
        if result['kind'] == 'production_failure':
            assert row.get('rejection') or any(s['result']['kind'] in ('complete_quantum_resource_plan', 'deterministic_correction_prior') for s in row['strategies'])
            continue
        if result['kind'] == 'deterministic_correction_prior':
            radius = A*(F(prior['probability_upper'])-F(prior['probability_lower']))/2
            assert result['T_count'] == 0 and F(result['total_price_error_upper']) == det+radius <= EPS
            zeros += 1
            valid.append(strategy)
            continue
        assert result['kind'] == 'complete_quantum_resource_plan' and det < EPS and len(result['flags']) == 1
        profile = profiles[strategy['architecture']]['centred']
        flag = result['flags'][0]
        assert flag['financial_profile_sha256'] == digest(profile)
        target, upper = (EPS-det)/A, F(prior['probability_upper'])
        assert F(flag['probability_tolerance']) == target and F(flag['prior_upper']) == upper
        statistical, failure = own_plan(flag['plan'], target, upper, flag['fourier_profile'], rotations, bounds)
        own_charge(flag['complete_charge'], profile, flag['fourier_profile'], flag['plan'])
        assert F(result['statistical_price_error_upper']) == A*statistical
        assert F(result['total_price_error_upper']) == det+A*statistical <= EPS
        assert F(result['failure_probability_upper']) == failure <= ALPHA
        assert result['T_count'] == flag['complete_charge']['metrics']['T']
        assert result['quantum_shots'] == flag['plan']['repetitions']
        assert result['logical_qubits_upper'] == flag['complete_charge']['peak_logical_qubits']
        plans += 1
        valid.append(strategy)
    if valid:
        assert row['selected_strategy'] == min(valid, key=lambda s: (s['result']['T_count'], s['architecture']))
        assert not row.get('rejection')
    else:
        assert row.get('rejection') and not row.get('selected_strategy')
        assert det >= EPS or all(s['result']['kind'] == 'production_failure' for s in row['strategies'])
    return plans, zeros, int(not valid)

def unary_selector_count(b):
    """Count visited internal prefixes, independently of the QROM emitter."""
    N = 1 << b
    def count(lo, hi):
        if lo >= N+1 or hi-lo == 1:
            return 0
        mid = (lo+hi)//2
        return 1+count(lo, mid)+count(mid, hi)
    return count(0, 1 << (b+1))

def high_width_bound(raw_rows, baseline_by_pair):
    """Conditional exclusion for b15..22 in this literal fixed architecture."""
    out = []
    for mid in dict.fromkeys(r['original_model_id'] for r in raw_rows):
        group = [r for r in raw_rows if r['original_model_id'] == mid]
        basis = baseline_by_pair[mid, group[0]['candidate_id']]
        budget, prior = basis['complete_deterministic_budget'], basis['coefficient_informed_prior']
        A = F(budget['external_factor_upper'])
        zero_eta_det = F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+F(budget['external_factor_radius'])
        target = (EPS-zero_eta_det)/A
        upper = F(prior['finite_ideal_mean_upper'])
        assert upper < F(1, 2)
        var = upper*(1-upper)
        M = 2
        while 2*PI*sqrt_cap(var)/M+PI**2/M**2+POST > target:
            M *= 2
        perfect_success = 8/PI**2
        reps = 1
        while binomial_failure(reps, perfect_success) > ALPHA:
            reps += 2
        assert reps == 25
        fixed = {}
        exp_components = []
        for row in group:
            sem = row['graph_semantic_nodes']
            exp_names = {sem['PG_exp'], sem['PA_exp']}
            assert len(exp_names) == 2
            for mode, enc in row['profiles'].items():
                p = enc['centred']
                blocks = [b for b in p['blocks'] if b['node'] in exp_names]
                assert len(blocks) == 2 and sum(b['lookup_invocations'] for b in blocks) == 4
                nonexp = p['logical_T_count']-2*sum(b['logical_T_count'] for b in blocks)
                key = row['normal_bits'], mode
                if key in fixed:
                    assert fixed[key] == nonexp
                fixed[key] = nonexp
                b = row['payoff_exp_bits']
                selectors = unary_selector_count(b)
                assert selectors == (1 << b)+b
                lookup_T = (112 if mode == 'exact7t' else 64)*selectors
                exp_T = 2*sum(z['logical_T_count'] for z in blocks)
                assert exp_T >= lookup_T
                exp_components.append({'normal_bits': row['normal_bits'], 'exp_bits': b,
                    'architecture': mode, 'nonexp_A_T': nonexp,
                    'exp_A_T': exp_T, 'eight_lookup_T': lookup_T,
                    'exp_nonlookup_A_T': exp_T-lookup_T,
                    'interpolation_remainder_bits': 22-b,
                    'A_T': p['logical_T_count']})
        assert set(fixed) == {(bn, mode) for bn in (8, 10, 12) for mode in ('exact7t', 'unitary8t')}
        floor = min(fixed.values())
        # Both modes exceed this common lower bound, ignoring every exp
        # interpolation operation, controlled reflection and Fourier T gate.
        selectors15 = unary_selector_count(15)
        assert selectors15 == 32768+15
        A_floor = floor+64*selectors15
        cap = reps*(2*M-1)*A_floor
        out.append({'model_id': mid, 'zero_eta_phase_length_lower': M,
          'perfect_qft_repetitions_lower': reps, 'minimum_fixed_nonexp_A_T': floor,
          'width15_selectors_per_lookup': selectors15,
          'width15_A_T_lower': A_floor, 'width15_complete_T_lower': cap,
          'excluded_widths': list(range(15, 23)), 'profile_decomposition': exp_components,
          'scope': 'Common PG/PA exp width; q28; span1/64; unchanged graph and outer widths; inner normal8/10/12; exact7t/unitary8t unary lookup; retained coefficient/mean certificates and canonical QPE only.'})
    return out

def main(args):
    start = time.monotonic()
    repo, newdir, oldraw, oldsmooth = map(Path, (args.repo, args.new_results, args.old_raw, args.old_smoothing))
    output = Path(args.output)
    assert not output.exists()
    frozen_paths = {'new_resource': newdir/'expanded_raw_payoff_resources.json.gz',
         'new_protocol': newdir/'expanded_raw_payoff_protocol.json',
         'new_archive': newdir/'expanded_raw_payoff_scalar_ordinates.json.gz',
         'new_producer': Path(args.producer_source),
         'old_resource': oldraw/'raw_payoff_resources.json.gz',
         'old_archive': oldraw/'raw_payoff_scalar_ordinates.json.gz',
         'old_smoothing': oldsmooth/'raw_payoff_smoothing_resources.json',
         'factor': repo/'results/factor_basis_quantum_resources.json',
         'numerical': repo/'results/strengthened_r1_resources.json',
         'threshold': repo/'results/shared_threshold_resources.json',
         'qft_rotations': repo/'results/qft_rotation_tapes.json',
         'qft_bounds': repo/'results/qft_rotation_certificates.json'}
    frozen = {k: sha(p) for k, p in frozen_paths.items()}
    data, protocol = read(frozen_paths['new_resource']), read(frozen_paths['new_protocol'])
    assert data['protocol_sha256'] == frozen['new_protocol']
    assert protocol['source_sha256'] == data['source_sha256'] == frozen['new_producer']
    assert protocol['new_narrow_exp_address_bits'] == [11, 12, 13, 14]
    assert protocol['combined_narrow_exp_address_bits'] == [4, 6, 8, 10, 11, 12, 13, 14]
    assert protocol['inner_normal_bits'] == [8, 10, 12]
    assert protocol['registered_before_new_scalar_and_profile_costs']
    assert data['new_configurations'] == 48 and data['combined_configurations'] == 96
    assert data['old_raw_resource_sha256'] == frozen['old_resource']
    assert data['old_raw_archive_sha256'] == frozen['old_archive']
    assert data['old_smoothing_resource_sha256'] == frozen['old_smoothing']
    assert data['scalar_archive_sha256'] == frozen['new_archive']
    dependencies = {}
    for name, h in data['dependencies_sha256'].items():
        assert sha(repo/name) == h
        dependencies[name] = h
    for name, h in data['inputs_sha256'].items():
        assert sha(repo/'results'/name) == h
        dependencies['results/'+name] = h
    assert data['raw_builder_source_sha256'] == sha(repo/'raw_payoff_oracle_tight.py')
    assert data['smoothing_source_sha256'] == sha(repo/'raw_payoff_smoothing_amendment.py')
    old_rows = read(frozen_paths['old_resource'])['records']
    old_amended = read(frozen_paths['old_smoothing'])['records']
    assert len(old_rows) == len(old_amended) == len(data['records']) == len(data['amended_records']) == 48
    factor = read(frozen_paths['factor'])
    baselines = {(r['original_model_id'], r['factor_basis_candidate_id']): r for r in factor['records'] if r['kind'] == 'positive_r1'}
    index, new_certs = archive_index(frozen_paths['new_archive'])
    old_index, old_certs = archive_index(frozen_paths['old_archive'])
    for k, value in old_index.items():
        if k in index:
            assert index[k] == value
        index[k] = value
    sys.path.insert(0, str(repo))
    import boundary_value_oracle as emitter
    import quantum_gates as gates
    rotations = read(frozen_paths['qft_rotations'])['rotations']
    bounds = {r['key']: F(r['operator_norm_upper_rational']) for r in read(frozen_paths['qft_bounds'])['certificates']}
    new_keys = {(r['original_model_id'], r['normal_bits'], r['payoff_exp_bits']) for r in data['records']}
    expected_new = {(mid, bn, be) for mid, _ in protocol['pairs'] for bn in (8, 10, 12) for be in (11, 12, 13, 14)}
    assert new_keys == expected_new and len(new_keys) == 48
    cache, graphs, profile_count = {}, {}, 0
    plans = zeros = rejects = 0
    # Re-emit all distinct blocks across every old/new circuit, then recompute
    # every new pointwise, new smoothing and old smoothing estimator plan.
    for generation, raw_rows, amended_rows in (('new48', data['records'], data['amended_records']), ('retained48', old_rows, old_amended)):
        for row, amended in zip(raw_rows, amended_rows):
            mid, bn, be = row['original_model_id'], row['normal_bits'], row['payoff_exp_bits']
            assert (amended['model_id'], amended['normal_bits'], amended['exp_bits']) == (mid, bn, be)
            assert amended['graph_sha256'] == digest(row['graph_nodes'])
            assert amended['configuration_sha256'] == digest(row['graph_configuration'])
            assert amended['profile_sha256'] == digest(row['profiles'])
            point_key = 'pointwise_row_sha256' if generation == 'new48' else 'retained_pointwise_row_sha256'
            assert amended[point_key] == digest(row)
            g, tables = own_graph(row, emitter, index)
            sem = g.semantic_nodes
            assert g.nodes[sem['upper_endpoint']].args == (sem['lower_endpoint'], sem['logT'])
            assert g.nodes[sem['PG_exp']].meta == g.nodes[sem['PA_exp']].meta
            assert g.nodes[sem['PG_exp']].meta.output_bounds == (F(1), F(64, 63))
            assert g.nodes[sem['PG_exp']].bits == g.nodes[sem['PA_exp']].bits == be
            graphs[generation, mid, bn, be] = g, tables, row
            for encoded in row['profiles'].values():
                assert set(encoded) == {'centred'}
                profile_check(g, tables, encoded['centred'], emitter, cache)
                profile_count += 1
            baseline = baselines[mid, row['candidate_id']]
            if generation == 'new48':
                c = budget_and_plan(row, row['profiles'], baseline, rotations, bounds)
                plans += c[0]; zeros += c[1]; rejects += c[2]
            c = budget_and_plan(amended, row['profiles'], baseline, rotations, bounds)
            plans += c[0]; zeros += c[1]; rejects += c[2]
            print('Independent profiles and full plans verified', generation, mid, bn, be, flush=True)
    all_amended = [(generation, row) for generation, rows in (('retained48', old_amended), ('new48', data['amended_records'])) for row in rows]
    expected_bindings = [{'model_id': r['model_id'], 'candidate_id': r['candidate_id'], 'normal_bits': r['normal_bits'], 'exp_bits': r['exp_bits'], 'evidence_generation': gen, 'record_sha256': digest(r), 'T_count': r.get('selected_strategy', {}).get('result', {}).get('T_count'), 'admitted': bool(r.get('selected_strategy'))} for gen, r in all_amended]
    assert data['combined_pool_bindings'] == expected_bindings
    numeric, threshold = read(frozen_paths['numerical']), read(frozen_paths['threshold'])
    bindings, minima = [], []
    for comp in data['selected_comparisons']:
        mid = comp['model_id']
        cohort = [(gen, r) for gen, r in all_amended if r['model_id'] == mid]
        valid = [(gen, r) for gen, r in cohort if r.get('selected_strategy')]
        assert len(cohort) == 24 and comp['admitted_configurations'] == len(valid)
        assert comp['rejected_configurations'] == 24-len(valid)
        gen, win = min(valid, key=lambda z: (z[1]['selected_strategy']['result']['T_count'], z[1]['normal_bits'], z[1]['exp_bits'], z[1]['selected_strategy']['architecture']))
        assert comp['selected_record_sha256'] == digest(win) and comp['selected_evidence_generation'] == gen
        assert comp['selected_raw_strategy'] == win['selected_strategy']
        assert (comp['selected_normal_bits'], comp['selected_exp_bits']) == (win['normal_bits'], win['exp_bits'])
        num = next(r['selected_numerical_strategy'] for r in numeric['selected_comparisons'] if r['original_model_id'] == mid)
        event = next(r['selected_strategy'] for r in threshold['rows'] if r['original_model_id'] == mid)
        assert comp['numerical_selected_strategy'] == num and comp['threshold_selected_strategy'] == event
        rawT, numT, eventT = win['selected_strategy']['result']['T_count'], num['result']['T_count'], event['result']['T_count']
        rival = min(rawT, numT)
        assert comp['raw_T'] == rawT and comp['numerical_T'] == numT and comp['shared_threshold_T'] == eventT
        assert comp['strongest_registered_rival_T'] == rival
        assert F(comp['event_saving_against_raw_exact']) == F(rawT-eventT, rawT)
        saving = F(rival-eventT, rival)
        assert F(comp['event_saving_against_strongest_exact']) == saving
        assert comp['exp_upper_boundary'] == (win['exp_bits'] == 14)
        g, tables, raw_row = graphs[gen, mid, win['normal_bits'], win['exp_bits']]
        bindings += whole_preparation(raw_row, g, tables, emitter, gates)
        flag = win['selected_strategy']['result']['flags'][0]
        minima.append({'model_id': mid, 'generation': gen, 'normal_bits': win['normal_bits'], 'exp_bits': win['exp_bits'], 'raw_T': rawT, 'numerical_T': numT, 'shared_threshold_T': eventT, 'strongest_rival_T': rival, 'saving_exact': str(saving), 'saving_per_cent_diagnostic': float(100*saving), 'M': flag['plan']['M'], 'repetitions': flag['plan']['repetitions'], 'upper_boundary': win['exp_bits'] == 14})
    exclusion = high_width_bound(data['records'], baselines)
    for bound in exclusion:
        rival = next(m['raw_T'] for m in minima if m['model_id'] == bound['model_id'])
        bound['selected_raw_T'] = rival
        assert bound['width15_complete_T_lower'] > rival
    assert {k: sha(p) for k, p in frozen_paths.items()} == frozen
    for name, h in dependencies.items():
        assert sha(repo/name) == h
    result = {'schema': 'qbp-expanded-precision-independent-resource-review-v1', 'status': 'passed',
      'audit_source_sha256': sha(__file__), 'frozen_input_sha256': frozen,
      'dependency_sha256': dependencies, 'new_configurations': 48,
      'combined_pool_configurations': 96, 'fresh_complete_profiles': profile_count,
      'fresh_distinct_literal_blocks': len(cache), 'own_formula_QPE_and_full_charge_plans': plans,
      'deterministic_prior_plans': zeros, 'rejected_estimator_configurations': rejects,
      'new_archive_certificates': new_certs, 'old_archive_certificates': old_certs,
      'winning_whole_preparations': bindings, 'independent_complete_minima': minima,
      'conditional_high_precision_exclusion': exclusion,
      'elapsed_seconds_diagnostic': time.monotonic()-start,
      'scope': 'Fresh distinct literal node emission for all 96 circuits; exact table-to-DAG/profile binding; monetary decode and unchanged mean-information prior; independent binary QPE minima, Fourier tapes/norm charges, exact binomial confidence, controlled phase/reflection, inverse, observations and full resets; 96-pool minimum and winning whole chronological preparations.',
      'limitations': ['Reuses established primitive emitter, with no new low-level implementation.', 'Mathematical scalar/finite-law/smoothing/outer transfer justification is separately audited.', 'Width15..22 exclusion is conditional on this precise common-width unary lookup, fixed graph and canonical QPE; widths1/2/3/5/7/9, separate PG/PA widths, other normal widths, state-preparation schemes and estimators remain outside this search.', 'Logical sufficient upper caps are compared. No physical outcome, quantum price estimate, runtime speedup, theoretical priority or global optimum is established.']}
    output.write_text(json.dumps(result, indent=2)+'\n')
    print('Expanded independent resource review passed', plans, 'full plans;', profile_count, 'fresh profiles;', len(bindings), 'whole winning preparations', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', required=True)
    parser.add_argument('--new-results', required=True)
    parser.add_argument('--old-raw', required=True)
    parser.add_argument('--old-smoothing', required=True)
    parser.add_argument('--producer-source', required=True)
    parser.add_argument('--output', required=True)
    main(parser.parse_args())
