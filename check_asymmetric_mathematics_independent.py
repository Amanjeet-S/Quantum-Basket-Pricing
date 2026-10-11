"""Independent retained-field/fresh-Arb asymmetric payoff certification.

No producer, retained graph constructor or smoothing mapping is imported.
Fresh transcendental node enclosures and exact rational propagation are
reconstructed from the retained mathematical contracts.
"""
import argparse, bisect, gzip, hashlib, json, sys, time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from flint import arb, ctx
sys.dont_write_bytecode = True
ctx.prec = 320

def read(p):
    p = Path(p)
    return json.load(gzip.open(p, 'rt')) if p.suffix == '.gz' else json.loads(p.read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def digest(v):
    return hashlib.sha256(json.dumps(v, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()

def ball(x):
    x = F(x)
    return arb(x.numerator)/x.denominator

def endpoint(a, upwards):
    m, e = map(int, (a.upper() if upwards else a.lower()).man_exp())
    v = F(m*(1 << e)) if e >= 0 else F(m, 1 << -e)
    scaled = v*(1 << 384)
    n = scaled.numerator//scaled.denominator
    if upwards and scaled.denominator != 1:
        n += 1
    return F(n, 1 << 384)

def interval(a):
    return endpoint(a, False), endpoint(a, True)

def nearest(x):
    x = F(x)+F(1, 2)
    return x.numerator//x.denominator

def phi(x):
    return (1+(x/arb(2).sqrt()).erf())/2

def coeff(x, q):
    x = F(x)
    if x.denominator & (x.denominator-1) == 0:
        return x.numerator, x.denominator.bit_length()-1
    return nearest(x*(1 << q)), q

def certificate_key(c):
    if 'L' in c and 'b_y' in c:
        return ('quantile', F(c['L']), c['b_y'], c['q'])
    return ('scalar', c['name'], F(c['lo']), F(c['span']),
            tuple(sorted(c['params'].items())), c['table_bits'], c['q'], F(c['output_offset']))

def node_key(n, q):
    if n['kind'] == 'quantile':
        return ('quantile', F(n['meta']), n['bits'], q)
    m = n['meta']
    return ('scalar', m['name'], F(m['lo']), F(m['span']),
            tuple(sorted(m['parameters'].items())), n['bits'], q, F(0))

def check_scalar(c):
    name, lo, span = c['name'], F(c['lo']), F(c['span'])
    assert name in ('exp', 'log') and not c['params'] and F(c['output_offset']) == 0
    q, b = c['q'], c['table_bits']
    s, N = 1 << q, 1 << b
    hi = lo+span
    fn = (lambda x: x.exp()) if name == 'exp' else (lambda x: x.log())
    if name == 'exp':
        low = endpoint(ball(lo).exp(), False)
        top = endpoint(ball(hi).exp(), True)
        first = second = top
    else:
        assert lo > 0
        low, top = endpoint(ball(lo).log(), False), endpoint(ball(hi).log(), True)
        first, second = 1/lo, 1/(lo*lo)
    assert list(map(F, c['exact_function_range'])) == [low, top]
    assert F(c['first_derivative_upper']) == first and F(c['second_derivative_upper']) == second
    assert len(c['nodes']) == N+1
    errors, vals = [], []
    for j, n in enumerate(c['nodes']):
        x = lo+span*F(j, N)
        left, right = interval(fn(ball(x)))
        val = max(0, nearest((left+right)*s/2))
        err = max(abs(F(val, s)-left), abs(F(val, s)-right))
        assert F(n['input']) == x
        assert F(n['lower']) == left and F(n['upper']) == right
        assert n['integer'] == val and F(n['error_upper']) == err
        errors.append(err); vals.append(val)
    nodal, interpolation, floor = max(errors), second*(span/N)**2/8, F(1, s)
    assert F(c['nodal_error_upper']) == nodal
    assert F(c['linear_interpolation_error_upper']) == interpolation
    assert F(c['interpolation_floor_error_upper']) == floor
    assert F(c['output_error_upper']) == nodal+interpolation+floor
    assert c['value_min_integer'] == min(vals) and c['value_max_integer'] == max(vals)
    assert c['integer_width'] == max(1, max(vals).bit_length())
    assert c['precision_bits'] == 320
    exponent = span.numerator.bit_length()-span.denominator.bit_length()
    assert span == F(2)**exponent
    assert c['fractional_remainder_bits'] == q+exponent-b
    return len(vals)

def check_quantile(c, nu0):
    L, b, q = F(c['L']), c['b_y'], c['q']
    N, s = 1 << b, 1 << q
    c0, a = phi(-ball(L)), phi(ball(L))-phi(-ball(L))
    al, au = interval(a)
    assert F(c['mass_interval']['lower_exact_fraction']) == al
    assert F(c['mass_interval']['upper_exact_fraction']) == au
    vals, errors = [], []
    assert len(c['nodes']) == N
    for j, n in enumerate(c['nodes']):
        midpoint = F(2*j+1, 2*N)
        prob = c0+ball(midpoint)*a
        yl, yu = interval((2*prob-1).erfinv()*arb(2).sqrt())
        value = nearest(((yl+yu)/2+L)*s)
        z = F(value, s)-L
        error = max(abs(z-yl), abs(z-yu))
        assert n['j'] == j and F(n['probability_midpoint']) == midpoint
        assert F(n['quantile_lower']) == yl and F(n['quantile_upper']) == yu
        assert n['offset_integer'] == value and F(n['actual_coordinate']) == z
        assert F(n['node_error_upper']) == error and -L <= z <= L
        vals.append(value); errors.append(error)
    assert vals == sorted(vals) and c['monotone'] is True
    e = max(errors)
    assert F(c['maximum_node_error_upper']) == e
    assert c['offset_integer_bits'] == max(1, max(vals).bit_length())
    D = (1-al)/2+al/F(2*N)+nu0*e
    # For a finite empirical CDF versus a continuous increasing CDF, its
    # largest absolute discrepancy occurs at a closed/left atom endpoint.
    count = Counter(vals); cumulative = 0; worst = F(0); atom_checks = 0
    for value, multiplicity in sorted(count.items()):
        z = F(value, s)-L
        pl, pu = interval(phi(ball(z)))
        for closed_count in (cumulative, cumulative+multiplicity):
            empirical = F(closed_count, N)
            worst = max(worst, abs(empirical-pl), abs(empirical-pu))
            assert max(abs(empirical-pl), abs(empirical-pu)) <= D
            atom_checks += 1
        cumulative += multiplicity
    assert cumulative == N
    return len(vals), atom_checks, worst

def body_replay(row, certificates):
    q = row['graph_configuration']['q']; s = 1 << q; ulp = F(1, s)
    graph = row['graph_nodes']; expected = row['probability_certificate']['body_arithmetic_certificate']
    assert expected['status'] == 'passed' and not expected['failures']
    records = {n['node']: n for n in expected['nodes']}
    assert len(records) == len(graph)
    actual, errors, ideals = {}, {}, {}
    width = graph[0]['output_bits']; assert all(n['output_bits'] == width for n in graph)
    limit = F(1 << width, s)
    def floor(x):
        y = x*s
        return F(y.numerator//y.denominator, s)
    for node in graph:
        name, kind = node['name'], node['kind']
        ideal = tuple(map(F, node['ideal_bounds'])); ideals[name] = ideal
        rr = records[name]; assert rr['kind'] == kind and tuple(map(F, rr['ideal_bounds'])) == ideal
        aa = [actual[x] for x in node['args']] if kind != 'quantile' else []
        ee = [errors[x] for x in node['args']] if kind != 'quantile' else []
        err = F(0)
        if kind == 'quantile':
            c = certificates[node_key(node, q)]
            vals = [n['offset_integer'] for n in c['nodes']]
            ran = F(min(vals), s), F(max(vals), s)
        elif kind == 'constant':
            val = F(nearest(F(node['meta'])*s), s)
            ran = val, val; err = abs(val-F(node['meta']))
        elif kind == 'add':
            ran = aa[0][0]+aa[1][0], aa[0][1]+aa[1][1]; err = sum(ee)
        elif kind == 'sub':
            ran = aa[0][0]-aa[1][1], aa[0][1]-aa[1][0]; err = sum(ee)
        elif kind == 'positive_sub':
            raw = aa[0][0]-aa[1][1], aa[0][1]-aa[1][0]
            assert max(map(abs, raw)) < limit
            ran = max(F(0), raw[0]), max(F(0), raw[1]); err = sum(ee)
        elif kind == 'scale':
            target = F(node['meta']); integer, drop = coeff(target, q)
            implemented = F(integer, 1 << drop)
            assert min(aa[0]) >= 0 and target >= 0
            ran = floor(aa[0][0]*implemented), floor(aa[0][1]*implemented)
            err = target*ee[0]+max(map(abs, aa[0]))*abs(implemented-target)+(ulp if drop else 0)
        elif kind == 'clamp':
            target = F(node['meta']); cap = F(nearest(target*s), s)
            ran = max(F(0), min(cap, aa[0][0])), max(F(0), min(cap, aa[0][1]))
            err = ee[0]+abs(cap-target)
        elif kind == 'function':
            c = certificates[node_key(node, q)]
            assert 0 <= aa[0][0] <= aa[0][1] <= F(node['meta']['span'])
            assert 0 <= c['value_min_integer'] <= c['value_max_integer'] < 1 << width
            ran = F(c['value_min_integer'], s), F(c['value_max_integer'], s)
            err = F(c['output_error_upper'])+F(c['first_derivative_upper'])*ee[0]
        else:
            raise AssertionError(('Unsupported unchecked graph kind', kind))
        ran = max(ran[0], ideal[0]-err), min(ran[1], ideal[1]+err)
        assert 0 <= ran[0] <= ran[1] < limit
        assert tuple(map(F, rr['actual_bounds'])) == ran
        assert F(rr['absolute_output_error_upper']) == err
        actual[name], errors[name] = ran, err
    out = row['graph_semantic_nodes']['probability']
    assert 0 <= actual[out][0] <= actual[out][1] <= 1
    assert F(expected['probability_error_upper']) == errors[out]
    assert list(map(F, expected['probability_actual_bounds'])) == list(actual[out])
    return {n['name']: n for n in graph}, errors, actual

def conditional_replay(row, amended, certificates, nu0):
    graph, errors, actual = body_replay(row, certificates)
    sem, config = row['graph_semantic_nodes'], row['graph_configuration']
    q, lam, L = config['q'], F(row['model']['lambda']), F(config['inner_normal_L'])
    assert q == 28 and L == 6 and lam > 0
    # Verify the exact signed target identification, while retaining the
    # separately guarded positive unsigned circuit construction.
    lower, upper = graph[sem['lower_endpoint']], graph[sem['upper_endpoint']]
    assert lower['kind'] == upper['kind'] == 'add'
    assert upper['args'] == [sem['lower_endpoint'], sem['logT']]
    scale = graph[lower['args'][1]]; shifted = graph[lower['args'][0]]
    assert scale['kind'] == 'scale' and F(scale['meta']) == lam
    assert scale['args'] == [sem['normal_plus_L']]
    assert shifted['kind'] == 'add' and shifted['args'][0] == sem['mG_offset']
    C, O = graph[shifted['args'][1]], graph[sem['common_offset']]
    assert C['kind'] == O['kind'] == 'constant'
    # Recover the ideal affine financial loading without the graph builder.
    # Quantile outputs are offset coordinates y_j+L; arithmetic constants
    # must cancel to the classical weighted geometric log drift exactly.
    weights = list(map(F, row['model']['weights']))
    loadings = [list(map(F, x)) for x in row['model']['loadings']]
    drifts = list(map(F, row['model']['log_drifts']))
    p = row['model']['outer_dimension']
    assert sum(weights) == 1 and min(weights) > 0
    rate = F(row['model']['risk_free_rate'])-F(row['model']['dividend_yield'])
    assert all(mu+(sum(x*x for x in rr)+lam*lam)/2 == rate for mu, rr in zip(drifts, loadings))
    rg = [sum(w*rr[j] for w, rr in zip(weights, loadings)) for j in range(p)]
    mg = sum(w*mu for w, mu in zip(weights, drifts))
    def affine(name):
        n = graph[name]
        if n['kind'] == 'constant':
            return F(n['meta']), [F(0)]*p
        if n['kind'] == 'quantile':
            assert n['args'][0].startswith('iy') and F(n['meta']) == L
            vector = [F(0)]*p; vector[int(n['args'][0][2:])] = F(1)
            return F(0), vector
        if n['kind'] == 'scale':
            c, vector = affine(n['args'][0]); factor = F(n['meta'])
            return c*factor, [v*factor for v in vector]
        assert n['kind'] in ('add', 'sub')
        ca, va = affine(n['args'][0]); cb, vb = affine(n['args'][1])
        sign = 1 if n['kind'] == 'add' else -1
        return ca+sign*cb, [a+sign*b for a, b in zip(va, vb)]
    affine_constant, affine_vector = affine(sem['mG_offset'])
    assert affine_vector == rg
    assert affine_constant+F(C['meta'])-F(O['meta'])+lam*L == mg-L*sum(rg)
    assert actual[sem['logT']][0] >= 0
    projection = graph[sem['probability']]
    assert projection['kind'] == 'clamp' and F(projection['meta']) == 1
    positive = graph[projection['args'][0]]
    assert positive['kind'] == 'positive_sub' and positive['args'] == [sem['PA'], sem['PG']]
    eexp = {}
    for side in ('PG', 'PA'):
        exp = graph[sem[side+'_exp']]
        assert exp['kind'] == 'function' and exp['meta']['name'] == 'exp'
        assert F(exp['meta']['lo']) == 0 and F(exp['meta']['span']) == F(1, 64)
        assert list(map(F, exp['meta']['output_bounds'])) == [1, F(64, 63)]
        assert exp['bits'] == row[side+'_exp_bits']
        cert = certificates[node_key(exp, q)]; eexp[side] = F(cert['output_error_upper'])
        assert all(a['integer'] <= b['integer'] for a, b in zip(cert['nodes'], cert['nodes'][1:]))
    quantile = certificates[node_key(graph[sem['normal_plus_L']], q)]
    a = F(quantile['mass_interval']['lower_exact_fraction']); eN = F(quantile['maximum_node_error_upper'])
    D = (1-a)/2+a/F(2*(1 << config['normal_bits']))+nu0*eN
    point = row['probability_certificate']
    assert F(point['inner_normal_CDF_error_upper']) == D
    assert F(point['inner_normal_payoff_difference_bias_upper']) == 2*D
    assert F(point['normal_quantile_rounding_error_upper']) == eN
    arith = errors[sem['probability']]
    assert F(point['arithmetic_probability_error_upper']) == arith
    assert F(point['probability_error_upper']) == min(F(1), arith+2*D)
    assert 'Identical monotone' not in point['positivity']
    im, drop = coeff(lam, q)
    difference = abs(F(im, 1 << drop)-lam); floor = F(1, 1 << q) if drop else F(0)
    Es = 2*L*difference+floor
    em, edelta, ec, eo = (errors[sem['mG_offset']], errors[sem['logT']], errors[C['name']], errors[O['name']])
    scalar_term = 100*(eexp['PG']+eexp['PA'])
    inner_term = 202*Es
    outer_term = nu0/lam*(2*(em+ec+eo)+edelta)
    smooth_arithmetic = scalar_term+inner_term+outer_term
    smooth = min(F(1), smooth_arithmetic+2*D)
    chosen = min(F(point['probability_error_upper']), smooth)
    ac = amended['probability_certificate']; t = ac['smoothing_derivation']
    expected = {'outer_m_error_upper': em, 'outer_logT_error_upper': edelta,
        'shared_shift_constant_error_upper': ec, 'shared_offset_constant_error_upper': eo,
        'implemented_inner_coefficient': F(im, 1 << drop), 'inner_coefficient_error_upper': difference,
        'inner_scale_floor_allowance': floor, 'inner_scale_error_upper': Es,
        'PG_intrinsic_exp_error_upper': eexp['PG'], 'PA_intrinsic_exp_error_upper': eexp['PA'],
        'scalar_payoff_error_upper': scalar_term, 'inner_scale_payoff_error_upper': inner_term,
        'outer_smoothed_payoff_error_upper': outer_term,
        'smoothing_arithmetic_error_upper': smooth_arithmetic,
        'inner_normal_CDF_error_upper': D, 'pointwise_probability_error_upper': F(point['probability_error_upper']),
        'smoothing_probability_error_upper': smooth, 'selected_probability_error_upper': chosen}
    # Field spelling is checked against the retained manifest below. No
    # producer output values are substituted into the recomputation.
    aliases = {'PG_intrinsic_exp_error_upper': 'PG_intrinsic_exp_approximation_upper',
        'PA_intrinsic_exp_error_upper': 'PA_intrinsic_exp_approximation_upper',
        'scalar_payoff_error_upper': 'raw_exp_output_contribution_upper',
        'inner_scale_payoff_error_upper': 'inner_scale_payoff_contribution_upper',
        'outer_smoothed_payoff_error_upper': 'smooth_outer_shift_contribution_upper'}
    for key, value in expected.items():
        actual_key = key if key in t else aliases.get(key)
        assert actual_key in t, ('Unknown reviewed component field', key, sorted(t))
        assert F(t[actual_key]) == value, (actual_key, row['original_model_id'], row['normal_bits'], row['PG_exp_bits'], row['PA_exp_bits'])
    assert t['inner_scale_node'] == scale['name'] and t['outer_constant_node'] == C['name']
    assert F(ac['probability_error_upper']) == chosen and F(ac['monetary_value_error_upper']) == chosen
    assert ac['retained_pointwise_certificate_sha256'] == digest(point)
    assert amended['pointwise_row_sha256'] == digest(row)
    budget = amended['complete_deterministic_budget']; old = row['complete_deterministic_budget']
    for key in ('coefficient_radius_upper', 'financial_transfer_upper', 'external_factor_upper', 'external_factor_radius'):
        assert budget[key] == old[key]
    monetary = F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+F(budget['external_factor_upper'])*chosen+F(budget['external_factor_radius'])
    assert F(budget['total_deterministic_undiscounted_price_error_upper']) == monetary
    assert F(budget['arithmetic_before_mass_upper']) == min(arith, smooth_arithmetic)
    assert F(budget['auxiliary_components']['bias_upper']) == 2*D
    prior = amended['coefficient_informed_prior']; oldprior = row['coefficient_informed_prior']
    assert F(prior['probability_lower']) == max(F(0), F(oldprior['finite_ideal_mean_lower'])-chosen)
    assert F(prior['probability_upper']) == min(F(1), F(oldprior['finite_ideal_mean_upper'])+chosen)
    assert prior['native_positive'] is True
    return {'model_id': row['original_model_id'], 'normal_bits': row['normal_bits'],
        'PG_exp_bits': row['PG_exp_bits'], 'PA_exp_bits': row['PA_exp_bits'],
        'checked_graph_nodes': len(graph), 'conditional_error_upper': str(chosen),
        'complete_deterministic_price_error_upper': str(monetary)}

def main(args):
    start = time.monotonic(); repo = Path(args.repo); directory = Path(args.results)
    protocol_path = directory/'asymmetric_raw_payoff_protocol.json'
    resource_path = directory/'asymmetric_raw_payoff_resources.json.gz'
    archive_path = directory/'asymmetric_raw_payoff_scalar_ordinates.json.gz'
    protocol, data, archive = read(protocol_path), read(resource_path), read(archive_path)
    assert data['protocol_sha256'] == sha(protocol_path)
    assert data['scalar_archive_sha256'] == sha(archive_path)
    assert protocol['configurations'] == 300 and protocol['novel_configurations'] == 252
    assert protocol['inner_normal_bits'] == [8, 10, 12]
    assert protocol['PG_exp_address_bits'] == protocol['PA_exp_address_bits'] == [8, 9, 10, 11, 12]
    norms = read(repo/'results/gaussian_derivative_norms.json')
    nu0 = F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order'] == 0 and r['passed']))
    # The stored Gaussian norm used 384-bit arithmetic; a 320-bit fresh
    # upper enclosure can be looser than that certified upper endpoint.
    # Use a separate 512-bit enclosure solely for this prerequisite check.
    ctx.prec = 512
    assert nu0 >= endpoint(1/(2*arb.pi()).sqrt(), True)
    ctx.prec = 320
    certificates, certificate_widths, fresh_nodes, atom_checks = {}, set(), 0, 0
    cdf_worst = F(0)
    for record in archive:
        c, width = record['certificate'], record['graph_word_width']
        assert record['certificate_sha256'] == digest(c)
        assert record['packed_table_sha256'] == digest(record['packed_table'])
        key = certificate_key(c)
        if key not in certificates:
            if key[0] == 'quantile':
                checked, checks, worst = check_quantile(c, nu0)
                atom_checks += checks; cdf_worst = max(cdf_worst, worst)
            else:
                checked = check_scalar(c)
            fresh_nodes += checked; certificates[key] = c
        else:
            assert certificates[key] == c
        if key[0] == 'quantile':
            packed = [n['offset_integer'] for n in c['nodes']]
        else:
            values = [n['integer'] for n in c['nodes']]; packed = []
            for i, v in enumerate(values):
                delta = values[i+1]-v if i+1 < len(values) else 0
                assert 0 <= v < 1 << width and abs(delta) < 1 << width
                packed.append(v | (abs(delta) << width) | (int(delta < 0) << (2*width)))
        assert list(record['packed_table']) == packed
        certificate_widths.add((key, width))
    assert len(data['records']) == len(data['amended_records']) == 252
    results, failures = [], []
    for row, amended in zip(data['records'], data['amended_records']):
        assert row['PG_exp_bits'] in range(8, 13) and row['PA_exp_bits'] in range(8, 13)
        if row['probability_certificate']['status'] != 'passed':
            failures.append({'model_id': row['original_model_id'], 'normal_bits': row['normal_bits'], 'PG_exp_bits': row['PG_exp_bits'], 'PA_exp_bits': row['PA_exp_bits'], 'rejection': row['rejection']})
            continue
        results.append(conditional_replay(row, amended, certificates, nu0))
    out = {'schema': 'qbp-asymmetric-independent-mathematics-v1', 'status': 'passed',
        'auditor_sha256': sha(__file__), 'protocol_sha256': sha(protocol_path),
        'resources_sha256': sha(resource_path), 'scalar_archive_sha256': sha(archive_path),
        'fresh_Arb_precision_bits': 320, 'fresh_scalar_and_quantile_nodes': fresh_nodes,
        'unique_mathematical_certificates': len(certificates),
        'checked_packed_certificate_width_pairs': len(certificate_widths),
        'checked_CDF_closed_and_left_atom_cases': atom_checks,
        'worst_fresh_atom_CDF_difference_upper': str(cdf_worst),
        'novel_conditional_certificates_recomputed': len(results),
        'checked_graph_nodes': sum(r['checked_graph_nodes'] for r in results),
        'retained_guard_failures': failures, 'records': results,
        'elapsed_seconds_diagnostic': time.monotonic()-start,
        'scope': 'Fresh Arb transcendental enclosures, exact scalar and packed-table fields, direct independent rational graph propagation and unsigned guards, atom-safe CDF discrepancy, separate asymmetric scalar errors, conservative smoothing conditional certificates, monetary budgets and prior enlargement. No producer, graph constructor or smoothing mapping imported. Full literal gate/estimator resource replay is a separate review.'}
    output = Path(args.output)
    if output.exists():
        raise FileExistsError('Preserve completed prior audit')
    output.write_text(json.dumps(out, indent=2)+'\n')
    print('Independent mathematics passed:', len(results), 'novel conditional certificates;', fresh_nodes, 'fresh nodes;', atom_checks, 'CDF atom checks')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', required=True); parser.add_argument('--results', required=True)
    parser.add_argument('--output', required=True)
    main(parser.parse_args())
