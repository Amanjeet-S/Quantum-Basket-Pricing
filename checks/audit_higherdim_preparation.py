"""Independent direct finite-law audit of separable preparation."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json, sys, time
from flint import arb, ctx
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import higherdim_preparation as h
from integral_value_oracle import model_parameters
from integral_value_certificate import quantile_certificate
from continuous_quadrature import frac_arb


def overlap(a, record):
    lo, hi = h.endpoints(a)
    rl, rh = F(record['lower_exact_fraction']), F(record['upper_exact_fraction'])
    if max(lo, rl) > min(hi, rh):
        raise AssertionError('Independent direct enclosure does not overlap retained bound')


def direct_moments(model, mgf):
    ctx.prec = 768
    w, d, s, mg, rg, lam = model_parameters(model)
    sums = [arb(0) for _ in range(6)]
    gap = arb(0)
    N = len(mgf.coordinates) ** len(rg)
    for y in product(mgf.coordinates, repeat=len(rg)):
        T = sum((frac_arb(wi) * frac_arb(di + sum(v * z for v, z in zip(si, y))).exp() for wi, di, si in zip(w, d, s)), arb(0))
        U = T - 1
        power = arb(1)
        for r in range(6):
            power *= U
            sums[r] += power
        level = frac_arb(mg + lam * lam / 2 + sum(v * z for v, z in zip(rg, y))).exp()
        gap += 100 * level * U
    return [a / N for a in sums], gap / N, N


def main():
    started = time.perf_counter()
    retained = json.loads((ROOT / 'results/higherdim_preparation.json').read_text())
    protocol = json.loads((ROOT / 'results/generalisation_protocol.json').read_text())
    models = protocol['models']
    by, q = retained['b_y'], retained['q']
    mgf = h.get_mgf(6, by, q)
    old = quantile_certificate(6, by, q)
    actual = tuple(F(row['actual_coordinate']) for row in old['nodes'])
    if actual != mgf.coordinates:
        raise AssertionError('Compressed symmetric construction changed an actual coordinate')
    domain_node_checks = 0
    vertex_checks = 0
    for record in retained['records']:
        domain = record['domain']
        dlo, dhi = map(F, domain['coordinate_bounds'])
        for node in old['nodes']:
            if not (dlo <= F(node['quantile_lower']) <= F(node['quantile_upper']) <= dhi and
                    dlo <= F(node['actual_coordinate']) <= dhi):
                raise AssertionError('Finite normaliser domain failed ideal/actual node coverage')
            domain_node_checks += 1
        model = next(x for x in models if x['model_id'] == record['model_id'])
        w, d, s, mg, rg, lam = model_parameters(model)
        ctx.prec = 768
        for vertex in domain['vertex_upper_enclosures']:
            y = tuple(dhi if sign else dlo for sign in vertex['signs'])
            value = sum((frac_arb(wi) * frac_arb(di + sum(v * z for v, z in zip(si, y))).exp() for wi, di, si in zip(w, d, s)), arb(0))
            if h.endpoints(value)[1] > F(vertex['T_upper']) or F(vertex['T_upper']) > F(domain['T_star']):
                raise AssertionError('Independent vertex enclosure exceeds the finite normaliser')
            vertex_checks += 1
    ctx.prec = 768
    rates = sorted(mgf.cache)
    # Shared construction regenerated just the nodes. Populate exact frozen
    # rate combinations independently from all retained moment models.
    from taylor_pricing import compositions
    allrates = set((F(0), F(1), F(-1), F(7, 5), F(-7, 5)))
    for model in models:
        if model['outer_dimension'] not in (3, 7):
            continue
        w, d, s, mg, rg, lam = model_parameters(model)
        for r in range(7):
            for alpha in compositions(r, len(w)):
                allrates.update(abs(sum(alpha[i] * s[i][j] for i in range(len(w)))) for j in range(len(rg)))
        allrates.update(abs(v) for v in rg)
        allrates.update(abs(a + b) for si in s for a, b in zip(rg, si))
    rates = sorted(allrates)
    selected = sorted(set(rates[::max(1, len(rates) // 48)] + rates[-4:]))
    rate_checks = []
    for rate in selected:
        direct = sum((frac_arb(rate * y).exp() for y in actual), arb(0)) / len(actual)
        enclosure = mgf(rate)
        dl, dh = h.endpoints(direct)
        el, eh = h.endpoints(enclosure)
        if not (el <= dl <= dh <= eh):
            # Two ball computations can share their exact answer while their
            # enclosure radii differ. Retained MGF bounds must nevertheless
            # enclose this much tighter high-precision independent direct sum.
            raise AssertionError('Series MGF did not contain the direct high-precision sum')
        rate_checks.append({'rate': str(rate), 'direct_interval': h.serial(direct), 'retained_interval': h.serial(enclosure)})
    weighted = h.EmpiricalMGF(actual, tuple(F(i % 5 - 2, len(actual)) for i in range(len(actual))))
    ctx.prec = 768
    weighted_checks = []
    for rate in (F(-3, 2), F(-1, 7), F(0), F(1, 3), F(3, 2)):
        direct = sum((frac_arb(F(i % 5 - 2, len(actual))) * frac_arb(rate * y).exp() for i, y in enumerate(actual)), arb(0))
        enclosure = weighted(rate)
        dl, dh = h.endpoints(direct)
        el, eh = h.endpoints(enclosure)
        if not el <= dl <= dh <= eh:
            raise AssertionError('Signed empirical MGF did not contain the independent direct sum')
        weighted_checks.append({'rate': str(rate), 'direct_interval': h.serial(direct), 'series_interval': h.serial(enclosure)})
    small = h.get_mgf(6, 1, 12)
    tensor_checks, labels = [], 0
    for model in models:
        if model['outer_dimension'] not in (1, 3, 7):
            continue
        factorised = h.product_grid_moments(model, 6, 1, 12, mgf=small)
        direct, gap, count = direct_moments(model, small)
        labels += count
        for r, value in enumerate(direct, 1):
            overlap(value, factorised['centred_U_moments'][str(r)])
        overlap(gap, factorised['correlated_underlying_gap_interval'])
        tensor_checks.append({'model_id': model['model_id'], 'cartesian_labels': count,
            'direct_U_moments': [h.serial(a) for a in direct], 'direct_underlying_gap': h.serial(gap)})
    previous = json.loads((ROOT / 'results/integral_finite_grid_moments.json').read_text())['records']
    old_moment_overlaps = 0
    for row in retained['records']:
        oldrow = next(x for x in previous if x['model_id'] == row['model_id'])
        for r in ('1', '2', '3'):
            a, b = row['moments']['centred_U_moments'][r], oldrow['centred_U_moments'][r]
            if max(F(a['lower_exact_fraction']), F(b['lower_exact_fraction'])) > min(F(a['upper_exact_fraction']), F(b['upper_exact_fraction'])):
                raise AssertionError('Historical exact finite moment failed overlap')
            old_moment_overlaps += 1
    oldcoeff = json.loads((ROOT / 'results/taylor_pricing.json').read_text())['configurations']
    coefficient_overlaps = 0
    for row in retained['records']:
        for case in row['cases']:
            if case['cap'] is None:
                continue
            match = next(x for x in oldcoeff if x['model_id'] == row['model_id'] and F(x['cap']) == F(case['cap']))
            summed = [F(0), F(0)]
            for order in match['orders'][:2]:
                summed[0] += F(order['coefficient_interval']['lower_exact_fraction'])
                summed[1] += F(order['coefficient_interval']['upper_exact_fraction'])
            new = case['coefficients']['coefficient']
            if max(summed[0], F(new['lower_exact_fraction'])) > min(summed[1], F(new['upper_exact_fraction'])):
                raise AssertionError('Geometric coefficient failed independent historical overlap')
            coefficient_overlaps += 1
    result = {'schema': 'qbp-higherdim-preparation-independent-audit-v1', 'status': 'passed',
        'exact_original_quantile_coordinates_checked': len(actual), 'direct_MGF_rate_checks': len(rate_checks),
        'ideal_and_actual_domain_node_checks': domain_node_checks, 'independent_domain_vertex_checks': vertex_checks,
        'signed_weighted_MGF_checks': len(weighted_checks), 'small_tensor_models': len(tensor_checks),
        'small_tensor_labels': labels, 'small_tensor_moment_checks': 6 * len(tensor_checks),
        'small_tensor_correlated_gap_checks': len(tensor_checks), 'historical_finite_moment_overlaps': old_moment_overlaps,
        'historical_coefficient_overlaps': coefficient_overlaps,
        'direct_precision_bits': 768, 'rate_checks': rate_checks, 'weighted_checks': weighted_checks,
        'tensor_checks': tensor_checks, 'source_sha256': h.sha(__file__),
        'implementation_sha256': h.sha(ROOT / 'higherdim_preparation.py'),
        'input_sha256': {p: h.sha(ROOT / p) for p in ['results/higherdim_preparation.json', 'results/generalisation_protocol.json',
            'results/integral_finite_grid_moments.json', 'results/taylor_pricing.json']},
        'seconds_diagnostic': time.perf_counter() - started}
    (ROOT / 'results/higherdim_preparation_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if isinstance(v, (str, int, float)) and
        (k.endswith('checks') or k.endswith('overlaps') or k in ('status', 'exact_original_quantile_coordinates_checked', 'small_tensor_labels'))}))


if __name__ == '__main__':
    main()
