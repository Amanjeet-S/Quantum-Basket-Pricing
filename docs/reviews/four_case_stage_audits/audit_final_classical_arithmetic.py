"""Independent exact trace, probability composition and synthetic input audit.

This auditor imports no production helper and replays no financial
transcendental scores. The separate full-score auditor supplies that check.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import hashlib, json, struct, time
from flint import arb, ctx

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE / 'classical'
RESULTS = ROOT / 'results'
REPO = Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def A(value):
    value = F(value)
    return arb(value.numerator) / value.denominator

def bounds(value):
    def end(v):
        m, e = map(int, v.man_exp())
        return F(m * (1 << e)) if e >= 0 else F(m, 1 << -e)
    return end(value.lower()), end(value.upper())

def pair(record):
    return F(record['lower_exact_fraction']), F(record['upper_exact_fraction'])

def contains(record, value):
    a, b = pair(record)
    lo, hi = bounds(value)
    assert a <= lo <= hi <= b

def binary(record, code):
    path = RESULTS / record['file']
    raw = path.read_bytes()
    assert len(raw) == record['bytes'] and sha(path) == record['sha256']
    values = tuple(v[0] for v in struct.iter_unpack('<' + code, raw))
    assert len(values) == record.get('count', record.get('coordinate_count'))
    return values

def completed(record, finite_lo, finite_hi, transfer):
    ml, mh = pair(transfer['retained_mass_interval'])
    beta = F(transfer['financial_grid_transfer_upper'])
    lo = max(F(0), ml * finite_lo - beta)
    hi = min(F(1), mh * finite_hi + beta + 1 - ml)
    assert pair(record['undiscounted_continuous_price_interval']) == (lo, hi)
    assert F(record['undiscounted_radius_upper']) == (hi - lo) / 2
    assert record['certifies_undiscounted_cent'] == ((hi - lo) / 2 <= F(1, 100))
    contains(record['discount_interval'], A(F(-3, 100)).exp())
    dl, dh = pair(record['discount_interval'])
    assert pair(record['discounted_continuous_price_interval']) == (dl * lo, dh * hi)
    assert F(record['discounted_radius_upper']) == (dh * hi - dl * lo) / 2

def main():
    start = time.perf_counter()
    ctx.prec = 1024
    protocol_path = RESULTS / 'four_case_classical_protocol.json'
    result_path = RESULTS / 'four_case_classical_pricing.json'
    protocol = json.loads(protocol_path.read_text())
    result = json.loads(result_path.read_text())
    assert sha(protocol_path) == result['protocol_sha256']
    assert sha(ROOT / 'four_case_classical_pricing.py') == protocol['producer_sha256'] == result['producer_sha256']
    assert result['input_sha256'] == protocol['input_sha256']
    for name, checksum in protocol['input_sha256'].items():
        assert sha(ROOT / name) == checksum
        assert sha(REPO / name) == checksum
    stress = json.loads((REPO / 'results/dispersion_separation_protocol.json').read_text())
    models = {m['model_id']: m for m in stress['models']}
    outcomes = json.loads((REPO / 'results/dispersion_separation_results.json').read_text())['records']
    expected = {r['model_id'] for r in outcomes if r['cap'] == '1' and r['separates_after_degree5']}
    assert set(protocol['case_ids']) == expected and len(expected) == 4
    assert [r['model_id'] for r in result['records']] == protocol['case_ids']
    assert protocol['cap'] == '1' and protocol['draws_per_case'] == 100000
    assert protocol['outward_score_bits'] == 40
    assert 'independent uniformly distributed' in protocol['confidence_premise']
    assert 'not physically proved' in protocol['confidence_premise']
    # Each low-ten-bit value has exactly 64 uint16 preimages. This checks
    # the deterministic push-forward, not the physical randomness premise.
    preimages = Counter(w & 1023 for w in range(1 << 16))
    assert set(preimages) == set(range(1024)) and set(preimages.values()) == {64}
    eps = F(protocol['statistical_epsilon'])
    assert eps == F(1, 150)
    failure = 2 * A(-2 * protocol['draws_per_case'] * eps * eps).exp()
    assert bounds(failure)[1] <= F(protocol['per_case_failure_allowance']) == F(1, 3600)
    assert 4 * F(protocol['per_case_failure_allowance']) == F(protocol['family_failure_allowance']) == F(1, 900)
    assert F(result['family_hoeffding_failure_allowance']) == F(1, 900)
    assert F(result['family_confidence_lower']) == F(899, 900)
    mass = (A(6) / A(2).sqrt()).erf()
    grid = result['gaussian_grid']
    contains(grid['mass_interval'], mass)
    mass_lo, mass_hi = pair(grid['mass_interval'])
    b1 = F(grid['rho1_mass_upper'])
    _, b2_hi = pair(grid['rho2_mass_interval'])
    assert b1 >= 0 and b2_hi >= 0 and len(grid['actual_cells']) == 1024
    norms = json.loads((RESULTS / 'boundary_derivative_norms.json').read_text())['records']
    rows = []
    counts = {'sample_scores': 0, 'captured_uint16_words': 0, 'exactsum_scores': 0,
              'n2_sample_score_lookup_links': 0, 'score_endpoints': 0,
              'word_mask_equalities': 0, 'retained_mass_intervals': 0,
              'continuous_price_compositions': 0, 'synthetic_cases': 0}
    for row in result['records']:
        model = models[row['model_id']]
        p = model['outer_dimension']
        assert row['assets'] == model['assets'] and row['outer_dimension'] == p
        assert F(row['lambda']) == F(model['lambda']) and F(row['dispersion_scale']) == F(model['dispersion_scale'])
        assert row['cap'] == '1'
        carry = F(model['risk_free_rate']) - F(model['dividend_yield'])
        assert model['spot'] == model['strike'] == 100 and F(model['maturity']) == 1
        for drift, loading in zip(model['log_drifts'], model['loadings']):
            assert F(drift) + (F(model['lambda']) ** 2 + sum(F(v) ** 2 for v in loading)) / 2 == carry
        counts['synthetic_cases'] += 1
        t = row['transfer_certificate']
        contains(t['retained_mass_interval'], mass ** p)
        counts['retained_mass_intervals'] += 1
        ml, mh = pair(t['retained_mass_interval'])
        assert F(t['bounded_tail_upper']) == 1 - ml
        assert F(t['rho1_mass_upper']) == b1 and t['rho2_mass_interval'] == grid['rho2_mass_interval']
        norm = next(n for n in norms if n['payoff'] == 'capped' and F(n['lambda']) == F(model['lambda']) and F(n['cap']) == 1)
        m1, m2 = F(norm['M1_upper']), F(norm['M2_upper'])
        assert F(t['global_M1_upper']) == m1 and F(t['global_M2_upper']) == m2
        beta = F(0)
        for j, derivative in enumerate(t['coordinate_derivative_bounds']):
            loading = max(abs(F(r[j])) for r in model['loadings'])
            d1, d2 = m1 * loading, (m1 + m2) * loading ** 2
            assert F(derivative['maximum_full_asset_loading']) == loading
            assert F(derivative['first_derivative_upper']) == d1 and F(derivative['second_derivative_upper']) == d2
            beta += d1 * b1 + d2 * b2_hi
        assert beta * mass_hi ** (p - 1) == F(t['financial_grid_transfer_upper'])
        det = row['complete_exactsum']
        det_lower = det_upper = None
        if det is not None:
            assert p == 1 and det['evaluations'] == 1024
            det_lower = binary(det['scores_lower'], 'Q')
            det_upper = binary(det['scores_upper'], 'Q')
            assert len(det_lower) == len(det_upper) == 1024
            assert all(0 <= l <= u <= 1 << 40 for l, u in zip(det_lower, det_upper))
            lo, hi = F(sum(det_lower), 1024 * (1 << 40)), F(sum(det_upper), 1024 * (1 << 40))
            assert pair(det['finite_grid_mean_interval']) == (lo, hi)
            completed(det, lo, hi, t)
            counts['exactsum_scores'] += 1024
            counts['continuous_price_compositions'] += 1
        mc = row['complete_monte_carlo']
        words = binary(mc['random_words'], 'H')
        lower, upper = binary(mc['scores_lower'], 'Q'), binary(mc['scores_upper'], 'Q')
        n = protocol['draws_per_case']
        assert len(words) == n * p and len(lower) == len(upper) == n
        assert all(0 <= l <= u <= 1 << 40 for l, u in zip(lower, upper))
        labels = tuple(w % 1024 for w in words)
        assert all(k == (w & 1023) for k, w in zip(labels, words))
        counts['word_mask_equalities'] += len(words)
        label_digest = hashlib.sha256(b''.join(struct.pack('<H', k) for k in labels)).hexdigest()
        if det is not None:
            assert all(l == det_lower[k] and u == det_upper[k] for l, u, k in zip(lower, upper, labels))
            counts['n2_sample_score_lookup_links'] += n
        sample_lo, sample_hi = F(sum(lower), n * (1 << 40)), F(sum(upper), n * (1 << 40))
        assert pair(mc['sample_mean_interval']) == (sample_lo, sample_hi)
        finite_lo, finite_hi = max(F(0), sample_lo - eps), min(F(1), sample_hi + eps)
        assert pair(mc['finite_grid_mean_confidence_interval']) == (finite_lo, finite_hi)
        assert bounds(failure)[1] <= F(mc['hoeffding_failure_upper']) <= F(1, 3600)
        assert F(mc['statistical_epsilon']) == eps and mc['draws'] == mc['evaluations'] == n
        assert F(mc['maximum_outward_score_width']) == F(max(u-l for l, u in zip(lower, upper)), 1 << 40)
        completed(mc, finite_lo, finite_hi, t)
        counts['sample_scores'] += n
        counts['captured_uint16_words'] += len(words)
        counts['score_endpoints'] += 2 * n
        counts['continuous_price_compositions'] += 1
        primary = det or mc
        rows.append({'model_id': row['model_id'], 'masked_labels_sha256': label_digest,
                     'exact_sample_lower_integer_sum': str(sum(lower)),
                     'exact_sample_upper_integer_sum': str(sum(upper)),
                     'primary_cent_admission': primary['certifies_undiscounted_cent'],
                     'primary_radius_exact_fraction': primary['undiscounted_radius_upper']})
    assert all(r['primary_cent_admission'] for r in rows)
    assert result['all_four_primary_prices_certify_cent'] and result['all_four_monte_carlo_prices_certify_cent']
    assert counts['sample_scores'] == 400000 and counts['captured_uint16_words'] == 1400000
    summary = json.loads((RESULTS / 'four_case_classical_summary.json').read_text())
    assert summary['full_result_sha256'] == sha(result_path)
    assert summary['records'] == result['records']
    audit = {'schema': 'qbp-four-case-independent-exact-arithmetic-audit-v1', 'status': 'passed',
             'date': '2026-10-09', 'auditor_sha256': sha(__file__),
             'producer_sha256': sha(ROOT / 'four_case_classical_pricing.py'),
             'protocol_sha256': sha(protocol_path), 'result_sha256': sha(result_path),
             'input_sha256': protocol['input_sha256'], 'counts': counts, 'records': rows,
             'independent_probability_and_mass_precision_bits': 1024,
             'uniform_word_preimage_counts': {'labels': 1024, 'uint16_preimages_each': 64},
             'confidence_scope': 'The simultaneous guarantee is conditional on independent uniform uint16 source words. Retention, masking checks and realised score arithmetic do not prove that physical premise.',
             'financial_transcendental_replays': 0,
             'limits': 'The separate all-score auditor verifies every financial score and fresh cell integrals. This audit independently checks retained binary arithmetic, exact sums, model provenance, the Gaussian retained mass, deterministic derivative/budget composition, Hoeffding, complete intervals and discounting.',
             'LSEG_inputs_found': False, 'licensed_data_acquired_or_deleted': False,
             'seconds_diagnostic': time.perf_counter() - start}
    (STAGE / 'licence_review/final_classical_arithmetic_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print('Passed', counts)

if __name__ == '__main__':
    main()
