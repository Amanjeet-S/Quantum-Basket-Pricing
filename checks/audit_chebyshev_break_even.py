"""Small exact-arithmetic audit of the added expected-charge derivation.

I preserve the full earlier mathematical execution record and identify its
proof version. This additive review checks the new section only. The fixed
mesh averages below are descriptive, not estimates at the unknown financial
probability and not a demonstration of expected quantum advantage.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=ROOT/'results/chebyshev_break_even_audit.json')
    args = parser.parse_args()
    own = sha(__file__)
    proof = ROOT/'chebyshev_estimator_proof.tex'
    full_audit = ROOT/'results/chebyshev_math_audit.json'
    result = ROOT/'results/chebyshev_resources.json'
    previous = json.loads(full_audit.read_text())
    assert previous['status'] == 'passed'
    text = proof.read_text()
    begin = '% BEGIN RETAINED CHEBYSHEV PROOF BODY\n'
    end_marker = '% END RETAINED CHEBYSHEV PROOF BODY\n'
    wrapped = begin in text or end_marker in text
    if wrapped:
        assert text.count(begin) == text.count(end_marker) == 1
        text = text.split(begin, 1)[1].split(end_marker, 1)[0]
    proof_body_sha = hashlib.sha256(text.encode()).hexdigest()
    start = text.index('\\subsection{An expected-gate break-even condition}')
    end = text.index('For an undiscounted financial decoding satisfying', start)
    old_text = text[:start]+text[end:]
    old_sha = hashlib.sha256(old_text.encode()).hexdigest()
    assert old_sha == previous['source_sha256']['chebyshev_estimator_proof.tex']
    assert sha(result) == previous['result_sha256']
    data = json.loads(result.read_text())
    rows = []
    histories = 0
    for case in data['cases']:
        t0 = case['retained_financial_information']['previous_complete_QPE_T']
        tf = case['half_alpha_QPE_fallback']['T_count']
        cap = case['pilot_T_cap']
        assert F(tf, t0) == F(27, 25)
        assert cap == tf//10
        threshold = 1-F(t0, tf)
        assert threshold == F(2, 27)
        assert threshold+F(cap, tf) <= F(47, 270)
        records = case['diagnostics']
        pilot_sum = successes = total_sum = 0
        for row in records:
            pilot = row['controller_result']['pilot_T_count']
            success = row['controller_result']['status'] == 'certified'
            assert row['actual_history_T_count'] == pilot+int(not success)*tf
            assert 0 <= pilot <= cap
            pilot_sum += pilot
            successes += int(success)
            total_sum += row['actual_history_T_count']
            histories += 1
        n = len(records)
        pi = F(successes, n)
        pilot_mean = F(pilot_sum, n)
        charge_mean = F(total_sum, n)
        break_even = threshold+pilot_mean/tf
        assert charge_mean == pilot_mean+(1-pi)*tf
        assert (charge_mean < t0) == (pi > break_even)
        if charge_mean < t0:
            assert pi > F(2, 27)
        if pi > threshold+F(cap, tf):
            assert charge_mean < t0
        rows.append({'model_id': case['model_id'],
                     'fallback_to_original_QPE_ratio': str(F(tf, t0)),
                     'necessary_strict_success_threshold': '2/27',
                     'sufficient_strict_budget_threshold': str(threshold+F(cap, tf)),
                     'uniform_conservative_sufficient_strict_threshold': '47/270',
                     'descriptive_mesh_success_fraction': str(pi),
                     'descriptive_mesh_mean_pilot_T': str(pilot_mean),
                     'descriptive_mesh_break_even_threshold': str(break_even),
                     'descriptive_mesh_mean_complete_T': str(charge_mean),
                     'descriptive_mesh_mean_to_original_QPE_ratio': str(charge_mean/t0),
                     'descriptive_mesh_mean_below_original_QPE': charge_mean < t0,
                     'scope': 'Equal-weight registered mesh and seeds only; not expectations at the unknown true financial probability.'})
    # Exact equality examples check the strictness of both criteria.
    tf, t0, cap = 270, 250, 27
    def charge(pi, mean_pilot):
        return F(mean_pilot)+(1-F(pi))*tf
    assert charge(F(2, 27), 0) == t0
    assert charge(F(47, 270), cap) == t0
    assert charge(F(48, 270), cap) < t0
    assert charge(F(3, 27), cap) > t0
    assert own == sha(__file__)
    out = {'schema': 'qbp-chebyshev-added-break-even-proof-audit-v1',
           'status': 'passed', 'auditor_sha256': own,
           'proof_sha256': sha(proof),
           'retained_proof_body_sha256': proof_body_sha,
           'standalone_presentation_wrapper': wrapped,
           'previous_full_math_audit_sha256': sha(full_audit),
           'previous_proof_sha256': old_sha,
           'proof_revision_scope': 'After extracting the explicitly marked mathematical body from any standalone presentation wrapper, removing only the added expected-gate subsection reconstructs the exact proof version recorded in the full earlier mathematical audit. Its checked mathematics and execution sources remain unchanged.',
           'resource_result_sha256': sha(result),
           'checked_history_charge_identities': histories,
           'checked_fixed_flag_resource_ratios_and_mean_identities': len(rows),
           'checked_strict_threshold_examples': 4,
           'cases': rows,
           'interpretation': 'The general break-even identity follows from boundedness and linear expectation without independence. Finite mesh success fractions and means are descriptive only. The first three mesh means are below the old full-alpha QPE count; the fourth is above it. These do not establish expected performance at the financial flag.'}
    output = args.output
    output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
