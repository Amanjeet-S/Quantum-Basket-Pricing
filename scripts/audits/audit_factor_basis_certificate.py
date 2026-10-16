"""Independent algebraic replay of the frozen exact factor-basis screen.

The candidate matrix, its multiplication and invariant checks are rebuilt
without importing the new producer. Existing reviewed Gaussian-cell and
direct-transfer implementations are used as independently retained checks.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse, hashlib, json, sys

ROOT = Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')
sys.path.insert(0, str(ROOT))
from flint import ctx
import tensor_transfer as tt
from boundary_grid_certificate import parameters
from boundary_derivative_certificate import payoff_norms
from four_survivor_quantum_resources import direct_first_transfer
from taylor_pricing import compositions

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dot(a, b):
    return sum((x*y for x,y in zip(a,b)), F(0))

def timeless(value):
    if isinstance(value, dict):
        return {k:timeless(v) for k,v in value.items() if not k.endswith('seconds_diagnostic')}
    if isinstance(value,list):
        return list(map(timeless,value))
    return value

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', type=Path, default=ROOT)
    ap.add_argument('--results', type=Path, default=Path(__file__).resolve().parent/'results')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    protocol_path = args.results/'factor_basis_protocol.json'
    result_path = args.results/'factor_basis_certificates.json'
    summary_path = args.results/'factor_basis_summary.json'
    protocol = json.loads(protocol_path.read_text())
    result = json.loads(result_path.read_text())
    summary = json.loads(summary_path.read_text())
    assert result['protocol_sha256'] == sha(protocol_path)
    assert summary['certificate_sha256'] == sha(result_path)
    assert result['source_sha256'] == protocol['source_sha256']
    for field in ('dependency_source_sha256', 'input_sha256'):
        assert result[field] == protocol[field]
        for name, value in protocol[field].items():
            assert value == sha(args.repo/name), name
    ctx.prec = 512
    transport = tt.TransportMGF(6, 10, 28)
    assert result['actual_financial_node_integers'] == [x['actual_node_integer'] for x in transport.cells]
    assert timeless(result['transport_profile']) == timeless(transport.profile())
    records = result['records']
    allkeys = {(row['original_model_id'], row['candidate_id']) for row in records}
    assert len(allkeys) == len(records) == protocol['candidate_count'] == 58
    n_orthogonality = n_grams = n_monomials = n_direct = 0
    identity_baselines = json.loads((args.repo/'results/four_survivor_quantum_resources.json').read_text())['preparations']
    for original in protocol['models']:
        n, p = original['assets'], original['outer_dimension']
        expected = [{'id': 'identity', 'dimension': p, 'plane': None, 'cosine': '1', 'sine': '0'}]
        for a, b in combinations(range(p), 2):
            for sign, word in ((-1, 'minus'), (1, 'plus')):
                expected.append({'id': f'plane_{a}_{b}_{word}', 'dimension': p,
                    'plane': [a,b], 'cosine': '3/5', 'sine': str(sign*F(4,5))})
        assert expected == protocol['candidate_pools'][original['model_id']]
        family = [r for r in records if r['original_model_id'] == original['model_id']]
        assert {r['candidate_id'] for r in family} == {c['id'] for c in expected}
        w,d,s,mg,rg,lam = parameters(original, False)
        R = tuple(tuple(map(F, x)) for x in original['loadings'])
        for candidate in expected:
            row = next(r for r in family if r['candidate_id'] == candidate['id'])
            assert row['candidate'] == candidate
            O = [[F(a == b) for b in range(p)] for a in range(p)]
            if candidate['plane'] is not None:
                a,b = candidate['plane']
                sine = F(candidate['sine'])
                O[a][a] = O[b][b] = F(3,5)
                O[a][b],O[b][a] = -sine,sine
            for a in range(p):
                for b in range(p):
                    assert dot(O[a], O[b]) == F(a == b)
                    assert sum(O[k][a]*O[k][b] for k in range(p)) == F(a == b)
                    n_orthogonality += 2
            actual_model = row['rotated_model']
            computed = tuple(tuple(sum(R[i][k]*O[k][j] for k in range(p)) for j in range(p)) for i in range(n))
            assert computed == tuple(tuple(map(F,x)) for x in actual_model['loadings'])
            expected_id = original['model_id'] if candidate['id'] == 'identity' else original['model_id']+'__'+candidate['id']
            assert actual_model['model_id'] == expected_id
            for name in original:
                if name not in ('model_id', 'loadings'):
                    assert original[name] == actual_model[name], name
            ww,dd,ss,mmg,rrg,llam = parameters(actual_model,False)
            assert (ww,dd,mmg,llam) == (w,d,mg,lam)
            for i in range(n):
                for k in range(n):
                    assert dot(computed[i],computed[k]) == dot(R[i],R[k])
                    assert dot(ss[i],ss[k]) == dot(s[i],s[k])
                    assert dot(computed[i],computed[k])+lam*lam == F(original['covariance'][i][k])
                    n_grams += 3
                assert dot(rrg,ss[i]) == dot(rg,s[i])
            assert dot(rrg,rrg) == dot(rg,rg)
            # This finite explicit check supplements the all-orders structural
            # proof. Every rate vector is independently rebuilt in both bases.
            for order in range(4):
                for alpha in compositions(order,n):
                    v = tuple(sum(alpha[i]*s[i][j] for i in range(n)) for j in range(p))
                    vv = tuple(sum(alpha[i]*ss[i][j] for i in range(n)) for j in range(p))
                    assert dot(v,v) == dot(vv,vv)
                    assert dot(rg,v) == dot(rrg,vv)
                    n_monomials += 2
            norm = payoff_norms(actual_model['lambda'],F(1))
            independent = direct_first_transfer(actual_model,transport,norm)
            claimed = row['direct_global_transfer']
            for name, value in independent.items():
                if name != 'scope':
                    assert value == claimed[name], (candidate['id'],name)
            B1 = sum((max(abs(computed[i][j]) for i in range(n))+abs(rrg[j]) for j in range(p)),F(0))
            B2 = sum((max(abs(computed[i][j]) for i in range(n))**2 for j in range(p)),F(0))+dot(rrg,rrg)
            assert F(claimed['geometric_B1']) == B1
            assert F(claimed['geometric_B2']) == B2
            assert F(claimed['financial_transfer_upper']) == sum(F(claimed[name]) for name in
                ('linear_transport_upper','quadratic_transport_upper','tail_upper'))
            assert row['invariants']['all_continuous_U_moments_invariant']
            assert not row['invariants']['finite_grid_invariance_claimed']
            n_direct += 1
            if candidate['id'] == 'identity':
                old = next(x for x in identity_baselines if x['model']['model_id'] == original['model_id'])['R1_direct_global_financial_transfer']
                for name,value in old.items():
                    if name != 'scope':
                        assert value == claimed[name], ('retained identity',name)
        winner = min(family,key=lambda r:(F(r['direct_global_transfer']['financial_transfer_upper']),r['candidate_id']))
        assert next(x for x in result['selected_models'] if x['original_model_id']==original['model_id']) == winner
        assert next(x for x in summary['records'] if x['model_id']==original['model_id'])['selected_candidate'] == winner['candidate_id']
    out = {'schema':'qbp-exact-factor-basis-independent-audit-v1','status':'passed',
        'protocol_sha256':sha(protocol_path),'certificates_sha256':sha(result_path),'summary_sha256':sha(summary_path),
        'audit_source_sha256':sha(__file__), 'candidate_count':len(records),
        'exact_orthogonality_scalar_identities':n_orthogonality,
        'exact_asset_and_centred_covariance_scalar_identities':n_grams,
        'explicit_gaussian_monomial_rate_and_cross_checks_orders_0_through_3':n_monomials,
        'full_direct_transfer_replays':n_direct,'identity_baseline_transfer_matches':4,
        'financial_price_answers_used':0,
        'scope':'Independent exact matrix, covariance and finite-pool selection replay. Fresh established Gaussian cell/direct-transfer implementation agrees. These checks supplement the all-order structural proof and do not claim invariance of the finite grid.'}
    output = args.output or args.results/'factor_basis_independent_audit.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__ == '__main__':
    main()
