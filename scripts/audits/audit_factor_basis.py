"""Independent exact-rational audit of latent Gaussian basis changes.

The transform producer is deliberately not imported. All pair mathematics
is rederived here from Fraction operations. `audit_pair` is reusable by the
registered experiment auditor. The CLI self-test also demonstrates why
continuous Gaussian invariance does not imply finite tensor-grid invariance.
No root repository is written by this script.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

DEFAULT_REPO = Path('/Users/amanjeetsingh/Projects/Quantum-Basket-Pricing')


def dot(a, b):
    assert len(a) == len(b)
    return sum((x*y for x, y in zip(a, b)), F(0))


def transpose(a):
    return [list(x) for x in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[dot(x, y) for y in bt] for x in a]


def identity(p):
    return [[F(int(i == j)) for j in range(p)] for i in range(p)]


def plane_rotation(p, a, b, sine=F(4, 5)):
    assert 0 <= a < b < p and sine in (F(4, 5), F(-4, 5))
    o = identity(p)
    o[a][a] = o[b][b] = F(3, 5)
    o[a][b] = -sine
    o[b][a] = sine
    return o


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for j in range(total+1):
            for rest in compositions(total-j, parts-1):
                yield (j,)+rest


def parameters(model):
    w = list(map(F, model['weights']))
    r = [list(map(F, row)) for row in model['loadings']]
    mu = list(map(F, model['log_drifts']))
    lam = F(model['lambda'])
    n, p = model['assets'], model['outer_dimension']
    assert len(w) == n == len(r) == len(mu)
    assert all(len(row) == p for row in r)
    assert sum(w) == 1 and min(w) > 0 and lam > 0
    rate = F(model['risk_free_rate']) - F(model['dividend_yield'])
    maturity = F(model['maturity'])
    # The retained parameters describe log returns at T=1. Do not silently
    # extrapolate this assertion to a different maturity representation.
    assert maturity == 1 and model['spot'] == model['strike'] == 100
    assert all(mu[i]+(lam*lam+dot(r[i], r[i]))/2 == rate for i in range(n))
    mg = dot(w, mu)
    rg = [sum((w[i]*r[i][j] for i in range(n)), F(0)) for j in range(p)]
    d = [x-mg for x in mu]
    s = [[x-y for x, y in zip(row, rg)] for row in r]
    assert dot(w, d) == 0
    assert all(sum((w[i]*s[i][j] for i in range(n)), F(0)) == 0 for j in range(p))
    return w, r, mu, lam, mg, rg, d, s


def covariance(r, lam):
    return [[lam*lam+dot(x, y) for y in r] for x in r]


def gram(r):
    return [[dot(x, y) for y in r] for x in r]


def strip_derived_model(original, o):
    # Test helper only. Registered candidates are supplied by the producer
    # and then independently compared in audit_pair.
    m = dict(original)
    m['loadings'] = [[str(x) for x in row] for row in matmul([list(map(F, row)) for row in original['loadings']], o)]
    return m


def audit_derivatives(r, rotated, o, rg, rotated_rg):
    """Exact chain rule and envelope checks for three positive softmax laws.

    Arbitrary positive pi laws can be realised by positive asset levels.
    P' and P''-P' stand in for real bounded smooth kernel derivatives.
    The checks below supplement, rather than replace, the universal proof.
    """
    n, p = len(r), len(r[0])
    distributions = [[F(1,n)]*n,
                     [F(i+1, n*(n+1)//2) for i in range(n)],
                     [F(n-i, n*(n+1)//2) for i in range(n)]]
    checks = 0
    for pi in distributions:
        lj = [sum((pi[i]*r[i][j] for i in range(n)), F(0)) for j in range(p)]
        jl = [sum((pi[i]*rotated[i][j] for i in range(n)), F(0)) for j in range(p)]
        assert [dot(lj, col) for col in transpose(o)] == jl
        h = [[sum((pi[i]*r[i][j]*r[i][k] for i in range(n)), F(0))-lj[j]*lj[k] for k in range(p)] for j in range(p)]
        hr = [[sum((pi[i]*rotated[i][j]*rotated[i][k] for i in range(n)), F(0))-jl[j]*jl[k] for k in range(p)] for j in range(p)]
        assert matmul(matmul(transpose(o), h), o) == hr
        for u, v, ug, vg in product((F(-2), F(0), F(2)), (F(-3), F(0), F(3)), repeat=2):
            m1, m2 = F(2), F(3)
            # Paired full derivative matrices are conjugated by O.
            f1 = [u*lj[j]-ug*rg[j] for j in range(p)]
            f1r = [u*jl[j]-ug*rotated_rg[j] for j in range(p)]
            assert [dot(f1, col) for col in transpose(o)] == f1r
            f2 = [[(u+v)*lj[j]*lj[k]+u*h[j][k]-(ug+vg)*rg[j]*rg[k] for k in range(p)] for j in range(p)]
            f2r = [[(u+v)*jl[j]*jl[k]+u*hr[j][k]-(ug+vg)*rotated_rg[j]*rotated_rg[k] for k in range(p)] for j in range(p)]
            assert matmul(matmul(transpose(o), f2), o) == f2r
            for j in range(p):
                a = max(abs(row[j]) for row in rotated)
                assert abs(f1r[j]) <= m1*(a+abs(rotated_rg[j]))
                expanded = (u+v)*jl[j]**2+u*hr[j][j]-(ug+vg)*rotated_rg[j]**2
                cancelled = v*jl[j]**2+u*sum((pi[i]*rotated[i][j]**2 for i in range(n)), F(0))-(ug+vg)*rotated_rg[j]**2
                assert expanded == cancelled
                assert abs(cancelled) <= (m1+m2)*(a*a+rotated_rg[j]**2)
                checks += 1
    return checks


def audit_pair(original, transformed, orthogonal_matrix, order=6):
    o = [list(map(F, row)) for row in orthogonal_matrix]
    w, r, mu, lam, mg, rg, d, s = parameters(original)
    wn, rn, mun, lamn, mgn, rgn, dn, sn = parameters(transformed)
    p, n = len(rg), len(w)
    assert len(o) == p and all(len(row) == p for row in o)
    assert matmul(transpose(o), o) == identity(p)
    assert matmul(o, transpose(o)) == identity(p)
    assert matmul(r, o) == rn
    assert (w, mu, lam, mg, d) == (wn, mun, lamn, mgn, dn)
    for key in ('assets','outer_dimension','weights','log_drifts','lambda','spot','strike','maturity','risk_free_rate','dividend_yield'):
        assert original[key] == transformed[key], key
    assert matmul([rg], o)[0] == rgn
    assert matmul(s, o) == sn
    assert covariance(r, lam) == covariance(rn, lamn)
    assert gram(r) == gram(rn)
    assert gram(s) == gram(sn)
    assert dot(rg, rg) == dot(rgn, rgn)
    dispersion = sum((w[i]*dot(r[i],r[i]) for i in range(n)), F(0))-dot(rg,rg)
    dispersion_n = sum((wn[i]*dot(rn[i],rn[i]) for i in range(n)), F(0))-dot(rgn,rgn)
    assert dispersion == dispersion_n
    # If derived fields are retained they must still be correct. Covariance
    # and scalar invariants may be reused because they were exactly proved.
    for model, rr, rgx in ((original,r,rg),(transformed,rn,rgn)):
        if 'covariance' in model:
            assert [list(map(F,row)) for row in model['covariance']] == covariance(rr, F(model['lambda']))
        if 'geometric_outer_variance' in model:
            assert F(model['geometric_outer_variance']) == dot(rgx,rgx)
        if 'weighted_residual_dispersion' in model:
            assert F(model['weighted_residual_dispersion']) == dispersion
    # Pointwise identity F_RO(y)=F_R(Oy) for the asset/geometric arguments,
    # established without floating exponentials or imported producer code.
    node_vectors = [[F((i+1)*(j+2)-7, 11) for j in range(p)] for i in range(7)]
    pointwise = 0
    for y in node_vectors:
        oy = [dot(row,y) for row in o]
        for i in range(n):
            assert mu[i]+dot(r[i],oy) == mun[i]+dot(rn[i],y)
            assert d[i]+dot(s[i],oy) == dn[i]+dot(sn[i],y)
            pointwise += 2
        assert mg+dot(rg,oy) == mgn+dot(rgn,y)
        pointwise += 1
    # Every continuous multinomial moment exponent through r=6 is identical,
    # although the transformed rates and finite products of MGFs differ.
    exponents = 0
    for total in range(order+1):
        for alpha in compositions(total,n):
            v = [sum((F(alpha[i])*s[i][j] for i in range(n)), F(0)) for j in range(p)]
            vn = [sum((F(alpha[i])*sn[i][j] for i in range(n)), F(0)) for j in range(p)]
            assert matmul([v],o)[0] == vn and dot(v,v) == dot(vn,vn)
            const = dot(list(map(F,alpha)),d)
            constn = dot(list(map(F,alpha)),dn)
            assert const+dot(v,v)/2 == constn+dot(vn,vn)/2
            assert dot(rg,v) == dot(rgn,vn)
            exponents += 1
    # Existing c1 formula uses each centred row norm and rg.s_i, both invariant.
    assert all(dot(x,x) == dot(y,y) and dot(rg,x) == dot(rgn,y) for x,y in zip(s,sn))
    derivative_checks = audit_derivatives(r,rn,o,rg,rgn)
    a0 = [max(abs(row[j]) for row in r) for j in range(p)]
    a1 = [max(abs(row[j]) for row in rn) for j in range(p)]
    s2 = sum((x*x for x in a1), F(0))
    assert max(dot(row,row) for row in rn) <= s2 <= sum((dot(row,row) for row in rn), F(0))
    return {'assets':n,'outer_dimension':p,'multinomial_exponent_checks':exponents,
            'pointwise_affine_checks':pointwise,'derivative_scalar_checks':derivative_checks,
            'original_full_loading_sum':str(sum(a0)),
            'transformed_full_loading_sum':str(sum(a1)),
            'original_full_loading_square_sum':str(sum(x*x for x in a0)),
            'transformed_full_loading_square_sum':str(s2),
            'continuous_invariants':'Exact no-arbitrage drifts, asset covariance/Gram matrix, centred Gram matrix, geometric variance and c0/c1 moment arguments.',
            'finite_invariants_not_asserted':'Product rounded grid, retained cube, finite moments, hybrid moments, tails other than cap-mass bound, domain/normaliser, arithmetic certificates, priors and resource charges must be recomputed.'}


def finite_counterexample():
    r = [[F(1), F(0)]]
    o = plane_rotation(2,0,1)
    rn = matmul(r,o)
    points = list(product((F(-1),F(1)),repeat=2))
    fourth = sum((dot(r[0],y)**4 for y in points),F(0))/4
    rotated_fourth = sum((dot(rn[0],y)**4 for y in points),F(0))/4
    assert dot(r[0],r[0]) == dot(rn[0],rn[0]) == 1
    assert fourth == 1 and rotated_fourth == F(1201,625)
    assert fourth != rotated_fourth
    return {'finite_uniform_sign_law_fourth_moment_original':str(fourth),
            'finite_uniform_sign_law_fourth_moment_rotated':str(rotated_fourth),
            'continuous_standard_gaussian_fourth_moment_both':'3',
            'conclusion':'Exact Gaussian covariance invariance does not justify finite tensor-grid invariance.'}


def find_pairs(obj):
    """Protocol-agnostic adapter; producer may call audit_pair directly."""
    if isinstance(obj,dict):
        keys=('original_model','rotated_model','orthogonal_matrix')
        if all(k in obj for k in keys):
            yield obj
        for value in obj.values():
            yield from find_pairs(value)
    elif isinstance(obj,list):
        for value in obj:
            yield from find_pairs(value)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo',type=Path,default=DEFAULT_REPO)
    ap.add_argument('--protocol',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    frozen = json.loads((args.repo/'results/four_survivor_quantum_protocol.json').read_text())
    rows=[]
    for model in frozen['models']:
        p=model['outer_dimension']
        matrices=[identity(p)]
        if p==1:
            matrices.append([[F(-1)]])
        else:
            for a,b in ((0,1),(p-2,p-1)):
                for sine in (F(4,5),F(-4,5)):
                    matrices.append(plane_rotation(p,a,b,sine))
        for o in matrices:
            rows.append(audit_pair(model,strip_derived_model(model,o),o))
    registered=[]
    if args.protocol:
        protocol=json.loads(args.protocol.read_text())
        for pair in find_pairs(protocol):
            registered.append(audit_pair(pair['original_model'],pair['rotated_model'],pair['orthogonal_matrix']))
        if not registered:
            raise ValueError('No adapted exact model pairs found; call audit_pair with the producer schema or adapt protocol keys.')
    result={'schema':'qbp-independent-factor-basis-mathematics-v1','status':'passed',
            'self_test_pair_count':len(rows),'self_test_pairs':rows,
            'registered_pair_count':len(registered),'registered_pairs':registered,
            'finite_grid_counterexample':finite_counterexample(),
            'source_sha256':sha(__file__),
            'frozen_four_case_protocol_sha256':sha(args.repo/'results/four_survivor_quantum_protocol.json'),
            'registered_protocol_sha256':None if not args.protocol else sha(args.protocol),
            'scope':'Exact rational mathematical checks and an explicit finite-law counterexample. No financial pricing run, resource reconstruction or originality certification.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'self_test_pairs':len(rows),'registered_pairs':len(registered),'output':str(args.output)}))


if __name__=='__main__':
    main()
