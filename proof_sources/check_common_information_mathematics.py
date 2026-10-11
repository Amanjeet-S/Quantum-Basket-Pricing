"""Independent Gaussian/common-shift algebra and retained-field error audit.

No producer, graph builder, smoothing mapper or confidence frontend is imported.
Fresh Arb enclosures supplement the general proofs; they are not substitutes
for those proofs. Scalar ordinate correctness is a retained audited premise.
"""
import argparse, gzip, hashlib, json, sys
from fractions import Fraction as F
from pathlib import Path
from flint import arb, ctx
sys.dont_write_bytecode = True
ctx.prec = 512
COUNTS = {}

def mark(key, amount=1):
    COUNTS[key] = COUNTS.get(key, 0)+amount

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    path = Path(path)
    return json.load(gzip.open(path, 'rt')) if path.suffix == '.gz' else json.loads(path.read_text())

def ball(x):
    x = F(x)
    return arb(x.numerator)/x.denominator

def endpoint(a, up):
    m, e = map(int, (a.upper() if up else a.lower()).man_exp())
    return F(m*(1 << e)) if e >= 0 else F(m, 1 << -e)

def Phi(x):
    return (1+(x/arb(2).sqrt()).erf())/2

def normal_density(x):
    return (-x*x/2).exp()/(2*arb.pi()).sqrt()

def intrinsic_index(path):
    """Stream individual archived certificates; retain no giant node arrays."""
    decoder, buffer, eof = json.JSONDecoder(), '', False
    result = []
    with gzip.open(path, 'rt') as stream:
        while True:
            buffer = buffer.lstrip()
            while buffer.startswith(('[', ',')):
                buffer = buffer[1:].lstrip()
            if buffer.startswith(']'):
                break
            if not buffer and eof:
                raise ValueError('Unterminated scalar array')
            try:
                value, stop = decoder.raw_decode(buffer)
            except json.JSONDecodeError:
                part = stream.read(1 << 20)
                if not part:
                    if eof:
                        raise
                    eof = True
                buffer += part
                continue
            buffer = buffer[stop:]
            c = value['certificate']
            if c.get('name') == 'cap':
                result.append({k: c[k] for k in ('name','lo','span','params','table_bits','q',
                        'output_error_upper','first_derivative_upper','second_derivative_upper',
                        'nodal_error_upper','linear_interpolation_error_upper',
                        'interpolation_floor_error_upper')})
            del value
    return result

def own_coeff(x, q):
    x = F(x)
    if x.denominator & (x.denominator-1) == 0:
        return x.numerator, x.denominator.bit_length()-1
    t = x*(1 << q)+F(1, 2)
    return t.numerator//t.denominator, q

def analytic_tests(nu0, nu1):
    # Three independently varied positive measures include atoms and mixtures.
    measures = [((F(0), F(1)),),
                ((F(-1, 5), F(1, 2)), (F(3, 10), F(1, 2))),
                ((F(-1, 50), F(1, 2)), (F(0), F(1, 3)), (F(1, 10), F(1, 6)))]
    lams = (F(1, 20), F(1, 5), F(2, 5))
    shifts = (F(-1, 100), F(-1, 10000), F(1, 10000), F(1, 100))
    for mu in measures:
        assert sum(w for _, w in mu) == 1 and min(w for _, w in mu) > 0
        for lam in lams:
            c0, c1 = nu0/lam, nu1/(lam*lam)
            def P(m):
                return sum((ball(w)*Phi(ball((m-a)/lam)) for a,w in mu), arb(0))
            def Rp(m, d):
                return P(m+d)-P(m)
            for m in (F(-1),F(-1,5),F(0),F(3,20),F(1,2),F(1)):
                for d in (F(0),F(1,10000),F(1,1000),F(1,100),F(1,10),F(1,2)):
                    derivative = sum((ball(w)*(normal_density(ball((m+d-a)/lam))
                        -normal_density(ball((m-a)/lam)))/ball(lam) for a,w in mu),arb(0))
                    if d:
                        assert endpoint(abs(derivative),True) <= min(c0,c1*d)
                    mark('fresh_atomic_mixture_common_derivative_checks')
                    for dm in shifts:
                        for dd in shifts:
                            bound = min(c0,c1*d)*abs(dm)+c0*abs(dd)
                            err = abs(Rp(m+dm,d+dd)-Rp(m,d))
                            assert endpoint(err,True) <= bound
                            # Projecting a negative represented gap fixes the
                            # true positive correction. Check it numerically.
                            if d+dd < 0:
                                projected_error = abs(Rp(m,d))
                                assert endpoint(projected_error,True) <= bound
                                mark('negative_represented_gap_projection_checks')
                            mark('fresh_Arb_common_shift_path_checks')
    # An inner-label-dependent perturbation is not a fixed Gaussian shift.
    lam, gap, e = F(1,5),F(1,1000),F(1,10000)
    # Move both outside boundary strips into (0,gap]; leave interior unchanged.
    increase = Phi(ball(gap/lam+e/lam))-Phi(ball(gap/lam))+Phi(ball(e/lam))-arb(1)/2
    false_allowance = min(nu0/lam,nu1*gap/(lam*lam))*e
    assert endpoint(increase,False) > false_allowance
    mark('inner_dependent_shift_counterexamples')
    counterexample = {'lambda':str(lam),'gap':str(gap),'error':str(e),
          'actual_increase_lower':str(endpoint(increase,False)),
          'invalid_common_shift_allowance':str(false_allowance),
          'valid_sandwich_allowance':str(2*nu0/lam*e),
          'scope':'A deterministic normal-label-dependent boundary perturbation. It refutes applying the fixed-shift coefficient to an arbitrary dependent scaling remainder; it is not asserted to be the retained emitted scaling remainder.'}
    # The payoff finite-difference Lipschitz bound handles kink crossings.
    ell = (1+arb(1)/100).log()
    ell_lo,ell_hi = endpoint(ell,False),endpoint(ell,True)
    def f(x):
        if x <= 0:
            return arb(0)
        if x >= ell_hi:
            return arb(1)
        y = 100*(ball(x).exp()-1)
        # Avoid undecidable tiny balls at the exact irrational cap boundary.
        assert x < ell_lo
        return y
    xs = (F(-1,100),F(-1,10000),F(0),F(1,10000),F(1,1000),F(9,1000),F(1,100),F(1,50))
    for d in (F(0),F(1,10000),F(1,1000),F(1,100),F(1,10)):
        for x in xs:
            for z in shifts:
                err = abs((f(x+d+z)-f(x+z))-(f(x+d)-f(x)))
                assert endpoint(err,True) <= 101*abs(z)
                mark('fresh_payoff_common_perturbation_checks')
    # Duplicate atomic labels exercise both strict and non-strict endpoints.
    nodes = (F(-1),F(0),F(0),F(1))
    D = F(0)
    for x in sorted(set(nodes)):
        left,right = F(sum(v<x for v in nodes),4),F(sum(v<=x for v in nodes),4)
        normal = Phi(ball(x))
        D = max(D,endpoint(abs(ball(left)-normal),True),endpoint(abs(ball(right)-normal),True))
        mark('duplicate_atom_CDF_left_right_checks',2)
    ends = (F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2))
    for lo in ends:
        for hi in ends:
            if lo>hi:
                continue
            for strict_lo, strict_hi in ((True,False),(False,True),(True,True),(False,False)):
                count = sum((x>lo if strict_lo else x>=lo) and
                            (x<hi if strict_hi else x<=hi) for x in nodes)
                probability = F(count,4)
                if lo==hi and not strict_lo and not strict_hi:
                    # A singleton atom is an interval using left and right
                    # limits at the same point, still subject to 2D.
                    continuous = arb(0)
                elif lo==hi:
                    continuous = arb(0)
                else:
                    continuous = Phi(ball(hi))-Phi(ball(lo))
                assert endpoint(abs(ball(probability)-continuous),True) <= 2*D
                mark('finite_atom_interval_transfer_checks')
    return counterexample

def metadata_match(capcerts,row,side,body):
    sem,cfg = row['graph_semantic_nodes'],row['graph_configuration']
    raw = sem['mA_offset'] if side=='PA' else sem['mG_offset']
    eraw,efull = F(body[raw]['absolute_output_error_upper']),F(body[sem[side]]['absolute_output_error_upper'])
    candidates = [c for c in capcerts if c['table_bits']==row[side+'_bits'] and c['q']==cfg['q']
      and F(c['params']['lambda'])==F(cfg['lambda']) and F(c['params']['K'])==100
      and F(c['params']['cap'])==1 and F(c['output_error_upper'])+F(c['first_derivative_upper'])*eraw==efull]
    assert len(candidates)==1,(row['original_model_id'],side,len(candidates))
    c = candidates[0]
    assert F(c['lo'])==-F(c['span'])/2
    actual = tuple(map(F,body[raw]['actual_bounds']))
    assert 0 <= actual[0] <= actual[1] <= F(c['span'])
    own_intrinsic = F(c['nodal_error_upper'])+F(c['second_derivative_upper'])*(F(c['span'])/(1 << c['table_bits']))**2/8+F(1,1 << c['q'])
    assert own_intrinsic==F(c['output_error_upper'])
    assert F(c['interpolation_floor_error_upper'])==F(1,1 << c['q'])
    mark('numerical_intrinsic_decomposition_and_inactive_clamp_checks')
    return c

def amended_checks(repo,path,nu0,nu1):
    d = read(path)
    capcerts = intrinsic_index(repo/'results/strengthened_r1_scalar_ordinates.json.gz')
    cache = {}
    for row in d['records']:
        source = row['source_path']
        if source not in cache:
            cache[source] = read(repo/'results'/source)
        sd = cache[source]
        old = sd['rows' if row['method']=='shared_threshold' else 'records'][row['source_row_index']]
        cfg,sem = old['graph_configuration'],old['graph_semantic_nodes']
        lam = F(cfg['lambda']);c0=nu0/lam;C=min(c0,nu1*F(cfg['U_star'])/(lam*lam))
        term = row['probability_certificate']['common_shift_derivation']
        assert F(term['c0'])==c0 and F(term['common_shift_coefficient'])==C and F(term['Delta_star'])==F(cfg['U_star'])
        if row['method']=='raw_payoff':
            body={x['node']:x for x in old['probability_certificate']['body_arithmetic_certificate']['nodes']}
            graph={x['name']:x for x in old['graph_nodes']}
            low=graph[sem['lower_endpoint']];add=graph[low['args'][0]];scale=graph[low['args'][1]]
            assert low['kind']==add['kind']=='add' and add['args'][0]==sem['mG_offset']
            assert scale['kind']=='scale' and scale['args']==[sem['normal_plus_L']]
            em=F(body[sem['mG_offset']]['absolute_output_error_upper'])+F(body[add['args'][1]]['absolute_output_error_upper'])+F(body[sem['common_offset']]['absolute_output_error_upper'])
            ed=F(body[sem['logT']]['absolute_output_error_upper'])
            I,drop=own_coeff(lam,cfg['q'])
            Es=2*F(cfg['inner_normal_L'])*abs(F(I,1 << drop)-lam)+(F(1,1 << cfg['q']) if drop else 0)
            amended_source=row['amended_source_path']
            if amended_source not in cache:
                cache[amended_source]=read(repo/'results'/amended_source)
            asrc=cache[amended_source]
            ar=asrc['amended_records' if 'amended_records' in asrc else 'records'][row['amended_source_row_index']]
            retained=ar['probability_certificate']['smoothing_derivation']
            exp=[F(retained[side+'_intrinsic_exp_approximation_upper']) for side in ('PG','PA')]
            assert list(map(F,term['exp_errors']))==exp and F(term['Em'])==em and F(term['E_s'])==Es and F(term['e_Delta'])==ed
            D=F(old['probability_certificate']['inner_normal_CDF_error_upper'])
            arithmetic=100*sum(exp)+101*Es+C*em+c0*ed
            bias=2*D
            mark('raw_common_information_certificate_checks')
        elif row['method']=='numerical_conditional_price':
            assert F(cfg['monetary_signal_normaliser'])==F(cfg['cap'])==1
            body={x['node']:x for x in old['base_graph_certificate']['nodes']}
            pa,pg=(metadata_match(capcerts,old,side,body) for side in ('PA','PG'))
            ba,bg=-F(pa['lo']),-F(pg['lo'])
            constants=[v for v in body.values() if v['kind']=='constant' and tuple(map(F,v['ideal_bounds']))==(ba-bg,ba-bg)]
            assert constants
            errors={F(v['absolute_output_error_upper']) for v in constants}
            assert len(errors)==1
            ec=next(iter(errors));em=F(body[sem['mG_offset']]['absolute_output_error_upper'])
            ed=F(body[sem['logT']]['absolute_output_error_upper'])+ec
            assert F(body[sem['mA_offset']]['absolute_output_error_upper'])==em+ed
            intrinsic=[F(pa['output_error_upper']),F(pg['output_error_upper'])]
            assert F(term['Em'])==em and F(term['e_Delta_including_offset'])==ed
            assert list(map(F,term['intrinsic_cap_errors']))==intrinsic and term['decoder_scale_rounding_upper']=='0'
            arithmetic=sum(intrinsic)+C*em+c0*ed;bias=F(0)
            mark('numerical_common_information_certificate_checks')
        else:
            body={v['node']:v for v in old['probability_certificate']['body_arithmetic_certificate']['nodes']}
            graph={v['name']:v for v in old['graph_nodes']}
            low=graph[sem['lower_endpoint']];upper=graph[sem['upper_endpoint']];thr=graph[sem['threshold']]
            assert low['kind']=='add' and low['args'][0]==sem['mG_offset']
            assert upper['kind']=='add' and upper['args']==[sem['lower_endpoint'],sem['logT']]
            assert thr['kind']=='sub'
            ta,scale=graph[thr['args'][0]],graph[thr['args'][1]]
            assert ta['kind']=='add' and ta['args'][0]==sem['strike_log']
            assert scale['kind']=='scale' and scale['args']==[sem['normal_plus_L']]
            errors=lambda n:F(body[n]['absolute_output_error_upper'])
            em=sum(errors(n) for n in (sem['mG_offset'],low['args'][1],ta['args'][1],sem['strike_log']))
            ed=errors(sem['logT']);I,drop=own_coeff(lam,cfg['q'])
            Es=2*F(cfg['inner_normal_L'])*abs(F(I,1 << drop)-lam)+(F(1,1 << cfg['q']) if drop else 0)
            assert F(term['common_error_including_strike'])==em and F(term['E_s'])==Es and F(term['e_Delta'])==ed
            arithmetic=2*c0*Es+C*em+c0*ed
            bias=2*F(old['probability_certificate']['inner_normal_CDF_error_upper'])+F(old['probability_certificate']['strike_midpoint_bias_upper'])
            mark('threshold_common_information_certificate_checks')
        assert F(term['arithmetic_upper'])==arithmetic and F(term['auxiliary_bias_upper'])==bias
        assert F(term['new_probability_error_upper'])==min(F(1),arithmetic+bias)
        previous=list(map(F,row['probability_certificate']['previous_probability_error_bounds']))
        eta=min(*previous,min(F(1),arithmetic+bias))
        assert F(row['probability_certificate']['probability_error_upper'])==eta
        budget=row['complete_deterministic_budget'];A=F(budget['external_factor_upper'])
        det=F(budget['coefficient_radius_upper'])+F(budget['financial_transfer_upper'])+A*eta+F(budget['external_factor_radius'])
        assert F(budget['total_deterministic_undiscounted_price_error_upper'])==det
        prior=row['coefficient_informed_prior']
        assert F(prior['probability_lower'])==max(F(0),F(prior['finite_ideal_mean_lower'])-eta)
        assert F(prior['probability_upper'])==min(F(1),F(prior['finite_ideal_mean_upper'])+eta)
        assert eta<=min(previous)
        mark('complete_deterministic_budget_and_prior_checks')
    assert COUNTS['raw_common_information_certificate_checks']==348
    assert COUNTS['numerical_common_information_certificate_checks']==64
    assert COUNTS['threshold_common_information_certificate_checks']==4
    factor=read(repo/'results/factor_basis_quantum_resources.json')
    moment_checks=[]
    for comp in d['selected_comparisons']:
        # The frozen protocol gives the same basis for all three competitors.
        key=comp['model_id'],comp['methods']['shared_threshold']['candidate_id'] if 'candidate_id' in comp['methods']['shared_threshold'] else None
        if key[1] is None:
            sr=next(v for v in d['records'] if v['model_id']==key[0] and v['method']=='shared_threshold')
            key=key[0],sr['candidate_id']
        prep=next(v for v in factor['preparations'] if (v['original_model_id'],v['candidate_id'])==key)
        EU=F(prep['finite_moments']['centred_U_moments']['1']['upper_exact_fraction'])
        value=min(F(1),101*EU/(1+EU))
        assert value==1
        moment_checks.append({'model_id':key[0],'candidate_id':key[1],'finite_U_first_moment_upper':str(EU),'Jensen_gap_variation_upper':str(value)})
        mark('moment_only_variation_no_gain_checks')
    return moment_checks

def dominance_algebra():
    # Exact affine identities and finite cases support the all-integer proof.
    for ac in range(1,9):
        for ab in range(1,9):
            for sc in range(0,8):
                for sb in range(0,8):
                    cond=ac<=ab and 2*ac+sc<=2*ab+sb
                    r=max(F(ac,ab),F(2*ac+sc,2*ab+sb))
                    for d in (1,2,3,4,5,6,7,8,31,32,2047,2048,100001,100002):
                        tc=d*ac+((d-1)//2)*sc;tb=d*ab+((d-1)//2)*sb
                        assert F(tc,tb)<=r
                        if cond:
                            assert tc<=tb
                        mark('exact_all_degree_affine_ratio_checks')
                    if not cond:
                        witness=1 if ac>ab else 2*(max(0,ab-ac)//((2*ac+sc)-(2*ab+sb))+1)+1
                        assert witness*ac+((witness-1)//2)*sc>witness*ab+((witness-1)//2)*sb
                        mark('exact_necessity_counterexample_checks')

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--resources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();norms=read(a.repo/'results/gaussian_derivative_norms.json')
    nu={k:F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==k and r['passed'])) for k in (0,1)}
    assert endpoint(normal_density(arb(0)),True)<=nu[0]
    assert endpoint(normal_density(arb(1)),True)<=nu[1]
    counterexample=analytic_tests(nu[0],nu[1]);dominance_algebra()
    moments=amended_checks(a.repo,a.resources,nu[0],nu[1])
    out={'schema':'qbp-common-information-independent-mathematics-v1','status':'passed',
      'checker_sha256':sha(__file__),'resources_sha256':sha(a.resources),
      'norms_sha256':sha(a.repo/'results/gaussian_derivative_norms.json'),
      'frozen_proof_sha256':sha(Path(__file__).with_name('common_information_stability_fragment.tex')),
      'presentation_proof_sha256':sha(Path(__file__).with_name('common_information_stability_presentation_fragment.tex')),
      'reflection_sampler_proof_sha256':sha(Path(__file__).with_name('reflection_sampler_dominance_fragment.tex')),
      'counts':COUNTS,'dependent_inner_shift_counterexample':counterexample,'moment_only_variation_screen':moments,
      'scope':'General full proofs plus fresh directed-Arb finite stress checks, exact affine sampler identities and independent retained-field reconstruction of all416 mathematical allowances, budgets and priors. No producer/graph-builder/frontend imports. Intrinsic scalar ordinate correctness and existing whole-graph guards remain bound to their earlier independent audits. No claim that finite stress checks replace proofs, no physical execution or theoretical priority.'}
    assert not a.output.exists();a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'counts':COUNTS},indent=2),flush=True)

if __name__=='__main__':
    main()
