"""Independent financial and resource inequalities for the complete records."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from audit_estimation_prior import independent_mean, density_reference, REF_ERROR

D=Decimal
F=Fraction


def rational(record):
    return F(record['numerator'],record['denominator'])


def main():
    path=ROOT/'results/complete_estimation_resources.json'
    data=json.loads(path.read_text())
    config_path=ROOT/'results/architecture_resource_candidates.json'
    assert data['configuration_sha256']==hashlib.sha256(config_path.read_bytes()).hexdigest()
    for name,expected in data['source_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected
    config=json.loads(config_path.read_text())
    models={m['model_id']:m for m in config['models']}
    key=lambda r:(r['case_id'],r['method'],r.get('architecture',r.get('resources',{}).get('architecture')),r['q'],r['exp_address_bits'],r['log_address_bits'],r['probability_address_bits'])
    source={key(r):r for r in config['configurations']}
    assert len(config['configurations'])==len(data['configurations'])==819
    assert len(source)==len(config['configurations'])
    assert Counter(map(key,data['configurations']))==Counter(map(key,config['configurations']))
    assert len({key(r) for r in data['configurations']})==len(data['configurations'])
    phi0=density_reference();records=[]
    with localcontext() as ctx:
        ctx.prec=300
        for row in data['configurations']:
            r=source[key(row)]
            rr=r['resources'];model=models[r['model_id']]
            ref=independent_mean(model,row['cap'],phi0)
            p=row['financial_prior_and_classical_single_price']
            numerical=D(p['decimal_numerical_allowance'])
            assert abs(ref['ag']-D(p['arithmetic_geometric_mean_currency_centre_decimal']))+REF_ERROR<=numerical
            assert abs(ref['slope_mean']-D(p['slope_Jensen_mean_currency_centre_decimal']))+REF_ERROR<=numerical
            expected_u=min(D(1),(ref['mean']+D.from_float(r['arithmetic_price_error_upper'])+D.from_float(r['distribution_price_bound']))/D(row['cap']))
            assert D(p['finite_mean_probability_upper_decimal'])+REF_ERROR>=expected_u
            u=D(p['finite_mean_probability_upper_decimal'])
            assert 0<=u<=1
            assert D(p['finite_variance_upper_decimal'])+REF_ERROR>=(u*(1-u) if u<=D('.5') else D('.25'))
            joint=row['classical_simultaneous_27_prices']
            for quantity in ('finite_mean_probability_upper_decimal','finite_variance_upper_decimal','statistical_probability_target_decimal'):
                assert joint[quantity]==p[quantity]
            for field in ('financial_prior_and_classical_single_price','classical_simultaneous_27_prices'):
                cp=row[field];v=D(cp['finite_variance_upper_decimal']);epsilon=D('.005')/D(row['cap']);n=cp['classical_certified_samples']
                assert D(cp['statistical_probability_target_decimal'])==epsilon
                failure=(2*(-D(n)*epsilon*epsilon/(2*v+2*epsilon/3)).exp() if cp['classical_selected_bound']=='Bernstein' else 2*(-2*D(n)*epsilon*epsilon).exp())+REF_ERROR
                assert failure<=D(cp['failure_probability'])
                if field.startswith('classical_simultaneous'):assert D(cp['failure_probability'])<=D(1)/2700
                else:assert D(cp['failure_probability'])==D('.01')
            κ={'exact7t':14,'unitary8t':8,'measured4t':4}[row['architecture']]
            assert rr['logical_T_count']==κ*rr['lookup_selector_internal_nodes']+7*rr['arithmetic_exact_Toffoli_count']
            for field in ('quantum_single_price_prior','quantum_single_price_full_range','quantum_simultaneous_27_prices_prior','quantum_simultaneous_27_prices_full_range'):
                budget=row[field];plan=budget['plan'];alpha=rational(plan['alpha']);ep=rational(plan['epsilon'])
                assert ep==F(1,200*row['cap'])
                assert alpha==(F(1,100) if field.startswith('quantum_single') else F(1,2700))
                prior=rational(plan['prior_upper'])
                assert prior==(F(1) if field.endswith('full_range') else F(p['finite_mean_probability_upper_decimal']))
                guard=min(F(1),prior+F(1,10**18))
                assert rational(plan['guarded_prior_upper'])==guard
                var=guard*(1-guard) if guard<=F(1,2) else F(1,4)
                gamma=rational(plan['amplitude_angle_lipschitz_upper'])
                assert gamma*gamma>=4*var
                R=plan['maximum_rounds'];N=plan['shots_per_round'];Q=budget['grover_invocations_upper'];shots=budget['flag_measurements_upper']
                if R:
                    assert 2*(-2*D(N)/625).exp()+REF_ERROR<=D(alpha.numerator)/D(alpha.denominator)/D(R)
                    m0=plan['initial_multiplier'];target=plan['target_multiplier'];M=plan['multiplier_upper']
                    assert m0>=1 and m0%2==1
                    t=gamma/(5*ep)
                    assert target==max(1,-((-t.numerator)//t.denominator))
                    assert m0*3**(R-1)>=target
                    if R>1:assert m0*3**(R-2)<target
                    assert M==(m0 if m0>=target else min(5*target,m0*5**(R-1)))
                    assert plan['grover_power_upper']==(M-1)//2
                    sum_max=(m0*((5**R-1)//4)-R)//2
                    K=((m0-1)//2 if R==1 else min((3*M)//4,sum_max))
                    assert Q==N*K==plan['total_grover_calls_upper']
                    assert rational(plan['amplitude_angle_lipschitz_upper'])/(3*target)+F(2,10**70)<2*ep
                else:assert shots==Q==0 and prior<=2*ep
                calls=shots+2*Q
                assert budget['clean_A_and_adjoint_invocations_upper']==calls
                assert budget['T_count_upper']==shots*rr['logical_T_count']+Q*(2*rr['logical_T_count']+rr['grover_zero_reflection']['T_count_exact_7T'])
                gates=budget['assigned_Clifford_T_and_instrument_operation_counts_upper']
                assert gates['T']+gates['TD']==budget['T_count_upper']
                assert budget['all_measurements_upper']==shots+calls*rr['measurement_count']
                assert gates.get('MZ',0)==budget['lookup_measurements_upper']==calls*rr['measurement_count']
                assert budget['conditional_CZ_feedback_upper']==calls*rr['classical_feedback_CZ_count']
                assert gates.get('CZ_IF',0)==budget['conditional_CZ_feedback_upper']
                assert budget['conditional_lookup_resets_upper']==calls*rr['classical_reset_X_count']
                assert gates.get('X_IF',0)==budget['conditional_lookup_resets_upper']
                assert shots==R*N
                assert shots==plan['total_shots_upper']
            records.append({'case_id':row['case_id'],'method':row['method'],'architecture':row['architecture'],'independent_moment_prior_and_classical_confidence':True,'quantum_confidence_and_integer_cap_inequalities':True,'all_T_measurement_feedback_charges_match':True})
    groups={}
    for row in data['configurations']:
        groups.setdefault((row['case_id'],row['architecture']),[]).append(row)
    comparisons={(r['case_id'],r['architecture']):r for r in data['method_comparison']}
    assert len(comparisons)==len(data['method_comparison'])==len(groups)==81
    assert set(comparisons)==set(groups)
    for group,rows in groups.items():
        c=comparisons[group]
        for field,minimum_key,tie_key in (
            ('quantum_single_price_prior','minimum_complete_T_upper','tied_methods'),
            ('quantum_simultaneous_27_prices_prior','minimum_simultaneous_27_complete_T_upper','simultaneous_27_tied_methods')):
            minimum=min(r[field]['T_count_upper'] for r in rows)
            assert c[minimum_key]==minimum
            assert c[tie_key]==sorted({r['method'] for r in rows if r[field]['T_count_upper']==minimum})
        arms={r['method']:r for r in c['method_minima']}
        assert len(arms)==len(c['method_minima'])==4
        for method,a in arms.items():
            pool=[r for r in rows if r['method']==method]
            for field,armfield in (('quantum_single_price_prior','single_price_prior'),('quantum_simultaneous_27_prices_prior','simultaneous_27_prior')):
                minimum=min(r[field]['T_count_upper'] for r in pool)
                assert a[armfield]['T_upper']==minimum
                pref=lambda r:tuple(r[k] for k in ('q','exp_address_bits','log_address_bits','probability_address_bits'))
                exact=Counter(pref(r) for r in pool if r[field]['T_count_upper']==minimum)
                assert Counter(pref(r) for r in a[armfield]['tied_prefixes'])==exact
                winning=[r for r in pool if r[field]['T_count_upper']==minimum]
                expected=min(winning,key=lambda r:(r[field]['all_measurements_upper'],r[field]['peak_logical_qubits_upper'],r['q'],r['exp_address_bits'],r['log_address_bits'],r['probability_address_bits']))
                assert pref(a[armfield]['representative_prefix'])==pref(expected)
    out={'status':'independent complete financial-confidence and resource inequalities passed',
        'resource_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'configurations':records,'checked_configurations':len(records),
        'scope':'mathematical and finite-logical record consistency; no hardware, timing or originality certification'}
    (ROOT/'results/complete_estimation_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Independent complete budgets passed:',len(records))


if __name__=='__main__':main()
