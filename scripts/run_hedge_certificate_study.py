"""Audit retained synthetic static hedges with exact-input interval certificates.

This is an informed reconstruction following interrupted development. Earlier
outcomes are known; the study is not described as untouched preregistration.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import mpmath as mp


ABSOLUTE_WIDTH='1e-50'
RELATIVE_WIDTH='1e-12'
INTEGRAL_RELATIVE_TOLERANCE='1e-70'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def binary(value):
    return F.from_float(float(value))


def bounds(record):
    return F(record['lower_exact_fraction']),F(record['upper_exact_fraction'])


def output_destination(output,root):
    output,root=output.expanduser().resolve(),root.resolve()
    if output.is_relative_to(root):
        raise ValueError('Choose an output directory outside the repository.')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError('Choose a new or empty output directory; evidence cannot be overwritten.')
    return output


def overall_status(manifest):
    return (not manifest['failures'] and manifest['sources_unchanged']
            and manifest['fixed_risk_rows']==168 and manifest['refined_rows']==126
            and manifest['slsqp_rows']==126 and manifest['positive_integral_rows']==60)


def run(args):
    root=args.project_root.resolve()
    output=output_destination(args.output,root)
    module_path=args.module.resolve()
    references_path=module_path.with_name('gbm_hedge_references.py')
    source_paths=[module_path,references_path,Path(__file__).resolve(),root/'src/gbm_price_certificate.py']
    data=root/'results/capped_gbm/hedging'
    source_paths.extend(data/name for name in ('study_manifest.json','constrained_mse.csv','unconstrained.csv','nohedge.csv'))
    before={str(path):sha(path) for path in source_paths}
    for directory in (str(module_path.parent),str(root/'src')):
        if directory not in sys.path:sys.path.insert(0,directory)
    engine=load(module_path,'gbm_hedge_certificate_audit_engine')
    reference=load(references_path,'gbm_hedge_audit_references')
    manifest_original=json.loads((data/'study_manifest.json').read_text())
    contracts={entry['case']:entry['parameters'] for entry in manifest_original['contracts']}
    tables={name:list(csv.DictReader((data/name).open())) for name in ('constrained_mse.csv','unconstrained.csv','nohedge.csv')}
    output.mkdir(parents=True,exist_ok=True)
    protocol={
        'schema':'qbp-validated-hedge-protocol-v1','reconstruction_date':'2026-10-20',
        'design_status':'Informed reconstruction; interrupted temporary development and its outcomes were known. Fresh numerical verification follows unchanged declared cases/tolerances.',
        'target':'Exact model risk and fair-price constraints at original binary64 primitive parameters, not exact printed decimals.',
        'retained_cases':7,'physical_laws_per_case':3,'budget_multipliers':[0,.25,.5,1,2,4],
        'expected_fixed_holding_rows':168,'expected_constrained_rows':126,'expected_slsqp_rows':126,
        'positive_integral_selection':'All21 no-hedge,21 unrestricted and9 archived primary-floor holdings; additionally9 refined floor holdings, total60.',
        'interval_absolute_width':ABSOLUTE_WIDTH,'interval_relative_width_and_refined_gap':RELATIVE_WIDTH,
        'reference_digits':160,'positive_integral_digits':100,'positive_integral_relative_tolerance':INTEGRAL_RELATIVE_TOLERANCE,
        'qp_reference_feasibility_allowance':'1e-130 scaled; diagnostic reference only, not a true feasibility certificate.',
        'archive_policy':'Preserve all original rows; report original, robustly repaired and refined portfolios separately. No negative-risk clipping.',
        'price_accounting':'New certificates use exact transcendental fair prices; archived integration used computed binary64 quotes, so archived values need not lie in new exact-model risk intervals.',
        'baseline_interpretation':'No-hedge and unrestricted are retained. The latter is a mathematical L2 relaxation, not necessarily an admissible trading strategy.',
        'numerical_reference_scope':'Independent high-precision piecewise moments and positive integration under the same stipulated GBM; no empirical market or Monte Carlo tail validation.',
        'source_hashes_before':before,
        'environment':{'python':platform.python_version(),'mpmath':mp.__version__,'python_flint':engine.flint.__version__},
    }
    write(output/'protocol.json',protocol)
    failures=[];fixed=[];refined=[];slsqp=[];integrals=[]
    model_cache={};contract_cache={}
    def target(row):
        name=row['case'];mu=binary(row['physical_total_return'])
        primitive=contracts[name]
        if name not in contract_cache:
            contract_cache[name]=engine.ExactCappedGBM(**{k:binary(v) for k,v in primitive.items() if k!='volatility'})
        c=contract_cache[name];sigma=binary(primitive['volatility'])
        key=(name,mu)
        if key not in model_cache:model_cache[key]=reference.parameters(c,sigma,mu)
        return c,sigma,mu,model_cache[key]
    def mp_fraction(value):return reference.number(F(value))
    def contains(interval,value):
        lo,hi=bounds(interval)
        return mp_fraction(lo)<=value<=mp_fraction(hi)
    def display(value):return mp.nstr(value,175)
    def integral_check(key,model,theta,reference_risk,interval):
        try:
            with mp.workdps(100):value=reference.positive_risk(model,theta)
            relative=abs(value-reference_risk)/abs(reference_risk)
            record={'key':key,'reference_risk':display(reference_risk),'positive_integral':display(value),
                    'relative_error':display(relative),'relative_agreement':relative<=mp.mpf(INTEGRAL_RELATIVE_TOLERANCE),
                    'integral_contained':contains(interval,value)}
            if not record['relative_agreement'] or not record['integral_contained']:failures.append(key+':positive_integral')
        except Exception as error:
            record={'key':key,'error':{'type':type(error).__name__,'message':str(error)}}
            failures.append(key+':positive_integral_error')
        integrals.append(record)
    with mp.workdps(180):
        for filename,rows in tables.items():
            for index,row in enumerate(rows):
                key=f'{filename}:{index}'
                try:
                    c,sigma,mu,model=target(row)
                    theta=tuple(binary(row[name]) for name in ('stock_holding','call_holding'))
                    certificate=engine.certified_hedge_risk(c,sigma,mu,theta,
                                     absolute_tolerance=ABSOLUTE_WIDTH,relative_tolerance=RELATIVE_WIDTH)
                    with mp.workdps(160):risk_reference=reference.fixed_risk(model,theta)
                    passed=contains(certificate['risk'],risk_reference)
                    if not passed:failures.append(key+':fixed_risk')
                    record={'source_file':filename,'source_row_zero_based':index,'case':row['case'],
                            'physical_total_return_binary64_exact':str(mu),
                            'budget_multiplier':row.get('budget_multiplier'),
                            'archived_risk_integral':row['population_loss_nonnegative_integral'],
                            'archived_risk_quadratic':row['population_loss_moment_quadratic'],
                            'archived_negative_quadratic':row['negative_moment_quadratic']=='True',
                            'archived_floor':row['risk_below_heuristic_roundoff_floor']=='True',
                            'certificate':certificate,'independent_risk_160digits':display(risk_reference),
                            'reference_contained':passed}
                    fixed.append(record)
                    if filename!='constrained_mse.csv' or record['archived_floor']:
                        integral_check(key+':fixed',model,theta,risk_reference,certificate['risk'])
                    if filename!='constrained_mse.csv':continue
                    kappa=binary(row['budget_multiplier'])
                    multipliers=tuple(binary(x) for x in json.loads(row['kkt_multipliers']))
                    cert=engine.certify_hedge(c,sigma,mu,kappa,theta,multipliers,
                              absolute_tolerance=ABSOLUTE_WIDTH,relative_tolerance=RELATIVE_WIDTH)
                    newtheta=tuple(F(x) for x in cert['holdings_exact'])
                    with mp.workdps(160):
                        qp=reference.optimum(model,kappa)
                        newrisk=reference.fixed_risk(model,newtheta)
                    lower=F(cert['optimal_loss_lower_bound_exact']);upper=F(cert['risk']['upper_exact_fraction'])
                    gap=F(cert['optimality_gap_upper_exact'])
                    eigen=F(cert['minimum_eigenvalue_lower_exact'])
                    distance=sum((mp_fraction(newtheta[i])-qp['holdings'][i])**2 for i in range(2))
                    distance_bound=F(cert['squared_distance_to_unique_optimum_upper_exact'])
                    checks={'robust_feasibility':bool(cert['robust_feasibility']),
                            'refined_risk_reference_contained':contains(cert['risk'],newrisk),
                            'optimum_reference_contained':mp_fraction(lower)<=qp['risk']<=mp_fraction(upper),
                            'gap_rational_identity':upper-lower==gap,
                            'relative_gap':gap<=F(RELATIVE_WIDTH)*F(cert['risk']['lower_exact_fraction']),
                            'holding_distance_bound':distance<=mp_fraction(distance_bound),
                            'distance_rational_identity':distance_bound==gap/eigen,
                            'feasible_rational_spend':F(cert['robust_spend_upper_exact'])<=F(cert['robust_budget_lower_exact'])}
                    if not all(checks.values()):failures.append(key+':refined_check')
                    refined.append({'source_row_zero_based':index,'case':row['case'],'physical_total_return':row['physical_total_return'],
                                    'budget_multiplier':row['budget_multiplier'],'archived_floor':record['archived_floor'],
                                    'certificate':cert,'independent_qp':{'risk':display(qp['risk']),
                                    'holdings':[display(x) for x in qp['holdings']],
                                    'distance_squared':display(distance)},'checks':checks})
                    if record['archived_floor']:
                        integral_check(key+':refined',model,newtheta,newrisk,cert['risk'])
                    s_theta=tuple(binary(row[name]) for name in ('slsqp_stock_holding','slsqp_call_holding'))
                    s_cert=engine.certify_hedge(c,sigma,mu,kappa,s_theta,refine=False,
                              absolute_tolerance=ABSOLUTE_WIDTH,relative_tolerance=RELATIVE_WIDTH)
                    repaired=tuple(F(x) for x in s_cert['repaired_original_holdings_exact'])
                    with mp.workdps(160):
                        s_raw_ref=reference.fixed_risk(model,s_theta)
                        s_repaired_ref=reference.fixed_risk(model,repaired)
                    s_checks={'raw_risk_reference_contained':contains(s_cert['original_risk'],s_raw_ref),
                              'repaired_risk_reference_contained':contains(s_cert['repaired_original_risk'],s_repaired_ref),
                              'repaired_feasible':s_cert['robust_feasibility']}
                    if not all(s_checks.values()):failures.append(key+':slsqp')
                    slsqp.append({'source_row_zero_based':index,'case':row['case'],'physical_total_return':row['physical_total_return'],
                                  'budget_multiplier':row['budget_multiplier'],
                                  **{name:s_cert[name] for name in ('original_holdings_exact','original_risk','original_budget_violation',
                                     'original_budget_status','original_box_violation_exact','original_feasibility_status',
                                     'repaired_original_holdings_exact','repaired_original_risk','robust_spend_upper_exact','robust_budget_lower_exact')},
                                  'checks':s_checks})
                except Exception as error:
                    failures.append(key+':exception')
                    fixed.append({'source_file':filename,'source_row_zero_based':index,
                                  'error':{'type':type(error).__name__,'message':str(error)}})
        write(output/'fixed_risk_audit.json',fixed)
        write(output/'refined_hedges.json',refined)
        write(output/'slsqp_audit.json',slsqp)
        write(output/'positive_integrals.json',integrals)
        after={str(path):sha(path) for path in source_paths}
        def maximum_fraction(values):return str(max(values,default=F(0)))
        primary=[row for row in fixed if row['source_file']=='constrained_mse.csv' and 'certificate' in row]
        summary={
            'schema':'qbp-validated-hedge-validation-v1','fixed_risk_rows':sum('certificate' in row for row in fixed),
            'refined_rows':len(refined),'slsqp_rows':len(slsqp),'positive_integral_rows':len(integrals),
            'failures':failures,'sources_unchanged':before==after,'source_hashes_before':before,'source_hashes_after':after,
            'archived_primary_floor_rows':sum(row['archived_floor'] for row in primary),
            'archived_primary_negative_quadratics':sum(row['archived_negative_quadratic'] for row in primary),
            'archived_all_negative_quadratics':sum(row.get('archived_negative_quadratic',False) for row in fixed),
            'original_primary_budget_statuses':{status:sum(row['certificate']['original_budget_status']==status for row in refined)
                                              for status in ('certified_feasible','certified_infeasible','unresolved')},
            'original_primary_joint_feasibility_statuses':{status:sum(row['certificate']['original_feasibility_status']==status for row in refined)
                                                          for status in ('certified_feasible','certified_infeasible','unresolved')},
            'maximum_original_budget_violation_upper_exact':maximum_fraction(F(row['certificate']['original_budget_violation']['upper_exact_fraction']) for row in refined),
            'maximum_original_box_violation_exact':maximum_fraction(F(row['certificate']['original_box_violation_exact']) for row in refined),
            'maximum_relative_optimality_gap_exact':maximum_fraction(F(row['certificate']['optimality_gap_upper_exact'])/F(row['certificate']['risk']['lower_exact_fraction']) for row in refined),
            'floor_results':[{'case':row['case'],'physical_total_return':row['physical_total_return'],'budget_multiplier':row['budget_multiplier'],
                             'risk':row['certificate']['risk'],'original_risk':row['certificate']['original_risk'],
                             'repaired_original_risk':row['certificate']['repaired_original_risk'],
                             'optimal_loss_lower_bound_exact':row['certificate']['optimal_loss_lower_bound_exact'],
                             'optimality_gap_upper_exact':row['certificate']['optimality_gap_upper_exact']} for row in refined if row['archived_floor']],
            'interpretation':'Model-conditional exact-input arithmetic/optimisation certificates; high-precision references remain diagnostics. No empirical, Monte Carlo rare-tail confirmation or quantum advantage claim.',
        }
        summary['overall_numerical_verification']=overall_status(summary)
        summary['files']={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(output.iterdir()) if p.is_file()}
        write(output/'validation.json',summary)
    print(json.dumps({k:summary[k] for k in ('overall_numerical_verification','fixed_risk_rows','refined_rows','slsqp_rows','positive_integral_rows','failures','archived_primary_floor_rows','archived_primary_negative_quadratics','archived_all_negative_quadratics','original_primary_budget_statuses','original_primary_joint_feasibility_statuses')},indent=2))
    return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--module',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.module is None:args.module=args.project_root/'src/gbm_hedge_certificate.py'
    return 0 if run(args)['overall_numerical_verification'] else 1


if __name__=='__main__':raise SystemExit(main())
