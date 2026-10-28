"""Bounded synthetic audit of compact-drift robust static hedging.

The protocol is retained before numerical outcomes. All output must be new and
outside the repository. Independent numerical references do not replace proofs.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import numpy as np
import scipy


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def binary(value):return F.from_float(float(value))


def output_destination(output,root):
    output,root=output.expanduser().resolve(),root.resolve()
    if output.is_relative_to(root):raise ValueError('Choose an output directory outside the repository.')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError('Choose a new or empty output directory; retained evidence cannot be overwritten.')
    return output


def overall_status(manifest):
    return (manifest.get('primary_rows')==12 and manifest.get('confidence_rows')==30
            and manifest.get('dense_checks')==12*5*257 and manifest.get('finite_grid_reference_rows')==12
            and manifest.get('bridge_rows')==3 and manifest.get('bridge_dense_checks')==3*4*257
            and manifest.get('sources_unchanged') is True and not manifest.get('failures'))


def run(args):
    root=args.project_root.resolve();output=output_destination(args.output,root)
    module_path=args.module.resolve();reference_path=module_path.with_name('gbm_drift_references.py')
    protocol_path=args.protocol.resolve()
    for directory in (str(root/'src'),str(module_path.parent)):
        if directory not in sys.path:sys.path.insert(0,directory)
    e=load(module_path,'drift_audit_engine');ref=load(reference_path,'drift_audit_reference')
    import gbm_hedge_certificate as fixed
    import gbm_hedge_references as static
    protocol=json.loads(protocol_path.read_text())
    source_paths=[module_path,reference_path,protocol_path,Path(__file__).resolve(),
                  root/'src/gbm_price_certificate.py',root/'src/gbm_hedge_certificate.py',
                  root/'src/gbm_hedge_references.py',root/'results/capped_gbm/hedging/study_manifest.json']
    before={str(path):sha(path) for path in source_paths}
    original=json.loads((root/'results/capped_gbm/hedging/study_manifest.json').read_text())
    chosen={entry['case']:entry['parameters'] for entry in original['contracts'] if entry['case'] in protocol['contracts']}
    output.mkdir(parents=True,exist_ok=True)
    registration={'protocol':protocol,'protocol_sha256':sha(protocol_path),'source_hashes_before':before,
                  'environment':{'python':platform.python_version(),'numpy':np.__version__,
                                 'scipy':scipy.__version__,'mpmath':mp.__version__,'python_flint':fixed.flint.__version__}}
    write(output/'registration.json',registration)
    failures=[];primary=[];confidence=[];bridge=[];dense_count=0;reference_count=0;bridge_dense_count=0
    display=lambda v:mp.nstr(v,175)
    mpf=lambda x:ref.number(F(x))
    centre=F(protocol['drift_centre']);true_mu=F('3/100')
    with mp.workdps(180):
        for name in protocol['contracts']:
            primitive=chosen[name]
            contract=fixed.ExactCappedGBM(**{k:binary(v) for k,v in primitive.items() if k!='volatility'})
            sigma=binary(primitive['volatility'])
            for width in map(F,protocol['drift_half_widths']):
                lower,upper=centre-width,centre+width
                dense_nodes=tuple(lower+(upper-lower)*i/(protocol['dense_drift_points']-1) for i in range(protocol['dense_drift_points']))
                dense_models=ref.models(contract,sigma,dense_nodes)
                true_model=static.parameters(contract,sigma,true_mu)
                for kappa in map(F,protocol['budget_multipliers']):
                    key=f'{name}:width={width}:budget={kappa}'
                    record={'key':key,'case':name,'drift_interval_exact':[str(lower),str(upper)],
                            'budget_multiplier_exact':str(kappa),'true_drift_exact':str(true_mu)}
                    try:
                        cert=e.certify_drift_robust_hedge(contract,sigma,lower,upper,kappa,
                            relative_gap_tolerance=protocol['relative_gap_tolerance'],
                            cell_inflation_tolerance=protocol['cell_inflation_tolerance'],
                            initial_cells=protocol['initial_cells'],max_nodes=protocol['max_nodes'],precision_bits=protocol['precision_bits'])
                        record['certificate']=cert
                        if cert['status']!='certified' or not cert['robust_feasibility']:
                            failures.append(key+':unresolved_certificate')
                        theta=tuple(F(v) for v in cert['holdings_exact'])
                        grid=tuple(F(v) for v in cert['node_grid_exact'])
                        qp=ref.finite_grid_optimum(contract,sigma,grid,kappa,digits=160)
                        reference_count+=1
                        lo=F(cert['optimal_loss_lower_bound_exact']);hi=F(cert['continuum_upper_exact'])
                        grid_models=ref.models(contract,sigma,grid)
                        grid_values=[ref.risk(model,theta) for model in grid_models]
                        grid_risk=max(grid_values)
                        node_checks=[mpf(interval['lower_exact_fraction'])<=value<=mpf(interval['upper_exact_fraction'])
                                     for interval,value in zip(cert['node_risks'],grid_values)]
                        weights=tuple(F(x) for x in cert['mixture_weights_exact'])
                        multipliers=tuple(F(x) for x in cert['financial_multipliers_exact'])
                        independent_dual=ref.mixture_dual(contract,sigma,grid,weights,kappa,multipliers)
                        raw_dual=cert['raw_mixture_dual']
                        dual_check=mpf(raw_dual['lower_exact_fraction'])<=independent_dual<=mpf(raw_dual['upper_exact_fraction'])
                        if len(node_checks)!=len(grid) or not all(node_checks):failures.append(key+':independent_node_risk_containment')
                        if not dual_check:failures.append(key+':independent_mixture_dual_containment')
                        if sum(weights,F(0))!=1 or min(weights)<0 or min(multipliers)<0:failures.append(key+':invalid_mixture_witness')
                        record['independent_certificate_checks']={'all_node_risks_contained':all(node_checks),'node_risks_checked':len(node_checks),
                            'mixture_dual_contained':dual_check,'mixture_dual_160digits':display(independent_dual)}
                        comparison=mpf(lo)<=qp['risk']<=grid_risk+mp.mpf('1e-100')
                        if not comparison:failures.append(key+':finite_grid_reference_bracketing')
                        record['independent_finite_grid_reference']={
                            'risk':display(qp['risk']),'holdings':[display(x) for x in qp['holdings']],
                            'kkt_residual':qp['kkt_residual'],'refinement':qp['refinement'],
                            'active_nodes':qp.get('active_nodes',[]),'active_faces':qp.get('active_faces',[]),
                            'certificate_lower_to_feasible_grid_risk_check':comparison}
                        centre_cert=fixed.certify_hedge(contract,sigma,centre,kappa,(0,0))
                        oracle_cert=fixed.certify_hedge(contract,sigma,true_mu,kappa,(0,0))
                        endpoint=ref.finite_grid_optimum(contract,sigma,(lower,upper),kappa,digits=160)
                        endpoint_theta=tuple(F(str(x)) for x in endpoint['holdings'])
                        endpoint_repaired=fixed.certify_hedge(contract,sigma,centre,kappa,endpoint_theta,refine=False)
                        portfolios={'robust':theta,'centre':tuple(F(v) for v in centre_cert['holdings_exact']),
                                    'oracle_true_drift':tuple(F(v) for v in oracle_cert['holdings_exact']),
                                    'no_hedge':(F(0),F(0)),
                                    'endpoint_only':tuple(F(v) for v in endpoint_repaired['holdings_exact'])}
                        comparisons={}
                        for label,holdings in portfolios.items():
                            values=[ref.risk(model,holdings) for model in dense_models]
                            dense_count+=len(values)
                            worst=max(values);imax=values.index(worst)
                            comparisons[label]={'holdings_exact':[str(v) for v in holdings],
                                'true_drift_risk':display(ref.risk(true_model,holdings)),
                                'dense_maximum':display(worst),'dense_argmax_drift_exact':str(dense_nodes[imax]),
                                'endpoint_maximum':display(max(values[0],values[-1])),
                                'interior_excess_over_endpoints':display(worst-max(values[0],values[-1])),
                                'grid_points':len(values)}
                            if label=='robust' and worst>mpf(hi):failures.append(key+':continuum_upper')
                        record['comparisons']=comparisons
                        record['robust_worse_at_true_drift_than_centre']=(mp.mpf(comparisons['robust']['true_drift_risk'])>mp.mpf(comparisons['centre']['true_drift_risk']))
                        record['robust_dense_worst_less_than_centre']=(mp.mpf(comparisons['robust']['dense_maximum'])<mp.mpf(comparisons['centre']['dense_maximum']))
                    except Exception as error:
                        record['error']={'type':type(error).__name__,'message':str(error)}
                        failures.append(key+':error')
                    primary.append(record)
        cfg=protocol['confidence_diagnostics'];rng=np.random.default_rng(cfg['seed'])
        n=cfg['observations_per_dataset'];delta=F(cfg['observation_spacing']);mu=F(cfg['true_drift'])
        sigma=F(cfg['known_volatility']);q=F(cfg['dividend']);radius=F(1,2*10**cfg['rounding_decimal_places'])
        for index in range(cfg['datasets']):
            record={'dataset_zero_based':index}
            try:
                latent=(float(mu-q-sigma*sigma/2)*float(delta)+float(sigma)*np.sqrt(float(delta))*rng.normal(size=n))
                observed=tuple(F(str(round(float(x),cfg['rounding_decimal_places']))) for x in latent)
                latent_mean=sum((F.from_float(float(x)) for x in latent),F(0))/n
                cert=e.certified_known_sigma_drift_interval(observed,sigma,delta,dividend=q,
                    rounding_radii=[radius]*n,alpha=F(cfg['alpha']),precision_bits=protocol['precision_bits'])
                record['certificate']=cert
                if cert['status']!='certified':failures.append(f'confidence:{index}:unresolved')
                lower,upper=(F(cert['drift_interval'][x]) for x in ('lower_exact_fraction','upper_exact_fraction'))
                observed_mean=F(cert['observed_log_return_mean_exact']);mean_radius=F(cert['mean_rounding_radius_exact'])
                sample_lower,sample_upper=observed_mean-mean_radius,observed_mean+mean_radius
                rounding_check=sample_lower<=latent_mean<=sample_upper
                if not rounding_check:failures.append(f'confidence:{index}:rounding_inclusion')
                record.update({'contains_true_drift':lower<=mu<=upper,'latent_sample_mean_contained':rounding_check,
                               'latent_sample_mean_exact':str(latent_mean),
                               'raw_increments_retained':False})
            except Exception as error:
                record['error']={'type':type(error).__name__,'message':str(error)};failures.append(f'confidence:{index}:error')
            confidence.append(record)
    # The amendment selects indices 0,1,2 before their bridge outcomes.
    bridge_cfg=protocol['inference_decision_bridge']
    primitive=chosen[bridge_cfg['contract']]
    bridge_contract=fixed.ExactCappedGBM(**{k:(F(cfg['dividend']) if k=='dividend' else binary(v))
        for k,v in primitive.items() if k!='volatility'})
    bridge_sigma=F(cfg['known_volatility']);bridge_mu=F(cfg['true_drift']);bridge_kappa=F(bridge_cfg['budget_multiplier'])
    with mp.workdps(180):
        for index in bridge_cfg['dataset_indices_zero_based']:
            record={'dataset_zero_based':index,'case':bridge_cfg['contract'],'true_drift_exact':str(bridge_mu)}
            try:
                inference=confidence[index]['certificate']
                if (bridge_contract.dividend!=F(inference['dividend_exact'])
                    or bridge_sigma!=F(inference['volatility_exact'])):
                    raise ValueError('Historical inference and future hedge must use identical known dividend and volatility inputs.')
                lower,upper=(F(inference['drift_interval'][x]) for x in ('lower_exact_fraction','upper_exact_fraction'))
                record['inference_interval_exact']=[str(lower),str(upper)]
                centre=F(inference['drift_centre_exact'])
                cert=e.certify_drift_robust_hedge(bridge_contract,bridge_sigma,lower,upper,bridge_kappa,
                    relative_gap_tolerance=protocol['relative_gap_tolerance'],
                    cell_inflation_tolerance=protocol['cell_inflation_tolerance'],
                    initial_cells=protocol['initial_cells'],max_nodes=protocol['max_nodes'],precision_bits=protocol['precision_bits'])
                record['certificate']=cert
                if cert['status']!='certified' or not cert['robust_feasibility']:
                    failures.append(f'bridge:{index}:unresolved_certificate')
                plug=fixed.certify_hedge(bridge_contract,bridge_sigma,centre,bridge_kappa,(0,0))
                oracle=fixed.certify_hedge(bridge_contract,bridge_sigma,bridge_mu,bridge_kappa,(0,0))
                portfolios={'robust':tuple(F(v) for v in cert['holdings_exact']),
                            'plug_in':tuple(F(v) for v in plug['holdings_exact']),
                            'oracle_true_drift':tuple(F(v) for v in oracle['holdings_exact']),
                            'no_hedge':(F(0),F(0))}
                dense_nodes=tuple(lower+(upper-lower)*F(i,256) for i in range(257))
                dense_models=ref.models(bridge_contract,bridge_sigma,dense_nodes)
                true_model=static.parameters(bridge_contract,bridge_sigma,bridge_mu)
                comparisons={}
                for label,theta in portfolios.items():
                    values=[ref.risk(model,theta) for model in dense_models];bridge_dense_count+=len(values)
                    true_risk=ref.risk(true_model,theta)
                    comparisons[label]={'holdings_exact':[str(x) for x in theta],
                        'true_drift_risk':display(true_risk),'dense_maximum':display(max(values)),
                        'dense_argmax_drift_exact':str(dense_nodes[values.index(max(values))]),'grid_points':len(values)}
                    if label=='robust':
                        if max(values)>mpf(cert['continuum_upper_exact']):failures.append(f'bridge:{index}:continuum_upper')
                        record['true_risk_exceeds_continuum_upper']=(true_risk>mpf(cert['continuum_upper_exact']))
                record['true_drift_in_interval']=(lower<=bridge_mu<=upper)
                if record['true_drift_in_interval'] and record['true_risk_exceeds_continuum_upper']:
                    failures.append(f'bridge:{index}:covered_true_drift_upper_bound')
                record['comparisons']=comparisons
                record['robust_worse_at_true_drift_than_plugin']=(mp.mpf(comparisons['robust']['true_drift_risk'])>mp.mpf(comparisons['plug_in']['true_drift_risk']))
            except Exception as error:
                record['error']={'type':type(error).__name__,'message':str(error)};failures.append(f'bridge:{index}:error')
            bridge.append(record)
    after={str(path):sha(path) for path in source_paths}
    manifest={'schema':'qbp-drift-hedge-audit-v1','study_date':'2026-10-24',
              'primary_rows':len(primary),'finite_grid_reference_rows':reference_count,
              'dense_checks':dense_count,'confidence_rows':len(confidence),
              'bridge_rows':len(bridge),'bridge_dense_checks':bridge_dense_count,
              'bridge_true_risk_understatements':sum(row.get('true_risk_exceeds_continuum_upper',False) for row in bridge),
              'bridge_robust_worse_at_true_drift_than_plugin':sum(row.get('robust_worse_at_true_drift_than_plugin',False) for row in bridge),
              'confidence_contains_true_drift':sum(row.get('contains_true_drift',False) for row in confidence),
              'robust_worse_at_true_drift_than_centre':sum(row.get('robust_worse_at_true_drift_than_centre',False) for row in primary),
              'robust_dense_worst_less_than_centre':sum(row.get('robust_dense_worst_less_than_centre',False) for row in primary),
              'sources_unchanged':before==after,'source_hashes_before':before,'source_hashes_after':after,
              'failures':failures,'scope':'Synthetic known-volatility GBM. Dense values and high-precision optimisers are diagnostics; continuum and primal-dual bounds use retained interval certificates. No market validation or claim of risk improvement at every drift.'}
    manifest['all_checks_pass']=overall_status(manifest)
    write(output/'primary.json',primary);write(output/'confidence.json',confidence);write(output/'bridge.json',bridge);write(output/'validation.json',manifest)
    write(output/'file_hashes.json',{p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()})
    return manifest


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root',type=Path,default=root)
    parser.add_argument('--module',type=Path,default=root/'src/gbm_drift_hedging.py')
    parser.add_argument('--protocol',type=Path,default=root/'configs/drift_hedge_protocol.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:manifest=run(args)
    except Exception as error:
        print(json.dumps({'all_checks_pass':False,'error':{'type':type(error).__name__,'message':str(error)}}));return 1
    print(json.dumps({k:manifest[k] for k in ('all_checks_pass','primary_rows','finite_grid_reference_rows','dense_checks','confidence_rows','bridge_rows','bridge_dense_checks','failures')},indent=2))
    return 0 if manifest['all_checks_pass'] else 1


if __name__=='__main__':raise SystemExit(main())
