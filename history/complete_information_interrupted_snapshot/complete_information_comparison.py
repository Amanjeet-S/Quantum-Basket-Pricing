"""Price-independent coefficient-informed recomparison of all 36 contracts.

The completed lower- and higher-dimensional scalar/circuit/precision pools
are retained unchanged. Only provable first-correction probability information
and analytical interval stopping are strengthened. Previous plans remain
admitted alongside the stronger plans. No pricing-answer file is imported.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,time
import correction_information_priors as information
import boundary_precision_selection as retained
from boundary_grid_certificate import interval,serial
ROOT=Path(__file__).resolve().parent
KINDS=retained.KINDS

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(name):
    p=ROOT/'results'/name;return json.load(gzip.open(p,'rt')) if p.suffix=='.gz' else json.loads(p.read_text())
def key(row):return row['model_id'],row['payoff'],None if row['cap'] is None else int(F(row['cap']))
def profiles(record,strategy):
    family=record['profiles']
    return family[strategy['implementation']][strategy['architecture']] if 'implementation' in strategy else family[strategy['architecture']]

def main():
    start=time.perf_counter();models={m['model_id']:m for m in load('generalisation_protocol.json')['models']}
    moment={mid:information.continuous_second_moment(m) for mid,m in models.items()}
    lower=load('boundary_financial_grid_transfer_summary.json');preparation=load('higherdim_preparation.json')
    coefficients={key(c):c for c in lower['coefficients']}
    coefficients.update({key(c):c['coefficients'] for r in preparation['records'] for c in r['cases']})
    financial={key(r):r for r in lower['records'] if r['order']==1}
    base_mass=preparation['shared_MGF_profile']['node_certificate']['mass_interval'];ml=F(base_mass['lower_exact_fraction']);mh=F(base_mass['upper_exact_fraction'])
    masses={p:serial(interval(ml,mh)**p) for p in (3,7)}
    for r in load('tensor_first_transfer.json')['records']:
        if models[r['model_id']]['outer_dimension'] in (3,7):
            tr=dict(r);tr['mass_interval']=masses[models[r['model_id']]['outer_dimension']];financial[key(r)]=tr
    names=['boundary_precision_comparison.json.gz','higherdim_precision_comparison.json.gz']
    families={name:load(name) for name in names};records=[];base_hashes={name:sha(ROOT/'results'/name) for name in names}
    for name,family in families.items():
        for index,record in enumerate(family['records']):
            row={'family_archive':name,'family_archive_sha256':base_hashes[name],'candidate_index':index,
                'model_id':record['model_id'],'payoff':record['payoff'],'cap':record['cap'],'kind':record['kind'],
                'base_graph_certificate_sha256':record['base_graph_certificate']['precision_search_certificate_sha256'],
                'node_address_bits':record['node_address_bits'],'allocation_limits':record['allocation_limits'],'strategies':[],'selected_strategy':None}
            for strategy in record['strategies']:
                original=dict(strategy);original.setdefault('implementation','original');original['prior_information']='previous_proved_prior'
                row['strategies'].append(original)
            if record['kind']=='positive_r1' and record['complete_deterministic_budget'] is not None:
                case=key(record);prior=information.first_probability_interval(models[record['model_id']],case[2],coefficients[case],financial[case],
                    record['base_graph_certificate'],record['correction_prior'],moment[record['model_id']])
                row['coefficient_informed_prior']=prior
                for strategy in record['strategies']:
                    new={k:v for k,v in strategy.items() if k!='result'};new.setdefault('implementation','original');new['prior_information']='coefficient_informed_interval'
                    new['result']=information.informed_complete_plan(strategy['strategy'],record['complete_deterministic_budget'],record['base_graph_certificate'],prior,profiles(record,strategy))
                    row['strategies'].append(new)
            eligible=[s for s in row['strategies'] if s['result']['kind']!='rejected']
            row['selected_strategy']=min(eligible,key=lambda s:(s['result']['T_count'],s['implementation'],s['architecture'],s['strategy'],s['prior_information'])) if eligible else None
            records.append(row)
    selected=[]
    for family in families.values():
        for price in family['selected_prices']:
            case=key(price);row={'model_id':case[0],'payoff':case[1],'cap':case[2],'methods':{}}
            for kind in KINDS:
                eligible=[(ix,r) for ix,r in enumerate(records) if key(r)==case and r['kind']==kind and r['selected_strategy']]
                winner=min(eligible,key=lambda pair:(pair[1]['selected_strategy']['result']['T_count'],pair[0])) if eligible else None
                row['methods'][kind]={'refined_candidate_index':winner[0],'family_archive':winner[1]['family_archive'],'candidate_index':winner[1]['candidate_index'],'selected':winner[1]['selected_strategy']} if winner else None
            available=[k for k in KINDS if row['methods'][k]]
            row['quantum_resource_winner']=min(available,key=lambda k:row['methods'][k]['selected']['result']['T_count']) if available else None
            screen=information.shared_continuous_screen(models[case[0]],case[2],coefficients[case],moment[case[0]])
            row['shared_continuous_screen']=screen
            if screen['stops_without_quantum']:
                for method in row['methods'].values():
                    if method is None:continue
                    method['conditional_estimator_selection']=method['selected']
                    method['selected']={'implementation':'shared_continuous_screen','architecture':None,'strategy':'analytical_price_interval','prior_information':'shared_low_order_control',
                        'result':{'kind':'deterministic_correction_prior','T_count':0,'quantum_shots':0,'total_price_error_upper':screen['price_radius_upper'],
                            'certified_price_interval':screen['price_interval'],'price_midpoint':screen['price_midpoint'],'reason':'Shared continuous curvature price certificate; no quantum gain.'}}
            row['strong_classical_policy']='All36 financial targets are certified by the retained capped and ordinary-call Taylor/validated quadrature methods; zero quantum estimation calls.'
            selected.append(row)
    assert len(selected)==36 and len({key(r) for r in selected})==36
    assert base_hashes=={name:sha(ROOT/'results'/name) for name in names}
    inputs=names+['generalisation_protocol.json','higherdim_preparation.json','boundary_financial_grid_transfer_summary.json','tensor_first_transfer.json','boundary_derivative_norms.json','gaussian_derivative_norms.json','qft_rotation_tapes.json','qft_rotation_certificates.json']
    dependencies=['correction_information_priors.py','boundary_precision_selection.py','boundary_grid_certificate.py','boundary_derivative_certificate.py','integral_pricing_budget.py','qpe_estimator.py','qft_synthesis.py','taylor_pricing.py','continuous_quadrature.py','gaussian_derivative_norm.py','certified_scalar.py','quantum_gates.py','value_oracle.py','architecture.py','reversible.py']
    result={'schema':'qbp-complete-information-pricing-comparison-v1','status':'complete36-case finite-pool comparison','records':records,'selected_prices':selected,
        'continuous_ratio_moment_certificates':list(moment.values()),'source_sha256':sha(__file__),'dependency_source_sha256':{n:sha(ROOT/n) for n in dependencies},
        'input_sha256':{n:sha(ROOT/'results'/n) for n in inputs},'monetary_radius':'1/100','alpha_per_selected_price':'1/3600','union_over_selected36_failure_upper':'1/100',
        'elapsed_seconds_diagnostic':time.perf_counter()-start,'scope':'Stronger coefficient-informed R1 probability intervals and analytical stops, with every previous strategy retained. Same finite financial law, base arithmetic certificates, circuits, precision maps and full costs. No financial price answer enters the prior refinement; classical pricing results inform only the separate practical policy.'}
    target=ROOT/'results/complete_information_comparison.json.gz'
    with target.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as z:z.write(json.dumps(result,default=str,separators=(',',':')).encode())
    summary={k:v for k,v in result.items() if k not in ('records','continuous_ratio_moment_certificates')};summary['archive_sha256']=sha(target);summary['candidate_programs']=len(records)
    (ROOT/'results/complete_information_summary.json').write_text(json.dumps(summary,default=str,indent=2)+'\n')
    for row in selected:print(row['model_id'],row['payoff'],row['cap'],{k:v['selected']['result']['T_count'] if v else None for k,v in row['methods'].items()},flush=True)
    print('Refined',len(records),'candidates',len(selected),'prices',time.perf_counter()-start,'seconds',flush=True)
if __name__=='__main__':main()
