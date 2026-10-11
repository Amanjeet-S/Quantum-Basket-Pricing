"""A registered independent PA/PG precision pool for the numerical R1 rival.

Only the two conditional-price address widths vary. Every other financial
map, finite law, architecture, confidence and monetary target is unchanged.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime,timezone
import argparse,copy,gzip,hashlib,json,sys,time
sys.dont_write_bytecode=True

ROOT=None;FACTOR=None;THRESHOLD=None;RESULTS=None
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def write(p,v):
    if p.exists():raise FileExistsError('Preserve frozen evidence')
    p.write_text(json.dumps(v,indent=2,default=str)+'\n')
def configure(repo,factor,threshold,output):
    global ROOT,FACTOR,THRESHOLD,RESULTS,o,c,informed
    ROOT=Path(repo).resolve();FACTOR=Path(factor).resolve();THRESHOLD=Path(threshold).resolve();RESULTS=Path(output).resolve()
    RESULTS.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(ROOT))
    import boundary_value_oracle as o
    import boundary_value_certificate as c
    import correction_information_priors as informed
DEPENDENCIES=('boundary_value_oracle.py','boundary_value_certificate.py','integral_value_oracle.py',
 'integral_value_certificate.py','boundary_derivative_certificate.py','value_oracle.py',
 'reversible.py','architecture.py','quantum_gates.py','boundary_precision_selection.py',
 'correction_information_priors.py','integral_pricing_budget.py','qpe_estimator.py',
 'qft_synthesis.py','certified_scalar.py','taylor_pricing.py','continuous_quadrature.py')
INPUTS=('gaussian_derivative_norms.json','boundary_derivative_norms.json','qft_rotation_tapes.json','qft_rotation_certificates.json')
def signatures():
    return {'source_sha256':sha(__file__),'factor_resources_sha256':sha(FACTOR),
      'threshold_resources_sha256':sha(THRESHOLD),
      'dependencies_sha256':{n:sha(ROOT/n) for n in DEPENDENCIES},
      'inputs_sha256':{n:sha(ROOT/'results'/n) for n in INPUTS}}
def register():
    threshold=json.loads(THRESHOLD.read_text())
    write(RESULTS/'strengthened_r1_protocol.json',{'schema':'qbp-strengthened-r1-protocol-v1',
      'registered_before_new_pool_execution':True,'registration_utc':datetime.now(timezone.utc).isoformat(),
      'research_record_date':'2026-10-10',**signatures(),
      'pairs':[[r['original_model_id'],r['candidate_id']] for r in threshold['rows']],
      'independent_PA_address_bits':[6,8,10,12],'independent_PG_address_bits':[6,8,10,12],
      'configurations_per_pair':16,'configurations_total':64,
      'unchanged_outer_policy':'Quantile labels10, q28, L6, outer exp and log addresses12, exact same graph headroom and data-independent finite normaliser.',
      'architectures':['exact7t','unitary8t'],'no_repairs':True,
      'confidence':'Same canonical QPE, alpha=1/3600 per price, monetary radius1/100.',
      'probability_information':'Identical price-independent finite R1 mean interval, enlarged using each fresh numerical graph probability-error certificate.',
      'selection':'Minimum complete T cap over every admitted configuration and both architectures; all budget/range rejects retained.',
      'known_predecessors':'All12 numerical and fixed threshold outcomes are already known. This is a new stronger-rival pool registered before its new lower-precision outcomes; the threshold source and policy are not tuned.',
      'financial_price_answers_excluded':True,'data_origin':'Synthetic exact rational models and scalar tables only.'})

def run():
    protocol=json.loads((RESULTS/'strengthened_r1_protocol.json').read_text())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Frozen inputs changed')
    factor=json.loads(FACTOR.read_text());threshold=json.loads(THRESHOLD.read_text())
    rows=[];selected=[];table_cache={};block_cache={};archive=[];seen=set();start=time.perf_counter()
    for original,candidate in protocol['pairs']:
        prep=next(p for p in factor['preparations'] if p['original_model_id']==original and p['candidate_id']==candidate)
        old=next(r for r in factor['records'] if r['original_model_id']==original and r['factor_basis_candidate_id']==candidate and r['kind']=='positive_r1')
        prior_old=old['coefficient_informed_prior'];meanlo=F(prior_old['finite_ideal_mean_lower']);meanhi=F(prior_old['finite_ideal_mean_upper'])
        for bpa,bpg in product(protocol['independent_PA_address_bits'],protocol['independent_PG_address_bits']):
            g=o.build_graph(prep['model'],F(1),q=28,b_y=10,table_bits=12,outer_bound=6,
              U_star=F(prep['domain']['U_star']),retained_bound=max(abs(F(v)) for v in prep['domain']['coordinate_bounds']),kind='positive_r1')
            PA,PG=g.semantic_nodes['PA'],g.semantic_nodes['PG'];g.nodes[PA].bits=bpa;g.nodes[PG].bits=bpg
            g.configuration['independent_price_table_address_bits']={'PA':bpa,'PG':bpg}
            tables,ords=o.make_tables(g,table_cache,True);cert=c.certify_graph(g,tables,ords)
            row={'original_model_id':original,'candidate_id':candidate,'PA_bits':bpa,'PG_bits':bpg,
              'graph_configuration':g.configuration,'graph_semantic_nodes':g.semantic_nodes,
              'base_graph_certificate':cert,'profiles':{},'strategies':[]}
            for name,cc in ords.items():
                key=(digest(cc),g.width)
                if key not in seen:
                    seen.add(key);table=tables[o.table_key(name,g.nodes[name])]
                    archive.append({'certificate_sha256':key[0],'certificate':cc,'packed_table':table,
                       'packed_table_sha256':digest(table),'graph_word_width':g.width})
            if cert['status']!='passed':row['rejection']={'reason':'Arithmetic guard failure'}
            else:
                eta=F(cert['probability_error_upper']);B=F(g.configuration['monetary_signal_normaliser'])
                bud=copy.deepcopy(old['complete_deterministic_budget'])
                mass=F(bud['external_factor_upper']);mr=F(bud['external_factor_radius']);ec=F(bud['coefficient_radius_upper']);beta=F(bud['financial_transfer_upper'])
                det=ec+beta+mass*B*eta+mr*B
                bud.update(arithmetic_before_mass_upper=str(B*eta),
                    total_deterministic_undiscounted_price_error_upper=str(det))
                prior={'probability_lower':str(max(F(0),meanlo/B-eta)),
                       'probability_upper':str(min(F(1),meanhi/B+eta)),
                       'native_positive':True,'finite_ideal_mean_lower':str(meanlo),'finite_ideal_mean_upper':str(meanhi),
                       'previous_mean_information_sha256':digest(prior_old),
                       'information_scope':'Same price-independent finite-mean interval; only fresh scalar arithmetic error changes.'}
                row.update(complete_deterministic_budget=bud,coefficient_informed_prior=prior)
                if det>=F(1,100):row['rejection']={'reason':'Deterministic monetary tolerance exhausted','deterministic_error_upper':str(det)}
                else:
                    for mode in protocol['architectures']:
                        profile=o.counted_composition(g,tables,block_cache,mode);row['profiles'][mode]={'centred':profile}
                        result=informed.informed_complete_plan('positive',bud,cert,prior,row['profiles'][mode])
                        row['strategies'].append({'architecture':mode,'result':result})
                    row['selected_strategy']=min(row['strategies'],key=lambda s:(s['result']['T_count'],s['architecture']))
            if bpa==bpg==12:
                if not row.get('selected_strategy'):raise AssertionError('All12 historical admitted baseline unexpectedly failed')
                assert row['selected_strategy']['result']['T_count']==old['selected_strategy']['result']['T_count']
                assert F(row['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper'])==F(old['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper'])
                row['all12_historical_T_and_deterministic_error_exactly_reproduced']=True
            rows.append(row)
            print(original,candidate,bpa,bpg,'T',row.get('selected_strategy',{}).get('result',{}).get('T_count'),flush=True)
        cohort=[r for r in rows if r['original_model_id']==original and r['candidate_id']==candidate and r.get('selected_strategy')]
        winner=min(cohort,key=lambda r:(r['selected_strategy']['result']['T_count'],r['PA_bits'],r['PG_bits']))
        new=next(r for r in threshold['rows'] if r['original_model_id']==original and r['candidate_id']==candidate)
        numerical=winner['selected_strategy']['result']['T_count'];binary=new['selected_strategy']['result']['T_count']
        selected.append({'original_model_id':original,'candidate_id':candidate,
          'selected_PA_bits':winner['PA_bits'],'selected_PG_bits':winner['PG_bits'],
          'selected_numerical_strategy':winner['selected_strategy'],
          'fixed_threshold_selected_strategy':new['selected_strategy'],
          'threshold_relative_T_saving_exact':str(F(numerical-binary,numerical)),
          'threshold_relative_T_saving_per_cent_diagnostic':100*(numerical-binary)/numerical,
          'admitted_configurations':len(cohort),'rejected_configurations':16-len(cohort)})
    archive_path=RESULTS/'strengthened_r1_scalar_ordinates.json.gz'
    if archive_path.exists():raise FileExistsError('Preserve archive')
    with archive_path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as stream:stream.write(json.dumps(archive,separators=(',',':'),default=str).encode())
    if any(protocol[k]!=v for k,v in signatures().items()):raise RuntimeError('Frozen sources changed during run')
    write(RESULTS/'strengthened_r1_resources.json',{'schema':'qbp-strengthened-r1-resources-v1',
      'status':'Completed retaining every rejection','protocol_sha256':sha(RESULTS/'strengthened_r1_protocol.json'),
      **signatures(),'records':rows,'selected_comparisons':selected,
      'scalar_archive_file':archive_path.name,'scalar_archive_sha256':sha(archive_path),
      'unique_scalar_certificates_and_wordwidths':len(archive),'emitted_block_profiles':len(block_cache),
      'generation_seconds_diagnostic':time.perf_counter()-start,
      'scope':'Registered finite independent PA/PG precision pool, matched complete monetary accuracy and confidence. Not global precision optimisation or physical runtime.'})
    print('COMPARISONS',json.dumps([{k:r[k] for k in ('original_model_id','selected_PA_bits','selected_PG_bits','threshold_relative_T_saving_per_cent_diagnostic')} for r in selected]),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run'));p.add_argument('--repo',required=True)
    p.add_argument('--factor-resources',required=True);p.add_argument('--threshold-resources',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();configure(a.repo,a.factor_resources,a.threshold_resources,a.output)
    register() if a.mode=='register' else run()
