"""Independent full replay of the frozen PA/PG precision rival pool.

The new rival producer is not imported. The retained numerical financial
builder supplies the unchanged graph, with its two registered price-address
widths assigned independently. Every new table and numeric certificate is
reconstructed. Exact priors, monetary targets, clean comparator profiles,
QPE confidence and complete charges use separate own-formula calculations.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import importlib.util
from itertools import product
import json
import sys
sys.dont_write_bytecode=True


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def lo(iv):return F(iv['lower_exact_fraction'])
def hi(iv):return F(iv['upper_exact_fraction'])

def import_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=Path,required=True);ap.add_argument('--repo',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    sys.path.insert(0,str(args.repo))
    import boundary_value_oracle as oracle
    import boundary_value_certificate as certifier
    qdir=args.step/'adversarial/quantum';rp=qdir/'strengthened_baseline_results';path=rp/'strengthened_r1_resources.json'
    producer=qdir/'strengthened_r1_baseline.py'
    data=json.loads(path.read_text());protocol=json.loads((rp/'strengthened_r1_protocol.json').read_text())
    factorpath=args.step/'resources/results/factor_basis_quantum_resources.json';factor=json.loads(factorpath.read_text())
    thresholdpath=qdir/'threshold_results/shared_threshold_resources.json';threshold=json.loads(thresholdpath.read_text())
    assert data['protocol_sha256']==sha(rp/'strengthened_r1_protocol.json')
    assert data['source_sha256']==protocol['source_sha256']==sha(producer)
    assert data['factor_resources_sha256']==protocol['factor_resources_sha256']==sha(factorpath)
    assert data['threshold_resources_sha256']==protocol['threshold_resources_sha256']==sha(thresholdpath)
    assert data['dependencies_sha256']==protocol['dependencies_sha256']
    assert data['inputs_sha256']==protocol['inputs_sha256']
    for name,value in protocol['dependencies_sha256'].items():assert sha(args.repo/name)==value,name
    for name,value in protocol['inputs_sha256'].items():assert sha(args.repo/'results'/name)==value,name
    assert protocol['registered_before_new_pool_execution'] is True and protocol['financial_price_answers_excluded'] is True
    pairs=[[r['original_model_id'],r['candidate_id']] for r in threshold['rows']]
    assert protocol['pairs']==pairs and len(pairs)==4
    assert protocol['independent_PA_address_bits']==protocol['independent_PG_address_bits']==[6,8,10,12]
    assert protocol['configurations_total']==64 and protocol['configurations_per_pair']==16 and protocol['no_repairs'] is True
    expectedkeys=[(original,candidate,bpa,bpg) for original,candidate in pairs for bpa,bpg in product((6,8,10,12),repeat=2)]
    assert [(r['original_model_id'],r['candidate_id'],r['PA_bits'],r['PG_bits']) for r in data['records']]==expectedkeys
    archivepath=rp/data['scalar_archive_file'];assert sha(archivepath)==data['scalar_archive_sha256']
    with gzip.open(archivepath,'rt') as f:archive=json.load(f)
    archived={(item['certificate_sha256'],item['graph_word_width']):item for item in archive}
    assert len(archived)==len(archive)==data['unique_scalar_certificates_and_wordwidths']
    for (h,width),item in archived.items():
        assert digest(item['certificate'])==h and digest(item['packed_table'])==item['packed_table_sha256']
    previous_argv=sys.argv[:]
    sys.argv=[str(args.step/'resources/audit_factor_basis_quantum_resources.py'),'--repo',str(args.repo),'--results',str(args.step/'resources/results'),'--producer',str(args.step/'resources/factor_basis_quantum_resources.py'),'--certificate',str(args.step/'certificate/results/factor_basis_certificates.json')]
    own=import_module('previous_r1_own_formula_audit',args.step/'resources/audit_factor_basis_quantum_resources.py')
    sys.argv=previous_argv
    tapes=json.loads((args.repo/'results/qft_rotation_tapes.json').read_text())
    qftcert=json.loads((args.repo/'results/qft_rotation_certificates.json').read_text())
    assert sha(args.repo/'results/qft_rotation_tapes.json')==qftcert['tapes_sha256']
    qftbounds={r['key']:F(r['operator_norm_upper_rational']) for r in qftcert['certificates']}
    tablecache={};blockcache={};profiles=plans=scalarchecks=budgetchecks=priorchecks=0;rejections=[];cohorts={tuple(pair):[] for pair in pairs};all12=0
    for row in data['records']:
        original,candidate,bpa,bpg=(row[k] for k in ('original_model_id','candidate_id','PA_bits','PG_bits'))
        prep=next(p for p in factor['preparations'] if p['original_model_id']==original and p['candidate_id']==candidate)
        old=next(r for r in factor['records'] if r['original_model_id']==original and r['factor_basis_candidate_id']==candidate and r['kind']=='positive_r1')
        domain=prep['domain']
        g=oracle.build_graph(prep['model'],F(1),q=28,b_y=10,table_bits=12,outer_bound=6,
            U_star=F(domain['U_star']),retained_bound=max(abs(F(v)) for v in domain['coordinate_bounds']),kind='positive_r1')
        PA,PG=g.semantic_nodes['PA'],g.semantic_nodes['PG'];g.nodes[PA].bits=bpa;g.nodes[PG].bits=bpg
        g.configuration['independent_price_table_address_bits']={'PA':bpa,'PG':bpg}
        assert g.configuration==row['graph_configuration'] and g.semantic_nodes==row['graph_semantic_nodes']
        assert all(node.bits==12 for name,node in g.nodes.items() if node.kind=='function' and name not in (PA,PG))
        tables,ordinates=oracle.make_tables(g,tablecache,True)
        cert=certifier.certify_graph(g,tables,ordinates)
        assert cert==row['base_graph_certificate']
        for name,cc in ordinates.items():
            key=(digest(cc),g.width);assert key in archived and archived[key]['certificate']==cc
            table=tables[oracle.table_key(name,g.nodes[name])]
            assert digest(table)==archived[key]['packed_table_sha256'] and list(table)==archived[key]['packed_table']
            scalarchecks+=len(cc['nodes'])
        if cert['status']!='passed':
            assert row['rejection']['reason']=='Arithmetic guard failure' and not row['profiles'] and not row['strategies']
            rejections.append({'original_model_id':original,'candidate_id':candidate,'PA_bits':bpa,'PG_bits':bpg,'reason':row['rejection']['reason']})
            continue
        eta=F(cert['probability_error_upper']);B=F(g.configuration['monetary_signal_normaliser'])
        assert B==1
        financial=prep['R1_selected_financial_transfer'];mass=financial['mass_interval'];mh=hi(mass);mr=(hi(mass)-lo(mass))/2
        coefficient=prep['coefficient_certificate']['c0'];cr=(hi(coefficient)-lo(coefficient))/2
        beta=F(financial['financial_transfer_upper']);det=cr+beta+mh*B*eta+mr*B
        budget=row['complete_deterministic_budget']
        for key,value in {'external_factor_upper':mh,'external_factor_radius':mr,'coefficient_radius_upper':cr,
            'financial_transfer_upper':beta,'arithmetic_before_mass_upper':B*eta,'total_deterministic_undiscounted_price_error_upper':det}.items():assert F(budget[key])==value,key
        assert budget['mass_interval']==mass and budget['coefficient_interval']==coefficient
        assert F(budget['auxiliary_components']['bias_upper'])==0
        budgetchecks+=7
        previous=old['coefficient_informed_prior'];meanlo=F(previous['finite_ideal_mean_lower']);meanhi=F(previous['finite_ideal_mean_upper'])
        prior=row['coefficient_informed_prior'];pl=max(F(0),meanlo/B-eta);pu=min(F(1),meanhi/B+eta)
        assert 0<=pl<=pu<=1
        assert F(prior['probability_lower'])==pl and F(prior['probability_upper'])==pu
        assert F(prior['finite_ideal_mean_lower'])==meanlo and F(prior['finite_ideal_mean_upper'])==meanhi
        assert prior['native_positive'] is True and prior['previous_mean_information_sha256']==digest(previous)
        priorchecks+=6
        if det>=F(1,100):
            assert row['rejection']['reason']=='Deterministic monetary tolerance exhausted'
            assert F(row['rejection']['deterministic_error_upper'])==det and not row['profiles'] and not row['strategies']
            rejections.append({'original_model_id':original,'candidate_id':candidate,'PA_bits':bpa,'PG_bits':bpg,'reason':row['rejection']['reason']})
            continue
        assert 'rejection' not in row and len(row['profiles'])==len(row['strategies'])==2
        for architecture,family in row['profiles'].items():
            assert set(family)=={'centred'}
            profile=family['centred']
            fresh=oracle.counted_composition(g,tables,blockcache,architecture)
            assert profile==fresh
            own.verify_profile(profile,set(item['packed_table_sha256'] for item in archive));profiles+=1
        for strategy in row['strategies']:
            out=strategy['result'];profile=row['profiles'][strategy['architecture']]['centred']
            if out['kind']=='deterministic_correction_prior':
                assert F(out['total_price_error_upper'])==det+mh*B*(pu-pl)/2<=F(1,100)
                continue
            assert out['kind']=='complete_quantum_resource_plan' and out['correction_prior']==prior
            assert F(out['deterministic_error_upper'])==det and F(out['monetary_signal_factor_upper'])==mh*B
            assert len(out['flags'])==1
            flag=out['flags'][0]
            assert F(flag['probability_tolerance'])==(F(1,100)-det)/(mh*B)
            assert F(flag['prior_upper'])==pu and F(flag['plan']['alpha'])==F(1,3600)
            T=own.verify_plan(flag,profile,tapes,qftbounds);plans+=1
            assert out['T_count']==T
            stat=mh*B*(F(flag['plan']['amplitude_error_upper'])+F(flag['plan']['postprocess_allowance']))
            assert F(out['statistical_price_error_upper'])==stat
            assert F(out['total_price_error_upper'])==det+stat<=F(1,100)
            assert F(out['failure_probability_upper'])==F(flag['plan']['binomial_failure_upper'])<=F(1,3600)
        selected=min(row['strategies'],key=lambda s:(s['result']['T_count'],s['architecture']))
        assert row['selected_strategy']==selected
        cohorts[(original,candidate)].append(row)
        if bpa==bpg==12:
            assert row['all12_historical_T_and_deterministic_error_exactly_reproduced'] is True
            assert selected['result']['T_count']==old['selected_strategy']['result']['T_count']
            assert det==F(old['complete_deterministic_budget']['total_deterministic_undiscounted_price_error_upper'])
            all12+=1
    assert len(data['selected_comparisons'])==4 and all12==4
    facts=[]
    for selected in data['selected_comparisons']:
        pair=(selected['original_model_id'],selected['candidate_id']);cohort=cohorts[pair]
        best=min(cohort,key=lambda r:(r['selected_strategy']['result']['T_count'],r['PA_bits'],r['PG_bits']))
        assert selected['selected_numerical_strategy']==best['selected_strategy']
        assert selected['selected_PA_bits']==best['PA_bits'] and selected['selected_PG_bits']==best['PG_bits']
        assert selected['admitted_configurations']==len(cohort) and selected['rejected_configurations']==16-len(cohort)
        shared=next(row for row in threshold['rows'] if (row['original_model_id'],row['candidate_id'])==pair)
        assert selected['fixed_threshold_selected_strategy']==shared['selected_strategy']
        numeric=best['selected_strategy']['result']['T_count'];binary=shared['selected_strategy']['result']['T_count']
        saving=F(numeric-binary,numeric)
        assert F(selected['threshold_relative_T_saving_exact'])==saving
        assert abs(selected['threshold_relative_T_saving_per_cent_diagnostic']-float(100*saving))<1e-12
        facts.append({'original_model_id':pair[0],'candidate_id':pair[1],
            'selected_PA_bits':best['PA_bits'],'selected_PG_bits':best['PG_bits'],
            'selected_numerical_T':numeric,'untuned_fixed_threshold_T':binary,
            'threshold_saving_exact':str(saving),'admitted':len(cohort),'rejected':16-len(cohort)})
    result={'schema':'qbp-independent-strengthened-r1-baseline-audit-v1','status':'passed',
        'registered_configurations':64,'retained_rejections':rejections,
        'fresh_scalar_ordinate_checks':scalarchecks,'fresh_complete_profiles':profiles,'complete_confidence_plans':plans,
        'financial_budget_scalar_checks':budgetchecks,'implemented_prior_scalar_checks':priorchecks,
        'all12_exact_historical_replays':all12,'selected_comparisons':facts,
        'resource_sha256':sha(path),'protocol_sha256':sha(rp/'strengthened_r1_protocol.json'),
        'producer_sha256':sha(producer),'audit_source_sha256':sha(__file__),
        'threshold_resource_sha256':sha(thresholdpath),'factor_resource_sha256':sha(factorpath),
        'own_formula_charge_auditor_sha256':sha(args.step/'resources/audit_factor_basis_quantum_resources.py'),
        'scalar_archive_sha256':sha(archivepath),
        'scope':'Independent full reconstruction of every finite PA/PG precision candidate, arithmetic/range reject, price-independent prior, emitted financial profile, confidence target and complete T cap. Exact matched all12 predecessor replay. Finite-pool comparison only; no global optimum, physical runtime or theoretical originality claim.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','registered_configurations','fresh_scalar_ordinate_checks','fresh_complete_profiles','complete_confidence_plans','all12_exact_historical_replays')}))
    print(json.dumps(facts))

if __name__=='__main__':main()
