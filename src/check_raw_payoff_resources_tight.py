"""Retained-DAG and full-preparation checks, separate from the new producer.

The new graph builder, certificate and confidence frontend are not imported.
Existing independently audited primitive emitters and own-formula QPE auditor
are reused. This is an interface/resource check, not independent novelty proof.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
import argparse,gzip,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def load_module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def reconstruct(row,protocol,o):
    g=o.IntegralGraph(protocol['q'],protocol['outer_scalar_table_bits'])
    for j in range(row['model']['outer_dimension']):g.inp('iy'+str(j),protocol['outer_quantile_bits'])
    g.inp('inner_normal_label',row['normal_bits'])
    for r in row['graph_nodes']:
        meta=r['meta'];kind=r['kind']
        if kind=='function':
            meta=o.ScalarSpec(meta['name'],F(meta['lo']),F(meta['span']),tuple(sorted((k,F(v)) for k,v in meta['parameters'].items())),tuple(map(F,meta['output_bounds'])))
        elif kind in ('constant','scale','clamp','clamp_signed','quantile'):meta=F(meta)
        elif kind in ('uniform','minimum'):meta=int(meta)
        elif meta=='None':meta=None
        n=g.op(kind,*r['args'],meta=meta,bits=r['bits'],output_bits=0,bounds=r['ideal_bounds']);assert n==r['name']
        if kind=='constant':g.constant_names[meta]=n
    g.configuration=dict(row['graph_configuration']);g.semantic_nodes=dict(row['graph_semantic_nodes'])
    g.finish(g.semantic_nodes['probability'])
    assert g.configuration==row['graph_configuration']
    assert all(g.nodes[r['name']].output_bits==r['output_bits'] for r in row['graph_nodes'])
    return g
def payoff_shape(g):
    s=g.semantic_nodes;assert g.nodes[s['upper_endpoint']].kind=='add'
    assert g.nodes[s['upper_endpoint']].args==(s['lower_endpoint'],s['logT'])
    for side in ('PG','PA'):
        endpoint=s['lower_endpoint' if side=='PG' else 'upper_endpoint']
        p=g.nodes[s[side+'_positive_log']];assert p.kind=='positive_sub' and p.args==(endpoint,s['common_offset'])
        b=g.nodes[s[side+'_bounded_log']];assert b.kind=='clamp' and b.meta==F(1,64) and b.args==(s[side+'_positive_log'],)
        e=g.nodes[s[side+'_exp']];assert e.kind=='function' and e.meta.name=='exp' and e.meta.lo==0 and e.meta.span==F(1,64)
        argument=g.nodes[e.args[0]];assert argument.kind=='clamp' and argument.meta==F(1,64) and argument.args==(s[side+'_bounded_log'],)
        cap=g.nodes[s[side]];assert cap.kind=='clamp' and cap.meta==1
        scale=g.nodes[cap.args[0]];assert scale.kind=='scale' and scale.meta==100
        sub=g.nodes[scale.args[0]];assert sub.kind=='positive_sub' and sub.args[0]==s[side+'_exp']
        assert g.nodes[sub.args[1]].kind=='constant' and g.nodes[sub.args[1]].meta==1
    assert g.nodes[s['PG_exp']].meta==g.nodes[s['PA_exp']].meta
    assert g.nodes[s['PG_exp']].meta.output_bounds==(F(1),F(64,63))
    output=g.nodes[g.output];assert output.kind=='clamp' and output.meta==1
    difference=g.nodes[output.args[0]];assert difference.kind=='positive_sub' and difference.args==(s['PA'],s['PG'])
def main():
    a=argparse.ArgumentParser();a.add_argument('--repo',required=True);a.add_argument('--results',required=True);a.add_argument('--output',required=True);args=a.parse_args()
    repo=Path(args.repo).resolve();results=Path(args.results).resolve();out=Path(args.output).resolve()
    if out.exists():raise FileExistsError('Preserve audit')
    sys.path.insert(0,str(repo));import boundary_value_oracle as o;import boundary_value_certificate as c
    import quantum_gates as qg
    confidence=load_module('retained_independent_charge',repo/'checks/audit_boundary_precision.py')
    data=json.load(gzip.open(results/'raw_payoff_resources.json.gz','rt'));protocol=json.loads((results/'raw_payoff_protocol.json').read_text())
    assert sha(results/'raw_payoff_protocol.json')==data['protocol_sha256']
    for n,h in data['dependencies_sha256'].items():assert sha(repo/n)==h
    for n,h in data['inputs_sha256'].items():assert sha(repo/'results'/n)==h
    assert sha(Path(__file__).parent/'raw_payoff_oracle_tight.py')==data['source_sha256']==protocol['source_sha256']
    factor=json.loads((repo/'results/factor_basis_quantum_resources.json').read_text())
    threshold=json.loads((repo/'results/shared_threshold_resources.json').read_text())
    numerical=json.loads((repo/'results/strengthened_r1_resources.json').read_text())
    for filename,key in [('factor_basis_quantum_resources.json','factor_resources_sha256'),('shared_threshold_resources.json','threshold_resources_sha256'),('strengthened_r1_resources.json','numerical_resources_sha256')]:assert sha(repo/'results'/filename)==data[key]==protocol[key]
    archive=results/data['scalar_archive_file'];assert sha(archive)==data['scalar_archive_sha256']
    with gzip.open(archive,'rt') as stream:ords=json.load(stream)
    archive_index={(r['certificate_sha256'],r['graph_word_width']):r for r in ords}
    assert len(archive_index)==len(ords)==data['unique_scalar_certificates_and_wordwidths']
    for r in ords:assert digest(r['certificate'])==r['certificate_sha256'] and digest(r['packed_table'])==r['packed_table_sha256']
    tape_data=json.loads((repo/'results/qft_rotation_tapes.json').read_text());qftcert=json.loads((repo/'results/qft_rotation_certificates.json').read_text())
    bounds={r['key']:F(r['operator_norm_upper_rational']) for r in qftcert['certificates']}
    norm=json.loads((repo/'results/gaussian_derivative_norms.json').read_text())
    nu0=F(next(r['norm_upper'] for r in norm['unit_norms'] if r['unit_derivative_order']==0 and r['passed']))
    expected={(m,bn,be) for m,_ in protocol['pairs'] for bn,be in product(protocol['inner_normal_bits'],protocol['narrow_exp_address_bits'])}
    assert {(r['original_model_id'],r['normal_bits'],r['payoff_exp_bits']) for r in data['records']}==expected
    assert len(data['records'])==48
    cache={};blockcache={};profiles=plans=integer_histories=scalar_certificates=0;bindings=[];start=time.perf_counter()
    winners={(r['original_model_id'],r['selected_normal_bits'],r['selected_payoff_exp_bits']) for r in data['selected_comparisons']}
    for row in data['records']:
        assert not row.get('rejection')
        assert [row['original_model_id'],row['candidate_id']] in protocol['pairs']
        prep=next(r for r in factor['preparations'] if r['original_model_id']==row['original_model_id'] and r['candidate_id']==row['candidate_id'])
        old=next(r for r in factor['records'] if r['original_model_id']==row['original_model_id'] and r['factor_basis_candidate_id']==row['candidate_id'] and r['kind']=='positive_r1')
        assert prep['model']==row['model'];g=reconstruct(row,protocol,o);payoff_shape(g)
        tables,certificates=o.make_tables(g,cache,True)
        for name,cc in certificates.items():
            archived=archive_index[(digest(cc),g.width)];assert tuple(archived['packed_table'])==tables[o.table_key(name,g.nodes[name])];scalar_certificates+=1
        body=c.certify_graph(g,tables,certificates);record=row['probability_certificate'];assert body==record['body_arithmetic_certificate'] and body['status']=='passed'
        normal=certificates[g.semantic_nodes['normal_plus_L']];a_lo=F(normal['mass_interval']['lower_exact_fraction']);eN=F(normal['maximum_node_error_upper'])
        D=(1-a_lo)/2+a_lo/F(2*(1<<row['normal_bits']))+nu0*eN;ear=F(body['probability_error_upper']);eta=min(F(1),ear+2*D)
        assert F(record['inner_normal_CDF_error_upper'])==D and F(record['arithmetic_probability_error_upper'])==ear and F(record['probability_error_upper'])==eta
        assert all(g.nodes[n].meta.span==F(1,64) for n in (g.semantic_nodes['PA_exp'],g.semantic_nodes['PG_exp']))
        narrow=certificates[g.semantic_nodes['PG_exp']];nodes=[r['integer'] for r in narrow['nodes']]
        assert all(a<=b for a,b in zip(nodes,nodes[1:]))
        # Strict integer order and range replay at fixed predeclared label patterns.
        for j in range(12):
            inputs={name:((j*j+17*j+3)*(k+1)*97)% (1<<bits) for k,(name,bits) in enumerate(g.inputs.items())}
            values=g.evaluate(inputs,tables,True);s=g.semantic_nodes
            assert values[s['upper_endpoint']]==values[s['lower_endpoint']]+values[s['logT']]
            assert 0<=values[s['PG']]<=values[s['PA']]<=1<<g.q
            assert values[g.output]==values[s['PA']]-values[s['PG']];integer_histories+=1
        bud=row['complete_deterministic_budget'];previous=old['complete_deterministic_budget']
        for k in ('external_factor_upper','external_factor_radius','coefficient_radius_upper','financial_transfer_upper'):assert bud[k]==previous[k]
        mass,mr,ec,beta=(F(previous[k]) for k in ('external_factor_upper','external_factor_radius','coefficient_radius_upper','financial_transfer_upper'))
        det=ec+beta+mass*eta+mr;assert det==F(bud['total_deterministic_undiscounted_price_error_upper'])<F(1,100)
        prior=row['coefficient_informed_prior'];op=old['coefficient_informed_prior'];lo,hi=F(op['finite_ideal_mean_lower']),F(op['finite_ideal_mean_upper'])
        assert F(prior['probability_lower'])==max(F(0),lo-eta) and F(prior['probability_upper'])==min(F(1),hi+eta)
        assert prior['previous_mean_information_sha256']==digest(op)
        for mode in protocol['architectures']:
            p=row['profiles'][mode]['centred'];fresh=o.counted_composition(g,tables,blockcache,mode);assert p==fresh;profiles+=1
            strategy=next(s for s in row['strategies'] if s['architecture']==mode)['result'];assert strategy['kind']=='complete_quantum_resource_plan'
            assert len(strategy['flags'])==1;flag=strategy['flags'][0];plan=flag['plan'];qft=flag['fourier_profile']
            target=(F(1,100)-det)/mass;upper=F(prior['probability_upper'])
            confidence.check_plan(plan,target,F(1,3600),upper,qft,tape_data['rotations'],bounds)
            confidence.charge_check(flag['complete_charge'],p,qft,plan)
            assert strategy['T_count']==flag['complete_charge']['metrics']['T']
            stat=mass*(F(plan['amplitude_error_upper'])+F(plan['postprocess_allowance']))
            assert stat==F(strategy['statistical_price_error_upper']) and det+stat==F(strategy['total_price_error_upper'])<=F(1,100)
            assert F(strategy['failure_probability_upper'])<=F(1,3600);plans+=1
            if (row['original_model_id'],row['normal_bits'],row['payoff_exp_bits']) in winners:
                # Emit the full clean chronological circuit without storing giant tapes.
                circuit,registers,comparison,flagwire=o.compile_flag(g,tables,mode,store_gates=False)
                metrics=circuit.resources();counts=Counter(metrics['clifford_T_expanded_gate_counts'])
                counts['H']+=len(circuit.uniform_preparation_gates)
                assert counts==Counter(p['primitive_counts'])
                core=tuple(w for r in registers.values() for w in r)+tuple(comparison)+(flagwire,)
                assert core==circuit.reflection_core and len(core)==p['reflection_core_qubits']
                highzero=flagwire-1;assert highzero not in core and highzero not in comparison
                assert metrics['wires']<=p['logical_qubits_upper']
                work=tuple(w for w in range(metrics['wires']) if w not in core and w!=highzero)[:len(core)-1]
                assert len(work)==len(core)-1
                reflection=qg.expand_toffoli(qg.grover([],core,flagwire,work,control=metrics['wires']))
                assert Counter(r[0] for r in reflection)==Counter(flag['complete_charge']['controlled_reflection']['gate_counts'])
                bindings.append({'original_model_id':row['original_model_id'],'architecture':mode,'actual_flag_wire':flagwire,
                  'excluded_endpoint_high_zero_wire':highzero,'reflection_core':list(core),'reflection_work':list(work),
                  'full_emitted_peak_wires':metrics['wires'],'serial_profile_peak_upper':p['logical_qubits_upper'],
                  'full_preparation_primitive_counts_verified':True,'controlled_reflection_phase_counts_verified':True})
        assert row['selected_strategy']==min(row['strategies'],key=lambda s:(s['result']['T_count'],s['architecture']))
    for comparison in data['selected_comparisons']:
        mid=comparison['original_model_id'];cohort=[r for r in data['records'] if r['original_model_id']==mid]
        winner=min(cohort,key=lambda r:(r['selected_strategy']['result']['T_count'],r['normal_bits'],r['payoff_exp_bits']))
        assert winner['selected_strategy']==comparison['selected_raw_payoff_strategy']
        tr=next(r for r in threshold['rows'] if r['original_model_id']==mid);nr=next(r for r in numerical['selected_comparisons'] if r['original_model_id']==mid)
        assert comparison['fixed_threshold_selected_strategy']==tr['selected_strategy'] and comparison['strengthened_numeric_selected_strategy']==nr['selected_numerical_strategy']
        raw=winner['selected_strategy']['result']['T_count'];thresholdT=tr['selected_strategy']['result']['T_count'];numericT=nr['selected_numerical_strategy']['result']['T_count']
        assert F(comparison['threshold_relative_T_saving_against_raw_exact'])==F(raw-thresholdT,raw)
        assert F(comparison['raw_relative_T_saving_against_numeric_exact'])==F(numericT-raw,numericT)
    out.write_text(json.dumps({'schema':'qbp-raw-payoff-interface-resource-check-tight-v2','status':'passed',
       'source_sha256':sha(__file__),'protocol_sha256':sha(results/'raw_payoff_protocol.json'),'resource_sha256':sha(results/'raw_payoff_resources.json.gz'),
       'previous_own_formula_charge_auditor_sha256':sha(repo/'checks/audit_boundary_precision.py'),
       'retained_DAGs':48,'regenerated_scalar_certificates':scalar_certificates,'whole_profiles':profiles,'independent_QPE_charge_plans':plans,
       'exact_integer_history_checks':integer_histories,'full_clean_preparation_bindings':bindings,
       'replay_seconds_diagnostic':time.perf_counter()-start,
       'scope':'Separate retained-DAG, primitive-emitter, complete clean preparation and independent own-formula QPE resource replay. Common-normal transfer derivation awaits independent reviewer; no originality or physical advantage claim.'},indent=2)+'\n')
    print('passed',profiles,'profiles',plans,'plans',len(bindings),'full preparations',flush=True)
if __name__=='__main__':main()
