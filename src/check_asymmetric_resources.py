"""Independent literal circuit and own-formula confidence/charge replay.

The asymmetric producer and confidence frontend are not imported. The prior
independent resource auditor supplies primitive reconstruction and exact QPE
formulae, while scalar and mathematical errors have their separate review.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',required=True)
    ap.add_argument('--results',required=True);ap.add_argument('--output',required=True)
    args=ap.parse_args();repo=Path(args.repo).resolve();results=Path(args.results).resolve();out=Path(args.output).resolve()
    assert not out.exists()
    sys.path.insert(0,str(repo))
    import boundary_value_oracle as e
    import quantum_gates as qg
    h=module('previous_independent_resource_formulae',repo/'check_expanded_precision_resources.py')
    data=h.read(results/'asymmetric_raw_payoff_resources.json.gz')
    protocol=h.read(results/'asymmetric_raw_payoff_protocol.json')
    assert h.sha(results/'asymmetric_raw_payoff_protocol.json')==data['protocol_sha256']
    assert len(data['records'])==len(data['amended_records'])==252
    expected={(m,n,g,a) for m,_ in protocol['pairs'] for n,g,a in product((8,10,12),range(8,13),range(8,13))}
    assert len(expected)==300 and len(data['combined_pool_bindings'])==300
    index,certs=h.archive_index(results/data['scalar_archive_file'])
    assert h.sha(results/data['scalar_archive_file'])==data['scalar_archive_sha256']
    factor=h.read(repo/'results/factor_basis_quantum_resources.json')
    baselines={(r['original_model_id'],r['factor_basis_candidate_id']):r for r in factor['records'] if r['kind']=='positive_r1'}
    tapes=h.read(repo/'results/qft_rotation_tapes.json')['rotations']
    bounds={r['key']:F(r['operator_norm_upper_rational']) for r in h.read(repo/'results/qft_rotation_certificates.json')['certificates']}
    cache={};profiles=plans=zeros=rejected=histories=0;whole=[];start=time.monotonic()
    winners={(r['model_id'],r['selected_normal_bits'],r['selected_PG_exp_bits'],r['selected_PA_exp_bits']) for r in data['selected_comparisons']}
    new={}
    for point,amend in zip(data['records'],data['amended_records']):
        key=(point['original_model_id'],point['normal_bits'],point['PG_exp_bits'],point['PA_exp_bits'])
        assert key in expected and key not in new
        assert (amend['model_id'],amend['normal_bits'],amend['PG_exp_bits'],amend['PA_exp_bits'])==key
        assert amend['pointwise_row_sha256']==h.digest(point)
        assert amend['graph_sha256']==h.digest(point['graph_nodes'])
        assert amend['configuration_sha256']==h.digest(point['graph_configuration'])
        assert amend['profile_sha256']==h.digest(point['profiles'])
        g,tables=h.own_graph(point,e,index);s=g.semantic_nodes
        assert g.nodes[s['PG_exp']].bits==key[2] and g.nodes[s['PA_exp']].bits==key[3]
        assert g.nodes[s['PG_exp']].meta==g.nodes[s['PA_exp']].meta
        output=g.nodes[g.output];assert output.kind=='clamp' and output.meta==1
        sub=g.nodes[output.args[0]];assert sub.kind=='positive_sub' and sub.args==(s['PA'],s['PG'])
        assert point['graph_configuration']['positive_projection_retained'] is True
        assert g.nodes[s['upper_endpoint']].args==(s['lower_endpoint'],s['logT'])
        # Integer replay deliberately checks projection rather than false ordering.
        for j in range(12):
            inp={name:((j*j+17*j+3)*(k+1)*97)%(1<<bits) for k,(name,bits) in enumerate(g.inputs.items())}
            v=g.evaluate(inp,tables,True)
            assert v[s['upper_endpoint']]==v[s['lower_endpoint']]+v[s['logT']]
            assert 0<=v[s['PG']]<=1<<28 and 0<=v[s['PA']]<=1<<28
            assert v[g.output]==min(1<<28,max(0,v[s['PA']]-v[s['PG']]))
            histories+=1
        for mode in ('exact7t','unitary8t'):
            h.profile_check(g,tables,point['profiles'][mode]['centred'],e,cache);profiles+=1
        baseline=baselines[point['original_model_id'],point['candidate_id']]
        for r in (point,amend):
            p,z,b=h.budget_and_plan(r,point['profiles'],baseline,tapes,bounds)
            plans+=p;zeros+=z;rejected+=b
        new[key]=(point,amend)
        if key in winners:
            # The existing helper uses this field for output annotation only.
            annotated=dict(point,payoff_exp_bits=[key[2],key[3]])
            whole.extend(h.whole_preparation(annotated,g,tables,e,qg))

    # Reconstruct exact old references independently of producer selection code.
    oldpoint=h.read(repo/'results/raw_payoff_tight/raw_payoff_resources.json.gz')
    oldamend=h.read(repo/'results/raw_payoff_smoothing/raw_payoff_smoothing_resources.json')
    expanded=h.read(repo/'results/expanded_raw_payoff/expanded_raw_payoff_resources.json.gz')
    old={}
    for p,a in list(zip(oldpoint['records'],oldamend['records']))+list(zip(expanded['records'],expanded['amended_records'])):
        key=(p['original_model_id'],p['normal_bits'],p['payoff_exp_bits'],p['payoff_exp_bits'])
        assert key not in old;old[key]=(p,a)
    assert len(old)==96
    bindings={};reuse=0
    for b in data['combined_pool_bindings']:
        key=(b['model_id'],b['normal_bits'],b['PG_exp_bits'],b['PA_exp_bits'])
        assert key in expected and key not in bindings
        p,a=new[key] if key in new else old[key]
        assert b['record_sha256']==h.digest(a)
        assert b['admitted']==bool(a.get('selected_strategy'))
        if a.get('selected_strategy'):assert b['selected_strategy']==a['selected_strategy']
        if key not in new:reuse+=1
        bindings[key]=a
    assert set(bindings)==expected and reuse==48
    numerical=h.read(repo/'results/strengthened_r1_resources.json')
    threshold=h.read(repo/'results/shared_threshold_resources.json')
    facts=[]
    for c in data['selected_comparisons']:
        mid=c['model_id'];group=[(k,a) for k,a in bindings.items() if k[0]==mid and a.get('selected_strategy')]
        key,a=min(group,key=lambda r:(r[1]['selected_strategy']['result']['T_count'],r[0][1],r[0][2],r[0][3],r[1]['selected_strategy']['architecture']))
        T=a['selected_strategy']['result']['T_count']
        assert key[1:]==(c['selected_normal_bits'],c['selected_PG_exp_bits'],c['selected_PA_exp_bits'])
        assert a['selected_strategy']==c['selected_raw_strategy'] and T==c['raw_T']
        previous=next(r for r in expanded['selected_comparisons'] if r['model_id']==mid)
        num=next(r['selected_numerical_strategy'] for r in numerical['selected_comparisons'] if r['original_model_id']==mid)
        event=next(r['selected_strategy'] for r in threshold['rows'] if r['original_model_id']==mid)
        strongest=min(T,previous['raw_T'],num['result']['T_count']);ET=event['result']['T_count']
        assert c['strongest_registered_rival_T']==strongest and c['shared_threshold_T']==ET
        assert F(c['event_saving_against_strongest_exact'])==F(strongest-ET,strongest)
        assert c['admitted_configurations']==75 and c['rejected_configurations']==0
        assert key[2:]==(9,9)
        facts.append({'model_id':mid,'normal_bits':key[1],'PG_PA_bits':list(key[2:]),'raw_T':T,'strongest_T':strongest,'event_T':ET,'saving_exact':str(F(strongest-ET,strongest))})
    assert profiles==504 and plans==1008 and zeros==rejected==0 and len(whole)==8 and histories==3024
    result={'schema':'qbp-asymmetric-independent-resource-audit-v1','status':'passed',
        'source_sha256':h.sha(__file__),'prior_resource_formula_source_sha256':h.sha(repo/'check_expanded_precision_resources.py'),
        'resources_sha256':h.sha(results/'asymmetric_raw_payoff_resources.json.gz'),
        'protocol_sha256':h.sha(results/'asymmetric_raw_payoff_protocol.json'),
        'fresh_configs':252,'old_diagonal_refs':48,'total_asymmetric_pool':300,
        'distinct_literal_blocks':len(cache),'complete_profiles':profiles,'own_formula_QPE_full_charges':plans,
        'integer_projection_histories':histories,'whole_winning_preparations':whole,'facts':facts,
        'independence':'New producer and confidence/charging frontend not imported. Prior independent exact formulae and primitive emitter reconstruction are explicitly reused. Mathematical scalar, CDF and projection proofs are separately checked.',
        'scope':'Matched finite pools and sufficient logical caps. Winners are symmetric9/9; improvement is due to newly tested diagonal9, not asymmetric allocation. No global optimum, theoretical priority or physical advantage.',
        'elapsed_seconds_diagnostic':time.monotonic()-start}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
    print('Passed independent resource audit:',profiles,'profiles,',plans,'full plans,',len(whole),'whole preparations.',flush=True)

if __name__=='__main__':main()
