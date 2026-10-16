"""Independent shared-threshold financial, graph and complete-charge replay.

No new producer transformation, graph builder, probability-certificate,
budget, prior or cost frontend is used. Graphs are reconstructed from their
retained DAG and independently checked symbolically against the financial
specification. Retained audited primitive libraries regenerate scalar tables,
node arithmetic certificates and actual body/predicate tapes. Complete QPE
charges use the previously retained own-formula tape/binomial auditor.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import gzip
import hashlib
import importlib.util
import json
from itertools import product
from pathlib import Path
import sys
sys.dont_write_bytecode=True


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def clean(v):return {k:x for k,x in v.items() if x}
def lo(iv):return F(iv['lower_exact_fraction'])
def hi(iv):return F(iv['upper_exact_fraction'])
def linear_add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return clean(out)
def linear_scale(a,c):return clean({k:v*F(c) for k,v in a.items()})
def linear_constant(c):return {} if not c else {'constant':F(c)}


def import_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def reconstruct_graph(row,protocol,oracle):
    g=oracle.IntegralGraph(protocol['q'],protocol['outer_scalar_table_bits'])
    for j in range(row['model']['outer_dimension']):g.inp(f'iy{j}',protocol['outer_quantile_bits'])
    g.inp('inner_normal_label',protocol['inner_normal_bits']);g.inp('uniform_strike_label',protocol['strike_uniform_bits'])
    zero=[]
    for item in row['graph_nodes']:
        kind=item['kind'];meta=item['meta']
        if kind=='function':
            meta=oracle.ScalarSpec(meta['name'],F(meta['lo']),F(meta['span']),tuple(sorted((k,F(v)) for k,v in meta['parameters'].items())),tuple(map(F,meta['output_bounds'])))
        elif kind in ('constant','scale','clamp','clamp_signed','quantile'):meta=F(meta)
        elif kind in ('uniform','minimum'):meta=int(meta)
        elif meta=='None':meta=None
        name=g.op(kind,*item['args'],meta=meta,bits=item['bits'],output_bits=0,bounds=tuple(map(F,item['ideal_bounds'])))
        assert name==item['name']
        if kind=='constant':
            g.constant_names[meta]=name
            if meta==0:zero.append(name)
    assert len(zero)==1
    g.configuration=dict(row['graph_configuration']);g.semantic_nodes=dict(row['graph_semantic_nodes'])
    g.finish(zero[0])
    assert all(g.nodes[item['name']].output_bits==item['output_bits'] for item in row['graph_nodes'])
    assert g.configuration==row['graph_configuration']
    return g


def symbolic_semantics(g,row,protocol,parameters):
    w,r,mu,lam,mg,rg,d,s=parameters(row['model']);p=len(rg)
    exps={};values={};log_count=strike_count=0
    node_rows={i['name']:i for i in row['graph_nodes']}
    for name,node in g.nodes.items():
        args=[values[a] for a in node.args if a in values]
        if node.kind=='constant':v=linear_constant(node.meta)
        elif node.kind=='quantile':
            label=node.args[0]
            if label.startswith('iy'):
                j=int(label[2:]);v={'y'+str(j):F(1),'constant':F(protocol['outer_L'])}
            else:
                assert label=='inner_normal_label' and node.bits==protocol['inner_normal_bits']
                v={'normal':F(1),'constant':F(protocol['inner_normal_L'])}
        elif node.kind=='uniform':
            assert node.args==('uniform_strike_label',) and node.meta==protocol['strike_uniform_bits']
            v={'uniform':F(1)}
        elif node.kind=='add':v=linear_add(args[0],args[1])
        elif node.kind=='sub':v=linear_add(args[0],linear_scale(args[1],-1))
        elif node.kind=='scale':v=linear_scale(args[0],node.meta)
        elif node.kind=='positive_sub':
            v=linear_add(args[0],linear_scale(args[1],-1))
            target={'exp'+str(i):wi for i,wi in enumerate(w)};target['constant']=F(-1)
            assert v==target  # AGM premise proved by exact parameter centring.
        elif node.kind=='clamp':
            v=args[0]
            if v==dict({'exp'+str(i):wi for i,wi in enumerate(w)},constant=F(-1)):
                assert F(node.meta)>=F(row['graph_configuration']['U_star'])
            else:
                bounds=g.bounds[node.args[0]]
                assert 0<=bounds[0]<=bounds[1]<=F(node.meta),(name,bounds,node.meta)
        elif node.kind=='function':
            argument=linear_add(args[0],linear_constant(node.meta.lo))
            if node.meta.name=='exp':
                matching=[i for i in range(len(w)) if argument==linear_add({'y'+str(j):s[i][j] for j in range(p) if s[i][j]},linear_constant(d[i]))]
                assert len(matching)==1
                i=matching[0];assert i not in exps;exps[i]=name;v={'exp'+str(i):F(1)}
            elif node.meta.name=='log':
                tsum={'exp'+str(i):wi for i,wi in enumerate(w)}
                if argument==tsum:
                    v={'logT':F(1)};log_count+=1
                else:
                    assert argument=={'constant':F(1),'uniform':F(1,100)}
                    assert node.bits==protocol['strike_log_table_bits']
                    v={'strike':F(1)};strike_count+=1
            else:raise ValueError(('Unexpected financial scalar in new event graph',node.meta.name))
        else:raise ValueError(('Unsupported event node',node.kind))
        values[name]=clean(v)
    assert len(exps)==len(w) and log_count==strike_count==1
    semantic=g.semantic_nodes;J=F(g.configuration['common_endpoint_offset'])
    m=linear_add({'y'+str(j):rg[j] for j in range(p) if rg[j]},linear_constant(mg+J))
    assert values[semantic['lower_endpoint']]==m
    assert values[semantic['upper_endpoint']]==linear_add(m,{'logT':F(1)})
    assert values[semantic['threshold']]=={'strike':F(1),'normal':-lam,'constant':J}
    assert values[g.output]=={}  # Numeric placeholder is not the event flag.
    return len(values)


def emit_predicate(c,lower,threshold,upper,flag):
    # Independent XOR-comparison construction of the same clean interval.
    work_lower=c.alloc(1)[0];work_upper=c.alloc(1)[0]
    c.compare_into(lower,threshold,work_lower)
    c.compare_into(upper,threshold,work_upper)
    c.x(work_upper);c.toffoli(work_lower,work_upper,flag);c.x(work_upper)
    c.compare_into(upper,threshold,work_upper)
    c.compare_into(lower,threshold,work_lower)
    c.release((work_lower,work_upper))


def check_profile(profile,g,tables,oracle,vo,bodycache):
    old=oracle.counted_composition(g,tables,bodycache,profile['architecture'])
    assert profile['blocks']==old['blocks']
    forward=Counter()
    for block in old['blocks']:
        forward.update(block['primitive_counts'])
        assert block['logical_T_count']==block['primitive_counts'].get('T',0)+block['primitive_counts'].get('Tdg',0)
    assert clean(forward)==clean(Counter(profile['primitive_counts_forward_body']))
    c=vo.ValueCircuit(store_gates=False,mode=profile['architecture'])
    lower,threshold,upper=(c.alloc(g.width) for _ in range(3));flag=c.alloc(1)[0]
    emit_predicate(c,lower,threshold,upper,flag);metric=c.resources()
    predicate=Counter(metric['clifford_T_expanded_gate_counts'])
    assert clean(predicate)==clean(Counter(profile['predicate_primitive_counts']))
    # Four complete clean comparison ladders, plus two X and one seven-T CCX.
    assert predicate['T']+predicate['Tdg']==4*42*g.width+7
    total=Counter({gate:2*count for gate,count in forward.items()});total.update(predicate)
    input_count=sum(g.inputs.values());total['H']+=input_count
    assert clean(total)==clean(Counter(profile['primitive_counts']))
    assert profile['logical_T_count']==total['T']+total['Tdg']
    assert profile['uniform_H']==input_count and profile['reflection_core_qubits']==input_count+1
    assert profile['active_qubits_after_clean_prep']==input_count+1
    assert profile['integer_fractional_bits']==g.q and profile['integer_word_width']==g.width
    # This alternate clean predicate can recycle scratch differently. Count
    # equality is exact; depth/peak use the producer's literal predicate and
    # are separately covered by its exact clean-phase audit.
    assert profile['explicit_T_depth_serial_upper']>=2*sum(x['explicit_T_depth_serial_upper'] for x in old['blocks'])
    assert profile['conservative_5layer_T_depth_serial_upper']>=2*sum(x['conservative_5layer_T_depth_serial_upper'] for x in old['blocks'])
    assert profile['scalar_table_entries']==old['scalar_table_entries']
    assert profile['scalar_table_program_bits']==old['scalar_table_program_bits']
    return total


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=Path,required=True);ap.add_argument('--repo',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    sys.path.insert(0,str(args.repo))
    import boundary_value_oracle as oracle
    import boundary_value_certificate as certifier
    import value_oracle as vo
    from flint import arb,ctx
    from boundary_grid_certificate import ball,endpoints
    independent=import_module('independent_factor_math',args.step/'review/audit_factor_basis.py')
    qdir=args.step/'adversarial/quantum';rp=qdir/'threshold_results'
    prototype=qdir/'shared_threshold_oracle.py'
    path=rp/'shared_threshold_resources.json';data=json.loads(path.read_text());protocol=json.loads((rp/'shared_threshold_protocol.json').read_text())
    factorpath=args.step/'resources/results/factor_basis_quantum_resources.json';factor=json.loads(factorpath.read_text())
    predicate_audit=json.loads((args.output.parent/'shared_threshold_predicate_audit.json').read_text())
    assert predicate_audit['status']=='passed' and predicate_audit['prototype_source_sha256']==sha(prototype)
    assert data['source_sha256']==protocol['source_sha256']==sha(prototype)
    assert data['factor_resources_sha256']==protocol['factor_resources_sha256']==sha(factorpath)
    assert data['protocol_sha256']==sha(rp/'shared_threshold_protocol.json')
    assert data['dependencies_sha256']==protocol['dependencies_sha256']
    assert data['inputs_sha256']==protocol['inputs_sha256']
    for name,value in protocol['dependencies_sha256'].items():assert sha(args.repo/name)==value,name
    for name,value in protocol['inputs_sha256'].items():assert sha(args.repo/'results'/name)==value,name
    expected_pairs=[[r['original_model_id'],r['selected_candidate_id']] for r in factor['paired_comparison']]
    assert protocol['pairs']==expected_pairs and len(expected_pairs)==len(data['rows'])==4
    assert protocol['cap']=='1' and protocol['inner_normal_bits']==10 and protocol['strike_uniform_bits']==4
    assert protocol['outer_scalar_table_bits']==12 and protocol['strike_log_table_bits']==6
    assert protocol['no_precision_repairs'] is True and protocol['financial_price_answers_excluded'] is True
    archivepath=rp/data['scalar_archive_file'];assert sha(archivepath)==data['scalar_archive_sha256']
    with gzip.open(archivepath,'rt') as f:archive=json.load(f)
    archived={item['certificate_sha256']:item for item in archive}
    assert len(archived)==len(archive)==data['unique_scalar_certificates']
    for key,item in archived.items():
        assert digest(item['certificate'])==key and digest(item['packed_table'])==item['packed_table_sha256']
    # Import the old own-formula complete-charge checker, without invoking it.
    previous_argv=sys.argv[:]
    sys.argv=[str(args.step/'resources/audit_factor_basis_quantum_resources.py'),'--repo',str(args.repo),'--results',str(args.step/'resources/results'),'--producer',str(args.step/'resources/factor_basis_quantum_resources.py'),'--certificate',str(args.step/'certificate/results/factor_basis_certificates.json')]
    cost=import_module('previous_independent_cost_audit',args.step/'resources/audit_factor_basis_quantum_resources.py')
    sys.argv=previous_argv
    tapes=json.loads((args.repo/'results/qft_rotation_tapes.json').read_text())
    qftcert=json.loads((args.repo/'results/qft_rotation_certificates.json').read_text())
    assert sha(args.repo/'results/qft_rotation_tapes.json')==qftcert['tapes_sha256']
    qftbounds={r['key']:F(r['operator_norm_upper_rational']) for r in qftcert['certificates']}
    norms=json.loads((args.repo/'results/gaussian_derivative_norms.json').read_text())
    nu={j:F(next(r['norm_upper'] for r in norms['unit_norms'] if r['unit_derivative_order']==j and r['passed'])) for j in (0,1)}
    ctx.prec=512;tablecache={};blockcache={};facts=[];packed_variants={};scalarchecks=profiles=plans=semanticchecks=domainvertices=0
    for row in data['rows']:
        original,candidate=row['original_model_id'],row['candidate_id']
        assert [original,candidate] in expected_pairs
        prep=next(p for p in factor['preparations'] if p['original_model_id']==original and p['candidate_id']==candidate)
        old=next(r for r in factor['records'] if r['original_model_id']==original and r['factor_basis_candidate_id']==candidate and r['kind']=='positive_r1')
        assert row['model']==prep['model'] and row['retained_direct_selected_strategy']==old['selected_strategy']
        assert F(old['base_graph_certificate']['configuration']['monetary_signal_normaliser'])==1
        g=reconstruct_graph(row,protocol,oracle)
        semanticchecks+=symbolic_semantics(g,row,protocol,independent.parameters)
        # Independently prove the finite U cap over the exact enclosing box.
        ctx.prec=768
        w,r,mu,lam,mg,rg,d,s=independent.parameters(row['model'])
        lower,upper=map(F,prep['domain']['coordinate_bounds']);tstar=F(prep['domain']['T_star'])
        for vertices in product((lower,upper),repeat=len(rg)):
            value=sum((ball(wi)*ball(di+independent.dot(si,list(vertices))).exp() for wi,di,si in zip(w,d,s)),arb(0))
            assert endpoints(value)[1]<=tstar
            domainvertices+=1
        ctx.prec=512
        tables,ordinates=oracle.make_tables(g,tablecache,True)
        freshbody=certifier.certify_graph(g,tables,ordinates)
        actual=row['probability_certificate'];assert actual['status']=='passed' and actual['body_arithmetic_certificate']==freshbody
        assert actual['nodes']==freshbody['nodes'] and freshbody['status']=='passed'
        assert freshbody['probability_error_upper']=='0'  # Placeholder only.
        for name,cc in ordinates.items():
            h=digest(cc);assert h in archived and cc==archived[h]['certificate']
            key=oracle.table_key(name,g.nodes[name]);packed=tables[key]
            item=archived[h];stored_width=item['graph_word_width']
            if cc['schema']=='qbp-integral-quantile-v1':
                expected=tuple(n['offset_integer'] for n in cc['nodes'])
            else:
                ordinates_int=[n['integer'] for n in cc['nodes']]
                expected=[]
                for j,v in enumerate(ordinates_int):
                    delta=ordinates_int[j+1]-v if j+1<len(ordinates_int) else 0
                    expected.append(v|(abs(delta)<<g.width)|(int(delta<0)<<(2*g.width)))
                expected=tuple(expected)
            assert packed==expected
            if g.width==stored_width:assert digest(packed)==item['packed_table_sha256']
            packed_variants[(h,g.width)]={'certificate_sha256':h,'graph_word_width':g.width,'packed_table_sha256':digest(packed),'packed_table':list(packed)}
            scalarchecks+=len(cc['nodes'])
        sem=g.semantic_nodes;nodeerrors={r['node']:F(r['absolute_output_error_upper']) for r in freshbody['nodes']}
        q=ordinates[sem['normal_plus_L']];a_lo=lo(q['mass_interval']);eN=F(q['maximum_node_error_upper'])
        assert nodeerrors[sem['normal_plus_L']]==0
        N=1<<protocol['inner_normal_bits'];Nt=1<<protocol['strike_uniform_bits']
        D=(1-a_lo)/2+a_lo/(2*N)+nu[0]*eN
        el,eu,ez=(nodeerrors[sem[k]] for k in ('lower_endpoint','upper_endpoint','threshold'))
        arithmetic=nu[0]*(el+eu+2*ez)/lam
        strike=F(1,100)**2*(nu[1]/lam**2+nu[0]/lam)/(12*Nt*Nt)
        eta=min(F(1),2*D+arithmetic+strike)
        for key,value in {'probability_error_upper':eta,'monetary_value_error_upper':eta,
            'inner_normal_CDF_error_upper':D,'inner_normal_node_rounding_error_upper':eN,
            'arithmetic_probability_sandwich_upper':arithmetic,'strike_midpoint_bias_upper':strike,
            'inner_normal_interval_bias_upper':2*D}.items():assert F(actual[key])==value,key
        assert actual['endpoint_arithmetic_errors']=={'lower':str(el),'upper':str(eu),'threshold':str(ez)}
        assert actual['inner_normal_mass_interval']==q['mass_interval']
        budget=row['complete_deterministic_budget'];oldbudget=old['complete_deterministic_budget']
        coefficient=prep['coefficient_certificate']['c0'];cr=(hi(coefficient)-lo(coefficient))/2
        financial=prep['R1_selected_financial_transfer'];mass=financial['mass_interval'];mh=hi(mass);mr=(hi(mass)-lo(mass))/2
        beta=F(financial['financial_transfer_upper']);det=cr+beta+mh*eta+mr
        assert budget['mass_interval']==mass and F(budget['external_factor_upper'])==mh
        assert F(budget['external_factor_radius'])==mr and F(budget['coefficient_radius_upper'])==cr
        assert F(budget['financial_transfer_upper'])==beta and F(budget['total_deterministic_undiscounted_price_error_upper'])==det
        assert F(budget['arithmetic_before_mass_upper'])==arithmetic
        assert F(budget['auxiliary_components']['bias_upper'])==2*D+strike
        assert det<F(1,100) and eta==2*D+arithmetic+strike
        previous=old['coefficient_informed_prior'];meanlo=F(previous['finite_ideal_mean_lower']);meanhi=F(previous['finite_ideal_mean_upper'])
        prior=row['coefficient_informed_prior'];pl=max(F(0),meanlo-eta);pu=min(F(1),meanhi+eta)
        assert 0<=pl<=pu<=1
        assert F(prior['probability_lower'])==pl and F(prior['probability_upper'])==pu
        assert F(prior['finite_ideal_mean_lower'])==meanlo and F(prior['finite_ideal_mean_upper'])==meanhi
        assert prior['native_positive'] is True and prior['previous_mean_information_sha256']==digest(previous)
        for architecture,family in row['profiles'].items():
            assert set(family)=={'centred'}
            profile=family['centred'];assert profile['architecture']==architecture and profile['configuration']==g.configuration
            assert profile['source_sha256']==sha(prototype)
            check_profile(profile,g,tables,oracle,vo,blockcache);profiles+=1
        for strategy in row['strategies']:
            out=strategy['result'];profile=row['profiles'][strategy['architecture']]['centred']
            if out['kind']=='deterministic_correction_prior':
                assert det+mh*(pu-pl)/2==F(out['total_price_error_upper'])<=F(1,100)
                continue
            assert out['kind']=='complete_quantum_resource_plan' and out['correction_prior']==prior
            assert F(out['deterministic_error_upper'])==det and F(out['monetary_signal_factor_upper'])==mh
            assert len(out['flags'])==1
            flag=out['flags'][0]
            assert F(flag['probability_tolerance'])==(F(1,100)-det)/mh
            assert F(flag['prior_upper'])==pu and F(flag['plan']['alpha'])==F(1,3600)
            T=cost.verify_plan(flag,profile,tapes,qftbounds);plans+=1
            assert out['T_count']==T
            stat=mh*(F(flag['plan']['amplitude_error_upper'])+F(flag['plan']['postprocess_allowance']))
            assert stat==F(out['statistical_price_error_upper'])
            assert det+stat==F(out['total_price_error_upper'])<=F(1,100)
            assert F(out['failure_probability_upper'])==F(flag['plan']['binomial_failure_upper'])<=F(1,3600)
        selected=min(row['strategies'],key=lambda x:(x['result']['T_count'],x['architecture']))
        assert selected==row['selected_strategy']
        out=selected['result'];oldout=old['selected_strategy']['result']
        assert [(f['plan']['M'],f['plan']['repetitions']) for f in out['flags']]==[(f['plan']['M'],f['plan']['repetitions']) for f in oldout['flags']]
        facts.append({'original_model_id':original,'candidate_id':candidate,'deterministic_error_upper':str(det),
            'conditional_probability_error_upper':str(eta),'inner_normal_CDF_error_upper':str(D),
            'implemented_prior_lower':str(pl),'implemented_prior_upper':str(pu),
            'normaliser':'1','no_extra_common_normal_mass':True,
            'retained_direct_T':oldout['T_count'],'new_selected_T':out['T_count'],
            'exact_saving_fraction':str(1-F(out['T_count'],oldout['T_count'])),
            'algorithm_phase_width_and_repetitions_unchanged':True})
    variants_path=args.output.parent/'shared_threshold_repacked_tables.json'
    variants={'schema':'qbp-shared-threshold-exact-table-repacking-v1','source_archive_sha256':sha(archivepath),
        'variants':list(packed_variants.values()),'scope':'Exact width-specific integer packing of retained synthetic scalar ordinates. Original archive and resource outcomes are unchanged.'}
    variants_path.write_text(json.dumps(variants,indent=2)+'\n')
    result={'schema':'qbp-independent-shared-threshold-financial-audit-v1','status':'passed','financial_cases':len(facts),
        'symbolic_graph_node_checks':semanticchecks,'independent_768bit_domain_vertices':domainvertices,
        'fresh_scalar_ordinate_checks':scalarchecks,'fresh_full_profiles':profiles,'complete_confidence_plans':plans,
        'width_specific_table_variants':len(packed_variants),'facts':facts,
        'resource_sha256':sha(path),'protocol_sha256':sha(rp/'shared_threshold_protocol.json'),
        'prototype_source_sha256':sha(prototype),'audit_source_sha256':sha(__file__),
        'predicate_audit_sha256':sha(args.output.parent/'shared_threshold_predicate_audit.json'),
        'previous_independent_charge_auditor_sha256':sha(args.step/'resources/audit_factor_basis_quantum_resources.py'),
        'factor_resources_sha256':sha(factorpath),'repacked_tables_sha256':sha(variants_path),
        'scope':'Independent exact symbolic financial DAG, high-precision finite-support proof, fresh scalar/node-arithmetic reconstruction, normalised finite-event CDF/arithmetic/strike errors, new implemented priors and complete canonical QPE charges. Exact toy clean phases are separately linked. Synthetic evidence only; no physical device or theoretical originality claim.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','financial_cases','symbolic_graph_node_checks','independent_768bit_domain_vertices','fresh_scalar_ordinate_checks','fresh_full_profiles','complete_confidence_plans')}))

if __name__=='__main__':main()
