"""Fresh production-node emission, independent metric and basis-path replay."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))

from integral_audit_support import packed_replay,register_values,independent_gate_accounting
from audit_integral_value_oracle import interpret,reference_table_key


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def normalise(counter):
    result=Counter()
    for gate,count in counter.items():result[{'TD':'Tdg','SD':'Sdg'}.get(gate,gate)]+=count
    return result


def build_from_profile(module,models,profile):
    c=profile['configuration']
    return module.build_graph(models[c['model_id']],F(c['cap']),q=c['q'],b_y=c['b_y'],
        beta_bits=c['beta_bits'],strike_bits=c['strike_bits'],table_bits=c['table_bits_requested'],
        outer_bound=F(c['L']),U_star=F(c['U_star']),D2=F(c['D2']),D1=F(c['D1']),
        kind=c['kind'],headroom_bits=c['headroom_bits'],distribution=c['distribution'],
        retained_bound=F(c['retained_coordinate_bound']),output_encoding=c.get('output_encoding','centred'))


def evaluate(root,profiles_path):
    root=Path(root);sys.path.insert(0,str(root));import integral_value_oracle as module
    profiles_path=Path(profiles_path);data=json.loads(profiles_path.read_text())
    assert data['source_sha256']==sha(root/'integral_value_oracle.py')
    for path,digest in data['input_hashes'].items():assert sha(root/path)==digest
    for path,digest in data['dependency_source_hashes'].items():assert sha(root/path)==digest
    models={m['model_id']:m for m in json.loads((root/'results/generalisation_protocol.json').read_text())['models']}
    tables_cache={};blocks_cache={};graph_records=[];node_basis_lanes=0
    for profile in data['profiles']:
        graph=build_from_profile(module,models,profile);tables=module.make_tables(graph,tables_cache)
        names=tuple(graph.inputs);samples=[]
        for seed in range(8):
            sample={name:([0,(1<<bits)-1,(1<<bits)//2,((1<<bits)-1)//2,
                ((1<<bits)//4),3*((1<<bits)//4),(1 if bits>0 else 0),max(0,(1<<bits)-2)][(seed+j)%8])
                for j,(name,bits) in enumerate(graph.inputs.items())}
            samples.append(sample)
        traces=[interpret(graph,sample,tables)[1] for sample in samples]
        expected_blocks={b['node']:b for b in profile['blocks']}
        assert set(expected_blocks)==set(graph.nodes)
        aggregate=Counter();raw=Counter();lookups=0;sum_depth=0;sum_five=0
        retained=sum(graph.inputs.values())+graph.q+2;peak=retained
        for name,node in graph.nodes.items():
            record=expected_blocks[name];distinct=tuple(dict.fromkeys(node.args));aliases=tuple(distinct.index(arg) for arg in node.args)
            widths=tuple(graph.inputs[arg] if arg in graph.inputs else graph.nodes[arg].output_bits for arg in node.args)
            assert record['argument_widths']==list(widths) and record['argument_alias_pattern']==list(aliases)
            table=tables.get(reference_table_key(name,node)) if node.kind in ('quantile','function') else None
            table_hash=hashlib.sha256(json.dumps(table,separators=(',',':')).encode()).hexdigest() if table is not None else None
            assert table_hash==record['table_sha256']
            meta=node.meta.as_dict() if hasattr(node.meta,'as_dict') else str(node.meta)
            cache_key=(profile['architecture'] if table is not None else 'shared_exact',node.kind,
                json.dumps(meta,sort_keys=True),node.bits,node.output_bits,graph.q,graph.width,widths,aliases,table_hash)
            if cache_key not in blocks_cache:
                circuit=module.vo.ValueCircuit(mode=profile['architecture'])
                registers=[circuit.alloc(widths[node.args.index(arg)]) for arg in distinct]
                args=[registers[index] for index in aliases];incoming=sum(map(len,registers))
                output=module.emit_node(circuit,node,args,graph.q,graph.width,table)
                accounting=independent_gate_accounting(circuit);metrics=circuit.resources()
                counts=normalise(accounting['expanded_literal_gate_counts'])
                assert counts==Counter(metrics['clifford_T_expanded_gate_counts'])
                assert accounting['actual_literal_T_dependency_depth']==metrics['explicit_7T_sequence_T_depth']
                cached={'counts':counts,'raw':Counter(accounting['primitive_gate_counts']),
                    'T':accounting['logical_T_count'],'depth':accounting['actual_literal_T_dependency_depth'],
                    'five':metrics['conservative_5layer_CCX_T_depth_upper'],'lookups':metrics['lookup_macro_invocations'],
                    'extra':circuit.wire_count-incoming,'output':len(output),
                    'tape_sha256':accounting['expanded_unitary_tape_sha256']}
                assignments=[{register:trace[arg] for register,arg in zip(registers,distinct)} for trace in traces]
                state,_=packed_replay(circuit,assignments)
                assert register_values(state,output,len(assignments))==[trace[name] for trace in traces],('Production node value discrepancy',profile['case_id'],name)
                for register,arg in zip(registers,distinct):assert register_values(state,register,len(assignments))==[trace[arg] for trace in traces]
                keep=set(output+tuple(w for register in registers for w in register))
                assert all(word==0 for wire,word in enumerate(state) if wire not in keep)
                cached['basis_paths_checked']=len(assignments);node_basis_lanes+=len(assignments)
                blocks_cache[cache_key]=cached
                del circuit,state
            block=blocks_cache[cache_key]
            assert block['counts']==Counter(record['primitive_counts'])
            assert block['T']==record['logical_T_count']
            assert block['depth']==record['explicit_T_depth_serial_upper']
            assert block['five']==record['conservative_5layer_T_depth_serial_upper']
            assert block['lookups']==record['lookup_invocations'] and block['output']==record['output_bits']
            assert block['extra']==record['extra_qubits_upper']
            aggregate.update(block['counts']);raw.update(block['raw']);lookups+=block['lookups']
            sum_depth+=block['depth'];sum_five+=block['five']
            peak=max(peak,retained+block['extra']);retained+=block['output']
        comparator_key=('comparator',graph.q)
        if comparator_key not in blocks_cache:
            comp=module.vo.ValueCircuit(mode='exact7t');a,b=comp.alloc(graph.q+1),comp.alloc(graph.q+1);flag=comp.alloc(1)[0]
            comp.compare_into(a,b,flag);accounting=independent_gate_accounting(comp);metrics=comp.resources()
            blocks_cache[comparator_key]={'counts':normalise(accounting['expanded_literal_gate_counts']),
                'depth':accounting['actual_literal_T_dependency_depth'],'five':metrics['conservative_5layer_CCX_T_depth_upper'],
                'wires':comp.wire_count,'tape_sha256':accounting['expanded_unitary_tape_sha256']}
        comp=blocks_cache[comparator_key]
        counts=Counter({gate:2*count for gate,count in aggregate.items()});counts.update(comp['counts'])
        uniform=sum(graph.inputs.values())+graph.q;counts['H']+=uniform
        assert counts==Counter(profile['primitive_counts'])
        assert raw==Counter(profile['pre_expansion_forward_counts'])
        assert counts['T']+counts['Tdg']==profile['logical_T_count']
        assert profile['lookup_invocations_forward_and_cleanup']==2*lookups
        assert profile['explicit_T_depth_serial_upper']==2*sum_depth+comp['depth']
        assert profile['conservative_5layer_T_depth_serial_upper']==2*sum_five+comp['five']
        peak=max(peak,retained+comp['wires']-2*(graph.q+1)-1)
        assert peak==profile['logical_qubits_upper']
        assert profile['reflection_core_qubits']==uniform+1
        assert profile['uniform_H']==uniform
        assert profile['endpoint_high_zero_bits']==1 and profile['maximum_probability_integer']==1<<graph.q
        graph_records.append({'case_id':profile['case_id'],'kind':graph.configuration['kind'],
            'output_encoding':graph.configuration['output_encoding'],'architecture':profile['architecture'],
            'node_count':len(graph.nodes),'logical_T_count':profile['logical_T_count'],
            'finite_integer_sample_outputs':[trace[graph.output] for trace in traces],
            'logical_qubits_upper':peak,'complete_primitive_count_verified':True})
    assert data['source_sha256']==sha(root/'integral_value_oracle.py'),'Oracle source changed during audit'
    return {'status':'independent production block and complete composition audit passed',
        'profile_count':len(graph_records),'unique_emitted_blocks':len(blocks_cache),
        'production_node_basis_paths_checked':node_basis_lanes,'records':graph_records,
        'profiles_sha256':sha(profiles_path),'source_sha256':data['source_sha256'],
        'actual_block_tape_hashes':[{'key_repr':repr(key),'actual_tape_sha256':value['tape_sha256']} for key,value in blocks_cache.items()],
        'scope':'Fresh production blocks, exact primitive and alias/count/depth replay plus retained sample paths. No exhaustive production tensor-grid or monetary-accuracy claim.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--profiles',default=None)
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_production_profiles_audit.json'))
    args=parser.parse_args();result=evaluate(args.root,args.profiles or Path(args.root)/'results/integral_production_profiles.json');result['audit_source_sha256']=sha(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','profile_count','unique_emitted_blocks','production_node_basis_paths_checked')}))

if __name__=='__main__':main()
