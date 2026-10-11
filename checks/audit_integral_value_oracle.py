"""Independent whole financial integer-oracle and primitive composition audit.

This module does not use graph.evaluate or Circuit.run to obtain reference
values. Small retained circuits are replayed exactly in packed basis lanes.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys

import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))

from integral_audit_support import packed_replay,register_values,independent_gate_accounting
from integral_phase_support import lookup_clean_isometry


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def exact_nearest_even(value):
    value=F(value);lower=value.numerator//value.denominator;rest=value-lower
    if rest>F(1,2) or rest==F(1,2) and lower%2:return lower+1
    return lower


def reference_scale(value,q):
    value=F(value)
    if value.denominator&(value.denominator-1)==0:
        return value.numerator,value.denominator.bit_length()-1
    return exact_nearest_even(value*(1<<q)),q


def reference_table_key(name,node):
    if node.kind=='quantile':return ('quantile',node.bits,str(node.meta))
    return ('scalar',json.dumps(node.meta.as_dict(),sort_keys=True),node.bits)


def interpret(graph,inputs,tables,*,domains=False):
    unit=1<<graph.q;modulus=1<<graph.width;values=dict(inputs);raw={};domain_records=[]
    for name,node in graph.nodes.items():
        arguments=[values[arg] for arg in node.args]
        if node.kind=='constant':value=exact_nearest_even(F(node.meta)*unit)
        elif node.kind=='add':value=arguments[0]+arguments[1]
        elif node.kind in ('sub','signed_sub'):value=arguments[0]-arguments[1]
        elif node.kind=='positive_sub':value=max(0,arguments[0]-arguments[1])
        elif node.kind=='mul':value=arguments[0]*arguments[1]//unit
        elif node.kind=='scale':
            integer,drop=reference_scale(node.meta,graph.q);value=arguments[0]*integer//(1<<drop)
        elif node.kind=='uniform':value=(2*arguments[0]+1)*(unit//(1<<(node.meta+1)))
        elif node.kind=='minimum':value=min(arguments)
        elif node.kind=='quantile':value=tables[reference_table_key(name,node)][arguments[0]]
        elif node.kind=='function':
            span=F(node.meta.span);exponent=span.numerator.bit_length()-span.denominator.bit_length()
            fractional=graph.q+exponent-node.bits
            index,remainder=divmod(arguments[0],1<<fractional)
            table=tables[reference_table_key(name,node)]
            assert 0<=index<len(table),('Function address outside retained table',name,index,len(table))
            packed=table[index];base=packed&(modulus-1);magnitude=(packed>>graph.width)&(modulus-1)
            decreasing=(packed>>(2*graph.width))&1
            adjustment=magnitude*remainder//(1<<fractional)
            value=base-adjustment if decreasing else base+adjustment
            assert value>=0,('Signed interpolant underflow',name,value)
            assert base>=magnitude if decreasing else base+magnitude<modulus
        elif node.kind in ('clamp','clamp_signed'):
            argument=arguments[0]
            if node.kind=='clamp_signed' and argument>=modulus//2:argument-=modulus
            value=max(0,min(argument,exact_nearest_even(F(node.meta)*unit)))
        else:raise ValueError(('Unknown independent node',node.kind))
        raw[name]=value
        if node.kind=='signed_sub':assert -modulus//2<=value<modulus//2,('Signed-domain wrap',name,value)
        else:assert 0<=value<1<<node.output_bits,('Unsigned-domain wrap',name,value,node.output_bits)
        values[name]=value% (1<<node.output_bits)
        if domains:domain_records.append((name,node.kind,value))
    return values[graph.output],values,raw,domain_records


def signed_interpolation_blocks(module):
    endpoints=(11,3,12,1,10);q=4;width=7;table_bits=2
    spec=module.ScalarSpec('audit_nonmonotone',F(0),F(1),(),(F(0),F(1)))
    table=[]
    for j,base in enumerate(endpoints):
        delta=endpoints[j+1]-base if j+1<len(endpoints) else 0
        table.append(base|(abs(delta)<<width)|((delta<0)<<(2*width)))
    records=[]
    for mode in ('exact7t','unitary8t'):
        circuit=module.arch.ArchitectureCircuit(mode=mode)
        inp=circuit.alloc(width);target=circuit.alloc(width);snapshot=circuit.snapshot()
        node=module.IntegralNode('function',('input',),spec,table_bits,width)
        output=module.emit_node(circuit,node,[inp],q,width,table)
        body_end=circuit.gate_count
        for src,dst in zip(output,target):circuit.cnot(src,dst)
        circuit.uncompute(snapshot,end_gate_count=body_end)
        assignments=[];expected=[]
        for argument in range(17):
            cell,rem=divmod(argument,4);delta=endpoints[cell+1]-endpoints[cell] if cell+1<len(endpoints) else 0
            adjustment=abs(delta)*rem//4;value=endpoints[cell]+(-adjustment if delta<0 else adjustment)
            for dirty in range(1<<width):assignments.append({inp:argument,target:dirty});expected.append(dirty^value)
        state,summary=packed_replay(circuit,assignments)
        assert register_values(state,target,len(assignments))==expected
        assert register_values(state,inp,len(assignments))==[a[inp] for a in assignments]
        retained=set(inp+target);assert all(word==0 for wire,word in enumerate(state) if wire not in retained)
        reversed_state,_=packed_replay(circuit,assignments,reverse=True,initial_state=state)
        assert register_values(reversed_state,target,len(assignments))==[a[target] for a in assignments]
        assert register_values(reversed_state,inp,len(assignments))==[a[inp] for a in assignments]
        accounting=independent_gate_accounting(circuit)
        actual=circuit.resources();assert accounting['logical_T_count']==actual['logical_T_count']
        assert Counter(accounting['primitive_gate_counts'])==Counter(actual['primitive_gate_counts'])
        assert summary['lookup_invocations']==actual['lookup_macro_invocations']==4
        records.append({'mode':mode,'basis_lanes':len(assignments),'negative_delta_cells':2,
            'upper_endpoint_address':4,'dirty_target_words':1<<width,'gate_accounting':accounting,
            'lookup_invocations':summary['lookup_invocations']})
    return records


def primitive_phase_checks(module):
    records=[]
    for table in ((0,7,1,5),(7,0,0,3),(0,0,0,0),(1,2,4,7)):
        for mode in ('exact7t','unitary8t'):
            for enabled in (False,True):
                plan=module.arch.lookup_plan((0,1),(2,3,4),table,mode,
                    (5,6) if enabled else (5,6,7),8 if enabled else None)
                records.append(lookup_clean_isometry(plan))
    return records


def saturated_comparison_and_split_checks(module,root):
    records=[]
    for q in range(2,8):
        circuit=module.vo.ValueCircuit(mode='exact7t')
        comparison=circuit.alloc(q);zero=circuit.alloc(1);value=circuit.alloc(q+1);flag=circuit.alloc(1)
        circuit.compare_into(comparison+zero,value,flag[0])
        assignments=[];expected=[]
        for a in range((1<<q)+1):
            positive=max(2*a-(1<<q),0);negative=max((1<<q)-2*a,0)
            assert positive-negative==2*a-(1<<q)
            assert 0<=positive<=1<<q and 0<=negative<=1<<q
            for uniform in range(1<<q):
                for dirty in (0,1):
                    assignments.append({comparison:uniform,value:a,flag:dirty});expected.append(dirty^(uniform<a))
        state,_=packed_replay(circuit,assignments)
        assert register_values(state,flag,len(assignments))==expected
        assert register_values(state,value,len(assignments))==[x[value] for x in assignments]
        assert register_values(state,comparison,len(assignments))==[x[comparison] for x in assignments]
        retained=set(comparison+value+flag);assert all(word==0 for wire,word in enumerate(state) if wire not in retained)
        records.append({'q':q,'basis_lanes':len(assignments),'saturated_endpoint_checked':True,
            'all_signed_split_thresholds_checked':True})
    model=next(m for m in json.loads((Path(root)/'results/generalisation_protocol.json').read_text())['models'] if m['model_id']=='n2_lambda0.2')
    for encoding,expected in (('centred',16),('positive_part',0),('negative_part',0)):
        graph=module.build_graph(model,25,q=5,b_y=1,beta_bits=1,strike_bits=1,table_bits=2,outer_bound=2,U_star=0,D2=16,kind='integral_r2',output_encoding=encoding)
        out,_,_,_=interpret(graph,{},{});assert out==expected,('Degenerate split output',encoding,out)
    return records


def full_financial_fixtures(module,root):
    models=json.loads((Path(root)/'results/generalisation_protocol.json').read_text())['models']
    model=next(m for m in models if m['model_id']=='n2_lambda0.2')
    records=[];cache={};total_value=total_flag=0;scalar_records={};resource_cache={}
    from audit_integral_scalar_certificate import audit_certificate
    for distribution in ('quantile','uniform_box'):
        cases=[(kind,encoding) for kind in ('integral_r2','direct_r2','positive_r1') for encoding in (('centred',) if kind=='positive_r1' else ('centred','positive_part','negative_part'))]
        for kind,encoding in cases:
            graph=module.build_graph(model,25,q=5,b_y=1,beta_bits=1,strike_bits=1,
                table_bits=2,outer_bound=2,U_star=F(1,8),D2=16,
                kind=kind,distribution=distribution,output_encoding=encoding)
            tables,certificates=module.make_tables(graph,cache,return_certificates=True)
            for name,certificate in certificates.items():
                node=graph.nodes[name];actual_words=list(tables[reference_table_key(name,node)])
                if node.kind=='quantile':
                    expected_words=[row['offset_integer'] for row in certificate['nodes']]
                else:
                    endpoints=[row['integer'] for row in certificate['nodes']];expected_words=[]
                    for j,left in enumerate(endpoints):
                        delta=endpoints[j+1]-left if j+1<len(endpoints) else 0
                        expected_words.append(left|(abs(delta)<<graph.width)|(int(delta<0)<<(2*graph.width)))
                assert actual_words==expected_words
                key=hashlib.sha256(json.dumps(certificate,sort_keys=True).encode()).hexdigest()
                if key not in scalar_records:scalar_records[key]=audit_certificate(certificate)
            input_names=tuple(graph.inputs)
            financial_inputs=[dict(zip(input_names,values)) for values in product(*(range(1<<graph.inputs[name]) for name in input_names))]
            target_values=[];traces=[];maximums={};minimums={}
            for sample in financial_inputs:
                target,values,raw,domains=interpret(graph,sample,tables,domains=True)
                assert 0<=target<=1<<graph.q
                target_values.append(target);traces.append(values)
                for name,_,value in domains:
                    maximums[name]=max(maximums.get(name,value),value)
                    minimums[name]=min(minimums.get(name,value),value)
            for mode in ('exact7t','unitary8t'):
                value_circuit,value_inputs,target=module.compile_value(graph,tables,mode)
                value_assignments=[];expected=[];retained_inputs=[]
                for sample,value in zip(financial_inputs,target_values):
                    for dirty in range(1<<len(target)):
                        assignment={value_inputs[name]:sample[name] for name in input_names};assignment[target]=dirty
                        value_assignments.append(assignment);expected.append(dirty^value)
                state,_=packed_replay(value_circuit,value_assignments)
                assert register_values(state,target,len(value_assignments))==expected
                for name in input_names:
                    assert register_values(state,value_inputs[name],len(value_assignments))==[sample[value_inputs[name]] for sample in value_assignments]
                retained=set(target+tuple(w for register in value_inputs.values() for w in register))
                assert all(word==0 for wire,word in enumerate(state) if wire not in retained)
                restored,_=packed_replay(value_circuit,value_assignments,reverse=True,initial_state=state)
                assert register_values(restored,target,len(value_assignments))==[a[target] for a in value_assignments]
                value_accounting=independent_gate_accounting(value_circuit)
                assert Counter(value_accounting['primitive_gate_counts'])==Counter(value_circuit.resources()['primitive_gate_counts'])
                assert value_accounting['logical_T_count']==value_circuit.resources()['logical_T_count']
                del state,restored
                flag_circuit,flag_inputs,comparison,flag=module.compile_flag(graph,tables,mode)
                flag_assignments=[];expected_flags=[]
                for sample,value in zip(financial_inputs,target_values):
                    for uniform in range(1<<graph.q):
                        for dirty in (0,1):
                            assignment={flag_inputs[name]:sample[name] for name in input_names}
                            assignment[comparison]=uniform;assignment[(flag,)]=dirty
                            flag_assignments.append(assignment);expected_flags.append(dirty^(uniform<value))
                state,flag_summary=packed_replay(flag_circuit,flag_assignments)
                assert register_values(state,(flag,),len(flag_assignments))==expected_flags
                for register in list(flag_inputs.values())+[comparison]:
                    assert register_values(state,register,len(flag_assignments))==[sample[register] for sample in flag_assignments]
                retained=set(comparison+(flag,)+tuple(w for register in flag_inputs.values() for w in register))
                assert all(word==0 for wire,word in enumerate(state) if wire not in retained)
                restored,_=packed_replay(flag_circuit,flag_assignments,reverse=True,initial_state=state)
                assert register_values(restored,(flag,),len(flag_assignments))==[sample[(flag,)] for sample in flag_assignments]
                flag_accounting=independent_gate_accounting(flag_circuit)
                actual=flag_circuit.resources()
                assert Counter(flag_accounting['primitive_gate_counts'])==Counter(actual['primitive_gate_counts'])
                assert flag_accounting['logical_T_count']==actual['logical_T_count']
                normalised=Counter({{'TD':'Tdg','SD':'Sdg'}.get(gate,gate):count for gate,count in flag_accounting['expanded_literal_gate_counts'].items()})
                assert normalised==Counter(actual['clifford_T_expanded_gate_counts'])
                assert flag_accounting['actual_literal_T_dependency_depth']==actual['explicit_7T_sequence_T_depth']
                complete_profile=module.counted_composition(graph,tables,resource_cache,mode=mode)
                normalised['H']+=len(flag_circuit.uniform_preparation_gates)
                assert normalised==Counter(complete_profile['primitive_counts'])
                assert complete_profile['logical_T_count']==flag_accounting['logical_T_count']
                assert complete_profile['lookup_invocations_forward_and_cleanup']==flag_summary['lookup_invocations']
                assert complete_profile['explicit_T_depth_serial_upper']>=flag_accounting['actual_literal_T_dependency_depth']
                assert complete_profile['conservative_5layer_T_depth_serial_upper']>=complete_profile['explicit_T_depth_serial_upper']
                assert complete_profile['logical_qubits_upper']>=flag_circuit.wire_count
                assert flag_summary['lookup_invocations']==actual['lookup_macro_invocations']
                expected_core=tuple(w for register in flag_inputs.values() for w in register)+comparison+(flag,)
                assert flag_circuit.reflection_core==expected_core
                assert flag_circuit.uniform_preparation_gates==[('H',w) for w in expected_core[:-1]]
                total_value+=len(value_assignments);total_flag+=len(flag_assignments)
                records.append({'kind':kind,'output_encoding':encoding,'distribution':distribution,'mode':mode,
                    'model_id':model['model_id'],'cap':25,'q':graph.q,'width':graph.width,
                    'financial_input_states':len(financial_inputs),'value_dirty_target_lanes':len(value_assignments),
                    'flag_comparison_dirty_target_lanes':len(flag_assignments),
                    'value_output_integer_min':min(target_values),'value_output_integer_max':max(target_values),
                    'node_count':len(graph.nodes),'value_gate_accounting':value_accounting,
                    'flag_gate_accounting':flag_accounting,
                    'complete_profile':complete_profile,
                    'node_integer_ranges':{name:[minimums[name],maximums[name]] for name in maximums},
                    'finite_outputs':target_values,
                    'qualification':'Coarse arithmetic and phase-composition fixture; no monetary accuracy claim is inferred'})
                del state,restored,value_circuit,flag_circuit
    for distribution in ('quantile','uniform_box'):
        for kind in ('integral_r2','direct_r2'):
            for mode in ('exact7t','unitary8t'):
                group={r['output_encoding']:r['finite_outputs'] for r in records if r['distribution']==distribution and r['kind']==kind and r['mode']==mode}
                for centre,plus,minus in zip(group['centred'],group['positive_part'],group['negative_part']):
                    assert plus==max(2*centre-(1<<5),0)
                    assert minus==max((1<<5)-2*centre,0)
                    assert plus-minus==2*centre-(1<<5)
    return {'status':'complete small financial value/flag/clean-inverse audit passed',
        'graph_mode_configurations':len(records),'value_dirty_target_basis_lanes':total_value,
        'flag_comparison_dirty_target_basis_lanes':total_flag,'records':records,
        'independent_scalar_node_certificates':list(scalar_records.values()),
        'independent_scalar_and_quantile_nodes':sum(r['nodes'] for r in scalar_records.values())}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_value_oracle_audit.json'))
    args=parser.parse_args();root=Path(args.root);sys.path.insert(0,str(root))
    import integral_value_oracle as module
    whole=full_financial_fixtures(module,root)
    result={'status':'whole finite financial integer-oracle audit passed',
        'whole_financial_graph_audit':whole,
        'saturated_comparison_and_split_predicate_audit':saturated_comparison_and_split_checks(module,root),
        'signed_interpolation_block_audit':signed_interpolation_blocks(module),
        'primitive_phase_audit':primitive_phase_checks(module),
        'source_sha256':digest(root/'integral_value_oracle.py'),'audit_source_sha256':digest(__file__),
        'dependency_sha256':{name:digest(root/name) for name in ['integral_value_certificate.py','architecture.py','value_oracle.py','reversible.py','quantum_gates.py']},
        'input_sha256':{name:digest(root/name) for name in ['results/generalisation_protocol.json','results/gaussian_derivative_norms.json']},
        'scope':'Whole financial finite-integer composition and unitary resource audit of coarse fixtures; production precision certificates and large graph profiles must be checked separately'}
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'signed_block_modes':2,'primitive_phase_cases':len(result['primitive_phase_audit']),'whole_graph_modes':whole['graph_mode_configurations']}))

if __name__=='__main__':main()
