"""Independent binding of frozen abstract samplers to the financial flag.

No producer, allocator, quantum charging API or estimator is imported.
Source inspection checks the allocation chain; own sequential allocation
formulas and macro expansion bind the high-zero omission, phase adjoint,
concrete reflected core, and full measurement/reset without altering costs.
"""
from pathlib import Path
from collections import Counter
import argparse, ast, hashlib, json

HERE=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def clean(v):return {k:n for k,n in v.items() if n}


def source_allocation_checks(repo):
    texts={n:(repo/n).read_text() for n in ['integral_value_oracle.py','value_oracle.py','architecture.py','reversible.py']}
    tree=ast.parse(texts['integral_value_oracle.py'])
    flag=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='compile_flag')
    src=ast.get_source_segment(texts['integral_value_oracle.py'],flag)
    ordered=["inputs={name:c.alloc(bits,name)","comparison=c.alloc(graph.q,'uniform comparator')", "zero=c.alloc(1,'comparator high zero')", "flag=c.alloc(1,'success flag')[0]", "snap=c.snapshot()", "value=emit_body"]
    positions=[src.index(x) for x in ordered];assert positions==sorted(positions)
    assert "c.reflection_core=tuple(w for reg in inputs.values() for w in reg)+comparison+(flag,)" in src
    for filename,class_name in [('value_oracle.py','ValueCircuit'),('architecture.py','ArchitectureCircuit')]:
        cls=next(n for n in ast.parse(texts[filename]).body if isinstance(n,ast.ClassDef) and n.name==class_name)
        assert not any(isinstance(n,ast.FunctionDef) and n.name=='alloc' for n in cls.body)
        base=ast.unparse(cls.bases[0]);assert base==('arch.ArchitectureCircuit' if class_name=='ValueCircuit' else 'base.Circuit')
    circuit=next(n for n in ast.parse(texts['reversible.py']).body if isinstance(n,ast.ClassDef) and n.name=='Circuit')
    init=next(n for n in circuit.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    init_src=ast.get_source_segment(texts['reversible.py'],init)
    assert 'self._next_wire = 0' in init_src and 'self._free: list[int] = []' in init_src
    alloc=next(n for n in circuit.body if isinstance(n,ast.FunctionDef) and n.name=='alloc')
    alloc_src=ast.get_source_segment(texts['reversible.py'],alloc)
    assert 'w = self._next_wire' in alloc_src and 'self._next_wire += 1' in alloc_src and 'reg.append(w)' in alloc_src
    return {name:sha(repo/name) for name in texts}


def main(repo,results):
    checked_sources=source_allocation_checks(repo)
    data=json.loads((results/'chebyshev_sampler_templates.json').read_text())
    protocol=json.loads((results/'chebyshev_sampler_template_protocol.json').read_text())
    frozen=json.loads((results/'chebyshev_resources.json').read_text())
    assert data['protocol_sha256']==sha(results/'chebyshev_sampler_template_protocol.json')
    assert data['frozen_resource_sha256']==protocol['frozen_resource_sha256']==sha(results/'chebyshev_resources.json')
    assert data['source_sha256']==protocol['source_sha256']==sha(HERE/'emit_chebyshev_sampler_templates.py')
    for name,h in protocol['financial_compiler_source_sha256'].items():assert sha(repo/name)==h
    count=0
    for case,old in zip(data['cases'],frozen['cases']):
        assert case['model_id']==old['model_id']
        p=old['retained_financial_information']['financial_preparation_profile']
        assert case['financial_preparation_sha256']==digest(p)
        q=p['comparison_bits'];inputs=p['uniform_H']-q;n=inputs+q+1
        # Circuit starts with no free wires and consecutive allocation.
        original=list(range(inputs));comparator=list(range(inputs,inputs+q))
        highzero=[inputs+q];flag=inputs+q+1
        core=original+comparator+[flag];work=list(range(flag+1,flag+n-1))
        assert case['original_uniform_input_wires']==original
        assert case['comparator_wires']==comparator
        assert case['excluded_endpoint_high_zero_wires']==highzero
        assert case['success_flag_wire']==flag and n==p['reflection_core_qubits']
        assert case['concrete_invariant_core_wires']==core and case['clean_reflection_scratch_wires']==work
        assert not set(core)&set(work) and not set(highzero)&set(core+work)
        assert max(core+work)<p['logical_qubits_upper']==case['full_financial_allocation_upper']
        assert case['good_phase_primitive_gate']==['Z',flag]
        mapping=dict(enumerate(core+work))
        concrete=[[g[0]]+[mapping[w] for w in g[1:]] for g in old['uncontrolled_zero_reflection']['chronological_primitive_gate_list']]
        assert case['concrete_zero_reflection_primitive_gate_list']==concrete
        assert case['concrete_zero_reflection_tape_sha256']==digest(concrete)
        adj=Counter();swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
        for op,num in p['primitive_counts'].items():adj[swap.get(op,op)]+=num
        assert case['full_adjoint_phase_opcode_map']==swap
        assert case['A_forward_primitive_counts']==p['primitive_counts']
        assert clean(adj)==clean(case['A_adjoint_primitive_counts'])
        for degree,template in case['templates'].items():
            d=int(degree);k=d//2;s=old['sampler_profiles'][degree]
            body=[{'primitive_gate':['Z',flag]},{'retained_financial_operation':'A_adjoint'},
                {'concrete_zero_reflection_tape_sha256':digest(concrete)},{'retained_financial_operation':'A'}]
            expected=[{'retained_financial_operation':'A'},{'repeat':k if d%2 else k-1,'chronological_body':body}]
            if not d%2:expected.extend([{'primitive_gate':['Z',flag]},{'retained_financial_operation':'A_adjoint'}])
            assert template['chronological_instructions']==expected
            observed={'wire':flag,'value':1} if d%2 else {'wires':core,'all_values':0}
            assert template['final_observable']==observed
            assert template['same_sampler_profile_sha256']==digest(s)
            for key in ['A_forward_calls','A_adjoint_calls','zero_reflection_calls','good_Z_calls','primitive_counts','per_shot_metrics']:assert template[key]==s[key]
            assert template['full_wire_measurement_and_conditional_reset_count']==s['peak_logical_qubits']==p['logical_qubits_upper']
            f,a,z,g=(k+1,k,k,k) if d%2 else (k,k,k-1,k)
            c=Counter({op:f*num for op,num in p['primitive_counts'].items()})
            c.update({op:a*num for op,num in adj.items()})
            c.update({op:z*num for op,num in Counter(x[0] for x in concrete).items()});c['Z']+=g
            assert clean(c)==s['primitive_counts'];count+=1
    result={'schema':'qbp-chebyshev-concrete-template-audit-v1','status':'passed',
        'source_sha256':sha(__file__),'template_result_sha256':sha(results/'chebyshev_sampler_templates.json'),
        'template_protocol_sha256':sha(results/'chebyshev_sampler_template_protocol.json'),
        'frozen_resource_sha256':sha(results/'chebyshev_resources.json'),
        'verified_allocation_chain_source_sha256':checked_sources,
        'financial_cases':4,'concrete_sampler_templates':count,
        'success_flag_is_after_excluded_high_zero':True,'core_and_clean_scratch_disjoint':True,
        'all_phase_adjoint_opcode_counts_checked':True,'retained_peak_allocations_suffice':True,
        'executed_source_or_resource_bytes_modified':False,
        'scope':'Source-checked concrete financial flag/core binding and independently reconstructed sampler macros/counts. Correctness of retained financial preparations and exact phase reflections remains in the separately retained certificates; no full million-gate pricing tape was recompiled.'}
    (results/'chebyshev_sampler_template_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=HERE)
    parser.add_argument('--results',type=Path,default=HERE/'results');args=parser.parse_args()
    main(args.repo.resolve(),args.results.resolve())
