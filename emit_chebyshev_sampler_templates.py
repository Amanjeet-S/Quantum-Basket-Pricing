"""Concrete core and good-flag binding for the frozen sampler resource record.

The financial compiler excludes a comparator high-zero wire from its invariant
core.  This supplement resolves its compressed-core gate labels explicitly;
it changes no executed estimator, profile, diagnostic or financial answer.
"""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import argparse, hashlib, json, sys

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),default=str).encode()).hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def main(repo,results):
    sys.path.insert(0,str(repo))
    from value_oracle import ValueCircuit
    source=json.loads((results/'chebyshev_resources.json').read_text())
    protocol_path=results/'chebyshev_sampler_template_protocol.json'
    if protocol_path.exists():raise FileExistsError('Preserve the already executed supplement')
    dependencies=['integral_value_oracle.py','boundary_value_oracle.py','value_oracle.py',
                  'lookup_architectures.py','reversible_arithmetic.py']
    dependencies=[name for name in dependencies if (repo/name).exists()]
    protocol={'schema':'qbp-chebyshev-concrete-sampler-protocol-v1',
        'registered_time_utc':datetime.now(timezone.utc).isoformat(),
        'registered_before_template_emission':True,'source_sha256':sha(__file__),
        'frozen_resource_sha256':sha(results/'chebyshev_resources.json'),
        'financial_compiler_source_sha256':{name:sha(repo/name) for name in dependencies},
        'scope':'Wire binding supplement to frozen compressed-core resource templates; no repeat diagnostic, changed outcome, new price or cost correction.'}
    write(protocol_path,protocol)
    cases=[]
    for row in source['cases']:
        p=row['retained_financial_information']['financial_preparation_profile']
        q=p['comparison_bits'];input_bits=p['uniform_H']-q
        c=ValueCircuit(store_gates=False,mode='unitary8t')
        original_inputs=c.alloc(input_bits,'original uniform inputs')
        comparator=c.alloc(q,'uniform comparator')
        highzero=c.alloc(1,'comparator high zero')
        flag=c.alloc(1,'success flag')[0]
        core=tuple(original_inputs)+tuple(comparator)+(flag,)
        n=len(core);assert n==p['reflection_core_qubits'] and highzero[0] not in core
        assert flag==n and highzero[0]==n-1
        work=tuple(range(flag+1,flag+n-1))
        assert len(work)==n-2 and max(work)<p['logical_qubits_upper']
        abstract=row['uncontrolled_zero_reflection']['chronological_primitive_gate_list']
        mapping=dict(enumerate(core+work))
        concrete=[[gate[0]]+[mapping[w] for w in gate[1:]] for gate in abstract]
        assert Counter(g[0] for g in concrete)==Counter(row['uncontrolled_zero_reflection']['primitive_counts'])
        forward=p['primitive_counts'];adj={}
        swap={'T':'Tdg','Tdg':'T','S':'Sdg','Sdg':'S'}
        for op,count in forward.items():adj[swap.get(op,op)]=adj.get(swap.get(op,op),0)+count
        templates={}
        for degree,old in row['sampler_profiles'].items():
            d=int(degree);k=d//2
            qbody=[{'primitive_gate':['Z',flag]}, {'retained_financial_operation':'A_adjoint'},
                   {'concrete_zero_reflection_tape_sha256':digest(concrete)}, {'retained_financial_operation':'A'}]
            instructions=[{'retained_financial_operation':'A'},
                {'repeat':k if d%2 else k-1,'chronological_body':qbody}]
            if not d%2:instructions.extend([{'primitive_gate':['Z',flag]}, {'retained_financial_operation':'A_adjoint'}])
            templates[degree]={'degree':d,'chronological_instructions':instructions,
                'final_observable':{'wire':flag,'value':1} if d%2 else {'wires':list(core),'all_values':0},
                'full_wire_measurement_and_conditional_reset_count':p['logical_qubits_upper'],
                'same_sampler_profile_sha256':digest(old),
                'A_forward_calls':old['A_forward_calls'],'A_adjoint_calls':old['A_adjoint_calls'],
                'zero_reflection_calls':old['zero_reflection_calls'],'good_Z_calls':old['good_Z_calls'],
                'primitive_counts':old['primitive_counts'],'per_shot_metrics':old['per_shot_metrics']}
        cases.append({'model_id':row['model_id'],'financial_preparation_sha256':digest(p),
            'original_uniform_input_wires':list(original_inputs),'comparator_wires':list(comparator),
            'excluded_endpoint_high_zero_wires':list(highzero),'success_flag_wire':flag,
            'concrete_invariant_core_wires':list(core),'clean_reflection_scratch_wires':list(work),
            'full_financial_allocation_upper':p['logical_qubits_upper'],
            'wire_binding_derivation':'compile_flag allocates original graph inputs, q comparator wires, one high-zero wire, then the success flag before financial-body allocation; ValueCircuit.alloc supplies those concrete initial indices. The high-zero wire is omitted from reflection_core.',
            'scratch_invariant':'Financial body and comparator work are clean after complete A. Reflection scratch uses this clean allocated private pool and is returned clean; these indices do not enter the core projector.',
            'global_minus_scope':'The unconditioned minus in Q contributes an irrelevant global phase. Conditional minus remains charged in the separate canonical QPE fallback.',
            'good_phase_primitive_gate':['Z',flag],
            'concrete_zero_reflection_primitive_gate_list':concrete,
            'concrete_zero_reflection_tape_sha256':digest(concrete),
            'A_forward_primitive_counts':forward,'A_adjoint_primitive_counts':adj,
            'full_adjoint_phase_opcode_map':swap,'templates':templates})
    output={'schema':'qbp-chebyshev-concrete-sampler-template-v1','status':'completed',
        'protocol_sha256':sha(protocol_path),'frozen_resource_sha256':sha(results/'chebyshev_resources.json'),
        'source_sha256':sha(__file__),'cases':cases,
        'scope':protocol['scope'],'compiler_version_scope':'Actual financial compiler allocation convention and retained full preparation, not a newly recompiled joint financial tape.'}
    write(results/'chebyshev_sampler_templates.json',output)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=HERE)
    parser.add_argument('--results',type=Path,default=HERE/'results');args=parser.parse_args()
    main(args.repo.resolve(),args.results.resolve())
