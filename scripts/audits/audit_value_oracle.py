"""Independent finite-value, bucket, phase and resource-composition checks.

State-vector arithmetic is a numerical diagnostic of emitted logical gates.
Finite Boolean replay and the algebraic corrected-instrument proof establish
the exact logical maps; sampled or floating state vectors are not substituted
for those arguments.
"""
from itertools import product
from collections import Counter
from pathlib import Path
import hashlib
import json
import math
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE if (HERE/'value_oracle.py').exists() else HERE.parent
sys.path.insert(0,str(ROOT))
import numpy as np
import architecture as arch
import encoding_resources as er
import generalisation as gen
import quantum_gates as qg
import value_oracle as vo


def _apply_primitive(state,operations,outcomes=None,normalise_branches=False):
    """Apply real primitive unitaries and specified measurement Kraus maps."""
    state=np.asarray(state,complex).copy();rows=np.arange(len(state));outcomes=outcomes or {}
    for op,wires,condition in operations:
        if op in ('CZ_IF','X_IF') and not outcomes[condition]:continue
        if op in ('X','X_IF','CX','CCX'):
            target=wires[-1];mask=np.ones(len(rows),bool)
            for control in wires[:-1]:mask&=(rows>>control&1).astype(bool)
            perm=rows^(mask.astype(np.int64)<<target);state=state[perm]
        elif op=='H':
            bit=1<<wires[0];low=rows[(rows&bit)==0];high=low|bit
            a,b=state[low].copy(),state[high].copy()
            state[low]=(a+b)/math.sqrt(2);state[high]=(a-b)/math.sqrt(2)
        elif op in ('T','TD','S','SD','Z','CZ','CZ_IF'):
            if op in ('CZ','CZ_IF'):
                active=np.ones(len(rows),bool)
                for w in wires:active&=(rows>>w&1).astype(bool)
                phase=-1
            else:
                active=(rows>>wires[0]&1).astype(bool)
                phase={'T':np.exp(1j*math.pi/4),'TD':np.exp(-1j*math.pi/4),'S':1j,'SD':-1j,'Z':-1}[op]
            state[active]*=phase
        elif op=='MZ':
            state[(rows>>wires[0]&1)!=outcomes[condition]]=0
            if normalise_branches:state*=math.sqrt(2)
        else:raise ValueError(op)
    return state


def apply_logical(gates,state,branch_bits=None,normalise_branches=False):
    """Replay logical tape with fresh fully expanded lookup instruments."""
    bits=iter(branch_bits or []);state=np.asarray(state,complex).copy()
    for gate in gates:
        if gate[0]=='LOOKUP':
            plan=gate[1];m=plan.summary()['measurement_count']
            outcomes={j:next(bits,0) for j in range(m)}
            state=_apply_primitive(state,plan.operations(),outcomes,normalise_branches)
        else:
            op={'Tdg':'TD','Sdg':'SD'}.get(gate[0],gate[0])
            state=_apply_primitive(state,[(op,tuple(gate[1:]),None)],normalise_branches=normalise_branches)
    return state


def toy_graph():
    g=er.Graph(2);g.width=3
    g.output=g.op('quantile',g.inp('x',2),bits=2)
    return g,{('quantile',2,2):[0,1,3,4]}


def predicate_audit():
    reports=[]
    for q in range(1,6):
        for b in range(q+1):
            c=vo.ValueCircuit(mode='exact7t');a=c.alloc(q+1)
            u=c.alloc(b+1 if b<q else 0);flag=c.alloc(1)[0]
            vo.bucket_predicate_into(c,a,u,b,flag);metrics=c.resources()
            count=0
            # Complete valid numerical-value domain, all comparison and dirty
            # flag states. Work is independent zero input, not dirty scratch.
            for value,U,dirty in product(range((1<<q)+1),range(1<<len(u)),(0,1)):
                good=value==(1<<q) if b==q else (1<<b)<=value<(1<<(b+1)) and U<value
                state=c.run({a:value,u:U,(flag,):dirty})
                assert state[flag]==dirty^bool(good)
                assert c.value(state,a)==value and c.value(state,u)==U
                assert all(x==0 for w,x in enumerate(state) if w not in a+u+(flag,))
                count+=1
            if b<q:
                assert metrics['toffoli']==2*q+4*b+7
                assert metrics['cnot']==8*(b+1)
                assert metrics['x']==2*(q+1)
                assert metrics['wires']==2*q+3*b+7
            else:
                assert metrics['toffoli']==metrics['x']==0 and metrics['cnot']==1
            reports.append({'q':q,'bucket':b,'valid_domain_complete_basis_checks':count,
                'predicate_Toffoli':metrics['toffoli'],'predicate_CX':metrics['cnot'],
                'predicate_X':metrics['x'],'all_private_work_zero':True})
    return reports


def exact_CCX_expansion_audit():
    identity=np.eye(8,dtype=complex)
    expected=identity.copy();expected[[3,7]]=expected[[7,3]]
    expanded=qg.expand_toffoli([('CCX',0,1,2)])
    actual=apply_logical(expanded,identity)
    maximum=float(np.max(abs(actual-expected)));assert maximum<5e-15
    metric=qg.resources(expanded,3)
    assert metric['gate_counts']=={'H':2,'CX':6,'Tdg':3,'T':4}
    assert metric['T_depth_assigned_3layer_Toffoli_or_explicit_T_upper']==4
    return {'complete_complex_matrix_maximum_error':maximum,'logical_basis_columns':8,
        'literal_expanded_gate_list':expanded,'seven_T_count':7,'literal_T_dependency_depth':4,
        'conservative_serial_T_depth_upper':5,'historical_assigned_depth':3}


def range_predicate_audit():
    reports=[]
    for q in range(1,5):
        for bits in range(q+1):
            for clipped in (False,True):
                c=vo.ValueCircuit(mode='exact7t');value=c.alloc(q+1)
                u=c.alloc(bits);flag=c.alloc(1)[0]
                function=vo.clipped_range_predicate_into if clipped else vo.range_predicate_into
                function(c,value,u,bits,flag);metric=c.resources();checks=0
                maximum=(1<<q) if clipped else (1<<bits)
                for a,U,dirty in product(range(maximum+1),range(1<<bits),(0,1)):
                    state=c.run({value:a,u:U,(flag,):dirty})
                    assert state[flag]==dirty^bool(U<a)
                    assert all(x==0 for w,x in enumerate(state) if w not in value+u+(flag,))
                    assert c.value(state,value)==a and c.value(state,u)==U;checks+=1
                expected_toffoli=6*(q+1) if clipped else 6*(bits+1) if bits else 0
                assert metric['toffoli']==expected_toffoli
                reports.append({'q':q,'range_bits':bits,'clipped':clipped,
                    'complete_valid_domain_basis_checks':checks,'predicate_Toffoli':metric['toffoli'],
                    'all_private_work_zero':True})
    return reports


def finite_value_audit():
    g,tables=toy_graph();reports=[]
    for mode in arch.MODES:
        c,inputs,z=vo.compile_value(g,tables,mode);metric=c.resources()
        composed=vo.counted_value_composition(g,tables,{},mode)
        assert metric['primitive_gate_counts']==composed['expanded_gate_counts']
        assert metric['clifford_T_expanded_gate_counts']==composed['clifford_T_expanded_gate_counts']
        assert composed['logical_qubits_serial_upper_bound']>=c.wire_count
        for x,dirty in product(range(4),range(8)):
            state=c.run({inputs['x']:x,z:dirty})
            assert c.value(state,z)==dirty^g.evaluate({'x':x},tables)
            assert c.value(state,inputs['x'])==x
            assert all(bit==0 for w,bit in enumerate(state) if w not in inputs['x']+z)
        # Every corrected measurement record is checked on a coherent complex
        # superposition of all logical dirty-output states. Extension to an
        # external reference follows from the algebraic instrument proof and
        # the separate existing complete lookup Kraus-matrix audit.
        rng=np.random.default_rng(81733);logical=inputs['x']+z
        state=np.zeros(1<<c.wire_count,complex)
        amp=rng.normal(size=1<<len(logical))+1j*rng.normal(size=1<<len(logical));amp/=np.linalg.norm(amp)
        expected=np.zeros_like(state)
        for j,a in enumerate(amp):
            index=sum(((j>>k)&1)<<w for k,w in enumerate(logical));state[index]=a
            expected[index^(tables[('quantile',2,2)][j&3]<<z[0])]=a
        m=metric['measurement_count'];maximum=0
        for branch in product((0,1),repeat=m):
            actual=apply_logical(list(vo.logical_gates(c)),state,branch)
            error=float(np.max(abs(actual-expected*2**(-m/2))));maximum=max(maximum,error)
            assert error<5e-13
        reports.append({'mode':mode,'dirty_value_basis_checks':32,'all_corrected_records_checked':1<<m,
            'measurement_bits':m,'coherent_complex_branch_error':maximum,
            'source_scope':'all actual primitive lookup gates; arbitrary coherent logical dirty outputs; private work clean',
            'complete_vs_composed_counts_match':True,'metrics':metric})
    return reports


def retained_block_adapter_audit():
    graph,tables=toy_graph();reports=[]
    for mode in arch.MODES:
        historical=arch.counted_composition(graph,tables,{},mode)
        body=vo.prepare_existing_body(historical,graph.q)
        checks=[]
        for b in range(graph.q+1):
            p=vo.profile_from_prepared_body(body,b)
            fresh=vo.counted_bucket_composition(graph,tables,{},b,mode)
            assert p['expanded_gate_counts']==fresh['expanded_gate_counts']
            assert p['clifford_T_expanded_gate_counts']==fresh['clifford_T_expanded_gate_counts']
            assert p['grover_iterate_T_count']==fresh['grover_iterate_T_count']
            c,inputs,u,f=vo.compile_bucket(graph,tables,b,mode)
            assert p['logical_qubits_serial_upper_bound']>=c.wire_count
            checks.append({'kind':'bucket','coordinate':b,'complete_counts_match_fresh_emission':True})
        for bits in range(graph.q+1):
            for clipped in (False,True):
                p=(vo.profile_from_prepared_clipped_range if clipped else vo.profile_from_prepared_range)(body,bits)
                c,inputs,u,f=(vo.compile_clipped_range if clipped else vo.compile_range)(graph,tables,bits,mode)
                metrics=c.resources()
                assert p['expanded_gate_counts']==metrics['primitive_gate_counts']
                assert p['clifford_T_expanded_gate_counts']==metrics['clifford_T_expanded_gate_counts']
                assert p['logical_qubits_serial_upper_bound']>=c.wire_count
                assert p['explicit_sequence_T_depth_serial_upper']>=metrics['explicit_7T_sequence_T_depth']
                assert p['conservative_5layer_CCX_T_depth_serial_upper']>=metrics['conservative_5layer_CCX_T_depth_upper']
                if clipped:
                    assert p['expanded_gate_counts']==historical['expanded_gate_counts']
                    good=sum(c.run({inputs['x']:x,u:U})[f] for x,U in product(range(4),range(1<<bits)))
                    assert good==sum(min(a,1<<bits) for a in tables[('quantile',2,2)])
                checks.append({'kind':'clipped_range' if clipped else 'certified_range','coordinate':bits,
                    'complete_counts_match_fresh_emission':True})
        reports.append({'mode':mode,'adapter_checks':checks})
    return reports


def bucket_grover_audit():
    g,tables=toy_graph();reports=[]
    for mode in arch.MODES:
        for b in range(g.q+1):
            c,inputs,u,flag=vo.compile_bucket(g,tables,b,mode);metric=c.resources()
            composed=vo.counted_bucket_composition(g,tables,{},b,mode)
            assert metric['primitive_gate_counts']==composed['expanded_gate_counts']
            assert metric['clifford_T_expanded_gate_counts']==composed['clifford_T_expanded_gate_counts']
            assert composed['logical_qubits_serial_upper_bound']>=c.wire_count
            assert composed['explicit_sequence_T_depth_serial_upper']>=metric['explicit_7T_sequence_T_depth']
            assert composed['conservative_5layer_CCX_T_depth_serial_upper']>=metric['conservative_5layer_CCX_T_depth_upper']
            core=inputs['x']+u+(flag,);work=tuple(w for w in range(c.wire_count) if w not in core)[:len(core)-2]
            oracle=list(vo.logical_gates(c));A=[('H',w) for w in inputs['x']+u]+oracle
            Q=qg.grover(A,core,flag,work)
            state=np.zeros(1<<c.wire_count,complex);state[0]=1
            state=apply_logical(A,state,normalise_branches=True)
            mean=sum((a/(1<<(b+1))) if (1<<b)<=a<(1<<(b+1)) else 0
                     for a in tables[('quantile',2,2)])/4 if b<g.q else 1/4
            actual=float(sum(abs(z)**2 for i,z in enumerate(state) if i>>flag&1))
            assert abs(actual-mean)<3e-12
            theta=math.asin(math.sqrt(mean));errors=[]
            for k in range(4):
                actual=float(sum(abs(z)**2 for i,z in enumerate(state) if i>>flag&1))
                expected=math.sin((2*k+1)*theta)**2
                errors.append(abs(actual-expected));assert errors[-1]<3e-11
                if k<3:state=apply_logical(Q,state,normalise_branches=True)
            # The controlled-Q operator has the correct conditional global
            # sign: control=0 is identity, control=1 is minus the unconstrained
            # global-sign convention of qg.grover's un-controlled output.
            control=c.wire_count;physical=c.wire_count+1
            controlled_work=tuple(w for w in range(c.wire_count) if w not in core)[:len(core)-1]
            controlled=qg.grover(A,core,flag,controlled_work,control)
            prepared=np.zeros(1<<c.wire_count,complex);prepared[0]=1
            prepared=apply_logical(A,prepared,normalise_branches=True)
            uncontrol=apply_logical(Q,prepared,normalise_branches=True)
            control_errors=[]
            for z in (0,1):
                big=np.zeros(1<<physical,complex);big[z*(1<<control):(z+1)*(1<<control)]=prepared
                actual_big=apply_logical(controlled,big,normalise_branches=True)
                expected_big=np.zeros_like(big)
                expected_big[z*(1<<control):(z+1)*(1<<control)]=prepared if z==0 else -uncontrol
                error=float(np.max(abs(actual_big-expected_big)));control_errors.append(error);assert error<3e-11
            # Direct dyadic reconstruction reproduces the unflagged mean.
            basis=0
            for x,U,dirty in product(range(4),range(1<<len(u)),(0,1)):
                a=tables[('quantile',2,2)][x]
                good=a==4 if b==2 else (1<<b)<=a<(1<<(b+1)) and U<a
                state_basis=c.run({inputs['x']:x,u:U,(flag,):dirty})
                assert state_basis[flag]==dirty^bool(good)
                assert all(z==0 for w,z in enumerate(state_basis) if w not in core);basis+=1
            reports.append({'mode':mode,'bucket':b,'complete_good_predicate_basis_checks':basis,
                'Born_probability':mean,'maximum_Grover_probability_error':max(errors),
                'Grover_powers_checked':[0,1,2,3],'controlled_Q_zero_one_errors':control_errors,
                'complete_vs_composed_counts_match':True,'statevector_qubits':c.wire_count,
                'instrument_branch_scope':'measured mode uses fresh zero records normalised by their proven input-independent factors; every record of the value oracle is separately checked',
                'profile':composed})
    reconstructed=sum((2**(r['bucket']+1-g.q) if r['bucket']<g.q else 1)*r['Born_probability']
                      for r in reports if r['mode']=='exact7t')
    assert reconstructed==sum(tables[('quantile',2,2)])/(4*(1<<g.q))
    return reports


def financial_audit(protocol):
    model=next(m for m in protocol['models'] if m['assets']==2 and m['lambda']=='0.2')
    reports=[]
    for mode in arch.MODES:
        for method in ('direct','interval','full_cap'):
            graph=gen.build_graph(model,method,9,1,1,1,3,4,4,10);tables={};gen.make_tables(graph,tables)
            cache={};c,inputs,z=vo.compile_value(graph,tables,mode)
            historical=arch.counted_composition(graph,tables,{},mode)
            prepared_body=vo.prepare_existing_body(historical,graph.q)
            cm=vo.counted_value_composition(graph,tables,cache,mode);metric=c.resources()
            assert metric['primitive_gate_counts']==cm['expanded_gate_counts']
            assert metric['clifford_T_expanded_gate_counts']==cm['clifford_T_expanded_gate_counts']
            assert cm['logical_qubits_serial_upper_bound']>=c.wire_count
            mappings=[dict(zip(inputs,values)) for values in product((0,1),repeat=len(inputs))]
            valuechecks=0
            for mapping in mappings:
                a=graph.evaluate(mapping,tables);assert 0<=a<=1<<graph.q
                for dirty in (0,1,(1<<(graph.q+1))-1):
                    initial={inputs[name]:value for name,value in mapping.items()};initial[z]=dirty
                    state=c.run(initial);assert c.value(state,z)==dirty^a
                    keep=set(z).union(*(set(reg) for reg in inputs.values()))
                    assert all(v==0 for w,v in enumerate(state) if w not in keep);valuechecks+=1
            rows=[]
            for b in (0,4,graph.q):
                c,inputs,u,flag=vo.compile_bucket(graph,tables,b,mode)
                cm=vo.counted_bucket_composition(graph,tables,cache,b,mode);metric=c.resources()
                adapted=vo.profile_from_prepared_body(prepared_body,b)
                assert adapted['expanded_gate_counts']==cm['expanded_gate_counts']
                assert adapted['clifford_T_expanded_gate_counts']==cm['clifford_T_expanded_gate_counts']
                assert metric['primitive_gate_counts']==cm['expanded_gate_counts']
                assert metric['clifford_T_expanded_gate_counts']==cm['clifford_T_expanded_gate_counts']
                assert cm['logical_qubits_serial_upper_bound']>=c.wire_count
                assert cm['explicit_sequence_T_depth_serial_upper']>=metric['explicit_7T_sequence_T_depth']
                checked=0
                for mapping in mappings:
                    a=graph.evaluate(mapping,tables)
                    levels=sorted({0,(1<<len(u))-1,min(a,(1<<len(u))-1),min(max(a-1,0),(1<<len(u))-1)})
                    for U,dirty in product(levels,(0,1)):
                        initial={inputs[name]:value for name,value in mapping.items()};initial[u]=U;initial[(flag,)]=dirty
                        state=c.run(initial)
                        good=a==(1<<graph.q) if b==graph.q else (1<<b)<=a<(1<<(b+1)) and U<a
                        assert state[flag]==dirty^bool(good)
                        keep={flag,*u,*[w for reg in inputs.values() for w in reg]}
                        assert all(v==0 for w,v in enumerate(state) if w not in keep)
                        checked+=1
                rows.append({'bucket':b,'boundary_basis_checks':checked,'all_work_clean':True,
                    'complete_vs_composed_counts_match':True,'T_count':cm['logical_T_count'],
                    'measurements':cm['measurement_count']})
            clip_checks=[]
            for bits in (0,4,graph.q):
                c,inputs,u,flag=vo.compile_clipped_range(graph,tables,bits,mode)
                adapted=vo.profile_from_prepared_clipped_range(prepared_body,bits);metric=c.resources()
                assert adapted['expanded_gate_counts']==metric['primitive_gate_counts']
                assert adapted['clifford_T_expanded_gate_counts']==metric['clifford_T_expanded_gate_counts']
                assert adapted['logical_qubits_serial_upper_bound']>=c.wire_count
                checked=0
                for mapping in mappings:
                    a=graph.evaluate(mapping,tables)
                    levels=sorted({0,(1<<bits)-1,min(a,(1<<bits)-1),min(max(a-1,0),(1<<bits)-1)})
                    for U,dirty in product(levels,(0,1)):
                        initial={inputs[name]:value for name,value in mapping.items()};initial[u]=U;initial[(flag,)]=dirty
                        state=c.run(initial);assert state[flag]==dirty^bool(U<a)
                        keep={flag,*u,*[w for reg in inputs.values() for w in reg]}
                        assert all(v==0 for w,v in enumerate(state) if w not in keep);checked+=1
                clip_checks.append({'range_bits':bits,'boundary_basis_checks':checked,
                    'complete_vs_retained_adapter_counts_match':True,'all_work_clean':True})
            reports.append({'mode':mode,'method':method,'model_id':model['model_id'],
                'q':graph.q,'value_dirty_basis_checks':valuechecks,'buckets':rows,
                'clipped_range_checks':clip_checks,
                'scope':'coarse primitive composition and finite semantics; no production price-accuracy assertion'})
    return reports


def main():
    protocol_path=ROOT/'results/generalisation_protocol.json'
    out={'status':'exact value XOR, finite dyadic predicates, emitted coherent gates and full block composition passed',
        'source_hashes':{name:hashlib.sha256((Path(__file__) if name=='audit_value_oracle.py' else ROOT/name).read_bytes()).hexdigest()
                         for name in ('value_oracle.py','audit_value_oracle.py','architecture.py','reversible.py','encoding_resources.py','quantum_gates.py','generalisation.py')},
        'exact_CCX_expansion_audit':exact_CCX_expansion_audit(),
        'predicate_audits':predicate_audit(),'range_predicate_audits':range_predicate_audit(),
        'value_instrument_audits':finite_value_audit(),
        'bucket_statevector_audits':bucket_grover_audit(),
        'retained_block_adapter_audits':retained_block_adapter_audit()}
    if protocol_path.exists():out['financial_audits']=financial_audit(json.loads(protocol_path.read_text()))
    out['summary']={'predicate_basis_checks':sum(x['valid_domain_complete_basis_checks'] for x in out['predicate_audits']),
        'corrected_value_instrument_records':sum(x['all_corrected_records_checked'] for x in out['value_instrument_audits']),
        'Grover_bucket_configurations':len(out['bucket_statevector_audits']),
        'complete_financial_graphs':len(out.get('financial_audits',[]))}
    (ROOT/'results').mkdir(exist_ok=True)
    output=ROOT/'results/value_oracle_audit.json';output.write_text(json.dumps(out,indent=2)+'\n');print(output)


if __name__=='__main__':main()
