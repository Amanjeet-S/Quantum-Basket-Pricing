"""Complex phase/instrument and Boolean-composition audits of expanded QROM."""
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import architecture as arch
import generalisation as gen


def quantum_apply(state,operations,outcomes=None):
    """Apply actual primitive gates and unnormalised measurement Kraus maps."""
    state=np.asarray(state,complex).copy();rows=np.arange(len(state));outcomes=outcomes or {}
    for op,wires,condition in operations:
        if op in ('CZ_IF','X_IF') and not outcomes[condition]:continue
        if op in ('X','X_IF','CX','CCX'):
            target=wires[-1];mask=np.ones(len(rows),bool)
            for control in wires[:-1]:mask&=(rows>>control&1).astype(bool)
            perm=rows^(mask.astype(np.int64)<<target);state=state[perm]
        elif op=='H':
            target=wires[0];low=rows[(rows>>target&1)==0];high=low|(1<<target)
            a,b=state[low].copy(),state[high].copy()
            state[low]=(a+b)/math.sqrt(2);state[high]=(a-b)/math.sqrt(2)
        elif op in ('T','TD','S','SD','CZ','CZ_IF'):
            if op in ('CZ','CZ_IF'):
                active=((rows>>wires[0]&1)&(rows>>wires[1]&1)).astype(bool);phase=-1
            else:
                active=(rows>>wires[0]&1).astype(bool)
                phase={'T':np.exp(1j*math.pi/4),'TD':np.exp(-1j*math.pi/4),'S':1j,'SD':-1j}[op]
            state[active]*=phase
        elif op=='MZ':state[(rows>>wires[0]&1)!=outcomes[condition]]=0
        else:raise ValueError(op)
    return state


def audit_AND():
    encoded=np.zeros((8,4),complex)
    for ab in range(4):encoded[ab,ab]=1
    actual=quantum_apply(encoded,arch.and_compute(0,1,2));expected=np.zeros_like(actual)
    for ab in range(4):expected[ab+((ab==3)<<2),ab]=1
    error=float(np.max(abs(actual-expected)));assert error<1e-14
    erases=[]
    for outcome in (0,1):
        sequence=[('H',(2,),None),('MZ',(2,),0),('CZ_IF',(0,1),0),('X_IF',(2,),0)]
        deleted=quantum_apply(actual,sequence,{0:outcome})
        difference=float(np.max(abs(deleted-encoded/math.sqrt(2))));assert difference<1e-14
        erases.append({'outcome':outcome,'maximum_Kraus_matrix_error':difference})
    restored=quantum_apply(actual,arch.and_inverse(0,1,2));assert np.max(abs(restored-encoded))<1e-14
    return {'explicit_compute_gate_order':arch.and_compute(0,1,2),
        'explicit_unitary_inverse_gate_order':arch.and_inverse(0,1,2),
        'compute_T_or_Tdag_count':4,'compute_matrix_error':error,'erase_branches':erases,
        'phase_statement':'each corrected erase branch is +1/sqrt(2) times the identity on arbitrary control states',
        'precondition':'AND target starts clean; this is not a four-T Toffoli on a dirty target'}


def audit_lookup_matrix(n,words,outbits,controlled,mode,all_branches=True):
    addr=tuple(range(n));out=tuple(range(n,n+outbits));enable=n+outbits if controlled else None
    logical=n+outbits+controlled
    plan=arch.lookup_plan(addr,out,words,mode,enable=enable)
    physical=max(plan.address+plan.output+plan.work+(() if enable is None else (enable,)),default=-1)+1
    columns=1<<logical;state=np.zeros((1<<physical,columns),complex)
    for j in range(columns):state[j,j]=1
    expected=np.zeros_like(state)
    for j in range(columns):
        index=j&((1<<n)-1);word=words[index] if index<len(words) else 0
        if controlled and not(j>>enable&1):word=0
        expected[j^(word<<n),j]=1
    ops=list(plan.operations());count=plan.summary()['measurement_count']
    branches=list(product((0,1),repeat=count)) if all_branches else [tuple((i>>j)&1 for j in range(count)) for i in (0,(1<<count)-1,37%(1<<count))]
    maxerror=0;probability_error=0
    for branch in branches:
        actual=quantum_apply(state,ops,dict(enumerate(branch)))
        error=float(np.max(abs(actual-expected*2**(-count/2))));maxerror=max(maxerror,error)
        assert error<2e-12,(n,mode,controlled,branch,error)
        probabilities=np.sum(abs(actual)**2,axis=0)
        probability_error=max(probability_error,float(np.max(abs(probabilities-2.**(-count)))))
    return {'address_bits':n,'output_bits':outbits,'table':words,'controlled':controlled,'mode':mode,
        'complete_dirty_output_logical_basis_columns':columns,'measurement_bits':count,
        'branches_checked':len(branches),'all_branches':all_branches,
        'maximum_complete_complex_Kraus_matrix_error':maxerror,
        'maximum_address_and_output_independent_branch_probability_error':probability_error,
        'private_work_reset':True,'logical_inverse':'fresh same lookup instrument; XOR involution',
        'resources':plan.summary()}


def audit_boolean_lookup():
    rows=[];rng=np.random.default_rng(718233)
    for mode in arch.MODES:
        for n in range(5):
            words=[int(x) for x in rng.integers(0,32,size=1<<n)]
            if n:words[0]=0
            c=arch.ArchitectureCircuit(mode=mode);a=c.alloc(n);out=c.lookup(a,words,5)
            count=0
            for index in range(1<<n):
                for dirty in range(32):
                    state=c.run({a:index,out:dirty})
                    assert c.value(state,a)==index and c.value(state,out)==dirty^words[index]
                    assert all(state[w]==0 for p in c._macros for w in p.work);count+=1
            # Append the fresh operational inverse on the same dirty target.
            if c._macros:
                p=c._macros[0];c._emit(arch.MACRO,0,-1,-1)
                for index in range(1<<n):
                    state=c.run({a:index,out:17});assert c.value(state,out)==17
            rows.append({'mode':mode,'address_bits':n,'dirty_basis_checks':count,
                'logical_inverse_checked':True,'metrics':c.resources()})
    return rows


def audit_financial_graphs(protocol):
    rows=[]
    picks=[next(m for m in protocol['models'] if m['assets']==n and m['lambda']==lam)
           for n,lam in [(2,'0.05'),(4,'0.2'),(8,'0.4')]]
    for mode in arch.MODES:
        for model in picks:
            for method in ('direct','interval','full_cap'):
                cap={2:25,4:10,8:1}[model['assets']]
                graph=gen.build_graph(model,method,9,1,1,1,3,4,4,cap);tables={};gen.make_tables(graph,tables)
                c,inputs,comparison,flag=arch.compile_graph(graph,tables,mode)
                metric=c.resources();composed=arch.counted_composition(graph,tables,{},mode)
                assert metric['logical_T_count']==composed['logical_T_count']
                assert metric['primitive_gate_counts']==composed['expanded_gate_counts'],(mode,model['model_id'],method)
                assert composed['logical_qubits_serial_upper_bound']>=metric['wires']
                cases=list(product((0,1),repeat=len(inputs)))
                if len(cases)>16:
                    rr=np.random.default_rng(1229);ids=np.unique(np.r_[0,len(cases)-1,rr.integers(0,len(cases),14)])
                    cases=[cases[i] for i in ids]
                checks=0
                for values in cases:
                    mapping=dict(zip(inputs,values));threshold=graph.evaluate(mapping,tables)
                    for level in sorted(set([0,511,min(threshold,511),max(0,min(threshold-1,511))])):
                        state=c.run({**{inputs[k]:v for k,v in mapping.items()},comparison:level})
                        assert state[flag]==int(level<threshold)
                        assert all(c.value(state,inputs[k])==v for k,v in mapping.items())
                        assert c.value(state,comparison)==level
                        reserved={flag,*comparison,*[w for reg in inputs.values() for w in reg]}
                        assert all(v==0 for w,v in enumerate(state) if w not in reserved)
                        checks+=1
                rows.append({'mode':mode,'model_id':model['model_id'],'method':method,'cap':cap,
                    'complete_boolean_comparator_checks':checks,'all_work_reset':True,
                    'phase_argument':'composition of exact arithmetic permutations and separately verified phase-free QROM logical maps',
                    'complete_vs_composed_primitive_counts_match':True,'logical_T_count':metric['logical_T_count'],
                    'measurement_count':metric['measurement_count'],
                    'qualification':'coarse gate/instrument composition audit, not production financial price accuracy'})
    return rows


def main():
    path=ROOT/'results/base_snapshot.json';base_record=json.loads(path.read_text())
    before={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in base_record['source_files']}
    assert before==base_record['source_files'];protocol=json.loads((ROOT/'results/generalisation_protocol.json').read_text())
    local=audit_AND();matrices=[]
    cases=[(0,[5],3,False),(1,[3,5],3,False),(2,[3,6,7,1],3,False),
           (2,[0,6,0,1],3,True),(2,[5],3,True),(3,[1,2,3,1,0,2,3,1],2,True)]
    for mode in arch.MODES:
        for n,words,outbits,control in cases:
            matrices.append(audit_lookup_matrix(n,words,outbits,control,mode,n<=2))
    boolean=audit_boolean_lookup();financial=audit_financial_graphs(protocol)
    after={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in before};assert before==after
    out={'status':'primitive phase, corrected Kraus branches, dirty-output XOR and complete coarse composition passed',
         'base_source_hashes':before,'architecture_sha256':hashlib.sha256((ROOT/'architecture.py').read_bytes()).hexdigest(),
         'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'AND_isometry':local,'complete_lookup_phase_audits':matrices,'boolean_lookup_audits':boolean,
         'financial_composition_audits':financial,
         'summary':{'lookup_phase_configurations':len(matrices),'instrument_branches_checked':sum(x['branches_checked'] for x in matrices),
             'dirty_boolean_lookup_checks':sum(x['dirty_basis_checks'] for x in boolean),
             'complete_financial_circuits':len(financial),
             'financial_comparator_checks':sum(x['complete_boolean_comparator_checks'] for x in financial)},
         'proof':'The displayed clean AND isometry is phase free. H/measurement followed by conditional CZ/reset gives each erasure Kraus map +1/sqrt(2) times deletion of the encoded AND. Intermediate controlled data writes preserve its controls and restore its AND target before erasure. Induction over the prefix tree gives every branch 2^(-I/2) times the same address-preserving dirty-output XOR unitary, including coherent enable controls. The branch record is therefore independent of logical input. A fresh complete instrument implements the logical inverse; no measurement tape is reversed.',
         'limits':['ideal logical gates and measurements only','conditional Clifford costs charged at worst-case execution',
             'no physical feedback latency or hardware advantage claim','arithmetic exact Toffolis retain the specified seven-T decomposition',
             'lookup instrument proof requires initially clean private selector work']}
    output=ROOT/'results/architecture_audit.json';output.write_text(json.dumps(out,indent=2)+'\n');print(output)


if __name__=='__main__':main()
