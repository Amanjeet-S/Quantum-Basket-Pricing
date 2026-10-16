"""Whole finite financial basis and literal preparation audit for boundaries."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).resolve().parent))
from boundary_audit_support import interpret
from integral_audit_support import packed_replay,register_values,independent_gate_accounting


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def normalise(counts):return Counter({{'TD':'Tdg','SD':'Sdg'}.get(g,g):v for g,v in counts.items()})


def audit(root):
    root=Path(root);sys.path.insert(0,str(root));import boundary_value_oracle as o
    import integral_value_oracle as retained
    models=json.loads((root/'results/generalisation_protocol.json').read_text())['models'];model=next(m for m in models if m['model_id']=='n2_lambda0.2')
    cache={};count_cache={};rows=[];value_lanes=flag_lanes=0
    for distribution in ('quantile','uniform_box'):
        for payoff,cap in [('capped',25),('vanilla',None)]:
            for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
                for encoding in (('centred','positive_part','negative_part') if payoff=='capped' and kind!='positive_r1' else ('centred',)):
                    graph=o.build_graph(model,cap,q=5,b_y=1,beta_bits=1,table_bits=2,outer_bound=2,U_star=F(1,8),
                        kind=kind,payoff=payoff,distribution=distribution,output_encoding=encoding)
                    tables=o.make_tables(graph,cache)
                    names=tuple(graph.inputs);samples=[dict(zip(names,values)) for values in product(*(range(1<<graph.inputs[n]) for n in names))]
                    outputs=[interpret(graph,s,tables)[0] for s in samples];assert all(0<=a<=1<<graph.q for a in outputs)
                    for mode in ('exact7t','unitary8t'):
                        c,inputs,out=o.compile_value(graph,tables,mode);assignments=[];expected=[]
                        for sample,a in zip(samples,outputs):
                            for dirty in range(1<<len(out)):
                                state={inputs[n]:sample[n] for n in names};state[out]=dirty;assignments.append(state);expected.append(dirty^a)
                        state,_=packed_replay(c,assignments);assert register_values(state,out,len(assignments))==expected
                        retained_wires=set(out+tuple(w for r in inputs.values() for w in r))
                        assert all(word==0 for wire,word in enumerate(state) if wire not in retained_wires)
                        for n in names:assert register_values(state,inputs[n],len(assignments))==[s[inputs[n]] for s in assignments]
                        reverse,_=packed_replay(c,assignments,reverse=True,initial_state=state)
                        assert register_values(reverse,out,len(assignments))==[s[out] for s in assignments]
                        value_lanes+=len(assignments);del c,state,reverse
                        c,inputs,comparison,flag=o.compile_flag(graph,tables,mode);assignments=[];expected=[]
                        for sample,a in zip(samples,outputs):
                            for word in range(1<<graph.q):
                                for dirty in (0,1):
                                    s={inputs[n]:sample[n] for n in names};s[comparison]=word;s[(flag,)]=dirty
                                    assignments.append(s);expected.append(dirty^(word<a))
                        state,_=packed_replay(c,assignments);assert register_values(state,(flag,),len(assignments))==expected
                        for r in list(inputs.values())+[comparison]:assert register_values(state,r,len(assignments))==[s[r] for s in assignments]
                        keep=set(comparison+(flag,)+tuple(w for r in inputs.values() for w in r))
                        assert all(word==0 for wire,word in enumerate(state) if wire not in keep)
                        reversed_state,_=packed_replay(c,assignments,reverse=True,initial_state=state)
                        assert register_values(reversed_state,(flag,),len(assignments))==[s[(flag,)] for s in assignments]
                        flag_lanes+=len(assignments)
                        actual=independent_gate_accounting(c);resources=c.resources()
                        assert actual['actual_literal_T_dependency_depth']==resources['explicit_7T_sequence_T_depth']
                        counts=normalise(actual['expanded_literal_gate_counts']);counts['H']+=len(c.uniform_preparation_gates)
                        profile=o.counted_composition(graph,tables,count_cache,mode)
                        assert counts==Counter(profile['primitive_counts']) and profile['logical_T_count']==actual['logical_T_count']
                        assert profile['explicit_T_depth_serial_upper']>=actual['actual_literal_T_dependency_depth']
                        assert profile['logical_qubits_upper']>=c.wire_count
                        assert profile['uniform_H']==sum(graph.inputs.values())+graph.q
                        assert graph.configuration['monetary_probability_decode_factor']==(1 if payoff=='vanilla' or kind=='positive_r1' or encoding!='centred' else 2)
                        rows.append({'payoff':payoff,'kind':kind,'encoding':encoding,'mode':mode,'distribution':distribution,
                            'finite_input_states':len(samples),'value_outputs':outputs,'T_count':actual['logical_T_count'],
                            'logical_wires':c.wire_count,'complete_counts_verified':True,
                            'qualification':'Coarse complete-circuit composition fixture, not an attained monetary accuracy claim'})
                        del c,state,reversed_state
    for distribution in ('quantile','uniform_box'):
        for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2'):
            for mode in ('exact7t','unitary8t'):
                parts={r['encoding']:r['value_outputs'] for r in rows if r['distribution']==distribution and r['payoff']=='capped' and r['kind']==kind and r['mode']==mode}
                for centre,plus,minus in zip(parts['centred'],parts['positive_part'],parts['negative_part']):assert plus-minus==2*centre-(1<<5)
    zero_checks=0
    for payoff,cap in [('capped',25),('vanilla',None)]:
        for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
            for encoding in (('centred','positive_part','negative_part') if payoff=='capped' and kind!='positive_r1' else ('centred',)):
                g=o.build_graph(model,cap,q=5,b_y=1,beta_bits=1,table_bits=2,outer_bound=2,U_star=0,kind=kind,payoff=payoff,output_encoding=encoding)
                value=interpret(g,{}, {})[0];expected=16 if payoff=='capped' and kind!='positive_r1' and encoding=='centred' else 0
                assert value==expected;zero_checks+=1
    return {'status':'whole boundary financial circuit audit passed','graph_mode_configurations':len(rows),
        'dirty_numerical_target_basis_lanes':value_lanes,'dirty_flag_comparison_basis_lanes':flag_lanes,
        'zero_support_encoding_checks':zero_checks,'records':rows,
        'source_sha256':sha(root/'boundary_value_oracle.py'),'scalar_source_sha256':sha(root/'boundary_value_certificate.py'),
        'norm_file_sha256':sha(root/'results/boundary_derivative_norms.json'),
        'retained_emitter_sha256':sha(root/'integral_value_oracle.py'),
        'scope':'Whole financial value/flag/clean-inverse and literal primitive-count checks. Existing phase-correct primitive proof is reused; no physical quantum price execution or monetary-accuracy assertion from coarse fixtures.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',default=str(ROOT));p.add_argument('--out',default=str(ROOT/'results/boundary_value_oracle_audit.json'))
    args=p.parse_args();result=audit(args.root);result['audit_source_sha256']=sha(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k in ('status','graph_mode_configurations','dirty_numerical_target_basis_lanes','dirty_flag_comparison_basis_lanes','zero_support_encoding_checks')}))

if __name__=='__main__':main()
