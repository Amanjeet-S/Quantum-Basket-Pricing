"""Literal basis replay of projection, comparator and clean inverse.

The negative branch is tested at the independently derived scalar witness.
This is a primitive financial-word test, not a sampled option-price claim.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,sys
sys.dont_write_bytecode=True

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',required=True)
    ap.add_argument('--resources',required=True);ap.add_argument('--output',required=True)
    a=ap.parse_args();repo=Path(a.repo).resolve();out=Path(a.output).resolve()
    assert not out.exists();sys.path.insert(0,str(repo))
    import boundary_value_oracle as e
    data=json.load(gzip.open(a.resources,'rt'))
    widths=sorted({p['profiles']['exact7t']['centred']['integer_word_width'] for p in data['records']})
    q=28;s=1<<q
    pairs=[(0,0),(4369400,4369500),(4369500,4369400),(s,0),(0,s),(s,s)]
    checked=[]
    for width in widths:
        for mode in ('exact7t','unitary8t'):
            c=e.old.vo.ValueCircuit(store_gates=True,mode=mode)
            A,B,U=c.alloc(width),c.alloc(width),c.alloc(q)
            zero=c.alloc(1);flagreg=c.alloc(1);flag=flagreg[0]
            snap=c.snapshot()
            difference=e.old.clean_selected_difference(c,A,B,positive=True)
            projected=e.old.clean_clamp(c,difference,s)
            end=c.gate_count
            c.compare_into(U+zero,projected[:q+1],flag)
            c.uncompute(snap,end_gate_count=end)
            retained=set(A+B+U+zero+flagreg)
            n=0
            for x,y in pairs:
                expected=min(s,max(0,x-y))
                for u in (0,99,100,s-1):
                    for dirty in (0,1):
                        state=c.run({A:x,B:y,U:u,flagreg:dirty})
                        assert c.value(state,A)==x and c.value(state,B)==y and c.value(state,U)==u
                        assert state[flag]==(dirty^(u<expected)) and state[zero[0]]==0
                        assert all(bit==0 for w,bit in enumerate(state) if w not in retained)
                        n+=1
            checked.append({'word_width':width,'architecture':mode,'basis_inputs':n,
                'negative_witness_payoff_words':[4369400,4369500],'negative_output_probability_integer':0,
                'dirty_flag_XOR_correct':True,'complete_private_workspace_clean':True,
                'literal_primitive_counts':c.resources()['clifford_T_expanded_gate_counts']})
    result={'schema':'qbp-projected-flag-negative-branch-v1','status':'passed',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'resources_sha256':hashlib.sha256(Path(a.resources).read_bytes()).hexdigest(),
        'quantum_fractional_bits':q,'checks':checked,
        'scope':'Literal Boolean permutation replay of the exact emitted positive-subtraction, clamp, q-bit comparator and inverse. Established exact Clifford/T decomposition is reused. No full-grid occurrence, hardware execution, price outcome or independent physical noise model is asserted.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print('Projected flag branch and clean inverse passed',sum(r['basis_inputs'] for r in checked),'basis assignments.')

if __name__=='__main__':main()
