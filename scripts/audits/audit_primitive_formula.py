"""Check selector/arithmetic resource formula against actual emitted blocks."""
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import architecture as arch
import encoding_resources as er
import generalisation as gen


def primitive(kind,bits,q,table,mode):
    width=q+7;span=er.function_span(kind);f=q+span.bit_length()-1-bits
    c=arch.ArchitectureCircuit(store_gates=False,mode=mode);a=c.alloc(width)
    node=er.Node('function',('a',),kind,bits)
    er.emit_node(c,node,[a],q,width,table)
    m=c.resources();internal=m['lookup_selector_internal_nodes'];k={'exact7t':14,'unitary8t':8,'measured4t':4}[mode]
    arithmetic=2*f*(4*width+3*f)+6*width
    expected=k*internal+7*arithmetic
    assert m['arithmetic_exact_Toffoli_count']==arithmetic,(kind,bits,q,m['arithmetic_exact_Toffoli_count'],arithmetic)
    assert m['logical_T_count']==expected
    return {'kernel':kind,'address_bits':bits,'q':q,'w':width,'f':f,'architecture':mode,
        'kappa':k,'internal_nodes':internal,'emitted_T_count':m['logical_T_count'],
        'arithmetic_Toffolis':arithmetic,'formula_T_count':expected,
        'measurement_count':m['measurement_count']}


def main():
    data=json.loads((ROOT/'results/generalisation_resources.json').read_text());rows=[];cache={};keys=set()
    # The retained winners fix the financial input profiles. Also check the
    # one-bit coarsening at fixed q, including a pruned cap/CDF support.
    for model_id,cap,method in [('n2_lambda0.05',25,'full_cap'),('n4_lambda0.05',25,'full_cap'),
                              ('n8_lambda0.05',25,'full_cap'),('n8_lambda0.4',25,'interval')]:
        cfg=next(r for r in data['configurations'] if r['model_id']==model_id and r['cap']==cap and r['method']==method)
        model=next(m for m in data['models'] if m['model_id']==model_id)
        g=gen.build_graph(model,method,*[cfg[k] for k in ('q','b_y','b_z','b_u','exp_address_bits','log_address_bits','probability_address_bits')],cap)
        gen.make_tables(g,cache)
        for node in g.nodes.values():
            if node.kind=='function':keys.add((node.meta,node.bits,g.q))
    for kind,bits,q in sorted(keys):
        table=cache[(kind,bits,q)]
        coarsekey=(kind,bits-1,q)
        if coarsekey not in cache:
            import general_kernels as kernels
            cache[coarsekey]=kernels.function_table(*coarsekey,q+7)
        for mode in arch.MODES:
            fine=primitive(kind,bits,q,table,mode);coarse=primitive(kind,bits-1,q,cache[coarsekey],mode)
            k=fine['kappa'];deltaI=fine['internal_nodes']-coarse['internal_nodes'];w=fine['w'];f=fine['f']
            saving=fine['emitted_T_count']-coarse['emitted_T_count']
            formula=k*deltaI-7*(8*w+12*f+6)
            assert saving==formula
            rows.append({'fine':fine,'coarse':coarse,'primitive_T_saving':saving,'primitive_T_saving_formula':formula,
                'four_call_flag_T_saving':4*saving,'four_call_flag_T_saving_formula':4*formula,
                'strict_saving_iff_kappa_deltaI_exceeds_arithmetic_increment':(saving>0)==(k*deltaI>7*(8*w+12*f+6))})
    out={'status':'primitive and fixed-q one-bit-coarsening formulas match emitted gates',
        'architecture_sha256':hashlib.sha256((ROOT/'architecture.py').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'historical_resource_sha256':hashlib.sha256((ROOT/'results/generalisation_resources.json').read_bytes()).hexdigest(),
        'primitive_formula':'T_b = kappa I_b + 7 [2 f(4w+3f) + 6w]',
        'four_call_coarsening_formula':'4 [kappa (I_b-I_(b-1)) - 7(8w+12f+6)]',
        'configurations':rows,
        'qualification':'fixed q/span/input profile only; coarsening is not asserted financially feasible by this resource-only check'}
    path=ROOT/'results/primitive_formula_audit.json';path.write_text(json.dumps(out,indent=2)+'\n');print(path)


if __name__=='__main__':main()
