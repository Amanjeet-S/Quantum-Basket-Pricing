"""Check the scoped matched-configuration Toffoli and certificate argument.

This maps each retained interval kernel to a full kernel at the same prefix
and probability precisions. It is not a circuit optimality assertion.
Correct endpoint rounding is a prerequisite verified in the retained
independent enclosure audit, whose configuration/source hashes are checked.
"""
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import encoding_resources as e
from resource_study import propagate


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def endpoints(packed,width):
    mask=(1<<width)-1
    return [row&mask for row in packed]+[(packed[-1]&mask)+(packed[-1]>>width)]


def main():
    config=ROOT/'results'/'conditioning_resources.json'
    audit=ROOT/'results'/'conditioning_independent_audit.json'
    data=json.loads(config.read_text());independent=json.loads(audit.read_text())
    assert independent['config_sha256']==sha(config)
    assert independent['root_encoding_sha256']==sha(ROOT/'encoding_resources.py')
    enclosed={(r['function'],r['address_bits'],r['fraction_bits'])
              for r in independent['function_tables'] if r['enclosed_endpoint_count']>0}
    tables={};cache={};rows=[]
    for cfg in data['configurations']:
        if cfg['method']!='interval_kernel':continue
        q=cfg['q'];cap=cfg['cap'];by=cfg['b_y'];bu=cfg['b_u']
        p=cfg['exp_address_bits'],cfg['log_address_bits'],cfg['cdf_address_bits']
        interval=e.build_graph('interval_kernel',q,by,0,bu,*p,cap)
        full=e.build_graph('full_kernel',q,by,0,0,*p,cap)
        e.make_tables(interval,tables);e.make_tables(full,tables)
        ikey=('scaled_cdf',p[2],q);fkey=(f'capped_{cap}',p[2],q)
        assert ikey in enclosed and fkey in enclosed
        iv,fv=tables[ikey],tables[fkey]
        assert all(a<=b for a,b in zip(endpoints(fv,full.width),endpoints(iv,interval.width)))
        assert all(b!=0 or a==0 for a,b in zip(fv,iv))
        ic,fc=propagate(interval,cap),propagate(full,cap)
        assert fc['arithmetic_price_error_upper']<=ic['arithmetic_price_error_upper']
        ir=e.counted_composition(interval,tables,cache)
        fr=e.counted_composition(full,tables,cache)
        assert ir['logical_T_count_exact_7T_Toffoli']==cfg['resources']['logical_T_count_exact_7T_Toffoli']
        it,ft=ir['expanded_gate_counts']['Toffoli'],fr['expanded_gate_counts']['Toffoli']
        def block_t(resources,predicate):
            return sum(b['counts']['Toffoli'] for b in resources['blocks'] if predicate(b))
        i_probability=block_t(ir,lambda b:b['function']=='scaled_cdf')
        f_probability=block_t(fr,lambda b:b['function']==f'capped_{cap}')
        strike=block_t(ir,lambda b:b['function']==f'strike_log_{cap}')
        assert f_probability<=i_probability
        assert it-ft==2*(strike+i_probability-f_probability)
        rows.append({'cap':cap,'price_tolerance':cfg['price_tolerance'],'q':q,
            'outer_gaussian_bits':by,'probability_address_bits':p[2],
            'full_endpoints_no_larger_than_scaled_cdf':True,
            'full_nonzero_row_support_contained_in_scaled_cdf':True,
            'full_probability_blocks_forward_Toffolis':f_probability,
            'interval_probability_blocks_forward_Toffolis':i_probability,
            'additional_strike_log_forward_Toffolis':strike,
            'flag_evaluator_Toffoli_gap':it-ft,'exact_gap_decomposition_verified':True,
            'interval_T':7*it,'matched_full_T':7*ft,
            'interval_arithmetic_price_bound':ic['arithmetic_price_error_upper'],
            'matched_full_arithmetic_price_bound':fc['arithmetic_price_error_upper'],
            'matched_full_distribution_bound':cfg['outer_coupling_or_BV_price_bound'],
            'interval_distribution_bound':cfg['distribution_price_bound'],
            'both_kernel_domains_verified':True})
    report={'status':'passed: six matched constructor mappings, not general resource optimality',
        'config_sha256':sha(config),'independent_enclosure_audit_sha256':sha(audit),
        'encoding_source_sha256':sha(ROOT/'encoding_resources.py'),
        'certificate_source_sha256':sha(ROOT/'resource_study.py'),
        'gate_library_sha256':sha(ROOT/'reversible.py'),'script_sha256':sha(Path(__file__)),
        'cases':rows,'limits':['assigned Toffoli/T dominance only within this exact prefix-tree and interpolation family',
            'does not order CNOT counts, depth, classical preprocessing or physical resources',
            'neither a quantum circuit lower bound nor an end-to-end estimator comparison']}
    target=ROOT/'results'/'kernel_dominance_audit.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'cases':len(rows),'output':str(target)}))


if __name__=='__main__':main()
