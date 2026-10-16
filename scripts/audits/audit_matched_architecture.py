"""Check matched conditional flag costs using retained emitted profiles only.

The cases and prefixes were fixed in the earlier matched-interval diagnostic.
This check does not launch or optimise new circuits and does not compare
retuned winners at unequal q or function address precision.
"""
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import encoding_resources as er


def internal_nodes(block,row):
    counts=block['primitive_gate_counts']
    if row['architecture']!='exact7t':
        return counts.get('SD',0)
    q=row['q'];w=q+7
    span=er.function_span(block['function'])
    f=q+span.bit_length()-1-block['address_bits']
    arithmetic=2*f*(4*w+3*f)+6*w
    selector=counts.get('CCX',0)-arithmetic
    assert selector>=0 and selector%2==0
    return selector//2


def matching(row,cfg,method,mode):
    fields=('q','b_y','exp_address_bits','log_address_bits','probability_address_bits')
    return row['case_id']==cfg['case_id'] and row['method']==method and row['architecture']==mode and all(row[k]==cfg[k] for k in fields)


def main():
    path=ROOT/'results/architecture_resource_candidates.json'
    data=json.loads(path.read_text())
    historical_path=ROOT/'results/generalisation_matched_interval.json'
    matched=json.loads(historical_path.read_text())
    rows=[]
    for cfg in matched['configurations']:
        for mode,kappa in (('exact7t',14),('unitary8t',8),('measured4t',4)):
            interval=[r for r in data['configurations'] if matching(r,cfg,'interval',mode)]
            full=[r for r in data['configurations'] if matching(r,cfg,'full_cap',mode)]
            assert len(interval)==len(full)==1,(cfg['case_id'],mode,len(interval),len(full))
            interval,full=interval[0],full[0]
            funcs=lambda r:[b for b in r['resources']['blocks'] if b['kind']=='function']
            ib=funcs(interval);fb=funcs(full)
            strike=[b for b in ib if b['function'].startswith('general_strike:')]
            cdf=[b for b in ib if b['function'].startswith('general_cdf:')]
            cap=[b for b in fb if b['function'].startswith('general_cap:')]
            assert len(strike)==1 and len(cdf)==len(cap)==2
            cdfi=[internal_nodes(b,interval) for b in cdf]
            capi=[internal_nodes(b,full) for b in cap]
            assert cdfi[0]==cdfi[1] and capi[0]==capi[1]
            assert capi[0]<=cdfi[0]
            strike_t=strike[0]['logical_T_count']
            actual=interval['resources']['logical_T_count']-full['resources']['logical_T_count']
            expected=2*strike_t+168*(cfg['q']+7)+4*kappa*(cdfi[0]-capi[0])
            assert actual==expected and actual>0,(cfg['case_id'],mode,actual,expected)
            # The two removed subtractions and their financial inverses
            # account together for 24w exact arithmetic Toffolis.
            arithmetic_gap=interval['resources']['arithmetic_exact_Toffoli_count']-full['resources']['arithmetic_exact_Toffoli_count']
            row={'case_id':cfg['case_id'],'architecture':mode,'kappa':kappa,
                 'q':cfg['q'],'w':cfg['q']+7,
                 'exp_address_bits':cfg['exp_address_bits'],
                 'log_address_bits':cfg['log_address_bits'],
                 'probability_address_bits':cfg['probability_address_bits'],
                 'cdf_selector_internal_nodes':cdfi[0],'cap_selector_internal_nodes':capi[0],
                 'strike_forward_block_T':strike_t,'emitted_flag_T_gap':actual,
                 'predicted_flag_T_gap':expected,'exact_arithmetic_Toffoli_gap':arithmetic_gap,
                 'interval_flag_T':interval['resources']['logical_T_count'],
                 'full_flag_T':full['resources']['logical_T_count'],
                 'full_to_interval_flag_T_ratio':full['resources']['logical_T_count']/interval['resources']['logical_T_count']}
            rows.append(row)
    out={'status':'all retained matched-profile architecture identities passed',
         'formula':'T_interval_flag - T_full_flag = 2 T_strike_block + 168 w + 4 kappa (I_cdf-I_cap)',
         'scope':'27 prespecified historical matched prefixes under three gate architectures; flag costs include financial inverse, not uniform preparation or reflection',
         'candidate_resources_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'matched_prefix_source_sha256':hashlib.sha256(historical_path.read_bytes()).hexdigest(),
         'architecture_sha256':hashlib.sha256((ROOT/'architecture.py').read_bytes()).hexdigest(),
         'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'configurations':rows,
         'summary':{'matched_cases':27,'architecture_checks':len(rows),
                    'minimum_positive_flag_T_gap':min(r['emitted_flag_T_gap'] for r in rows),
                    'maximum_full_to_interval_flag_T_ratio':max(r['full_to_interval_flag_T_ratio'] for r in rows)}}
    target=ROOT/'results/matched_architecture_audit.json'
    target.write_text(json.dumps(out,indent=2)+'\n');print(target)


if __name__=='__main__':main()
