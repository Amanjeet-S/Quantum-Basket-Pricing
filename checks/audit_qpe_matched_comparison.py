"""Independent archive-link and matched-resource audit, without repricing."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    old=json.loads(gzip.decompress((ROOT/'results/variance_oracle_comparison.json.gz').read_bytes()))
    new=json.loads(gzip.decompress((ROOT/'results/qpe_benchmark.json.gz').read_bytes()))
    result=json.loads((ROOT/'results/qpe_matched_comparison.json').read_text())
    assert result['source_sha256']==sha(ROOT/'qpe_comparison.py')
    for name,key in [('variance_oracle_comparison.json.gz','old_archive'),('qpe_benchmark.json.gz','qpe_benchmark.json.gz'),('architecture_resource_candidates.json','architecture_resource_candidates.json')]:
        assert result['input_hashes'][key]==sha(ROOT/'results'/name)
    old_rows={(r['candidate_id'],r['architecture']):r for r in old['configurations']}
    new_rows={(r['candidate_id'],r['architecture']):r for r in new['configurations']}
    old_winners={r['case_id']:r for r in old['selected_cases'] if r['family']=='best_single_flag'}
    new_winners={r['case_id']:r for r in new['selected_cases'] if r['family']=='best_single_flag'}
    ratios=[];records=[]
    assert len(result['configurations'])==27
    for r in result['configurations']:
        case=r['case_id'];ow=old_winners[case];key=(r['matched_candidate_id'],r['matched_architecture'])
        assert key==(ow['candidate_id'],ow['architecture'])
        previous=old_rows[key]['controllers'][ow['controller']]['best_single_flag']
        assert r['old_controller']==ow['controller'] and r['matched_method']==ow['method']
        assert r['matched_flag_kind']==previous['kind'] and r['matched_flag_c']==previous['c']
        assert r['old_complete_T_cap']==previous['T'] and r['old_peak_qubits_upper']==previous['qubits']
        q=r['matched_QPE']
        assert q['prior']==r['shared_probability_prior']==previous['prior']
        assert q['epsilon']==r['normalised_probability_tolerance']==previous['epsilon']
        assert F(q['alpha'])==F(1,100)
        choices=new_rows[key]['clipping_candidates' if previous['kind']=='clipped' else 'range_candidates']
        candidate=next(x for x in choices if x['c']==previous['c'])
        # The complete scalar is identical to a record checked by the independent
        # full-archive auditor; candidate-only tail metadata is excluded here.
        assert all(candidate[k]==v for k,v in q.items()),case
        ratio=F(q['T'],previous['T']);assert ratio==F(r['matched_T_ratio_exact_fraction'])
        assert r['optimised_QPE_selection']==new_winners[case]
        assert F(r['optimised_T_ratio_exact_fraction'])==F(new_winners[case]['T'],previous['T'])
        ratios.append(ratio);records.append({'case_id':case,'exact_ratio':str(ratio),'all_matched_inputs_and_scalar_metrics_identical':True})
    assert len({r['case_id'] for r in records})==27
    assert all(x<1 for x in ratios)
    summary=result['summary'];assert summary['matched_cases']==27 and summary['matched_QPE_lower_caps']==27
    assert summary['matched_minimum_ratio']==float(min(ratios)) and summary['matched_maximum_ratio']==float(max(ratios))
    report={'schema':'qbp-independent-matched-qpe-audit-v1','status':'passed','cases':27,
            'scope':'Archive identities and every matched scalar field linked to independently audited full QPE records; exact upper-budget ratios. No new execution or hardware claim.',
            'records':records,'minimum_ratio_exact_fraction':str(min(ratios)),'maximum_ratio_exact_fraction':str(max(ratios)),
            'source_sha256':{n:sha(ROOT/n) for n in ['qpe_comparison.py','checks/audit_qpe_matched_comparison.py']},
            'input_sha256':{n:sha(ROOT/n) for n in ['results/variance_oracle_comparison.json.gz','results/qpe_benchmark.json.gz','results/qpe_matched_comparison.json']}}
    (ROOT/'results/qpe_matched_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','matched_cases':27,'minimum_ratio':float(min(ratios)),'maximum_ratio':float(max(ratios))}))

if __name__=='__main__':main()
