"""Matched old-controller flags and optimised canonical-QPE price caps."""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib
import value_oracle as vo
from qpe_benchmark import scalar
ROOT=Path(__file__).resolve().parent
DEFAULT_OLD=ROOT/'results/variance_oracle_comparison.json.gz'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--old-archive',type=Path,default=DEFAULT_OLD);args=a.parse_args()
    prior=json.loads(gzip.decompress(args.old_archive.read_bytes()));new=json.loads(gzip.decompress((ROOT/'results/qpe_benchmark.json.gz').read_bytes()));original=json.loads((ROOT/'results/architecture_resource_candidates.json').read_text())
    oldrows={(r['candidate_id'],r['architecture']):r for r in prior['configurations']};src={(r['candidate_id'],r['architecture']):r for r in original['configurations']};selected={(r['case_id'],r['family']):r for r in new['selected_cases']};records=[]
    for winner in prior['selected_cases']:
        if winner['family']!='best_single_flag':continue
        row=oldrows[(winner['candidate_id'],winner['architecture'])];prev=row['controllers'][winner['controller']]['best_single_flag'];s=src[(winner['candidate_id'],winner['architecture'])];body=vo.prepare_existing_body(s['resources'],s['q']);c=prev['c']
        flag=vo.profile_from_prepared_clipped_range(body,c) if prev['kind']=='clipped' else vo.profile_from_prepared_range(body,c)
        matched=scalar(F(prev['epsilon']),F(prev['prior']),F(1,100),flag);best=selected[(winner['case_id'],'best_single_flag')]
        records.append({'case_id':winner['case_id'],'assets':int(winner['case_id'].split('_')[0][1:]),'cap':row['cap'],'lambda':winner['case_id'].split('lambda')[1].split('_')[0],
            'matched_candidate_id':winner['candidate_id'],'matched_architecture':winner['architecture'],'matched_method':winner['method'],'matched_flag_kind':prev['kind'],'matched_flag_c':c,'normalised_probability_tolerance':prev['epsilon'],'shared_probability_prior':prev['prior'],
            'old_controller':winner['controller'],'old_complete_T_cap':prev['T'],'old_peak_qubits_upper':prev['qubits'],
            'matched_QPE':matched,'matched_T_ratio_exact_fraction':str(F(matched['T'],prev['T'])),
            'optimised_QPE_selection':best,'optimised_T_ratio_exact_fraction':str(F(best['T'],prev['T']))})
    matched=[F(r['matched_T_ratio_exact_fraction']) for r in records];optim=[F(r['optimised_T_ratio_exact_fraction']) for r in records]
    result={'status':'matched quantum confidence budgets reconstructed','source_sha256':sha(Path(__file__)),
        'input_hashes':{'old_archive':sha(args.old_archive),'qpe_benchmark.json.gz':sha(ROOT/'results/qpe_benchmark.json.gz'),'architecture_resource_candidates.json':sha(ROOT/'results/architecture_resource_candidates.json')},
        'comparison_scope':'same financial circuit, fixed comparison flag, proved probability prior, probability tolerance and99% per-price confidence for matched rows;optimisedQPE additionalfinitepoolsearch separatelylabelled',
        'configurations':records,'summary':{'matched_cases':len(records),'matched_QPE_lower_caps':sum(x<1 for x in matched),'matched_minimum_ratio':float(min(matched)),'matched_maximum_ratio':float(max(matched)),'optimised_QPE_lower_caps':sum(x<1 for x in optim),'optimised_minimum_ratio':float(min(optim)),'optimised_maximum_ratio':float(max(optim))},
        'interpretation':'smaller proved logical resource upper caps, not measured stochastic cost dominance, a complexity lower bound or physical quantum advantage; classical financial tolerance already resolved onthese27contracts'}
    def pack(x):
        if isinstance(x,F):return str(x)
        if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
        if isinstance(x,(tuple,list)):return [pack(v) for v in x]
        return x
    (ROOT/'results/qpe_matched_comparison.json').write_text(json.dumps(pack(result),indent=2)+'\n');print(json.dumps(result['summary'],indent=2))
if __name__=='__main__':main()
