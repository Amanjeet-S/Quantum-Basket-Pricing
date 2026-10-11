"""Exact descriptive charges and break-even identity for the frozen diagnostics.

The reported finite mesh averages are not expectations for an unknown
financial flag probability. This calculation uses every retained history;
it does not select a favourable seed, infer a device runtime or alter the
registered estimator policy.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, statistics
ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(results, output):
    path = results/'chebyshev_resources.json'
    data = json.loads(path.read_text())
    rows = []
    for case in data['cases']:
        baseline = case['retained_financial_information']['previous_complete_QPE_T']
        fallback = case['half_alpha_QPE_fallback']['T_count']
        budget = case['pilot_T_cap']
        assert F(fallback, baseline) == F(27, 25)
        trials = case['diagnostics']
        count = len(trials)
        successes = sum(x['controller_result']['status'] == 'certified' for x in trials)
        pilot = F(sum(x['controller_result']['pilot_T_count'] for x in trials), count)
        charge = F(sum(x['actual_history_T_count'] for x in trials), count)
        rate = F(successes, count)
        assert charge == pilot+(1-rate)*fallback
        threshold = 1-F(baseline, fallback)+pilot/fallback
        assert (charge < baseline) == (rate > threshold)
        median = F(str(statistics.median([x['actual_history_T_count'] for x in trials])))
        assert max(x['actual_history_T_count'] for x in trials) <= fallback+budget
        assert 1-F(baseline, fallback) == F(2, 27)
        assert 1-F(baseline, fallback)+F(1, 10) == F(47, 270)
        rows.append({'model_id':case['model_id'], 'histories':count,
            'pilot_certified':successes, 'pilot_success_fraction':str(rate),
            'previous_complete_QPE_T':baseline,'half_alpha_fallback_T':fallback,
            'pilot_T_budget':budget,'complete_worst_history_T_cap':fallback+budget,
            'diagnostic_mean_pilot_T':str(pilot),
            'diagnostic_mean_complete_T':str(charge),
            'diagnostic_median_complete_T':str(median),
            'diagnostic_mean_over_previous_QPE':str(charge/baseline),
            'diagnostic_median_over_previous_QPE':str(median/baseline),
            'diagnostic_exact_break_even_threshold':str(threshold),
            'diagnostic_mean_below_previous_QPE':charge<baseline,
            'diagnostic_median_below_previous_QPE':median<baseline,
            'general_necessary_success_fraction':str(F(2,27)),
            'general_sufficient_success_fraction':str(F(2,27)+F(budget,fallback)),
            'uniform_conservative_sufficient_fraction':str(F(47,270))})
    result={'schema':'qbp-chebyshev-exact-interpretation-v1','status':'passed',
        'producer_sha256':sha(__file__),'resource_result_sha256':sha(path),
        'analysis_timing':'Exact descriptive interpretation after the registered diagnostic execution; no additional pre-outcome registration is claimed.',
        'rows':rows,'total_histories':sum(r['histories'] for r in rows),
        'total_pilot_certified':sum(r['pilot_certified'] for r in rows),
        'identity':'Complete charge = pilot charge + fallback charge times unsuccessful-pilot indicator. No independence of success and cost is assumed.',
        'scope':'Finite mesh fractions, means and medians are descriptive. They do not identify the unknown financial flag, certify expected runtime or prove quantum advantage. The general break-even identity uses expectations under an explicitly specified shot law.'}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('Exact descriptive and break-even identities passed for four cases.')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--results',type=Path,default=ROOT/'results')
    p.add_argument('--output',type=Path,default=ROOT/'results/chebyshev_interpretation.json')
    a=p.parse_args();run(a.results,a.output)
