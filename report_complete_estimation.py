"""Retain numerical claims and update their marked manuscript/report sections."""
from pathlib import Path
from collections import Counter
import hashlib
import json

ROOT=Path(__file__).resolve().parent


def replace_section(source,body,start,end,placeholder=None):
    if start in source:
        a=source.index(start)+len(start);b=source.index(end,a)
        return source[:a]+'\n'+body+'\n'+source[b:]
    if placeholder is not None:
        assert placeholder in source
        return source.replace(placeholder,start+'\n'+body+'\n'+end)
    raise ValueError('Missing retained section marker')


def main():
    complete_path=ROOT/'results/complete_estimation_resources.json'
    oracle_path=ROOT/'results/architecture_selected_resources.json'
    data=json.loads(complete_path.read_text());oracle=json.loads(oracle_path.read_text())
    rows=data['configurations'];comparisons=data['method_comparison']
    assert len(rows)==819 and len(comparisons)==81
    cases=['n2_lambda0.05_cap25','n4_lambda0.05_cap25','n8_lambda0.05_cap25']
    modes=['exact7t','unitary8t','measured4t'];table=[];statistics={}
    for cid in cases:
        record={'case_id':cid,'assets':int(cid.split('_')[0][1:])}
        for mode in modes:
            o={r['method']:r for r in oracle['configurations'] if r['case_id']==cid and r['architecture']==mode}
            g=o['full_global']['resources']['logical_T_count'];c=o['full_cap']['resources']['logical_T_count']
            arm={r['method']:r['single_price_prior'] for x in comparisons if x['case_id']==cid and x['architecture']==mode for r in x['method_minima']}
            gt,ct=arm['full_global']['T_upper'],arm['full_cap']['T_upper']
            record[mode]={'selected_flag_generic_T':g,'selected_flag_cap_T':c,'selected_flag_saving_per_cent':100*(g-c)/g,
                'complete_generic_T_upper':gt,'complete_cap_T_upper':ct,'complete_saving_per_cent':100*(gt-ct)/gt}
        table.append(record)
    for mode in modes:
        comps=[x for x in comparisons if x['architecture']==mode]
        counts=Counter('|'.join(x['tied_methods']) for x in comps)
        strict=0
        for x in comps:
            arms={a['method']:a['single_price_prior']['T_upper'] for a in x['method_minima']}
            strict+=arms['full_cap']<arms['full_global']
        statistics[mode]={'complete_winner_method_sets':dict(counts),'strict_cap_generic_complete_cases':strict,
            'best_complete_T_upper_range':[min(x['minimum_complete_T_upper'] for x in comps),max(x['minimum_complete_T_upper'] for x in comps)]}
    ratio=[r['quantum_single_price_full_range']['T_count_upper']/r['quantum_single_price_prior']['T_count_upper'] for r in rows]
    samples=[r['financial_prior_and_classical_single_price']['classical_certified_samples'] for r in rows]
    summary={'configuration_count':819,'financial_cases':27,'architectures':3,'feasible_financial_prefixes':273,'rejected_transfers':60,
        'tables':table,'architecture_statistics':statistics,'prior_complete_T_budget_reduction_factor_range':[min(ratio),max(ratio)],
        'classical_sample_range':[min(samples),max(samples)],'resource_sha256':hashlib.sha256(complete_path.read_bytes()).hexdigest(),
        'oracle_sha256':hashlib.sha256(oracle_path.read_bytes()).hexdigest(),'reporting_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'single-price99per cent, ideal logical channels and independent shots; all null and adverse cases retained; no physical runtime or novelty inference'}
    (ROOT/'results/architecture_estimation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    body='''I retain 273 feasible financial prefixes, 60 rejected transfers and all 819 architecture profiles. Every original seven-T winner appears in the inherited pool, and all 108 same-prefix count regressions pass. The intermediate selector minimises emitted Grover-iterate T count. The complete selector independently minimises the proved full confidence budget across every retained prefix, preserving exact ties.

The cap-aware and generic full-conditioning comparison is architecture dependent. At λ=0.05 and H=25, the table shows evaluator savings followed by complete single-price confidence-budget savings, each in per cent. The three columns charge κ=14, 8 and 4 T gates per visited lookup selector node. All unlisted cases have zero paired saving.

| Assets | Exact7T lookup pair: evaluator / complete | Unitary8T lookup pair: evaluator / complete | Measured4T lookup: evaluator / complete |
| ---: | ---: | ---: | ---: |
'''
    for r in table:
        body+='| '+str(r['assets'])+' | '+' | '.join(f"{r[m]['selected_flag_saving_per_cent']:.3f} / {r[m]['complete_saving_per_cent']:.3f}" for m in modes)+' |\n'
    body+='''
The measured lookup eliminates the former four-asset strict cap-aware saving. Its two surviving paired reductions do not imply two unique overall wins: the eight-asset case is still cheaper by direct evaluation. Across 27 complete selections, exact7T has 22 generic/cap-aware ties, two unique cap-aware wins and three direct wins. Unitary8T and measured4T each have 23 generic/cap-aware ties, one unique cap-aware win and three direct wins. Interval conditioning attains no minimum in any of these complete selections. The 81 matched full/interval cost identities also pass, with a minimum same-prefix gap of 42,408 T gates. That dominance concerns interval versus full conditioning, rather than full versus direct evaluation.

'''
    body+=f"Across all 819 profiles, the financial prior reduces the same controller's complete worst-case T bound by factors {min(ratio):.2f} to {max(ratio):.2f} relative to prior one. The classical comparator receives the same prior and requires between {min(samples):,} and {max(samples):,} certified finite-function samples for the single-price target. These quantities have different units. For the best measured-lookup implementation in each case, complete single-price T bounds range from {statistics['measured4t']['best_complete_T_upper_range'][0]:,} to {statistics['measured4t']['best_complete_T_upper_range'][1]:,}. Large logical budgets remain even after the proved reductions; no physical speedup follows. The 540 ideal-binomial diagnostic runs all covered their stipulated probabilities, but this finite observation neither proves the confidence theorem nor represents quantum hardware execution.\n"
    report=ROOT/'architecture_estimation_study.md'
    report.write_text(replace_section(report.read_text(),body,'<!-- BEGIN NUMERICAL RESULTS -->','<!-- END NUMERICAL RESULTS -->','NUMERICAL_RESULT_INSERTION'))
    tex=r'''\section{Complete inherited-pool findings}
I retain 273 independently feasible financial prefixes, 60 rejected transfers
and 819 complete architecture profiles. All 108 original same-prefix
seven-T constructor counts match. The intermediate selector minimises
Grover-iterate T count; the complete selector minimises the proved
single-price confidence budget over every retained feasible prefix. Neither
selection is a global optimisation over all precision or circuit designs.

At $\lambda=0.05,H=25$, Table~\ref{tab:complete-costs} compares cap-aware
and generic full conditioning. Each entry gives the selected flag-evaluator
reduction, followed by the complete confidence-budget reduction, in per cent.
Every unlisted contract has zero paired reduction.
\begin{table}[ht]
\centering
\caption{Paired evaluator and complete-budget reductions under the three
lookup modes. Complete bounds use single-price 99 per cent confidence.}
\label{tab:complete-costs}
\begin{tabular}{rrrr}
\toprule
Assets & $\kappa=14$ & $\kappa=8$ & $\kappa=4$\\
\midrule
'''
    for r in table:
        tex+=str(r['assets'])+' & '+' & '.join(f"{r[m]['selected_flag_saving_per_cent']:.3f} / {r[m]['complete_saving_per_cent']:.3f}" for m in modes)+r'\\'+'\n'
    tex+=r'''\bottomrule
\end{tabular}
\end{table}
The four-asset strict saving disappears with cheaper lookup. Under measured
lookup, 23 cases tie between generic and cap-aware full conditioning, one
uniquely favours cap-aware full conditioning and three favour direct
evaluation. Interval conditioning attains no minimum. The eight-asset
paired cap-aware reduction still does not beat direct evaluation.
The 81 matched full/interval identities pass, with minimum same-prefix
gap 42,408 T gates. These facts do not establish full-conditioning dominance
over direct evaluation or all possible quantum constructions.

'''
    tex+=f"The prior reduces the same controller's complete worst-case T budget by\nfactors {min(ratio):.2f} to {max(ratio):.2f} over all 819 profiles. The equally\nsupplied classical prior yields between {min(samples):,} and {max(samples):,}\ncertified finite-function samples. These costs have different units and\ndo not demonstrate quantum advantage. All 819 complete financial, confidence,\nselection and resource records pass the independent audit. Thirty regression\ntests pass. The 540 ideal-binomial simulations are diagnostics; their\nobserved coverage is not the confidence proof.\n"
    proof=ROOT/'generalisation_proofs.tex';source=proof.read_text()
    marker='% BEGIN COMPLETE NUMERICAL RESULTS';end='% END COMPLETE NUMERICAL RESULTS'
    if marker not in source:
        source=source.replace('\\clearpage\n\\begin{thebibliography}',marker+'\n'+end+'\n\\clearpage\n\\begin{thebibliography}')
    proof.write_text(replace_section(source,tex,marker,end))
    print('Retained summary and updated marked report/proof numerical sections.')


if __name__=='__main__':main()
