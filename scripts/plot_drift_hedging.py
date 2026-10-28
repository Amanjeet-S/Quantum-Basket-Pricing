#!/usr/bin/env python3
"""Reproduce two illustrative drift-risk panels from audited synthetic inputs.

The independent piecewise lognormal reference evaluates each curve at 160
decimal digits before the separate binary64 display conversion. These plotted
values are numerical diagnostics, not interval proofs or additional samples.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def generate(project_root, results, output):
    root, results, output = (Path(p).resolve() for p in (project_root, results, output))
    if output.exists() or output == root or root in output.parents:
        raise ValueError('Output must be a fresh directory outside the project.')
    if output == results or results in output.parents:
        raise ValueError('Output must not modify or be inside the input study.')
    sys.path.insert(0, str(root/'src'))
    import mpmath as mp
    import gbm_hedge_references as reference

    manifest = json.loads((results/'file_hashes.json').read_text())
    for name, expected in manifest.items():
        if Path(name).name != name or digest(results/name) != expected:
            raise ValueError('Input study manifest does not match its retained files.')
    validation = json.loads((results/'validation.json').read_text())
    if not validation.get('all_checks_pass', False):
        raise ValueError('The input study has not passed its retained validation.')
    records = json.loads((results/'primary.json').read_text())
    key = 'inside_spread:width=3/20:budget=1'
    matches = [r for r in records if r['key'] == key]
    if len(matches) != 1:
        raise ValueError('Exactly one declared illustrative primary record is required.')
    record = matches[0]
    registered = json.loads((results/'registration.json').read_text())
    paths = [Path(__file__).resolve(), Path(reference.__file__).resolve(),
             root/'results/capped_gbm/hedging/study_manifest.json',
             results/'file_hashes.json']+[results/name for name in manifest]
    paths += [Path(p) for p in registered.get('source_hashes_before', {}) if Path(p).is_file()]
    paths = sorted(set(paths))
    before = {str(p): digest(p) for p in paths}
    output.mkdir(parents=True, exist_ok=False)

    fixed_contract = SimpleNamespace(spot=F(100), strike=F(100), cap=F(20),
                                    maturity=F(1), rate=F(3, 100), dividend=F(1, 100))
    fixed_sigma, fixed_theta = F(1, 100), (F(0), F(1, 2))
    fixed_grid = [F(21, 200)+F(17, 100)*F(i, 256) for i in range(257)]
    contract = SimpleNamespace(**{k: F(v) for k, v in record['certificate']['contract_exact'].items()})
    sigma = F(record['certificate']['volatility_exact'])
    lo, hi = map(F, record['drift_interval_exact'])
    hedge_grid = [lo+(hi-lo)*F(i, 256) for i in range(257)]
    names = ('robust', 'centre', 'oracle_true_drift', 'endpoint_only')
    holdings = {name: tuple(map(F, record['comparisons'][name]['holdings_exact'])) for name in names}
    rows = []
    with mp.workdps(160):
        fixed_risks = [reference.fixed_risk(reference.parameters(fixed_contract, fixed_sigma, mu), fixed_theta)
                       for mu in fixed_grid]
        endpoint_max = max(fixed_risks[0], fixed_risks[-1])
        curves = {name: [] for name in names}
        for mu in hedge_grid:
            model = reference.parameters(contract, sigma, mu)
            for name in names:
                curves[name].append(reference.fixed_risk(model, holdings[name]))
        for i, (first_mu, second_mu) in enumerate(zip(fixed_grid, hedge_grid)):
            row = {'counterexample_drift_exact': str(first_mu),
                   'counterexample_drift': format(float(first_mu), '.17g'),
                   'counterexample_risk_160digits': mp.nstr(fixed_risks[i], 160),
                   'counterexample_risk': format(float(fixed_risks[i]), '.17g'),
                   'endpoint_maximum': format(float(endpoint_max), '.17g'),
                   'hedge_drift_exact': str(second_mu),
                   'hedge_drift': format(float(second_mu), '.17g')}
            for name in names:
                row[name+'_risk_160digits'] = mp.nstr(curves[name][i], 160)
                row[name+'_risk'] = format(float(curves[name][i]), '.17g')
            if not all(math.isfinite(float(row[column])) for column in
                       ('counterexample_risk', 'endpoint_maximum')+tuple(name+'_risk' for name in names)):
                raise ArithmeticError('A display curve is not finite.')
            rows.append(row)
        if not fixed_risks[128] > endpoint_max:
            raise ArithmeticError('The declared interior counterexample was not reproduced.')
        checks = {}
        for name in names:
            calculated = max(curves[name])
            retained = mp.mpf(record['comparisons'][name]['dense_maximum'])
            checks[name] = abs(calculated-retained) <= mp.mpf('1e-145')*max(1, abs(retained))
        if not all(checks.values()):
            raise ArithmeticError('Independent display curves differ from the retained study.')
        numerical = {'counterexample_markers': [
            {'drift_exact': str(fixed_grid[i]), 'risk_160digits': mp.nstr(fixed_risks[i], 160)}
            for i in (0, 128, 256)],
            'counterexample_endpoint_maximum_160digits': mp.nstr(endpoint_max, 160),
            'independent_curve_replay': checks}

    curve_path = output/'drift_robustness_curves.csv'
    with curve_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    fragment = r'''% Generated by scripts/plot_drift_hedging.py; see the curve manifest.
\definecolor{driftblue}{HTML}{0072B2}
\definecolor{driftorange}{HTML}{E69F00}
\definecolor{driftgreen}{HTML}{009E73}
\begin{tikzpicture}
\begin{groupplot}[group style={group size=2 by 1,horizontal sep=14mm},
 width=.44\textwidth,height=5.1cm,
 xlabel={Physical total return $\mu$},ylabel={Mean-square loss},
 tick label style={font=\scriptsize},label style={font=\small},
 title style={font=\small},grid=major,grid style={gray!15},
 scaled ticks=false,legend style={font=\scriptsize,draw=none,fill=white,
 at={(.03,.97)},anchor=north west,inner sep=1.5pt},
 legend cell align=left]
\nextgroupplot[title={(a) An interior worst case},xmin=.105,xmax=.275,
 ymin=0,ymax=82,xtick={.105,.19,.275},xticklabels={0.105,0.190,0.275},
 legend style={font=\scriptsize,draw=none,fill=white,at={(.03,.03)},anchor=south west,inner sep=1.5pt}]
\addplot[driftblue,thick,mark=*,mark size=1.5pt,mark indices={0,128,256}]
 table[x=counterexample_drift,y=counterexample_risk,col sep=comma]
 {figures/drift_robustness_curves.csv};
\addplot[black,dashed,thick] table[x=counterexample_drift,y=endpoint_maximum,col sep=comma]
 {figures/drift_robustness_curves.csv};
\nextgroupplot[title={(b) Funded hedging trade-off},xmin=-.09,xmax=.21,
 ymin=10,ymax=56,xtick={-.09,.03,.09,.21},xticklabels={$-0.09$,0.03,0.09,0.21}]
\addplot[driftblue,very thick] table[x=hedge_drift,y=robust_risk,col sep=comma]
 {figures/drift_robustness_curves.csv};
\addlegendentry{Robust}
\addplot[driftorange,thick,dashed] table[x=hedge_drift,y=centre_risk,col sep=comma]
 {figures/drift_robustness_curves.csv};
\addlegendentry{Centre}
\addplot[driftgreen,thick,dashdotted] table[x=hedge_drift,y=oracle_true_drift_risk,col sep=comma]
 {figures/drift_robustness_curves.csv};
\addlegendentry{True-drift oracle}
\addplot[black,densely dotted,forget plot] coordinates {(.03,10) (.03,56)};
\end{groupplot}
\end{tikzpicture}
'''
    (output/'drift_robustness.tex').write_text(fragment)
    after = {str(p): digest(p) for p in paths}
    if before != after:
        raise ArithmeticError('A source or input changed during figure generation.')
    save_json(output/'plot_manifest.json', {
        'schema': 'qbp-drift-risk-figure-v1', 'study_key': key,
        'selection_status': 'Informed illustrative selection after the study outcomes; all twelve primary records remain retained.',
        'panel_a_inputs_exact': {k: str(v) for k, v in vars(fixed_contract).items()},
        'panel_a_sigma_exact': str(fixed_sigma), 'panel_a_holdings_exact': ['0', '1/2'],
        'panel_b_holdings_exact': {k: [str(x) for x in v] for k, v in holdings.items()},
        'panel_b_contract_exact': record['certificate']['contract_exact'],
        'panel_b_sigma_exact': str(sigma), 'panel_b_true_drift_exact': record['true_drift_exact'],
        'panel_b_budget_multiplier_exact': record['budget_multiplier_exact'],
        'points_per_curve': 257, 'calculation_decimal_digits': 160,
        'display_rule': 'Evaluate independent piecewise moments at 160 digits, retain their decimal diagnostic, then convert separately to 17-digit binary64 display coordinates.',
        'omitted_curve_reason': 'The endpoint-only line is close to the robust line at this panel scale; its full curve is retained in the CSV. Three plotted policies make the nominal-risk versus worst-risk trade-off legible.',
        'source_hashes_before': before, 'source_hashes_after': after,
        'sources_unchanged': before == after,
        'registered_study_source_hashes': registered.get('source_hashes_before', {}),
        'numerical_checks': numerical,
        'output_sha256': {name: digest(output/name) for name in
                          ('drift_robustness_curves.csv', 'drift_robustness.tex')},
        'scope': 'Synthetic risk curves and illustrative comparisons. High-precision diagnostics are not certified continuum extrema, additional data or a general dominance result.'})
    return {'output': str(output), 'rows': len(rows), 'sources_unchanged': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(args.project_root, args.results, args.output), sort_keys=True))


if __name__ == '__main__':
    main()
