"""Colour budget-risk frontier and its normalised shadow-price diagnostic.

The dense grid illustrates three preregistered contracts. It is supplemental
numerical evidence, not a replacement for the full registered comparison.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from capped_gbm import synthetic_contracts
from gbm_hedging import hedge_moments, solve_constrained_hedge, residual_square_integral

CASES = [('lower_kink', 'Lower kink', '#0072B2'),
         ('negative_rate_dividend', 'Negative rate/dividend', '#D55E00'),
         ('high_volatility', 'High volatility', '#009E73')]


def run(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    contracts = dict(synthetic_contracts())
    rows, profiles = [], []
    for name, label, colour in CASES:
        c = contracts[name]
        moments = hedge_moments(c, .06)
        cash_risk, _ = residual_square_integral(c, .06, np.zeros(2))
        unrestricted_risk, _ = residual_square_integral(c, .06,
                                            np.linalg.solve(moments.matrix, moments.cross))
        values, shadows = [], []
        grid = np.linspace(0., 4., 121)
        for kappa in grid:
            solution = solve_constrained_hedge(moments.matrix, moments.cross,
                          moments.constant, moments.costs, float(kappa*moments.premium))
            risk, error = residual_square_integral(c, .06, solution.theta)
            shadow = float(solution.multipliers[-1]*moments.premium/cash_risk)
            if kappa == 0:
                # The zero-budget multipliers are non-unique; plot the proved
                # right derivative rather than an arbitrary representative.
                shadow = float(2*max(0., *(moments.cross/moments.costs))
                                     *moments.premium/cash_risk)
            values.append(risk/cash_risk)
            shadows.append(shadow)
            rows.append({'case': name, 'physical_total_return': .06,
                         'budget_multiplier': float(kappa), 'population_loss_integral': risk,
                         'quadrature_error_estimate': error, 'cash_risk': cash_risk,
                         'normalised_risk': risk/cash_risk, 'normalised_shadow': shadow,
                         'unrestricted_risk_ratio': unrestricted_risk/cash_risk,
                         'arithmetic_status': 'floating-point illustration, not a certificate'})
        profiles.append((label, colour, grid, values, shadows, unrestricted_risk/cash_risk))
    with (output/'budget_frontier.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(9., 3.7), layout='constrained')
    for label, colour, grid, values, shadows, relaxed in profiles:
        axes[0].plot(grid, values, color=colour, linewidth=1.8, label=label)
        axes[0].axhline(relaxed, color=colour, linestyle='--', linewidth=1.)
        axes[1].plot(grid, shadows, color=colour, linewidth=1.8)
    axes[0].axhline(1., color='#666666', linestyle=':', linewidth=.9)
    axes[0].set(ylabel=r'$v(\kappa V_0)/L(0)$', title='(a) Risk under funding limits', ylim=(0, 1.05))
    axes[1].set(ylabel=r'$\lambda_B V_0/L(0)$', title='(b) Value of additional budget')
    axes[0].legend(frameon=False, fontsize=8.5)
    for ax in axes:
        ax.set(xlabel=r'Instrument budget/premium, $\kappa$', xlim=(0, 4.))
        ax.grid(color='#dddddd', linewidth=.6)
    fig.savefig(output/'budget_frontier.png', dpi=240)
    fig.savefig(output/'budget_frontier.svg')
    plt.close(fig)
    tex = [r'\begin{figure}[tb]', r'\centering']
    tex += [r'\definecolor{hedge%d}{HTML}{%s}'%(i,item[1][1:]) for i,item in enumerate(profiles)]
    tex += [r'\begin{tikzpicture}',
           r'\begin{groupplot}[group style={group size=2 by 1,horizontal sep=13mm},',
           r'width=0.45\textwidth,height=5.6cm,xmin=0,xmax=4,grid=major,',
           r'grid style={gray!25},tick label style={font=\scriptsize},',
           r'label style={font=\small},title style={font=\small},',
           r'xlabel={Budget/premium $\kappa$}]',
           r'\nextgroupplot[title={(a) Risk under funding limits},ylabel={$v(\kappa V_0)/L(0)$},',
           r'ymin=0,ymax=1.05,legend style={font=\scriptsize,draw=none,fill=none,at={(.98,.98)},anchor=north east}]']
    def coordinates(grid, values):
        return ' '.join(f'({x:.8g},{y:.10g})' for x,y in zip(grid,values))
    for i, (label, colour, grid, values, shadows, relaxed) in enumerate(profiles):
        tex.append(r'\addplot[color=hedge%d,thick] coordinates {%s};'%(i,coordinates(grid,values)))
        tex.append(r'\addlegendentry{%s}'%label)
        tex.append(r'\addplot[color=hedge%d,dashed,forget plot] coordinates {(0,%.10g)(4,%.10g)};'%(i,relaxed,relaxed))
    tex.append(r'\addplot[gray,dotted,forget plot] coordinates {(0,1)(4,1)};')
    tex.append(r'\nextgroupplot[title={(b) Value of additional budget},ylabel={$\lambda_B V_0/L(0)$},ymin=0]')
    for i, (label, colour, grid, values, shadows, relaxed) in enumerate(profiles):
        tex.append(r'\addplot[color=hedge%d,thick] coordinates {%s};'%(i,coordinates(grid,shadows)))
    tex += [r'\end{groupplot}', r'\end{tikzpicture}',
        r'\caption{Funding and residual risk for three registered contractual scenarios at $\mu=0.06$. Solid curves use the constrained mean-square hedge; dashed horizontal levels give each two-instrument unrestricted relaxation. The dotted level is the cash-only baseline. The right panel shows the normalised budget multiplier, equal to the marginal reduction in normalised risk on regular active-set intervals. At zero I use the proved right derivative. A supplemental 121-point grid illustrates curvature and saturation; all seven contracts remain in the primary comparison. These are synthetic floating-point diagnostics.}',
        r'\label{fig:h-budget-frontier}', r'\end{figure}']
    (output/'hedging_figure.tex').write_text('\n'.join(tex)+'\n')
    (output/'figure_manifest.json').write_text(json.dumps({'cases':[x[0] for x in CASES],
          'physical_total_return':.06,'grid_points_per_case':121,'rows':len(rows),
          'purpose':'Illustrate risk/funding trade-off and dual marginal value; no 3D or decorative figure.',
          'scope':'Supplemental to preregistered126rows; not an additional comparative outcome claim.'},indent=2)+'\n')

if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output)
