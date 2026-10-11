"""Colour display of every registered numerical-rival precision configuration.

Integer costs and Fraction comparisons determine every state and minimum.
Floating point is used solely for drawing a shared logarithmic colour scale.
No cost confidence interval, continuous interpolation or hardware time is used.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, platform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

WIDTHS=(6,8,10,12)
COLOUR='viridis'
CAPTION = 'Complete T caps over the 64 registered numerical-rival configurations. Panels give asset count, common volatility $\\lambda$ and loading multiplier $c$. Columns and rows give the address widths $b_A$ and $b_G$ of $P(m_A)$ and $P(m_G)$; other widths and inputs are fixed. Colours and labels show each admitted cap/its panel minimum, using a shared logarithmic colour scale. Labels have three decimals; exact costs and ratios are retained. Grey crosses mark twelve rejected monetary budgets. Outlines identify the numerical minima, $b_A=b_G=10$. Orange footers compare the fixed shared-threshold encoding with these minima. Every alternative uses radius $0.01$ and per-price failure allowance $1/3600$. These are synthetic deterministic logical resource caps, not price observations or device timings.'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(path):
    d=json.loads(path.read_text());rows=[]
    assert len(d['records'])==64
    for selected in d['selected_comparisons']:
        mid=selected['original_model_id'];family=[r for r in d['records'] if r['original_model_id']==mid]
        assert len(family)==16
        base=selected['selected_numerical_strategy']['result']['T_count']
        assert base==min(r['selected_strategy']['result']['T_count'] for r in family if r.get('selected_strategy'))
        cells=[]
        for r in family:
            ok=r.get('selected_strategy') is not None
            cost=r['selected_strategy']['result']['T_count'] if ok else None
            ratio=F(cost,base) if ok else None
            winner=ok and (r['PA_bits'],r['PG_bits'])==(selected['selected_PA_bits'],selected['selected_PG_bits'])
            assert winner==(ok and cost==base)
            cells.append({'PA_bits':r['PA_bits'],'PG_bits':r['PG_bits'],
                'admitted':ok,'T_count':cost,'ratio_exact_fraction':str(ratio) if ok else None,
                'ratio_display':f'{float(ratio):.3f}' if ok else 'X','exact_numerical_winner':winner,
                'rejection':r.get('rejection')})
        assert {(r['PA_bits'],r['PG_bits']) for r in cells}==set((a,b) for a in WIDTHS for b in WIDTHS)
        event=selected['fixed_threshold_selected_strategy']['result']['T_count']
        rows.append({'model_id':mid,'candidate_id':selected['candidate_id'],
            'assets':2 if mid.startswith('n2') else 4 if mid.startswith('n4') else 8,
            'lambda':'0.4' if mid.startswith('n2') else '0.05' if '0.05' in mid else '0.2',
            'scale':'1.5' if 'over2' in mid else '3','base_T_count':base,
            'interval_T_count':event,'interval_ratio_exact_fraction':str(F(event,base)),
            'interval_ratio_display':f'{float(F(event,base)):.3f}','cells':cells})
    assert sum(not c['admitted'] for r in rows for c in r['cells'])==12
    top=max(float(F(c['ratio_exact_fraction'])) for r in rows for c in r['cells'] if c['admitted'])
    return {'schema':'qbp-pricing-precision-figure-data-v1','source_sha256':sha(path),
        'source_protocol_sha256':d['protocol_sha256'],'source_producer_sha256':d['source_sha256'],
        'widths':WIDTHS,'rows':rows,'colour_transform':'log2(T/numerical_minimum)',
        'colour_upper_ratio':max(2.25,math.ceil(top*4)/4),'admitted_cells':52,'rejected_cells':12,
        'figure_choice':'Post-result explanatory display of the entire registered pool, not a preregistered performance hypothesis.',
        'claims':'Exact deterministic logical caps with pricing confidence charged. No observed price or time, cost uncertainty or interpolation across widths.'}

def rgb(ratio,upper):
    return plt.get_cmap(COLOUR)(math.log2(ratio)/math.log2(upper))[:3]

def standard(d,out):
    plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],
        'mathtext.fontset':'dejavuserif','font.size':9,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,4,figsize=(10,3.75),layout='constrained')
    norm=matplotlib.colors.Normalize(0,math.log2(d['colour_upper_ratio']))
    cmap=plt.get_cmap(COLOUR).copy();cmap.set_bad('#E2E2E2')
    for ax,r in zip(axs,d['rows']):
        a=np.full((4,4),np.nan)
        for c in r['cells']:
            x=WIDTHS.index(c['PA_bits']);y=WIDTHS.index(c['PG_bits'])
            if c['admitted']:a[y,x]=math.log2(float(F(c['ratio_exact_fraction'])))
        ax.imshow(a,cmap=cmap,norm=norm,origin='lower',interpolation='none')
        for c in r['cells']:
            x=WIDTHS.index(c['PA_bits']);y=WIDTHS.index(c['PG_bits'])
            colour='black' if not c['admitted'] or a[y,x]>.62*norm.vmax else 'white'
            ax.text(x,y,c['ratio_display'],ha='center',va='center',color=colour,fontsize=8)
            if c['exact_numerical_winner']:
                ax.add_patch(matplotlib.patches.Rectangle((x-.5,y-.5),1,1,fill=False,edgecolor='white',linewidth=3.4))
                ax.add_patch(matplotlib.patches.Rectangle((x-.5,y-.5),1,1,fill=False,edgecolor='black',linewidth=2))
        ax.set_xticks(range(4),WIDTHS);ax.set_yticks(range(4),WIDTHS)
        ax.set_xlabel(r'$b_A$');ax.set_ylabel(r'$b_G$')
        ax.set_title(f"{r['assets']} assets\n"+rf"$\lambda={r['lambda']},c={r['scale']}$",fontsize=10)
        ax.text(.5,-.29,'Event/min='+r['interval_ratio_display'],transform=ax.transAxes,
            ha='center',va='top',fontsize=9,color='#D55E00')
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),ax=list(axs),
        location='top',fraction=.06,pad=.075,shrink=.8)
    ticks=[1,1.25,1.5,2];cb.set_ticks([math.log2(x) for x in ticks],labels=[str(x) for x in ticks])
    cb.set_label('Complete numerical T cap/numerical minimum',fontsize=10)
    fig.savefig(out/'pricing_precision.png',dpi=240,bbox_inches='tight',pad_inches=.1)
    fig.savefig(out/'pricing_precision.svg',bbox_inches='tight',pad_inches=.1)
    plt.close(fig)

def native(d,out):
    lines=[r'% Self-contained vector figure; xcolor and graphicx only.',r'\begin{figure}[!htb]',
        r'\centering\begingroup',r'\definecolor{QBPRejected}{HTML}{E2E2E2}',
        r'\definecolor{QBPInterval}{HTML}{D55E00}',r'\setlength{\unitlength}{1mm}',
        r'\noindent\makebox[\linewidth][c]{',r'\resizebox{0.96\linewidth}{!}{',r'\begin{picture}(162,84)(-2,4)',r'\fontsize{7}{9}\selectfont',
        r'\put(80,85){\makebox(0,0){Complete numerical T cap/numerical minimum}}']
    for j in range(60):
        ratio=d['colour_upper_ratio']**((j+.5)/60);c=rgb(ratio,d['colour_upper_ratio'])
        lines.append(rf'\definecolor{{QBPBar{j}}}{{rgb}}{{{c[0]:.6f},{c[1]:.6f},{c[2]:.6f}}}')
        lines.append(rf'\color{{QBPBar{j}}}\put({25+110*j/60:.5f},77){{\rule{{{110/60+.02:.5f}mm}}{{3mm}}}}')
    for z in (1,1.25,1.5,2):
        x=25+110*math.log2(z)/math.log2(d['colour_upper_ratio'])
        lines.append(rf'\color{{black}}\put({x:.5f},74){{\makebox(0,0){{{z}}}}}')
    coords=[]
    for i,r in enumerate(d['rows']):
        x0=8+39*i;y0=29;cell=7;centre=x0+14
        lines += [rf'\color{{black}}\put({centre},67){{\makebox(0,0){{{r["assets"]} assets}}}}',
            rf'\put({centre},62){{\makebox(0,0){{$\lambda={r["lambda"]},c={r["scale"]}$}}}}']
        for c in r['cells']:
            x=x0+cell*WIDTHS.index(c['PA_bits']);y=y0+cell*WIDTHS.index(c['PG_bits'])
            name=f'QBPCell{i}_{c["PA_bits"]}_{c["PG_bits"]}'
            if c['admitted']:
                col=rgb(float(F(c['ratio_exact_fraction'])),d['colour_upper_ratio'])
                lines.append(rf'\definecolor{{{name}}}{{rgb}}{{{col[0]:.6f},{col[1]:.6f},{col[2]:.6f}}}')
            else:name='QBPRejected'
            lines.append(rf'\color{{{name}}}\put({x},{y}){{\rule{{{cell}mm}}{{{cell}mm}}}}')
            text_colour='black' if not c['admitted'] or math.log2(float(F(c['ratio_exact_fraction'])))>.62*math.log2(d['colour_upper_ratio']) else 'white'
            lines.append(rf'\color{{{text_colour}}}\put({x+cell/2},{y+cell/2}){{\makebox(0,0){{{c["ratio_display"]}}}}}')
            if c['exact_numerical_winner']:
                lines += [r'\color{white}\linethickness{1.6pt}',
                    rf'\put({x},{y}){{\framebox({cell},{cell}){{}}}}',
                    r'\color{black}\linethickness{0.8pt}',
                    rf'\put({x},{y}){{\framebox({cell},{cell}){{}}}}']
            coords.append({'model_id':r['model_id'],'PA_bits':c['PA_bits'],'PG_bits':c['PG_bits'],
                'x_mm':x,'y_mm':y,'colour':name,'ratio_exact_fraction':c['ratio_exact_fraction'],
                'admitted':c['admitted'],'winner':c['exact_numerical_winner']})
        for j,w in enumerate(WIDTHS):
            lines += [rf'\color{{black}}\put({x0+cell*(j+.5)},26){{\makebox(0,0){{{w}}}}}',
                rf'\put({x0-1.5},{y0+cell*(j+.5)}){{\makebox(0,0)[r]{{{w}}}}}']
        lines += [rf'\put({centre},22){{\makebox(0,0){{$b_A$}}}}',
            rf'\put({x0-7},43){{\rotatebox{{90}}{{\makebox(0,0){{$b_G$}}}}}}',
            r'\color{QBPInterval}\put('+str(centre)+r',17){\makebox(0,0){Event/min='+r['interval_ratio_display']+'}}']
    lines += [r'\color{QBPRejected}\put(20,7){\rule{3mm}{3mm}}',
        r'\color{black}\put(24,8.5){\makebox(0,0)[l]{X: rejected budget}}',
        r'\linethickness{0.8pt}\put(73,7){\framebox(3,3){}}',
        r'\put(77,8.5){\makebox(0,0)[l]{Exact numerical minimum}}',r'\end{picture}',r'}',r'}\par',r'\endgroup',r'\small\setlength{\abovecaptionskip}{6pt}',
        r'\caption{'+CAPTION+'}',r'\label{fig:pricing-precision}',r'\end{figure}']
    (out/'pricing_precision_figure.tex').write_text('\n'.join(line if line.startswith('%') or line.startswith(r'\caption') or line.startswith(r'\label') or line.startswith(r'\end{figure}') else line+'%' for line in lines)+'\n')
    (out/'pricing_precision_coordinates.json').write_text(json.dumps(coords,indent=2)+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);d=read(a.source)
    (a.out/'pricing_precision_plot_data.json').write_text(json.dumps(d,indent=2)+'\n')
    standard(d,a.out);native(d,a.out)
    (a.out/'pricing_precision_render_record.json').write_text(json.dumps({
        'source_sha256':sha(__file__),'data_sha256':sha(a.out/'pricing_precision_plot_data.json'),
        'matplotlib_version':matplotlib.__version__,'python_version':platform.python_version(),
        'all64_registered_configurations_displayed':True,'colour':True,'3D_used':False,
        'native_figure_file':'pricing_precision_figure.tex','native_external_file_dependency':False,
        'caption':CAPTION},indent=2)+'\n')
    print('Rendered64 exact cost states,12 rejections and4 fixed interval comparisons.')

if __name__=='__main__':main()
