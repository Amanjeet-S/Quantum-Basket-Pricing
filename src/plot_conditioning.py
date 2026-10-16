"""Plot the synthetic retained-variance diagnostics and standalone TeX drawing."""
from pathlib import Path
import argparse,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES={'lower_kink':('Lower kink','#0072B2'),'deep_cap':('Deep cap','#D55E00'),'low_vol_short_dated':('Short dated','#009E73'),'high_volatility':('High volatility','#CC79A7')}
def run(source,output):
    output.mkdir(parents=True,exist_ok=True);rows=list(csv.DictReader(source.open()));curves=[]
    for case,(label,colour) in CASES.items():
        points=[(0.,0.)]+[(float(r['conditioning_fraction']),float(r['conditional_variance'])/float(r['raw_variance'])) for r in rows if r['case']==case]+[(1.,1.)]
        curves.append((label,colour,sorted(points)))
    fig,ax=plt.subplots(figsize=(6.4,3.55),layout='constrained')
    for label,col,points in curves:ax.plot(*zip(*points),marker='o',markersize=4,linewidth=1.6,color=col,label=label)
    ax.set(xlabel='Observation fraction t/T',ylabel='Retained price variance Var(m_t)/Var(Y)',xlim=(0,1),ylim=(0,1));ax.set_xticks([0,.25,.5,.75,1]);ax.grid(color='#ddd',linewidth=.6);ax.legend(fontsize=9,frameon=False,ncols=2)
    fig.savefig(output/'conditioning_variance.png',dpi=220);fig.savefig(output/'conditioning_variance.svg');plt.close(fig)
    lines=[r'\begin{figure}[tb]',r'\centering',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(160,66)%']
    x0,y0,w,h=20,21,115,36
    for tick in [0,.25,.5,.75,1]:
        x=x0+w*tick
        lines.append(r'{\color[HTML]{E5E5E5}\put(%.3f,%s){\line(0,1){%s}}}\put(%.3f,16){\makebox(0,3){\scriptsize %g}}%%'%(x,y0,h,x,tick))
        y=y0+h*tick
        lines.append(r'{\color[HTML]{E5E5E5}\put(%s,%.3f){\line(1,0){%s}}}\put(17,%.3f){\makebox(0,3)[r]{\scriptsize %g}}%%'%(x0,y,w,y-1.5,tick))
    lines.extend([r'\put(20,60){\makebox(115,4){\small Retained variance $\operatorname{Var}(m_t)/\operatorname{Var}(Y)$}}%',r'\put(20,12){\makebox(115,3){\scriptsize Observation fraction $t/T$}}%'])
    for j,(label,col,points) in enumerate(curves):
        points=[(x0+w*x,y0+h*y) for x,y in points]
        for (a,b),(c,d) in zip(points,points[1:]):lines.append(r'{\color[HTML]{%s}\linethickness{0.5pt}\qbezier(%.3f,%.3f)(%.3f,%.3f)(%.3f,%.3f)}%%'%(col[1:],a,b,(a+c)/2,(b+d)/2,c,d))
        for x,y in points:lines.append(r'{\color[HTML]{%s}\put(%.3f,%.3f){\circle*{.8}}}%%'%(col[1:],x,y))
        lx=x0+60*(j%2);ly=7-4*(j//2)
        lines.append(r'{\color[HTML]{%s}\put(%s,%s){\line(1,0){5}}}\put(%s,%.3f){\makebox(0,3)[l]{\scriptsize %s}}%%'%(col[1:],lx,ly,lx+6,ly-1.5,label))
    lines.extend([r'\end{picture}',r'\caption{Retained conditional-price variance for selected synthetic contracts. Interior markers use the one-dimensional moment diagnostics; endpoints follow the analytic model. Lines join evaluated points. Increasing the observation time retains more terminal-price randomness: estimator variance rises while payoff-prediction loss falls. These are distinct quantities. All seven contracts remain in the comparison files; the figure reports no computational speed or quantum advantage.}',r'\label{fig:cgbm-conditioning}',r'\end{figure}'])
    (output/'conditioning_figure.tex').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
