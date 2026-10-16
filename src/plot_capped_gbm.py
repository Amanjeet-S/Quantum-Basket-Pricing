"""Draw synthetic convergence evidence and self-contained LaTeX vector panels."""
from pathlib import Path
import argparse,csv,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLOURS={'lower_kink':'#0072B2','upper_kink':'#D55E00','low_vol_short_dated':'#009E73','high_volatility':'#CC79A7'}
LABELS={'lower_kink':'Lower kink','upper_kink':'Upper kink','low_vol_short_dated':'Short dated','high_volatility':'High volatility'}
def run(source,output):
    output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader(source.open()))
    price=[]
    for case in COLOURS:
        r=[x for x in rows if x['case']==case and x['quantity']=='price']
        price.append((LABELS[case],COLOURS[case],[(int(x['space_intervals']),abs(float(x['total_error']))) for x in r]))
    r=[x for x in rows if x['case']=='low_vol_short_dated' and x['quantity']=='gamma']
    gamma=[('Total error','#0072B2',[(int(x['space_intervals']),abs(float(x['total_error']))) for x in r]),('Bump bias','#D55E00',[(int(x['space_intervals']),abs(float(x['analytic_bump_bias']))) for x in r]),('PDE difference','#009E73',[(int(x['space_intervals']),abs(float(x['pde_bumped_minus_analytic_bumped']))) for x in r])]
    fig,ax=plt.subplots(1,2,figsize=(8.7,3.45),layout='constrained')
    for axis,series,title,ylabel in zip(ax,[price,gamma],['(a) Price convergence','(b) Short-dated gamma decomposition'],['Absolute price error','Absolute gamma error']):
        for label,col,points in series:
            axis.loglog(*zip(*points),color=col,marker='o',markersize=3.4,label=label,linewidth=1.5)
        axis.set_title(title,fontsize=10);axis.set_xlabel('Space intervals M');axis.set_ylabel(ylabel)
        axis.grid(True,which='major',color='#ddd',linewidth=.6);axis.legend(fontsize=8,frameon=False)
        axis.tick_params(labelsize=8)
    fig.savefig(output/'capped_gbm_convergence.png',dpi=220)
    fig.savefig(output/'capped_gbm_convergence.svg')
    fig.savefig(output/'capped_gbm_convergence.pdf')
    plt.close(fig)
    # The document compiler accepts a standalone source: no external figure input.
    lines=[r'\begin{figure}[tb]',r'\centering',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(164,71)']
    def panel(series,x0,title,ylow,yhigh):
        xleft,xright=7,14 # log2 interval count; the registered meshes span these.
        bottom,width,height=22,61,38
        lines.append(r'\put(%s,62){\makebox(61,5)[c]{\small %s}}'%(x0,title))
        for power in range(7,15):
            xx=x0+(power-xleft)*width/(xright-xleft)
            lines.append(r'\put(%.3f,%s){\line(0,1){%s}}'%(xx,bottom,height) if power in [7,14] else r'{\color[HTML]{E5E5E5}\put(%.3f,%s){\line(0,1){%s}}}'%(xx,bottom,height))
            if power in [7,9,11,13,14]:lines.append(r'\put(%.3f,17){\makebox(0,3){\scriptsize %s}}'%(xx,2**power))
        for power in range(ylow,yhigh+1):
            yy=bottom+(power-ylow)*height/(yhigh-ylow)
            lines.append(r'{\color[HTML]{E5E5E5}\put(%.3f,%.3f){\line(1,0){%s}}}'%(x0,yy,width))
            lines.append(r'\put(%.3f,%.3f){\makebox(0,3)[r]{\scriptsize $10^{%s}$}}'%(x0-2,yy-1.5,power))
        lines.append(r'\put(%s,13){\makebox(61,3)[c]{\scriptsize Space intervals $M$}}'%x0)
        for index,(label,col,points) in enumerate(series):
            transformed=[(x0+(math.log2(x)-xleft)*width/(xright-xleft),bottom+(math.log10(max(y,1e-20))-ylow)*height/(yhigh-ylow)) for x,y in points]
            for (x1,y1),(x2,y2) in zip(transformed,transformed[1:]):
                lines.append(r'{\color[HTML]{%s}\linethickness{0.5pt}\qbezier(%.3f,%.3f)(%.3f,%.3f)(%.3f,%.3f)}'%(col[1:],x1,y1,(x1+x2)/2,(y1+y2)/2,x2,y2))
            for xx,yy in transformed:lines.append(r'{\color[HTML]{%s}\put(%.3f,%.3f){\circle*{0.8}}}'%(col[1:],xx,yy))
            # A compact two-row legend, outside the plotting region.
            lx=x0+32*(index%2);ly=8-4*(index//2)
            lines.append(r'{\color[HTML]{%s}\put(%.3f,%s){\line(1,0){4}}}\put(%.3f,%s){\makebox(0,3)[l]{\scriptsize %s}}'%(col[1:],lx,ly,lx+5,ly-1.5,label))
    pmin=math.floor(math.log10(min(p[1] for _,_,xs in price for p in xs)));pmax=math.ceil(math.log10(max(p[1] for _,_,xs in price for p in xs)))
    gmin=math.floor(math.log10(min(p[1] for _,_,xs in gamma for p in xs)));gmax=math.ceil(math.log10(max(p[1] for _,_,xs in gamma for p in xs)))
    panel(price,14,r'(a) Absolute price error',pmin,pmax)
    panel(gamma,100,r'(b) Absolute gamma error',gmin,gmax)
    lines += [r'\end{picture}',r'\caption{Synthetic convergence evidence. Panel (a) contrasts the lower and upper payoff kinks, a short-dated low-volatility claim and a high-volatility claim. Panel (b) separates the short-dated gamma total error into analytic bump bias and the PDE difference at the same bump. Absolute component errors need not add because their signed contributions can cancel. Both axes are logarithmic. All seven contracts and all five Greeks are retained in the comparison files.}',r'\label{fig:cgbm-convergence}',r'\end{figure}']
    # Picture mode preserves inter-command glue: suppress every line ending
    # inside the drawing to prevent cumulative horizontal origin shifts.
    start=lines.index(r'\begin{picture}(164,71)');end=lines.index(r'\end{picture}')
    for j in range(start,end+1):lines[j]+='%'
    (output/'convergence_figure.tex').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
