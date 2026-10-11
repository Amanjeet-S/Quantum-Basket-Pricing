"""Retained exact radii, matching SVG/PNG and self-contained LaTeX picture."""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import argparse,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
DEFAULT=ROOT.parent/'results/taylor_pricing.json'
COLOURS={'.05':'#0072B2','.2':'#D55E00','.4':'#009E73'}
MARKERS={'.05':'o','.2':'s','.4':'^'}
LABELS={'.05':r'$\lambda=0.05$','.2':r'$\lambda=0.2$','.4':r'$\lambda=0.4$'}
YMIN,YMAX=-12,1
TAU=Fraction(1,100)

def read_exact(source):
    blob=source.read_bytes();raw=json.loads(blob)
    curves=[];aux=[]
    for c in raw['configurations']:
        lam=c['model_id'].split('_lambda')[1]
        lk={'.05':'.05','.2':'.2','.4':'.4','0.05':'.05','0.2':'.2','0.4':'.4'}[lam]
        points=[]
        for o in c['orders']:
            it=o['undiscounted_price_interval']
            lo=Fraction(it['lower_exact_fraction']);hi=Fraction(it['upper_exact_fraction']);radius=(hi-lo)/2
            assert radius>0
            certified=radius<=TAU
            assert certified==o['meets_undiscounted_radius_001']
            points.append({'degree':o['order'],'lower_exact_fraction':str(lo),'upper_exact_fraction':str(hi),'radius_exact_fraction':str(radius),'radius_display_float':float(radius),'log10_radius_display_float':math.log10(float(radius)),'certified':certified})
        first=next((p['degree']for p in points if p['certified']),None)
        assert first==c['first_certifying_order']
        if c['cap']==25:
            assert [p['degree']for p in points]==list(range(9))
            curves.append({'case_id':c['case_id'],'assets':c['assets'],'lambda_exact_decimal':lam,'lambda_key':lk,'cap':25,'first_certifying_degree':first,'points':points})
        if c['assets']==2 and lk=='.05' and c['cap']==10:
            aux.append({'case_id':c['case_id'],'first_certifying_degree':first,'minimum_radius_exact_fraction':str(min(Fraction(p['radius_exact_fraction'])for p in points))})
    assert len(curves)==9 and aux[0]['first_certifying_degree'] is None
    curves.sort(key=lambda c:(c['assets'],float(c['lambda_exact_decimal'])))
    for c in curves:
        assert all(YMIN<p['log10_radius_display_float']<YMAX for p in c['points'])
    return {'schema':'qbp-strategic-taylor-figure-data-v1','source_file':str(source),'source_sha256':hashlib.sha256(blob).hexdigest(),'threshold_exact_fraction':str(TAU),'radius_definition':'Half the retained undiscounted enclosure width, computed from exact rational interval endpoints','price_unit':raw['price_unit'],'data_scope':raw['data'],'plot_log10_range':[YMIN,YMAX],'cap':25,'degree_range':[0,8],'curves':curves,'caption_auxiliary_evidence':aux,'plotting_precision':'Exact rational radii and threshold determine certificate status. Floating point logarithms and coordinates are used only to draw the figure.'}

CAPTION='Finite-order Taylor price-enclosure radii for cap H=25 in the frozen synthetic two-, four- and eight-asset models. I calculate each radius from the exact retained undiscounted interval endpoints. The horizontal line is the certification threshold tau=0.01 in unspecified price currency units. The black ring marks the first integer degree that certifies this threshold. At lambda=0.05, no degree from zero to eight certifies the two-asset case; the corresponding H=10 case also fails this Taylor criterion and is resolved by the separately validated quadrature benchmark. Joined points are visual guides between tested integer degrees. These are enclosure radii, rather than observed pricing errors; the figure asserts neither convergence of an infinite Taylor series nor a practical quantum advantage.'

def draw_matplotlib(data,out):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'svg.fonttype':'none','figure.dpi':180})
    fig,axs=plt.subplots(1,3,figsize=(160/25.4,89/25.4),sharey=True)
    fig.subplots_adjust(left=.09,right=.99,bottom=.19,top=.89,wspace=.16)
    ticks=[-12,-10,-8,-6,-4,-2,0]
    for ax,n in zip(axs,(2,4,8)):
        ax.set_title(f'{n} assets',pad=6)
        ax.set_yscale('log');ax.set_ylim(10**YMIN,10**YMAX);ax.set_xlim(-.18,8.18)
        ax.set_xticks(range(9));ax.set_xlabel('Taylor degree',labelpad=3)
        ax.set_yticks([10**t for t in ticks]);ax.set_yticklabels([rf'$10^{{{t}}}$'for t in ticks])
        ax.minorticks_off();ax.grid(axis='y',color='#e6e6e6',lw=.5,zorder=0)
        ax.axhline(float(TAU),color='#666666',ls=(0,(4,2)),lw=.85,zorder=1)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.tick_params(length=2.4,width=.6)
        for c in [c for c in data['curves']if c['assets']==n]:
            lk=c['lambda_key'];p=c['points'];ax.plot([x['degree']for x in p],[x['radius_display_float']for x in p],color=COLOURS[lk],marker=MARKERS[lk],markersize=3,linewidth=1.15,zorder=3)
            d=c['first_certifying_degree']
            if d is not None:
                r=p[d]['radius_display_float'];ax.plot([d],[r],marker='o',markersize=7.2,mec='#222222',mfc='none',mew=.8,zorder=5)
        if n==2:
            ax.text(.965,.22,r'$\lambda=0.05$:'+ '\n'+'not certified',transform=ax.transAxes,ha='right',va='top',fontsize=7,color=COLOURS['.05'],bbox={'facecolor':'white','edgecolor':'none','pad':1.5})
        if n==8:ax.text(.98,(math.log10(float(TAU))-YMIN)/(YMAX-YMIN)+.012,r'$\tau=0.01$',transform=ax.transAxes,ha='right',va='bottom',fontsize=7,color='#555555')
    axs[0].set_ylabel('Enclosure radius (price units)',labelpad=2)
    legends=[Line2D([0],[0],color=COLOURS[k],marker=MARKERS[k],markersize=3,lw=1.15,label=LABELS[k])for k in COLOURS]
    legends.append(Line2D([0],[0],ls='none',marker='o',mfc='none',mec='#222222',markersize=6,label='First certifying degree'))
    fig.legend(handles=legends,loc='lower center',bbox_to_anchor=(.5,.025),ncol=4,frameon=False,columnspacing=1.6,handlelength=1.7,handletextpad=.5,fontsize=7.4)
    fig.savefig(out/'taylor_certificate_radii.svg',metadata={'Description':CAPTION,'Creator':'plot_taylor_certificate.py','Date':'2026-10-07'})
    fig.savefig(out/'taylor_certificate_radii.png',dpi=300,metadata={'Description':CAPTION})
    plt.close(fig)


def qline(x1,y1,x2,y2):
    return rf'\qbezier({x1:.4f},{y1:.4f})({(x1+x2)/2:.4f},{(y1+y2)/2:.4f})({x2:.4f},{y2:.4f})'
def marker(x,y,kind,colour,size=1.05):
    out=[rf'\color{{{colour}}}']
    if kind=='o':out.append(rf'\put({x:.4f},{y:.4f}){{\circle*{{{size:.3f}}}}}')
    elif kind=='s':out.append(rf'\put({x-size/2:.4f},{y-size/2:.4f}){{\rule{{{size:.3f}mm}}{{{size:.3f}mm}}}}')
    else:
        h=size*.85
        for i in range(9):
            yy=y-h/2+i*h/9;width=size*(1-i/9)
            out.append(rf'\put({x-width/2:.4f},{yy:.4f}){{\rule{{{width:.3f}mm}}{{{h/9+.01:.3f}mm}}}}')
    return '\n'.join(out)

def draw_latex(data,out):
    # The 160 by 89 mm picture is self-contained; no PGF, external graphics or input files.
    names={'.05':'QBPTaylorBlue','.2':'QBPTaylorOrange','.4':'QBPTaylorGreen'}
    W,H,Y0=43,60,17;lefts={2:11,4:62,8:113}
    coords=[]
    lines=[r'% Self-contained figure. Requires the standard xcolor and graphicx packages.',r'\begin{figure}[htbp]',r'\centering',r'\begingroup',r'\definecolor{QBPTaylorBlue}{HTML}{0072B2}',r'\definecolor{QBPTaylorOrange}{HTML}{D55E00}',r'\definecolor{QBPTaylorGreen}{HTML}{009E73}',r'\definecolor{QBPTaylorGrey}{HTML}{666666}',r'\definecolor{QBPTaylorGrid}{HTML}{E6E6E6}',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(160,89)',r'\fontsize{7.5}{9}\selectfont',r'\put(0,17){\rotatebox{90}{\makebox(60,0){Enclosure radius (price units)}}}']
    yticks=[-12,-10,-8,-6,-4,-2,0]
    for n,x0 in lefts.items():
        lines.append(rf'\color{{black}}\put({x0+W/2:.3f},81){{\makebox(0,0){{\fontsize{{9}}{{11}}\selectfont {n} assets}}}}')
        for t in yticks:
            y=Y0+H*(t-YMIN)/(YMAX-YMIN)
            lines.extend([r'\color{QBPTaylorGrid}\linethickness{0.18pt}',qline(x0,y,x0+W,y)])
            if n==2:lines.append(rf'\color{{black}}\put({x0-1:.4f},{y:.4f}){{\makebox(0,0)[r]{{$10^{{{t}}}$}}}}')
        lines.extend([r'\color{black}\linethickness{0.35pt}',qline(x0,Y0,x0+W,Y0),qline(x0,Y0,x0,Y0+H)])
        for d in range(9):
            x=x0+W*d/8
            lines.extend([qline(x,Y0,x,Y0-.65),rf'\put({x:.4f},{Y0-2:.4f}){{\makebox(0,0){{{d}}}}}'])
        lines.append(rf'\put({x0+W/2:.4f},10){{\makebox(0,0){{Taylor degree}}}}')
        ty=Y0+H*(-2-YMIN)/(YMAX-YMIN)
        lines.append(r'\color{QBPTaylorGrey}\linethickness{0.5pt}')
        xx=x0
        while xx<x0+W:
            lines.append(qline(xx,ty,min(xx+2.2,x0+W),ty));xx+=3.2
        if n==8:lines.append(rf'\put({x0+W-1:.4f},{ty+2:.4f}){{\makebox(0,0)[r]{{$\tau=0.01$}}}}')
        for c in [c for c in data['curves']if c['assets']==n]:
            lk=c['lambda_key'];colour=names[lk];points=[]
            for p in c['points']:
                x=x0+W*p['degree']/8;y=Y0+H*(p['log10_radius_display_float']-YMIN)/(YMAX-YMIN)
                points.append((x,y));coords.append({'case_id':c['case_id'],'degree':p['degree'],'radius_exact_fraction':p['radius_exact_fraction'],'x_mm':x,'y_mm':y})
            lines.append(rf'\color{{{colour}}}\linethickness{{0.65pt}}')
            for a,b in zip(points,points[1:]):lines.append(qline(*a,*b))
            for x,y in points:lines.append(marker(x,y,MARKERS[lk],colour))
            d=c['first_certifying_degree']
            if d is not None:
                x,y=points[d];lines.append(rf'\color{{black}}\linethickness{{0.45pt}}\put({x:.4f},{y:.4f}){{\circle{{2.3}}}}')
        if n==2:
            lines.extend([rf'\color{{QBPTaylorBlue}}\put({x0+W-1:.4f},34){{\makebox(0,0)[r]{{$\lambda=0.05$:}}}}',rf'\put({x0+W-1:.4f},30.5){{\makebox(0,0)[r]{{not certified}}}}'])
    for x,k in [(17,'.05'),(53,'.2'),(88,'.4')]:
        colour=names[k];lines.extend([rf'\color{{{colour}}}\linethickness{{0.65pt}}',qline(x,2,x+5,2),marker(x+2.5,2,MARKERS[k],colour),rf'\color{{black}}\put({x+7:.3f},2){{\makebox(0,0)[l]{{$\lambda={data["curves"][list(COLOURS).index(k)]["lambda_exact_decimal"]}$}}}}'])
    lines.extend([r'\color{black}\linethickness{0.45pt}\put(120,2){\circle{2.3}}',r'\put(123,2){\makebox(0,0)[l]{First certifying degree}}',r'\end{picture}',r'\endgroup',r'\caption{Finite-order Taylor price-enclosure radii for $H=25$ in the frozen synthetic two-, four- and eight-asset models. I calculate each radius from exact retained undiscounted interval endpoints. The horizontal line is $\tau=0.01$ in unspecified price currency units. The black ring marks the first certifying integer degree. At $\lambda=0.05$, no degree from zero to eight certifies the two-asset case; the corresponding $H=10$ case also fails this Taylor criterion and is resolved by separately validated quadrature. Joined points are visual guides between tested integer degrees. Enclosure radii are not observed pricing errors. The figure asserts neither infinite-series convergence nor practical quantum advantage.}',r'\label{fig:taylor-certificate-radii}',r'\end{figure}'])
    (out/'taylor_certificate_figure.tex').write_text('\n'.join(lines)+'\n')
    (out/'taylor_certificate_coordinates.json').write_text(json.dumps({'schema':'qbp-taylor-latex-coordinates-v1','mapping':{'picture_mm':[160,89],'plot_width_mm':W,'plot_height_mm':H,'y_origin_mm':Y0,'log10_y_min':YMIN,'log10_y_max':YMAX,'panel_x_origins_mm':lefts},'points':coords},indent=2)+'\n')


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=DEFAULT);p.add_argument('--out',type=Path,default=ROOT);args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    data=read_exact(args.source);(args.out/'taylor_certificate_plot_data.json').write_text(json.dumps(data,indent=2)+'\n')
    draw_matplotlib(data,args.out);draw_latex(data,args.out)
    (args.out/'taylor_certificate_caption.txt').write_text(CAPTION+'\n')
    print('Nine exact-radius curves, 81 integer-degree points, eight first-degree markers.')
if __name__=='__main__':main()
