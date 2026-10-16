"""Matched estimator resource caps, exact rows and self-contained native picture."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import argparse,csv,hashlib,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_taylor_certificate import qline
ROOT=Path(__file__).resolve().parent
DEFAULT=ROOT.parent/'results/qpe_matched_comparison.json'
COLOURS={2:'#0072B2',4:'#D55E00',8:'#009E73'}
NAMES={2:'QBPQPEBlue',4:'QBPQPEOrange',8:'QBPQPEGreen'}
MARKERS={1:'o',10:'s',25:'^'}
SIZES={2:3.4,4:5.0,8:6.6}
XL=[9,13.5];YL=[8,13.5];XR=[4.7,13.3];YR=[-1.1,.04]
CAPTION=('Matched complete logical T caps for all 27 frozen synthetic contracts. I hold the financial circuit, comparison flag, probability prior, probability tolerance and 99 per cent per-price confidence fixed while changing the amplitude estimator. Left, on logarithmic axes: the canonical QPE amplitude-estimation cap against the retained iterative-controller cap, with equal-cap reference line. Right, with a logarithmic ratio axis: their exact T-cap ratio against the QPE phase-register width, with ratio-one reference. Colours identify asset count and marker shapes identify payoff cap H. Concentric marker sizes distinguish nearly coincident cases without moving any coordinate; sizes do not encode resource amounts. The exact case-labelled rows are retained separately. The QPE caps include the emitted controlled oracle, actual certified Clifford and T inverse Fourier tapes, confidence amplification and the compiled norm-error allowance rounded upwards on a 2^-80 grid. Phase-register width is not an exact increase in peak qubits: both peak-qubit figures are conservative bounds. These are proved resource upper caps rather than measured running times or a practical quantum advantage; the classical financial tolerance is already certified for these contracts.')

def read_exact(source):
    blob=source.read_bytes();raw=json.loads(blob);rows=[]
    for r in raw['configurations']:
        q=r['matched_QPE'];old=int(r['old_complete_T_cap']);new=int(q['T']);ratio=Fraction(new,old)
        assert ratio==Fraction(r['matched_T_ratio_exact_fraction'])
        assert q['prior']==r['shared_probability_prior'] and q['epsilon']==r['normalised_probability_tolerance']
        assert Fraction(q['alpha'])==Fraction(1,100)
        allowance=Fraction(q['plan']['qft_norm_error_allowance']);actual=Fraction(q['qft_norm_error_upper'])
        assert allowance>=actual and allowance-actual<Fraction(1,2**80)
        assert Fraction(q['plan']['binomial_failure_upper'])<=Fraction(1,100)
        row={'case_id':r['case_id'],'assets':r['assets'],'lambda_exact_decimal':r['lambda'],'cap':r['cap'],'matched_candidate_id':r['matched_candidate_id'],'matched_architecture':r['matched_architecture'],'matched_method':r['matched_method'],'matched_flag_kind':r['matched_flag_kind'],'matched_flag_c':r['matched_flag_c'],'shared_probability_prior':r['shared_probability_prior'],'shared_probability_tolerance':r['normalised_probability_tolerance'],'shared_alpha':'1/100','old_controller':r['old_controller'],'old_complete_T_cap':old,'matched_QPE_T_cap':new,'ratio_exact_fraction':str(ratio),'ratio_display_float':float(ratio),'phase_register_bits':q['phase_register_bits'],'old_peak_qubits_upper':r['old_peak_qubits_upper'],'QPE_peak_qubits_upper':q['peak_logical_qubits'],'difference_between_qubit_upper_bounds':q['peak_logical_qubits']-r['old_peak_qubits_upper'],'qft_norm_error_upper_exact_fraction':str(actual),'qft_norm_error_allowance_exact_fraction':str(allowance),'qft_T_per_run':q['qft_T_per_run'],'confidence_repetitions':q['repetitions'],'binomial_failure_upper_exact_fraction':q['plan']['binomial_failure_upper']}
        assert XL[0]<math.log10(old)<XL[1] and YL[0]<math.log10(new)<YL[1]
        assert XR[0]<q['phase_register_bits']<XR[1] and YR[0]<math.log10(float(ratio))<YR[1]
        rows.append(row)
    assert len(rows)==27 and len({r['case_id']for r in rows})==27
    rows.sort(key=lambda r:(r['assets'],float(r['lambda_exact_decimal']),r['cap']))
    return {'schema':'qbp-strategic-matched-qpe-figure-data-v1','source_file':str(source),'source_sha256':hashlib.sha256(blob).hexdigest(),'comparison_scope':'Each row holds the financial circuit, flag, proved prior, probability tolerance and per-price confidence fixed. Only the amplitude estimator is changed. Optimised flag searches are not plotted.','reference_exact_fraction':'1','colours_by_assets':COLOURS,'markers_by_cap':MARKERS,'markers_visibility_sizes_pt':SIZES,'jitter_applied':False,'phase_register_interpretation':'Width of the QPE phase register, not the exact increment in peak qubits; reported old and new peaks are conservative upper bounds.','axis_mapping':{'left_log10_x_range':XL,'left_log10_y_range':YL,'right_x_range':XR,'right_log10_y_range':YR},'rows':rows,'plotting_precision':'Integers and exact rational ratios establish all resource quantities. Float logarithms determine drawing coordinates only.'}

def draw_matplotlib(data,out):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'svg.fonttype':'none','figure.dpi':180})
    fig,axs=plt.subplots(1,2,figsize=(160/25.4,89/25.4))
    fig.subplots_adjust(left=.09,right=.985,bottom=.255,top=.89,wspace=.33)
    a,b=axs
    a.set_xscale('log');a.set_yscale('log');a.set_xlim(10**XL[0],10**XL[1]);a.set_ylim(10**YL[0],10**YL[1])
    a.set_xticks([10**i for i in range(9,14)]);a.set_yticks([10**i for i in range(8,14)])
    a.plot([10**XL[0],10**XL[1]],[10**XL[0],10**XL[1]],color='#666666',ls=(0,(4,2)),lw=.85,zorder=1)
    a.set_title('Complete T caps',pad=6);a.set_xlabel('Iterative T cap',labelpad=3);a.set_ylabel('Matched QPE T cap',labelpad=2)
    a.text(.54,.90,'Equal T cap',transform=a.transAxes,ha='left',va='top',fontsize=7,color='#555555')
    b.set_yscale('log');b.set_ylim(10**YR[0],10**YR[1]);b.set_xlim(*XR);b.set_xticks(range(5,14))
    b.set_yticks([.1,.2,.4,1]);b.set_yticklabels(['0.1','0.2','0.4','1'])
    b.axhline(1,color='#666666',ls=(0,(4,2)),lw=.85,zorder=1)
    b.set_title('Ratio versus phase width',pad=6);b.set_xlabel('QPE phase-register width',labelpad=3);b.set_ylabel('Matched QPE / iterative T cap',labelpad=2)
    for ax in axs:
        ax.minorticks_off();ax.grid(axis='y',color='#e6e6e6',lw=.5,zorder=0)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.tick_params(length=2.4,width=.6)
    for n in (8,4,2):
        for r in [r for r in data['rows']if r['assets']==n]:
            style={'marker':MARKERS[r['cap']],'markersize':SIZES[n],'mec':COLOURS[n],'mfc':'white','mew':.9,'ls':'none','zorder':10-n}
            a.plot([r['old_complete_T_cap']],[r['matched_QPE_T_cap']],**style)
            b.plot([r['phase_register_bits']],[r['ratio_display_float']],**style)
    assets=[Line2D([0],[0],ls='none',marker='o',markersize=SIZES[n],mec=COLOURS[n],mfc='white',mew=.9,label=f'{n} assets')for n in (2,4,8)]
    caps=[Line2D([0],[0],ls='none',marker=MARKERS[h],markersize=5,mec='#222222',mfc='white',mew=.9,label=rf'$H={h}$')for h in (1,10,25)]
    fig.legend(handles=assets,loc='lower center',bbox_to_anchor=(.29,.092),ncol=3,frameon=False,columnspacing=1.25,handlelength=1.0,handletextpad=.45,fontsize=7.4)
    fig.legend(handles=caps,loc='lower center',bbox_to_anchor=(.76,.092),ncol=3,frameon=False,columnspacing=1.25,handlelength=1.0,handletextpad=.45,fontsize=7.4)
    fig.text(.5,.027,'Exact coordinates; marker sizes separate near-coincident cases',ha='center',va='bottom',fontsize=7,color='#555555')
    fig.savefig(out/'qpe_matched_resources.svg',metadata={'Description':CAPTION,'Creator':'plot_qpe_matched_resources.py','Date':'2026-10-07'})
    fig.savefig(out/'qpe_matched_resources.png',dpi=300,metadata={'Description':CAPTION})
    plt.close(fig)

def hollow_marker(x,y,kind,colour,size):
    half=size/2
    if kind=='o':return rf'\color{{{colour}}}\linethickness{{0.5pt}}\put({x:.4f},{y:.4f}){{\circle{{{size:.4f}}}}}'
    if kind=='s':points=[(x-half,y-half),(x+half,y-half),(x+half,y+half),(x-half,y+half),(x-half,y-half)]
    else:points=[(x-half,y-half*.85),(x+half,y-half*.85),(x,y+half*.85),(x-half,y-half*.85)]
    return rf'\color{{{colour}}}\linethickness{{0.5pt}}'+'\n'+'\n'.join(qline(*a,*b)for a,b in zip(points,points[1:]))

def dashed(x1,y1,x2,y2):
    length=math.hypot(x2-x1,y2-y1);out=[];s=0
    while s<length:
        end=min(s+2.2,length);a=s/length;b=end/length
        out.append(qline(x1+a*(x2-x1),y1+a*(y2-y1),x1+b*(x2-x1),y1+b*(y2-y1)));s+=3.2
    return out

def draw_latex(data,out):
    Y0,H=23,54;panels=[{'x0':13,'width':62,'xs':XL,'ys':YL,'title':'Complete T caps'},{'x0':99,'width':59,'xs':XR,'ys':YR,'title':'Ratio versus phase width'}]
    lines=[r'% Self-contained picture. Requires only standard xcolor and graphicx.',r'\begin{figure}[htbp]',r'\centering',r'\begingroup',r'\definecolor{QBPQPEBlue}{HTML}{0072B2}',r'\definecolor{QBPQPEOrange}{HTML}{D55E00}',r'\definecolor{QBPQPEGreen}{HTML}{009E73}',r'\definecolor{QBPQPEGrey}{HTML}{666666}',r'\definecolor{QBPQPEGrid}{HTML}{E6E6E6}',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(160,89)',r'\fontsize{7.5}{9}\selectfont']
    coords=[]
    for ix,p in enumerate(panels):
        x0,W=p['x0'],p['width'];xs,ys=p['xs'],p['ys']
        def xy(x,y):return x0+W*(x-xs[0])/(xs[1]-xs[0]),Y0+H*(y-ys[0])/(ys[1]-ys[0])
        lines.append(rf'\color{{black}}\put({x0+W/2:.4f},81){{\makebox(0,0){{\fontsize{{9}}{{11}}\selectfont {p["title"]}}}}}')
        xticks=range(9,14)if ix==0 else range(5,14)
        yticks=list(range(8,14))if ix==0 else [math.log10(y)for y in [.1,.2,.4,1]]
        ylabs=[rf'$10^{{{y}}}$'for y in yticks]if ix==0 else ['0.1','0.2','0.4','1']
        for y,label in zip(yticks,ylabs):
            _,yy=xy(xs[0],y);lines.extend([r'\color{QBPQPEGrid}\linethickness{0.18pt}',qline(x0,yy,x0+W,yy),rf'\color{{black}}\put({x0-1:.4f},{yy:.4f}){{\makebox(0,0)[r]{{{label}}}}}'])
        lines.extend([r'\color{black}\linethickness{0.35pt}',qline(x0,Y0,x0+W,Y0),qline(x0,Y0,x0,Y0+H)])
        for tick in xticks:
            xx,_=xy(tick,ys[0]);label=rf'$10^{{{tick}}}$'if ix==0 else str(tick)
            lines.extend([qline(xx,Y0,xx,Y0-.65),rf'\put({xx:.4f},{Y0-2:.4f}){{\makebox(0,0){{{label}}}}}'])
        xlabel='Iterative T cap'if ix==0 else 'QPE phase-register width'
        ylabel='Matched QPE T cap'if ix==0 else 'Matched QPE / iterative T cap'
        lines.extend([rf'\put({x0+W/2:.4f},15){{\makebox(0,0){{{xlabel}}}}}',rf'\put({x0-12:.4f},{Y0:.4f}){{\rotatebox{{90}}{{\makebox({H},0){{{ylabel}}}}}}}'])
        lines.append(r'\color{QBPQPEGrey}\linethickness{0.5pt}')
        if ix==0:
            lines.extend(dashed(*xy(XL[0],XL[0]),*xy(XL[1],XL[1])))
            lines.append(rf'\put({x0+W*.54:.4f},{Y0+H*.90:.4f}){{\makebox(0,0)[l]{{Equal T cap}}}}')
        else:lines.extend(dashed(*xy(XR[0],0),*xy(XR[1],0)))
        for n in (8,4,2):
            for r in [r for r in data['rows']if r['assets']==n]:
                xx,yy=xy(math.log10(r['old_complete_T_cap']),math.log10(r['matched_QPE_T_cap']))if ix==0 else xy(r['phase_register_bits'],math.log10(r['ratio_display_float']))
                size=SIZES[n]*25.4/72
                lines.append(hollow_marker(xx,yy,MARKERS[r['cap']],NAMES[n],size))
                coords.append({'panel':'left'if ix==0 else 'right','case_id':r['case_id'],'x_mm':xx,'y_mm':yy,'marker_size_mm':size,'marker_kind':MARKERS[r['cap']],'colour_assets':n,'exact_x':str(r['old_complete_T_cap']if ix==0 else r['phase_register_bits']),'exact_y':str(r['matched_QPE_T_cap']if ix==0 else r['ratio_exact_fraction'])})
    for x,n in [(24,2),(48,4),(72,8)]:
        lines.extend([hollow_marker(x,7,'o',NAMES[n],SIZES[n]*25.4/72),rf'\color{{black}}\put({x+2.5:.4f},7){{\makebox(0,0)[l]{{{n} assets}}}}'])
    for x,h in [(106,1),(124,10),(143,25)]:
        lines.extend([hollow_marker(x,7,MARKERS[h],'black',1.76),rf'\color{{black}}\put({x+2.5:.4f},7){{\makebox(0,0)[l]{{$H={h}$}}}}'])
    lines.extend([r'\color{QBPQPEGrey}\put(80,1.5){\makebox(0,0){\fontsize{7}{8}\selectfont Exact coordinates; marker sizes separate near-coincident cases}}',r'\end{picture}',r'\endgroup',r'\caption{Matched complete logical T caps for all 27 frozen synthetic contracts. I hold the financial circuit, flag, proved probability prior, probability tolerance and 99 per cent per-price confidence fixed while changing the amplitude estimator. Left, on logarithmic axes: canonical QPE amplitude-estimation cap against the retained iterative cap, with equal-cap reference. Right, with a logarithmic ratio axis: their exact T-cap ratio against QPE phase-register width, with ratio-one reference. Colours identify asset count; shapes identify payoff cap $H$. Concentric marker sizes distinguish near-coincident cases without displacing coordinates and do not encode resource amounts. The exact case-labelled rows are retained separately. QPE caps include emitted controlled oracles, actual certified Clifford and T Fourier tapes, confidence amplification and the compiled norm allowance rounded upwards on a $2^{-80}$ grid. Register width is not an exact increment in peak qubits; both peak figures are conservative bounds. These are proved resource upper caps, not measured running times or practical quantum advantage. The classical financial tolerance is already certified for these contracts.}',r'\label{fig:matched-qpe-resources}',r'\end{figure}'])
    (out/'qpe_matched_resources_figure.tex').write_text('\n'.join(lines)+'\n')
    (out/'qpe_matched_resources_coordinates.json').write_text(json.dumps({'schema':'qbp-qpe-resource-native-coordinates-v1','picture_mm':[160,89],'plot_y_origin_mm':Y0,'plot_height_mm':H,'panels':panels,'points':coords,'jitter_applied':False},indent=2)+'\n')

def main():
    a=argparse.ArgumentParser();a.add_argument('--source',type=Path,default=DEFAULT);a.add_argument('--out',type=Path,default=ROOT);args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    data=read_exact(args.source);(args.out/'qpe_matched_resources_plot_data.json').write_text(json.dumps(data,indent=2)+'\n')
    fields=['case_id','assets','lambda_exact_decimal','cap','old_controller','old_complete_T_cap','matched_QPE_T_cap','ratio_exact_fraction','phase_register_bits','old_peak_qubits_upper','QPE_peak_qubits_upper','difference_between_qubit_upper_bounds']
    with(args.out/'qpe_matched_resources_case_rows.csv').open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(data['rows'])
    draw_matplotlib(data,args.out);draw_latex(data,args.out);(args.out/'qpe_matched_resources_caption.txt').write_text(CAPTION+'\n')
    print('All 27 matched cases retained at exact coordinates; 54 native-panel points.')
if __name__=='__main__':main()
