"""Exact retained allowances, Matplotlib output and self-contained LaTeX figure.

This is a display of two differently scoped certificates, not pricing error
or runtime. Exact fractions determine every admission; float logs draw points.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,math,platform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
COLOURS={'transfer':'#0072B2','classical':'#D55E00'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(source,selection):
 s=json.loads(source.read_text());p=json.loads(selection.read_text());rows=[]
 for scale in p['scales']:
  r=next(r for r in s['records'] if r['original_model_id']==p['model_family'] and r['cap'] is None and r['dispersion_scale']==scale)
  a=F(r['transfer_upper'])/F(1,200);b=F(r['degree3_intersection_radius_upper'])/F(1,100)
  assert a>0 and b>0
  assert (a<1)==r['transport_pass']
  assert (b<=1)==(r['degree3_intersection_first_successful_order'] is not None)
  rows.append({'model_id':r['model_id'],'dispersion_scale':scale,'x_float':float(F(scale)),
    'transfer_ratio_exact_fraction':str(a),'classical_ratio_exact_fraction':str(b),
    'transfer_pass':a<1,'classical_price_pass':b<=1})
 counts=[]
 for z in s['summary_by_scale']:
  q=[r for r in s['records'] if r['dispersion_scale']==z['dispersion_scale']]
  assert len(q)==36
  before=sum(r['separates_after_degree3'] for r in q);after=sum(r['separates_after_degree5'] for r in q)
  assert before==z['degree3_separators'] and after==z['remaining_degree5_separators'] and after<=before
  counts.append({'dispersion_scale':z['dispersion_scale'],'contracts':36,'through_degree3':before,'after_targeted_degree5':after})
 logs=[math.log10(float(F(r[k]))) for r in rows for k in ['transfer_ratio_exact_fraction','classical_ratio_exact_fraction']]
 limits=[min(-1,math.floor(min(logs))),max(1,math.ceil(max(logs)))]
 return {'schema':'qbp-dispersion-figure-data-v1','source_summary_sha256':sha(source),'source_full_result_sha256':s['full_result_sha256'],
  'experiment_protocol_sha256':s['protocol_sha256'],'selection_sha256':sha(selection),'representative_model':p['model_family'],
  'representative_payoff':'ordinary call','scale_values':p['scales'],'ratio_axis_log10_limits':limits,
  'threshold_exact_fraction':'1','transfer_allocation_exact_fraction':'1/200','price_allocation_exact_fraction':'1/100',
  'curves':rows,'counts':counts,'claims':'Bounds divided by their different allocations. Degree5 applied only to provisional separators. Lines are guides between tested scales. Not observed errors, complete quantum prices or speed ratios.',
  'rendering_arithmetic':'Only logarithms/plot coordinates are floating point. Scientific bounds, counts and classifications use retained exact rational endpoints.'}

def line(x1,y1,x2,y2):
 return rf'\qbezier({x1:.4f},{y1:.4f})({(x1+x2)/2:.4f},{(y1+y2)/2:.4f})({x2:.4f},{y2:.4f})'

def native(data,out):
 xmin,xmax=.85,3.15;ymin,ymax=data['ratio_axis_log10_limits'];height=53;bottom=23;left=17;width=63;right=103;barwidth=51
 countmax=max(2,max(r['through_degree3'] for r in data['counts'])+1)
 lines=[r'% Self-contained scientific figure: standard xcolor and graphicx only.',r'\begin{figure}[htbp]',r'\centering',r'\begingroup',r'\definecolor{QBPStressBlue}{HTML}{0072B2}',r'\definecolor{QBPStressOrange}{HTML}{D55E00}',r'\definecolor{QBPStressGrid}{HTML}{DDDDDD}',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(160,85)',r'\fontsize{8}{10}\selectfont',r'\color{black}',r'\put(48.5,82){\makebox(0,0){(a) Fixed eight-asset ordinary call}}',r'\put(128.5,82){\makebox(0,0){(b) Complete stress panel}}']
 coords=[]
 def xy(x,z):return left+width*(x-xmin)/(xmax-xmin),bottom+height*(z-ymin)/(ymax-ymin)
 for y in range(ymin,ymax+1):
  xx,yy=xy(xmin,y);lines += [r'\color{QBPStressGrid}\linethickness{0.18pt}',line(left,yy,left+width,yy),rf'\color{{black}}\put({left-2:.3f},{yy:.4f}){{\makebox(0,0)[r]{{$10^{{{y}}}$}}}}']
 lines += [r'\color{black}\linethickness{0.4pt}',line(left,bottom,left+width,bottom),line(left,bottom,left,bottom+height),r'\put(3,23){\rotatebox{90}{\makebox(53,0){Allowance / its allocation}}}']
 for scale in data['scale_values']:
  x=xy(float(F(scale)),ymin)[0];lab='1.5' if scale=='3/2' else scale
  lines += [line(x,bottom,x,bottom-1),rf'\put({x:.4f},20){{\makebox(0,0){{{lab}}}}}']
 lines += [r'\put(48.5,15){\makebox(0,0){Loading scale $c$}}',r'\color{black}\linethickness{0.5pt}']
 ty=xy(xmin,0)[1];x=left
 while x<left+width:lines.append(line(x,ty,min(x+1.8,left+width),ty));x+=2.8
 for key,colour,kind in [('transfer_ratio_exact_fraction','QBPStressBlue','o'),('classical_ratio_exact_fraction','QBPStressOrange','s')]:
  pts=[]
  for row in data['curves']:
   x,y=xy(row['x_float'],math.log10(float(F(row[key]))));pts.append((x,y));coords.append({'panel':'a','model_id':row['model_id'],'series':key,'x_mm':x,'y_mm':y,'exact_fraction':row[key]})
  lines.append(rf'\color{{{colour}}}\linethickness{{0.65pt}}')
  for a,b in zip(pts,pts[1:]):lines.append(line(*a,*b))
  for x,y in pts:
   if kind=='o':lines.append(rf'\put({x:.4f},{y:.4f}){{\circle*{{1.2}}}}')
   else:lines.append(rf'\put({x-.55:.4f},{y-.55:.4f}){{\rule{{1.1mm}}{{1.1mm}}}}')
 # Compact legends stay inside the empty upper-left portion of the panel.
 lines += [r'\color{QBPStressBlue}\put(21,70){\circle*{1.2}}',r'\color{black}\put(23,70){\makebox(0,0)[l]{Grid transfer}}',r'\color{QBPStressOrange}\put(20.45,65.45){\rule{1.1mm}{1.1mm}}',r'\color{black}\put(23,66){\makebox(0,0)[l]{Taylor $\leq3$ price interval}}']
 for v in range(countmax+1):
  y=bottom+height*v/countmax;lines += [r'\color{QBPStressGrid}\linethickness{0.18pt}',line(right,y,right+barwidth,y),rf'\color{{black}}\put({right-2:.3f},{y:.4f}){{\makebox(0,0)[r]{{{v}}}}}']
 lines += [r'\color{black}\linethickness{0.4pt}',line(right,bottom,right+barwidth,bottom),line(right,bottom,right,bottom+height),r'\put(89,23){\rotatebox{90}{\makebox(53,0){Separating contracts}}}']
 for i,row in enumerate(data['counts']):
  x=right+barwidth*(i+.5)/4;lab='1.5' if row['dispersion_scale']=='3/2' else row['dispersion_scale']
  for key,colour,offset in [('through_degree3','QBPStressBlue',-3.8),('after_targeted_degree5','QBPStressOrange',.3)]:
   h=height*row[key]/countmax;xx=x+offset
   if h:lines.append(rf'\color{{{colour}}}\put({xx:.4f},{bottom:.4f}){{\rule{{3.5mm}}{{{h:.4f}mm}}}}')
   coords.append({'panel':'b','dispersion_scale':row['dispersion_scale'],'series':key,'count':row[key],'x_mm':xx,'y_mm':bottom,'height_mm':h})
  lines.append(rf'\color{{black}}\put({x:.4f},20){{\makebox(0,0){{{lab}}}}}')
 lines += [r'\put(128.5,15){\makebox(0,0){Loading scale $c$}}',r'\color{QBPStressBlue}\put(106,71){\rule{2mm}{1mm}}',r'\color{black}\put(110,71.5){\makebox(0,0)[l]{After Taylor $\leq3$}}',r'\color{QBPStressOrange}\put(106,67){\rule{2mm}{1mm}}',r'\color{black}\put(110,67.5){\makebox(0,0)[l]{After targeted $\leq5$}}',r'\end{picture}',r'\endgroup',r'\caption{Different certificates under a controlled dispersion stress. Panel (a) uses the eight-asset ordinary call at $\lambda=0.2$, selected before outcomes. I divide the actual-grid financial-transfer allowance by $0.005$ and the complete Taylor-price intersection radius by $0.01$; the dashed line is admission at one. These are different error targets. Joined points guide the eye between tested scales. Panel (b) counts transfer-admitted contracts whose strengthened Taylor-price intervals fail, among all 36 contracts per scale. Degree five is evaluated only for provisional degree-three separators. All price intervals include lower-order intersections and intrinsic bounds. The figure reports certified allowances, not observed errors, complete quantum prices or a speed-up.}',r'\label{fig:dispersion-separation}',r'\end{figure}']
 (out/'dispersion_separation_figure.tex').write_text('\n'.join(lines)+'\n')
 (out/'dispersion_separation_coordinates.json').write_text(json.dumps({'mapping':{'xmin':xmin,'xmax':xmax,'log10_ymin':ymin,'log10_ymax':ymax,'left':left,'width':width,'bottom':bottom,'height':height,'right':right,'bar_width':barwidth,'count_ymax':countmax},'points':coords},indent=2)+'\n')

def standard(data,out):
 plt.rcParams.update({'font.family':'serif','font.serif':['DejaVu Serif'],'mathtext.fontset':'dejavuserif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
 fig,axes=plt.subplots(1,2,figsize=(8.7,3.9),gridspec_kw={'width_ratios':[1.2,1]})
 a,b=axes;x=[r['x_float'] for r in data['curves']]
 for key,c,mark,label in [('transfer_ratio_exact_fraction',COLOURS['transfer'],'o','Grid transfer'),('classical_ratio_exact_fraction',COLOURS['classical'],'s',r'Taylor $\leq3$ price interval')]:
  a.plot(x,[float(F(r[key])) for r in data['curves']],color=c,marker=mark,lw=1.5,ms=5,label=label)
 a.axhline(1,color='black',ls='--',lw=1);a.set_yscale('log');a.set_ylim(10**data['ratio_axis_log10_limits'][0],10**data['ratio_axis_log10_limits'][1]);a.set_xlim(.85,3.15);a.set_xticks(x,['1','1.5','2','3']);a.set_ylabel('Allowance / its allocation');a.set_xlabel(r'Loading scale $c$');a.grid(axis='y',color='#DDDDDD',lw=.5);a.set_axisbelow(True);a.legend(loc='upper left',fontsize=9,frameon=False);a.set_title('(a) Fixed eight-asset ordinary call',fontsize=11)
 xs=list(range(4));bar=.32
 b.bar([z-bar/2 for z in xs],[r['through_degree3'] for r in data['counts']],width=bar,color=COLOURS['transfer'],label=r'After Taylor $\leq3$')
 b.bar([z+bar/2 for z in xs],[r['after_targeted_degree5'] for r in data['counts']],width=bar,color=COLOURS['classical'],label=r'After targeted $\leq5$')
 ymax=max(2,max(r['through_degree3'] for r in data['counts'])+1);b.set_ylim(0,ymax);b.set_yticks(range(ymax+1));b.set_xticks(xs,['1','1.5','2','3']);b.set_xlabel(r'Loading scale $c$');b.set_ylabel('Separating contracts');b.grid(axis='y',color='#DDDDDD',lw=.5);b.set_axisbelow(True);b.legend(loc='upper left',fontsize=9,frameon=False);b.set_title('(b) Complete stress panel',fontsize=11)
 fig.tight_layout(w_pad=2.2);fig.savefig(out/'dispersion_separation.png',dpi=220);fig.savefig(out/'dispersion_separation.svg')
 plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=ROOT.parent/'results/dispersion_separation_summary.json');ap.add_argument('--selection',type=Path,default=ROOT.parent/'figure1_protocol.json');ap.add_argument('--out',type=Path,default=ROOT);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 data=read(args.source,args.selection);(args.out/'dispersion_separation_plot_data.json').write_text(json.dumps(data,indent=2)+'\n');standard(data,args.out);native(data,args.out)
 (args.out/'dispersion_separation_render_record.json').write_text(json.dumps({'source_sha256':sha(__file__),'data_sha256':sha(args.out/'dispersion_separation_plot_data.json'),'matplotlib_version':matplotlib.__version__,'python_version':platform.python_version(),'format':'PNG/SVG plus self-contained standard-LaTeX vector picture; no external graphic dependency for native compilation'},indent=2)+'\n')
 print('Saved eight exact-fraction curve points and eight complete-panel count bars.')
if __name__=='__main__':main()
