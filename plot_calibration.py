"""Illustrate the deterministic non-monotone volatility image."""
import argparse,json
from pathlib import Path
from dataclasses import replace
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from capped_gbm import Contract,analytic

def run(source,output):
 output.mkdir(parents=True,exist_ok=True);d=json.loads(source.read_text())['nonmonotone_price'];c=Contract(**d['contract']);lo,hi=d['volatility_lower'],d['volatility_upper'];peak=d['volatility_peak'];vol=np.linspace(lo,hi,81);value=np.array([analytic(replace(c,volatility=float(x)))['price'] for x in vol]);ymin=float(value.min())-.002;ymax=float(value.max())+.002
 fig,ax=plt.subplots(figsize=(6.4,3.45),layout='constrained');ax.plot(vol,value,color='#0072B2',linewidth=1.9,label='Capped price');ax.axhline(d['endpoint_only_upper'],color='#D55E00',linestyle='--',label='Endpoint-only upper');ax.scatter([peak],[d['true_price_at_peak']],color='#009E73',s=36,zorder=4,label='Interior maximum');ax.axvline(peak,color='#aaa',linewidth=.7);ax.set(xlabel='Volatility',ylabel='Model price',ylim=(ymin,ymax));ax.grid(color='#ddd',linewidth=.6);ax.legend(fontsize=9,frameon=False);fig.savefig(output/'volatility_price_image.png',dpi=220);fig.savefig(output/'volatility_price_image.svg');plt.close(fig)
 x0,y0,w,h=20,19,115,36;lines=[r'\begin{figure}[tb]',r'\centering',r'\setlength{\unitlength}{1mm}',r'\begin{picture}(160,66)%']
 def point(x,y):return(x0+(x-lo)*w/(hi-lo),y0+(y-ymin)*h/(ymax-ymin))
 for x in [lo,peak,hi]:
  xx,_=point(x,ymin);lines.append(r'{\color[HTML]{DDDDDD}\put(%.3f,%s){\line(0,1){%s}}}\put(%.3f,15){\makebox(0,3){\scriptsize %.4f}}%%'%(xx,y0,h,xx,x))
 for y in np.linspace(ymin,ymax,5):
  _,yy=point(lo,y);lines.append(r'{\color[HTML]{DDDDDD}\put(%s,%.3f){\line(1,0){%s}}}\put(18,%.3f){\makebox(0,3)[r]{\scriptsize %.3f}}%%'%(x0,yy,w,yy-1.5,y))
 pp=[point(x,y) for x,y in zip(vol,value)]
 for (a,b),(cc,dd) in zip(pp,pp[1:]):lines.append(r'{\color[HTML]{0072B2}\linethickness{0.6pt}\qbezier(%.3f,%.3f)(%.3f,%.3f)(%.3f,%.3f)}%%'%(a,b,(a+cc)/2,(b+dd)/2,cc,dd))
 _,yy=point(lo,d['endpoint_only_upper']);lines.append(r'{\color[HTML]{D55E00}\linethickness{0.6pt}\put(%s,%.3f){\line(1,0){%s}}}%%'%(x0,yy,w));xx,yy=point(peak,d['true_price_at_peak']);lines.append(r'{\color[HTML]{009E73}\put(%.3f,%.3f){\circle*{1.2}}}%%'%(xx,yy))
 lines.extend([r'\put(4,37){\makebox(0,0){\rotatebox{90}{\scriptsize Model price}}}%',r'\put(20,60){\makebox(115,4){\small A deterministic volatility-image counterexample}}%',r'\put(20,11){\makebox(115,3){\scriptsize Volatility $\sigma$}}%',r'{\color[HTML]{0072B2}\put(20,7){\line(1,0){5}}}\put(26,5.5){\makebox(0,3)[l]{\scriptsize Capped price}}%',r'{\color[HTML]{D55E00}\put(80,7){\line(1,0){5}}}\put(86,5.5){\makebox(0,3)[l]{\scriptsize Endpoint-only upper}}%',r'{\color[HTML]{009E73}\put(20,3){\circle*{.8}}}\put(26,1.5){\makebox(0,3)[l]{\scriptsize Interior maximum}}%',r'\end{picture}',r'\caption{A deterministic interval $[0.9\sigma_*,1.1\sigma_*]$ for the synthetic volatility-maximum contract. No confidence level is assigned to this illustration. Its maximum price is attained inside the interval, so the endpoint-only upper value misses a legitimate price. The repeated-data confidence experiment uses the chi-squared intervals separately.}',r'\label{fig:cgbm-volatility-image}',r'\end{figure}'])
 (output/'calibration_figure.tex').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
