"""Native colour display of the deferred-cleanup numerical precision pool.

The complete resource replay determines every integer cost, rejection and
minimum. Fraction arithmetic determines all ratios. Floating point is used
only for the shared logarithmic colour mapping. The retained sixty viridis
samples are RGB knots, linearly interpolated in logarithmic colour position;
the extreme half-bin positions use the nearest endpoint knot. This keeps
the existing colour scale without a new external plotting dependency.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import argparse, gzip, hashlib, json, math, platform, re

WIDTHS=(6,8,10,12)
METHOD='numerical_conditional_price'
VARIANT='deferred_B'
CAPTION = (
    'Complete T caps for the deferred-cleanup preparation $B$ over the 64 '
    'registered numerical-rival configurations. Panels give asset count, '
    'common volatility $\\lambda$ and loading multiplier $c$. Columns and '
    'rows give the address widths $b_A$ and $b_G$ of $P(m_A)$ and $P(m_G)$; '
    'other widths and inputs are fixed. Each cell uses the cheaper complete '
    '$B$ cap across the two retained coherent architectures. Colours and '
    'labels show its cap/the numerical $B$ minimum in that panel, using the '
    'same shared logarithmic colour scale. Labels have three decimals; exact '
    'costs and ratios are retained. Grey crosses mark the same twelve '
    'rejected monetary budgets. Outlines identify the numerical minima, '
    '$b_A=b_G=10$. Orange footers compare the shared-threshold $B$ encoding '
    'with these numerical $B$ minima; Table~\\ref{tab:shared-threshold} also '
    'includes the raw-payoff pool. The same established cancellation is '
    'applied to every representation, with each configuration retaining its '
    'financial certificate, information prior and phase plan. All alternatives '
    'use monetary radius $0.01$ and per-price failure allowance $1/3600$. '
    'The figure displays deterministic logical resource caps for synthetic '
    'models.'
)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def display(value):
    with localcontext() as ctx:
        ctx.prec=60
        x=Decimal(value.numerator)/Decimal(value.denominator)
        return format(x.quantize(Decimal('0.001'),rounding=ROUND_HALF_EVEN),'f')

def parse_palette(path):
    text=Path(path).read_text()
    entries=re.findall(r'\\definecolor\{QBPBar(\d+)\}\{rgb\}\{([^}]+)\}',text)
    palette={int(i):tuple(float(x) for x in rgb.split(',')) for i,rgb in entries}
    assert set(palette)==set(range(60))
    return [palette[i] for i in range(60)]

def colour_at(position,palette):
    if position<=.5/60:
        return palette[0]
    if position>=59.5/60:
        return palette[-1]
    q=position*60-.5
    j=math.floor(q)
    w=q-j
    return tuple((1-w)*palette[j][k]+w*palette[j+1][k] for k in range(3))

def rgb(value,upper,palette):
    return colour_at(math.log2(float(value))/math.log2(upper),palette)

def read(path,old_data_path,palette_path):
    path=Path(path)
    d=json.loads(gzip.decompress(path.read_bytes()))
    old=json.loads(Path(old_data_path).read_text())
    assert d['status']=='completed retaining every registered outcome'
    assert d['registered_configurations']==416 and d['rejected_configurations']==12
    assert d['admitted_architecture_plans']==808
    assert len(d['selected_comparisons'])==4
    records=[r for r in d['records'] if r['method']==METHOD]
    rejected=[r for r in d['rejections'] if r['method']==METHOD]
    assert len(records)==104 and len(rejected)==12
    palette=parse_palette(palette_path)
    rows=[]
    for selected in d['selected_comparisons']:
        mid=selected['model_id']
        family=[r for r in records if r['model_id']==mid]
        rejected_family=[r for r in rejected if r['model_id']==mid]
        numerical=selected['variants'][VARIANT]['methods'][METHOD]
        event=selected['variants'][VARIANT]['methods']['shared_threshold']
        base=numerical['T_full_alpha']
        assert base==min(r[VARIANT]['T_full_alpha'] for r in family)
        assert numerical['canonical_configuration']=={'PA_bits':10,'PG_bits':10}
        cells=[]
        for a in WIDTHS:
            for g in WIDTHS:
                pair={'PA_bits':a,'PG_bits':g}
                alternatives=[r for r in family if r['configuration']==pair]
                rejection=[r for r in rejected_family if r['configuration']==pair]
                assert (len(alternatives),len(rejection)) in ((2,0),(0,1))
                if alternatives:
                    assert {r['architecture'] for r in alternatives}=={'exact7t','unitary8t'}
                    best=min(alternatives,key=lambda r:(r[VARIANT]['T_full_alpha'],r['architecture']))
                    cost=best[VARIANT]['T_full_alpha']
                    ratio=F(cost,base)
                    result=best['deferred_B_full_result']
                    assert result['T_count']==cost
                    assert best['financial_target_prior_budget_phaseplan_unchanged'] is True
                    assert F(result['total_price_error_upper'])<=F(1,100)
                    assert F(result['failure_probability_upper'])<=F(1,3600)
                    assert all(F(f['plan']['alpha'])==F(1,3600) for f in result['flags'])
                    winner=cost==base
                    assert winner==((a,g)==(10,10))
                    cells.append({'PA_bits':a,'PG_bits':g,'admitted':True,
                        'T_count':cost,'ratio_exact_fraction':str(ratio),
                        'ratio_display':display(ratio),'exact_numerical_winner':winner,
                        'selected_architecture':best['architecture'],
                        'selected_record_index':best['record_index'],
                        'source_row_index':best['source_row_index'],
                        'selected_record_sha256':best['record_sha256'],
                        'deferred_profile_sha256':best['deferred_profile_sha256'],
                        'architecture_caps':{r['architecture']:r[VARIANT]['T_full_alpha'] for r in alternatives},
                        'rejection':None})
                else:
                    r=rejection[0]
                    cells.append({'PA_bits':a,'PG_bits':g,'admitted':False,
                        'T_count':None,'ratio_exact_fraction':None,'ratio_display':'X',
                        'exact_numerical_winner':False,'selected_architecture':None,
                        'selected_record_index':r['record_index'],
                        'source_row_index':r['source_row_index'],
                        'selected_record_sha256':r['record_sha256'],
                        'deferred_profile_sha256':None,'architecture_caps':{},
                        'rejection':r['retained_rejection']})
        event_cost=event['T_full_alpha']
        event_alternatives=[r for r in d['records'] if r['method']=='shared_threshold' and r['model_id']==mid]
        assert len(event_alternatives)==2
        assert event_cost==min(r[VARIANT]['T_full_alpha'] for r in event_alternatives)
        rows.append({'model_id':mid,'assets':2 if mid.startswith('n2') else 4 if mid.startswith('n4') else 8,
            'lambda':'0.4' if mid.startswith('n2') else '0.05' if '0.05' in mid else '0.2',
            'scale':'1.5' if 'over2' in mid else '3','base_T_count':base,
            'numerical_minimum_architecture':numerical['canonical_architecture'],
            'numerical_minimum_record_index':numerical['canonical_record_index'],
            'interval_T_count':event_cost,'interval_ratio_exact_fraction':str(F(event_cost,base)),
            'interval_ratio_display':display(F(event_cost,base)),
            'event_architecture':event['canonical_architecture'],
            'event_record_index':event['canonical_record_index'],'cells':cells})
    actual_rejected={(r['model_id'],c['PA_bits'],c['PG_bits']) for r in rows for c in r['cells'] if not c['admitted']}
    old_rejected={(r['model_id'],c['PA_bits'],c['PG_bits']) for r in old['rows'] for c in r['cells'] if not c['admitted']}
    assert len(actual_rejected)==12 and actual_rejected==old_rejected
    top=max(F(c['ratio_exact_fraction']) for r in rows for c in r['cells'] if c['admitted'])
    upper=old['colour_upper_ratio']
    assert top<=F(str(upper))
    return {'schema':'qbp-pricing-precision-figure-data-v2','source_sha256':sha(path),
        'source_protocol_sha256':d['protocol_sha256'],'source_producer_sha256':d['producer_source_sha256'],
        'checked_proof_sha256':d['checked_proof_sha256'],'source_variant':VARIANT,
        'source_resource_filename':'deferred_cleanup_resources.json.gz',
        'widths':WIDTHS,'rows':rows,'colour_transform':'log2(T_B/numerical_B_minimum)',
        'colour_upper_ratio':upper,'maximum_ratio_exact_fraction':str(top),
        'colour_ticks':[1,1.25,1.5,2],'colour_clipped_cells':0,
        'colour_palette':'Retained sixty-sample viridis RGB knots, linearly interpolated in logarithmic colour position; endpoint half-bins use the nearest knot.',
        'colour_palette_RGB_knots':palette,'colour_palette_knot_positions':[(j+.5)/60 for j in range(60)],
        'colour_palette_source_sha256':sha(palette_path),'old_plot_data_sha256':sha(old_data_path),
        'admitted_cells':52,'rejected_cells':12,'architecture_selection':'Minimum complete full-alpha deferred-B cap at each admitted configuration across exact7t and unitary8t.',
        'footer_denominator':'Numerical deferred-B pool minimum only; excludes the raw-payoff minimum.',
        'figure_choice':'Post-result explanatory display of the complete registered numerical pool, preserving the existing two-dimensional layout.',
        'claims':'Exact deterministic logical caps with pricing confidence charged. Financial certificates, priors and phase plans are retained per configuration; no market data or device timing.'}

def native(d,out):
    palette=d['colour_palette_RGB_knots']
    upper=d['colour_upper_ratio']
    lines=[r'% Self-contained native vector figure; xcolor and graphicx only.',
        r'% Source variant: deferred_B; all integer caps and exact ratios are retained in the backing JSON.',
        r'% Retained sixty-sample viridis RGB knots; linear interpolation in logarithmic colour position.',
        r'\begin{figure}[!htb]',r'\centering\begingroup',
        r'\definecolor{QBPRejected}{HTML}{E2E2E2}',r'\definecolor{QBPInterval}{HTML}{D55E00}',
        r'\setlength{\unitlength}{1mm}',r'\noindent\makebox[\linewidth][c]{',
        r'\resizebox{0.96\linewidth}{!}{',r'\begin{picture}(162,84)(-2,4)',
        r'\fontsize{7}{9}\selectfont',
        r'\put(80,85){\makebox(0,0){Complete numerical $B$ T cap/numerical $B$ minimum}}']
    for j,col in enumerate(palette):
        lines.append(rf'\definecolor{{QBPBar{j}}}{{rgb}}{{{col[0]:.6f},{col[1]:.6f},{col[2]:.6f}}}')
        lines.append(rf'\color{{QBPBar{j}}}\put({25+110*j/60:.5f},77){{\rule{{{110/60+.02:.5f}mm}}{{3mm}}}}')
    for z in d['colour_ticks']:
        x=25+110*math.log2(z)/math.log2(upper)
        lines.append(rf'\color{{black}}\put({x:.5f},74){{\makebox(0,0){{{z}}}}}')
    coordinates=[]
    for i,row in enumerate(d['rows']):
        x0=8+39*i;y0=29;cell=7;centre=x0+14
        lines += [rf'\color{{black}}\put({centre},67){{\makebox(0,0){{{row["assets"]} assets}}}}',
            rf'\put({centre},62){{\makebox(0,0){{$\lambda={row["lambda"]},c={row["scale"]}$}}}}']
        for c in row['cells']:
            x=x0+cell*WIDTHS.index(c['PA_bits']);y=y0+cell*WIDTHS.index(c['PG_bits'])
            name=f'QBPCell{i}_{c["PA_bits"]}_{c["PG_bits"]}'
            col=rgb(F(c['ratio_exact_fraction']),upper,palette) if c['admitted'] else None
            if col is not None:
                lines.append(rf'\definecolor{{{name}}}{{rgb}}{{{col[0]:.6f},{col[1]:.6f},{col[2]:.6f}}}')
            else:
                name='QBPRejected'
            lines.append(rf'\color{{{name}}}\put({x},{y}){{\rule{{{cell}mm}}{{{cell}mm}}}}')
            text_colour='black' if not c['admitted'] or math.log2(float(F(c['ratio_exact_fraction'])))>.62*math.log2(upper) else 'white'
            lines.append(rf'\color{{{text_colour}}}\put({x+cell/2},{y+cell/2}){{\makebox(0,0){{{c["ratio_display"]}}}}}')
            if c['exact_numerical_winner']:
                lines += [r'\color{white}\linethickness{1.6pt}',rf'\put({x},{y}){{\framebox({cell},{cell}){{}}}}',
                    r'\color{black}\linethickness{0.8pt}',rf'\put({x},{y}){{\framebox({cell},{cell}){{}}}}']
            coordinates.append({'model_id':row['model_id'],'PA_bits':c['PA_bits'],'PG_bits':c['PG_bits'],
                'x_mm':x,'y_mm':y,'colour':name,'RGB':None if col is None else [f'{v:.6f}' for v in col],
                'ratio_exact_fraction':c['ratio_exact_fraction'],'admitted':c['admitted'],
                'winner':c['exact_numerical_winner'],'selected_record_index':c['selected_record_index']})
        for j,width in enumerate(WIDTHS):
            lines += [rf'\color{{black}}\put({x0+cell*(j+.5)},26){{\makebox(0,0){{{width}}}}}',
                rf'\put({x0-1.5},{y0+cell*(j+.5)}){{\makebox(0,0)[r]{{{width}}}}}']
        lines += [rf'\put({centre},22){{\makebox(0,0){{$b_A$}}}}',
            rf'\put({x0-7},43){{\rotatebox{{90}}{{\makebox(0,0){{$b_G$}}}}}}',
            r'\color{QBPInterval}\put('+str(centre)+r',17){\makebox(0,0){Event $B$/min $B$='+row['interval_ratio_display']+'}}']
    lines += [r'\color{QBPRejected}\put(20,7){\rule{3mm}{3mm}}',
        r'\color{black}\put(24,8.5){\makebox(0,0)[l]{X: rejected budget}}',
        r'\linethickness{0.8pt}\put(73,7){\framebox(3,3){}}',
        r'\put(77,8.5){\makebox(0,0)[l]{Exact numerical minimum}}',
        r'\end{picture}',r'}',r'}\par',r'\endgroup',r'\small\setlength{\abovecaptionskip}{6pt}',
        r'\caption{'+CAPTION+'}',r'\label{fig:pricing-precision}',r'\end{figure}']
    (out/'pricing_precision_figure.tex').write_text('\n'.join(line if line.startswith('%') or line.startswith(r'\caption') or line.startswith(r'\label') or line.startswith(r'\end{figure}') else line+'%' for line in lines)+'\n')
    (out/'pricing_precision_coordinates.json').write_text(json.dumps(coordinates,indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--previous-data',type=Path,required=True)
    parser.add_argument('--previous-native',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    d=read(args.source,args.previous_data,args.previous_native)
    (args.out/'pricing_precision_plot_data.json').write_text(json.dumps(d,indent=2)+'\n')
    native(d,args.out)
    (args.out/'pricing_precision_render_record.json').write_text(json.dumps({
        'generator_sha256':sha(__file__),'data_sha256':sha(args.out/'pricing_precision_plot_data.json'),
        'source_resource_sha256':d['source_sha256'],'native_figure_sha256':sha(args.out/'pricing_precision_figure.tex'),
        'python_version':platform.python_version(),'all64_registered_configurations_displayed':True,
        'colour':True,'3D_used':False,'native_figure_file':'pricing_precision_figure.tex',
        'native_external_file_dependency':False,'retained_canvas':'picture(162,84)(-2,4)',
        'retained_scale':'0.96 linewidth inside a centred linewidth makebox',
        'source_mapping_audited_separately':False,'native_pdf_compiled':False,
        'native_pdf_pixel_reviewed':False,'matplotlib_outputs_generated':False,
        'status':'Native TeX generated; compilation and native-PDF visual inspection are separate.',
        'caption':CAPTION},indent=2)+'\n')
    print(json.dumps({'displayed':64,'admitted':52,'rejected':12,
        'footer_ratios':[r['interval_ratio_display'] for r in d['rows']],
        'source_sha256':d['source_sha256'],'native_pdf_pixel_reviewed':False}))

if __name__=='__main__':
    main()
