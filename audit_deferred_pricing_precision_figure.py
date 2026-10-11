"""Independent source-to-native mapping audit for Figure 2.

This audit reads the immutable resource replay directly, without importing
the figure producer. It verifies data and native source, not rendered pixels.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from itertools import product
import argparse, gzip, hashlib, json, math, re

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def three_places(value):
    with localcontext() as context:
        context.prec=80
        return str((Decimal(value.numerator)/Decimal(value.denominator)).quantize(
            Decimal('0.001'),rounding=ROUND_HALF_EVEN))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--figure-dir',type=Path,required=True)
    args=parser.parse_args()
    source=json.loads(gzip.decompress(args.source.read_bytes()))
    base=args.figure_dir
    data=json.loads((base/'pricing_precision_plot_data.json').read_text())
    coordinates=json.loads((base/'pricing_precision_coordinates.json').read_text())
    native=(base/'pricing_precision_figure.tex').read_text()
    old=json.loads((base/'snapshot/pricing_precision_plot_data.json').read_text())
    old_native=(base/'snapshot/pricing_precision_figure.tex.txt').read_text()
    assert data['source_sha256']==sha(args.source)
    assert data['source_protocol_sha256']==source['protocol_sha256']
    assert data['source_producer_sha256']==source['producer_source_sha256']
    assert data['checked_proof_sha256']==source['checked_proof_sha256']
    assert data['source_variant']=='deferred_B'
    assert len(data['rows'])==4 and len(coordinates)==64
    assert data['colour_upper_ratio']==old['colour_upper_ratio']
    assert data['colour_clipped_cells']==0
    numeric=[x for x in source['records'] if x['method']=='numerical_conditional_price']
    rejected=[x for x in source['rejections'] if x['method']=='numerical_conditional_price']
    assert len(numeric)==104 and len(rejected)==12
    keys={(x['model_id'],x['configuration']['PA_bits'],x['configuration']['PG_bits']) for x in numeric+rejected}
    expected_keys=set(product([x['model_id'] for x in data['rows']],(6,8,10,12),(6,8,10,12)))
    assert len(keys)==64 and keys==expected_keys
    old_rejected={(row['model_id'],x['PA_bits'],x['PG_bits']) for row in old['rows'] for x in row['cells'] if not x['admitted']}
    current_rejected={(x['model_id'],x['configuration']['PA_bits'],x['configuration']['PG_bits']) for x in rejected}
    assert current_rejected==old_rejected
    palette={int(i):tuple(float(s) for s in rgb.split(',')) for i,rgb in re.findall(
        r'\\definecolor\{QBPBar(\d+)\}\{rgb\}\{([^}]+)\}',old_native)}
    assert set(palette)==set(range(60))
    assert [list(palette[i]) for i in range(60)]==data['colour_palette_RGB_knots']
    palette_native={int(i):tuple(float(s) for s in rgb.split(',')) for i,rgb in re.findall(
        r'\\definecolor\{QBPBar(\d+)\}\{rgb\}\{([^}]+)\}',native)}
    assert palette_native==palette
    native_colours={key:tuple(rgb.split(',')) for key,rgb in re.findall(
        r'\\definecolor\{(QBPCell\d+_\d+_\d+)\}\{rgb\}\{([^}]+)\}',native)}
    admitted=0;rejection_count=0;minimum_count=0;rows=[]
    for panel,row in enumerate(data['rows']):
        model=row['model_id']
        family=[x for x in numeric if x['model_id']==model]
        numerical_min=min(x['deferred_B']['T_full_alpha'] for x in family)
        assert numerical_min==row['base_T_count']
        assert {(x['PA_bits'],x['PG_bits']) for x in row['cells']}==set(product((6,8,10,12),repeat=2))
        outlined=[]
        for cell in row['cells']:
            a,g=cell['PA_bits'],cell['PG_bits']
            alternatives=[x for x in family if x['configuration']=={'PA_bits':a,'PG_bits':g}]
            x_mm=8+39*panel+7*((a-6)//2)
            y_mm=29+7*((g-6)//2)
            coordinate=next(x for x in coordinates if (x['model_id'],x['PA_bits'],x['PG_bits'])==(model,a,g))
            assert coordinate['x_mm']==x_mm and coordinate['y_mm']==y_mm
            if alternatives:
                assert len(alternatives)==2
                assert {x['architecture'] for x in alternatives}=={'exact7t','unitary8t'}
                cost=min(x['deferred_B']['T_full_alpha'] for x in alternatives)
                assert cell['admitted'] is True and cell['T_count']==cost
                exact=Fraction(cost,numerical_min)
                assert Fraction(cell['ratio_exact_fraction'])==exact
                assert cell['ratio_display']==three_places(exact)
                assert cell['architecture_caps']=={x['architecture']:x['deferred_B']['T_full_alpha'] for x in alternatives}
                chosen=next(x for x in alternatives if x['architecture']==cell['selected_architecture'])
                assert chosen['deferred_B']['T_full_alpha']==cost
                assert cell['selected_record_index']==chosen['record_index']
                assert cell['selected_record_sha256']==chosen['record_sha256']
                assert cell['deferred_profile_sha256']==chosen['deferred_profile_sha256']
                assert chosen['financial_target_prior_budget_phaseplan_unchanged'] is True
                result=chosen['deferred_B_full_result']
                assert Fraction(result['total_price_error_upper'])<=Fraction(1,100)
                assert Fraction(result['failure_probability_upper'])<=Fraction(1,3600)
                assert all(Fraction(flag['plan']['alpha'])==Fraction(1,3600) for flag in result['flags'])
                pos=math.log(float(exact))/math.log(data['colour_upper_ratio'])
                assert 0<=pos<=1
                scaled=60*pos-.5
                if scaled<=0:
                    colour=palette[0]
                elif scaled>=59:
                    colour=palette[59]
                else:
                    lo=int(scaled);weight=scaled-lo
                    colour=tuple(palette[lo][k]+weight*(palette[lo+1][k]-palette[lo][k]) for k in range(3))
                name=f'QBPCell{panel}_{a}_{g}'
                assert native_colours[name]==tuple(f'{v:.6f}' for v in colour)
                assert coordinate['RGB']==list(native_colours[name])
                is_minimum=cost==numerical_min
                assert cell['exact_numerical_winner']==is_minimum
                assert coordinate['winner']==is_minimum
                if is_minimum:
                    outlined.append((a,g))
                    assert native.count(rf'\put({x_mm},{y_mm}){{\framebox(7,7){{}}}}')==2
                    minimum_count+=1
                admitted+=1
            else:
                rejection=next(x for x in rejected if (x['model_id'],x['configuration'])==(model,{'PA_bits':a,'PG_bits':g}))
                assert cell['admitted'] is False and cell['ratio_display']=='X'
                assert cell['T_count'] is None and cell['ratio_exact_fraction'] is None
                assert cell['selected_record_sha256']==rejection['record_sha256']
                assert cell['rejection']==rejection['retained_rejection']
                assert coordinate['colour']=='QBPRejected' and coordinate['winner'] is False
                rejection_count+=1
            assert rf'\color{{{coordinate["colour"]}}}\put({x_mm},{y_mm}){{\rule{{7mm}}{{7mm}}}}' in native
            assert rf'\put({x_mm+3.5},{y_mm+3.5}){{\makebox(0,0){{{cell["ratio_display"]}}}}}' in native
        assert outlined==[(10,10)]
        event_records=[x for x in source['records'] if x['method']=='shared_threshold' and x['model_id']==model]
        assert len(event_records)==2
        event_min=min(x['deferred_B']['T_full_alpha'] for x in event_records)
        footer=Fraction(event_min,numerical_min)
        assert row['interval_T_count']==event_min
        assert Fraction(row['interval_ratio_exact_fraction'])==footer
        assert row['interval_ratio_display']==three_places(footer)
        assert 'Event $B$/min $B$='+three_places(footer) in native
        rows.append({'model_id':model,'numerical_B_minimum_T':numerical_min,
            'event_B_minimum_T':event_min,'event_B_over_numerical_B_exact':str(footer),
            'event_B_over_numerical_B_display':three_places(footer),'minimum':{'PA_bits':10,'PG_bits':10}})
    assert (admitted,rejection_count,minimum_count)==(52,12,4)
    assert len(native_colours)==52
    assert native.count(r'\framebox(7,7){}')==8
    assert r'\noindent\makebox[\linewidth][c]{%' in native
    assert r'\resizebox{0.96\linewidth}{!}{%' in native
    assert r'\begin{picture}(162,84)(-2,4)%' in native
    assert all(line.endswith('%') for line in native.splitlines()
        if line and not line.startswith(('%',r'\caption',r'\label',r'\end{figure}')))
    assert r'\includegraphics' not in native and r'\input{' not in native
    assert chr(8212) not in native
    report={'schema':'qbp-deferred-precision-native-mapping-audit-v1','status':'passed',
        'source_resource_sha256':sha(args.source),'plot_data_sha256':sha(base/'pricing_precision_plot_data.json'),
        'native_figure_sha256':sha(base/'pricing_precision_figure.tex'),'audit_source_sha256':sha(__file__),
        'configurations':64,'admitted_cells':52,'rejected_cells':12,'admitted_architecture_plans_checked':104,
        'exact_numerical_minima':4,'unchanged_rejection_identities':True,'ratio_clipping':False,
        'palette_reused_from_native_source':True,'colour_method':'Piecewise linear RGB interpolation at log-ratio positions over the retained sixty-sample viridis palette.',
        'exact_ratio_and_three_decimal_label_mapping':True,'all_native_cell_colours_and_coordinates_checked':True,
        'centring_and_existing_canvas_preserved':True,'native_external_file_dependency':False,
        'rows':rows,'native_pdf_compiled':False,'native_pdf_pixel_reviewed':False,
        'scope':'Independent direct source-to-data/native-source check. This verifies mathematical mapping and source geometry; it does not claim native PDF compilation or visual pixel inspection.'}
    (base/'pricing_precision_mapping_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    record_path=base/'pricing_precision_render_record.json'
    record=json.loads(record_path.read_text())
    assert record['data_sha256']==report['plot_data_sha256']
    assert record['native_figure_sha256']==report['native_figure_sha256']
    assert record['source_resource_sha256']==report['source_resource_sha256']
    record.update(source_mapping_audited_separately=True,
        independent_mapping_audit_sha256=sha(base/'pricing_precision_mapping_audit.json'),
        status='Native TeX and all source/data/colour mappings independently checked; compilation and native-PDF visual inspection are separate.')
    record_path.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'status':'passed','cells':64,'admitted':52,'rejected':12,'minima':4,
        'footer_ratios':[x['event_B_over_numerical_B_display'] for x in rows],
        'native_pdf_pixel_reviewed':False}))

if __name__=='__main__':
    main()
