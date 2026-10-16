"""Exact scientific data/state audit of the discrete 2D precision figure.

No plot-source function is imported. Every cell and winner is derived from
full registered resource records and compared to retained plot data. Display
rounding/colour coordinates are rendering properties, not statistical data.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--step',type=Path,required=True);ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();q=args.step/'adversarial/quantum';figure=args.step/'figure'
    source=q/'strengthened_baseline_results/strengthened_r1_resources.json'
    baseline=json.loads(source.read_text());plotfile=figure/'pricing_precision_plot_data.json';plot=json.loads(plotfile.read_text())
    cert=json.loads((args.step/'certificate/results/factor_basis_protocol.json').read_text())
    models={m['model_id']:m for m in cert['models']}
    assert plot['source_sha256']==sha(source)
    assert plot['source_protocol_sha256']==baseline['protocol_sha256']
    assert plot['source_producer_sha256']==baseline['source_sha256']
    assert plot['widths']==[6,8,10,12] and len(plot['rows'])==4
    assert plot['colour_transform']=='log2(T/numerical_minimum)'
    checked=admitted=rejected=winners=0;ratios=[];facts=[]
    for panel in plot['rows']:
        mid=panel['model_id'];selected=next(s for s in baseline['selected_comparisons'] if s['original_model_id']==mid)
        model=models[mid]
        assert panel['candidate_id']==selected['candidate_id']
        assert panel['assets']==model['assets'] and F(panel['lambda'])==F(model['lambda'])
        assert F(panel['scale'])==F(model['dispersion_scale'])
        cohort=[r for r in baseline['records'] if r['original_model_id']==mid]
        assert len(cohort)==16 and len(panel['cells'])==16
        live=[r for r in cohort if r.get('selected_strategy')]
        truebest=min(live,key=lambda r:(r['selected_strategy']['result']['T_count'],r['PA_bits'],r['PG_bits']))
        base=truebest['selected_strategy']['result']['T_count'];event=selected['fixed_threshold_selected_strategy']['result']['T_count']
        assert panel['base_T_count']==base and panel['interval_T_count']==event
        assert F(panel['interval_ratio_exact_fraction'])==F(event,base)
        assert panel['interval_ratio_display']==f'{float(F(event,base)):.3f}'
        assert {(c['PA_bits'],c['PG_bits']) for c in panel['cells']}==set(product((6,8,10,12),repeat=2))
        nwin=0
        for cell in panel['cells']:
            row=next(r for r in cohort if (r['PA_bits'],r['PG_bits'])==(cell['PA_bits'],cell['PG_bits']))
            ok=row.get('selected_strategy') is not None
            assert cell['admitted']==ok and cell['rejection']==row.get('rejection')
            if ok:
                cost=row['selected_strategy']['result']['T_count'];ratio=F(cost,base)
                assert cell['T_count']==cost and F(cell['ratio_exact_fraction'])==ratio>=1
                assert cell['ratio_display']==f'{float(ratio):.3f}'
                win=(row['PA_bits'],row['PG_bits'])==(truebest['PA_bits'],truebest['PG_bits'])
                assert cell['exact_numerical_winner']==win and win==(cost==base)
                admitted+=1;winners+=win;nwin+=win;ratios.append(ratio)
            else:
                assert cell['T_count'] is cell['ratio_exact_fraction'] is None
                assert cell['ratio_display']=='X' and cell['exact_numerical_winner'] is False
                assert row['rejection']['reason']=='Deterministic monetary tolerance exhausted'
                assert F(row['rejection']['deterministic_error_upper'])>=F(1,100)
                rejected+=1
            checked+=1
        assert nwin==1 and (truebest['PA_bits'],truebest['PG_bits'])==(10,10)
        facts.append({'model_id':mid,'minimum_T':base,'event_T':event,'event_to_minimum_exact':str(F(event,base)),
                      'admitted_cells':len(live),'rejected_cells':16-len(live),'exact_winner':[10,10]})
    assert checked==64 and admitted==plot['admitted_cells']==52 and rejected==plot['rejected_cells']==12 and winners==4
    maximum=max(ratios);expected_upper=max(2.25,math.ceil(float(maximum)*4)/4)
    assert plot['colour_upper_ratio']==expected_upper and F(str(expected_upper))>=maximum
    coordinate_file=figure/'pricing_precision_coordinates.json';coordinates_checked=0
    if coordinate_file.exists():
        coordinates=json.loads(coordinate_file.read_text());assert len(coordinates)==64
        for c in coordinates:
            i=next(i for i,p in enumerate(plot['rows']) if p['model_id']==c['model_id'])
            panel=plot['rows'][i];cell=next(k for k in panel['cells'] if (k['PA_bits'],k['PG_bits'])==(c['PA_bits'],c['PG_bits']))
            assert c['admitted']==cell['admitted'] and c['winner']==cell['exact_numerical_winner']
            assert c['ratio_exact_fraction']==cell['ratio_exact_fraction']
            assert c['x_mm']==8+39*i+7*[6,8,10,12].index(c['PA_bits'])
            assert c['y_mm']==29+7*[6,8,10,12].index(c['PG_bits'])
            assert c['colour']=='QBPRejected' if not c['admitted'] else c['colour'].startswith('QBPCell')
            coordinates_checked+=1
    render_file=figure/'pricing_precision_render_record.json'
    if render_file.exists():
        render=json.loads(render_file.read_text())
        assert render['data_sha256']==sha(plotfile) and render['source_sha256']==sha(args.step/'plot_pricing_precision.py')
        assert render['all64_registered_configurations_displayed'] and render['colour'] and render['3D_used'] is False
        assert render['native_external_file_dependency'] is False
    old_figure_data=args.repo/'figures/dispersion_separation_plot_data.json'
    assert json.loads((figure/'dispersion_separation_plot_data.json').read_text())==json.loads(old_figure_data.read_text())
    figure1render=json.loads((figure/'dispersion_separation_render_record.json').read_text())
    assert figure1render['source_sha256']==sha(args.step/'plot_dispersion_separation_v2.py')
    assert figure1render['data_sha256']==sha(figure/'dispersion_separation_plot_data.json')
    final_artifacts={name:sha(figure/name) for prefix in ('dispersion_separation','pricing_precision')
        for name in (prefix+'.png',prefix+'.svg',prefix+'_figure.tex',prefix+'_coordinates.json',prefix+'_plot_data.json',prefix+'_render_record.json')}
    result={'schema':'qbp-independent-pricing-precision-figure-audit-v1','status':'passed',
        'exact_cells_checked':checked,'admitted_cells':admitted,'rejected_cells':rejected,'exact_winners':winners,
        'native_coordinates_checked':coordinates_checked,'shared_colour_upper_ratio':plot['colour_upper_ratio'],
        'facts':facts,'plot_data_sha256':sha(plotfile),'baseline_resource_sha256':sha(source),
        'source_sha256':sha(__file__),
        'precision_plot_source_sha256':sha(args.step/'plot_pricing_precision.py'),
        'dispersion_plot_source_sha256':sha(args.step/'plot_dispersion_separation_v2.py'),
        'final_figure_artifact_sha256':final_artifacts,
        'figure1_numerical_data_exactly_preserved':True,
        'visual_review':{'surface':'Final exported PNGs inspected directly; native PDF layout remains outside this check.',
            'pricing_precision':'All Event/min footers now fully visible; all 64 cells and grey rejects readable; contrasting winner outlines clear; axes and colour legend legible.',
            'dispersion_separation':'Colour/marker/legend distinctions are clear; both panels have a specific certificate-comparison purpose; axes and titles fit without clipping.',
            '3D_disposition':'Rejected for this discrete precision pool because a surface would occlude cells or interpolate untested widths.'},
        'figure_choice':'A colour 2D discrete four-panel heatmap is justified by the complete registered rival pool and cost cliffs. A 3D surface would hide discrete cells and imply untested interpolation. Rounded labels accompany exact retained ratios; grey denotes failed monetary certificates, not missing cost observations. No cost confidence intervals are appropriate for deterministic ideal-logical caps.',
        'scope':'Exact data, admission, axis-coordinate and comparator mapping. Final exported PNGs have been visually reviewed, but this does not assert visual review of a native compiled PDF or infer theoretical priority/physical speed-up.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','exact_cells_checked','admitted_cells','rejected_cells','native_coordinates_checked')}))

if __name__=='__main__':main()
