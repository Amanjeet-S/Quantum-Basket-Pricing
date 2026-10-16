"""Independent closed-count verification of retained production block reuse."""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE if (HERE/'value_oracle.py').exists() else HERE.parent
sys.path.insert(0,str(ROOT))
import value_oracle as vo


def closed_predicate(q,kind,coordinate):
    if kind=='bucket':
        b=coordinate
        return Counter({'CX':1}) if b==q else Counter({'CCX':2*q+4*b+7,'CX':8*(b+1),'X':2*(q+1)})
    if kind=='range' and coordinate==0:return Counter({'CX':1})
    n=q+1 if kind=='clip' else coordinate+1
    return Counter({'CCX':6*n,'CX':8*n+1,'X':2*n})


def run(source):
    raw=source.read_bytes();data=json.loads(raw);rows=[];checks=0
    for row in data['configurations']:
        q=row['q'];resources=row['resources'];body=Counter()
        for block in resources['blocks']:body.update(block['primitive_gate_counts'])
        original=Counter({k:2*n for k,n in body.items()});original.update(closed_predicate(q,'clip',q))
        assert original==Counter(resources['expanded_gate_counts'])
        t_body=body['T']+body['TD']+7*body['CCX']
        assert 2*t_body+42*(q+1)==resources['logical_T_count']
        source_row_hash=hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        prepared=vo.prepare_existing_body(resources,q,provenance_sha256=source_row_hash)
        coordinates={'bucket':sorted({0,q-1,q}),'range':sorted({0,q//2,q}),'clip':sorted({0,q//2,q})}
        for kind,coords in coordinates.items():
            for coordinate in coords:
                function={'bucket':vo.profile_from_prepared_body,'range':vo.profile_from_prepared_range,
                          'clip':vo.profile_from_prepared_clipped_range}[kind]
                profile=function(prepared,coordinate)
                predicate=closed_predicate(q,kind,coordinate)
                expected=Counter({k:2*n for k,n in body.items()});expected.update(predicate)
                assert expected==Counter(profile['expanded_gate_counts'])
                t_predicate=7*predicate['CCX']
                assert profile['logical_T_count']==2*t_body+t_predicate
                comparison=coordinate+1 if kind=='bucket' and coordinate<q else 0 if kind=='bucket' else coordinate
                reflection=14*(resources['input_qubits']+comparison-1)
                assert profile['grover_iterate_T_count']==2*profile['logical_T_count']+reflection
                black=profile['blackbox_value_oracle_alternative']
                assert black['logical_T_count']==4*t_body+t_predicate
                expected_black=Counter({k:4*n for k,n in body.items()});expected_black.update(predicate)
                expected_black['CX']+=2*(q+1)
                assert expected_black==Counter(black['expanded_gate_counts'])
                assert black['grover_iterate_T_count']==2*black['logical_T_count']+reflection
                assert profile['conservative_5layer_CCX_T_depth_serial_upper']>=profile['explicit_sequence_T_depth_serial_upper']
                assert black['conservative_5layer_CCX_T_depth_serial_upper']>=black['explicit_sequence_T_depth_serial_upper']
                api=vo.resource_api(profile);assert api['A']['T']==profile['logical_T_count']
                assert api['Q']['T']==profile['grover_iterate_T_count']
                checks+=1
        rows.append({'candidate_id':row.get('candidate_id'),'architecture':resources['architecture'],
            'case_id':row['case_id'],'method':row['method'],'q':q,
            'source_row_sha256':source_row_hash,'forward_body_T':t_body,
            'historical_complete_primitive_and_T_reconstruction':True,
            'bucket_coordinates_checked':coordinates['bucket'],
            'range_coordinates_checked':coordinates['range'],
            'clip_coordinates_checked':coordinates['clip']})
    return {'status':'retained production forward-block reuse and fresh-predicate closed counts passed',
        'source_file':source.name,'source_sha256':hashlib.sha256(raw).hexdigest(),
        'value_oracle_sha256':hashlib.sha256((ROOT/'value_oracle.py').read_bytes()).hexdigest(),
        'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'production_profiles':len(rows),'checked_flag_profiles':checks,
        'architecture_modes':sorted({r['architecture'] for r in rows}),
        'rows':rows,'scope':'independent closed Boolean predicate counts and original emitted-body reconstruction; no rerun of production financial prices or physical hardware claim'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=ROOT/'results/architecture_resource_candidates.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/retained_value_profiles_audit.json');args=parser.parse_args()
    report=run(args.source);args.output.parent.mkdir(exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(args.output)


if __name__=='__main__':main()
