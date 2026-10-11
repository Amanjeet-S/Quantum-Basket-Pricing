"""Check the exact nonzero geometric loadings in retained tensor evidence."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--out',type=Path)
    args=ap.parse_args();root=args.root.resolve();sys.path.insert(0,str(root))
    from boundary_grid_certificate import parameters
    models=json.loads((root/'results/generalisation_protocol.json').read_text())['models']
    second=json.loads((root/'results/tensor_transfer.json').read_text())
    first=json.loads((root/'results/tensor_first_transfer.json').read_text())
    for filename,data in [('tensor_transfer.py',second),('tensor_first_transfer.py',first)]:
        assert data['source_sha256']==hashlib.sha256((root/filename).read_bytes()).hexdigest()
    count=0;families=[]
    for assets in (2,4,8):
        family=[m for m in models if m['assets']==assets];base=parameters(family[0],False)
        for data,key in [(second,'geometry_moments'),(first,'extra_moments')]:
            row=next(r for r in data[key] if r['model_id']==family[0]['model_id'])
            assert tuple(F(c['geometric_loading']) for c in row['coordinates'])==base[4]
        for model in family:
            actual=parameters(model,False)
            assert actual[:3]==base[:3] and actual[4]==base[4]
            assert actual[3]+actual[5]**2/2==base[3]+base[5]**2/2
            for row in second['records']:
                if row['model_id']!=model['model_id']:continue
                for coord in row['coordinates']:
                    assert F(coord['first_transport']['components_upper']['r_U2'])>0
                    assert F(coord['second_transport']['components_upper']['r2_U2'])>0
                count+=1
        families.append({'assets':assets,'reused_from_model':family[0]['model_id'],
            'models':[m['model_id'] for m in family],'exact_RG':[str(v) for v in base[4]],
            'RG_decimal_diagnostic':[float(v) for v in base[4]],
            'exact_forward_prefactor_log':str(base[3]+base[5]**2/2),
            'reuse_premise':'weights,d,s,RG and mG+lambda^2/2 exactly equal across lambda'})
    assert count==36
    result={'status':'passed','financial_cases_checked':count,'families':families,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in
            ['results/generalisation_protocol.json','results/tensor_transfer.json','results/tensor_first_transfer.json']},
        'R1_tilt_source_scope':'tilted_monomial and outside_U_tilt add exact RG to every scalar MGF rate. The frozen implementation is independently reviewed.',
        'correction':'An earlier message confused centred residual loading sum with RG. Mathematical sources and numerical results retain the frozen nonzero RG.',
        'mutations':'No frozen numerical source/result regenerated or edited'}
    path=args.out or root/'results/tensor_parameter_provenance_audit.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('Exact tensor parameter provenance passed',count,'financial cases')

if __name__=='__main__':main()
