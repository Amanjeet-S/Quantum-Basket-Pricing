"""Independent quantile-map and direct finite-sum review of variance MGFs."""
from decimal import Decimal,localcontext
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import importlib
import json
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import encoding_resources as er


def direct_mgf(values,q,coefficient):
    """300-digit direct finite exponentials, without moments or seed transfer."""
    with localcontext() as ctx:
        ctx.prec = 300
        t = Decimal(coefficient.numerator)/Decimal(coefficient.denominator)
        unit = Decimal(1 << q)
        total = Decimal(0)
        for y in values:
            total += (t*Decimal(y)/unit).exp()
        value = total/Decimal(len(values))
        # Fewer than 2**18 exponentials of magnitude below exp(5), plus
        # corresponding additions at 300 digits, are contained in this pad.
        return value,Decimal('1e-260')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--module-dir',type=Path,default=ROOT)
    p.add_argument('--resources',type=Path,default=ROOT/'results/architecture_resource_candidates.json')
    p.add_argument('--output',type=Path,default=ROOT/'results/variance_review_audit.json')
    p.add_argument('--maps-only',action='store_true')
    args = p.parse_args()
    sys.path.insert(0,str(args.module_dir.resolve()))
    vc = importlib.import_module('variance_certificate')
    data = json.loads(args.resources.read_text());rows = data['configurations']
    pairs = {(r['b_y'],r['q']) for r in rows}|{(r['b_z'],r['q']) for r in rows if r['b_z']}
    widths = sorted({b for b,q in pairs})
    started = time.perf_counter();references = {};maps = []
    for b in widths:
        actual,_ = er.quantile_table(b,vc.REFERENCE_Q)
        actual = [x-5*(1 << vc.REFERENCE_Q) for x in actual]
        reference = vc.quantile_integers(b,vc.REFERENCE_Q)
        assert actual == reference,b
        references[b] = np.asarray(reference,dtype=np.int64)
    for b,q in sorted(pairs):
        actual,_ = er.quantile_table(b,q)
        actual = [x-5*(1 << q) for x in actual]
        proposed = vc.quantile_integers(b,q)
        assert actual == proposed,(b,q)
        maximum_support = max(map(abs,actual))
        assert maximum_support < 19*(1 << (q-2)),(b,q,maximum_support)
        diff = np.abs(np.asarray(actual,dtype=np.int64)*(1 << (vc.REFERENCE_Q-q))-references[b])
        numerator = int(diff.max())
        assert 2*numerator <= (1 << (vc.REFERENCE_Q-q))+1,(b,q,numerator)
        maps.append({'bits':b,'q':q,'entries':len(actual),
            'actual_constructor_map_equal':True,
            'maximum_displacement_in_reference_units':numerator,
            'pointwise_rounding_bound_checked':True,
            'maximum_absolute_support_exact_fraction':str(Fraction(maximum_support,1 << q)),
            'inside_body_quantile_domain_minus4_75_to4_75':True})
    profiles = [(9,18,Fraction(1,1000)),(10,26,Fraction(4,5)),
                (12,30,Fraction(3,10)),(17,30,Fraction(1,100)),
                (17,18,Fraction(4,5))]
    # Add a long exact project coefficient from the largest asset model.
    model = max(data['models'],key=lambda m:m['assets'])
    long_coefficient = Fraction(model['loadings'][0][0])+Fraction(model['loadings'][-1][0])
    profiles.append((14,27,long_coefficient))
    if args.maps_only:profiles = []
    direct = []
    for b,q,t in profiles:
        values = vc.quantile_integers(b,q)
        reference,error = direct_mgf(values,q,t)
        lower,upper = vc.mgf_interval(b,q,t)
        assert lower < reference-error and reference+error < upper,(b,q,t,str(lower),str(reference),str(upper))
        direct.append({'bits':b,'q':q,'coefficient_exact_fraction':str(t),
            'direct_finite_sum_300_digits':str(reference),
            'reference_absolute_enclosure':str(error),
            'certificate_lower':str(lower),'certificate_upper':str(upper),
            'entire_direct_enclosure_strictly_inside_certificate':True})
        print('Reviewed finite MGF',b,q,str(t),flush=True)
    out = {'status':'passed','quantile_maps':maps,
        'quantile_map_configurations':len(maps),'quantile_states':sum(m['entries'] for m in maps),
        'reference_table_latent_widths':widths,'reference_q':vc.REFERENCE_Q,
        'direct_mgf_checks':direct,'elapsed_seconds_diagnostic':time.perf_counter()-started,
        'variance_source_sha256':hashlib.sha256((args.module_dir/'variance_certificate.py').read_bytes()).hexdigest(),
        'resources_sha256':hashlib.sha256(args.resources.read_bytes()).hexdigest(),
        'claim_limit':'all retained quantile map transfers checked; direct MGF validation covers selected coefficient profiles; the analytical Taylor proof covers all admissible coefficients'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('quantile_maps','direct_mgf_checks')},indent=2))


if __name__ == '__main__':main()
