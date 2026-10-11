"""Audit production transport MGFs against literal Gaussian-cell sums."""
import sys,json,hashlib,argparse
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
sys.dont_write_bytecode=True
STAGE=Path(sys.argv[sys.argv.index('--root')+1]) if '--root' in sys.argv else Path(__file__).resolve().parents[1]
if not (STAGE/'results/generalisation_protocol.json').exists():STAGE=Path('/private/tmp/qbp-higherdim-curvature-step')
sys.path.insert(0,str(STAGE))
import tensor_transfer as prod
from boundary_grid_certificate import ball,gaussian_cells,Phi,endpoints,serial,parameters,abs_upper

def direct_mgf(cells,mass,rate,kind):
    v=ball(rate)
    if kind==1:return sum((abs(c['_first'])*ball(F(rate)*F(c['actual_node'])).exp() for c in cells),arb(0))
    if not rate:return sum((c['_second']/2 for c in cells),arb(0))
    rounded=sum((mass*ball(F(1,len(cells)))*ball(F(rate)*F(c['actual_node'])).exp() for c in cells),arb(0))
    linear=sum((c['_first']*ball(F(rate)*F(c['actual_node'])).exp() for c in cells),arb(0))
    cont=(v*v/2).exp()*(Phi(ball(6)-v)-Phi(-ball(6)-v))
    return (cont-rounded-v*linear)/(v*v)

def audit():
    transport=prod.TransportMGF(6,10,28);ctx.prec=768;mass,cells=gaussian_cells(F(6),10,28,768)
    rates=[F(0),F(1,10**40),F(-1,10**40),F(1,1000),F(-1,1000),F(1,20),F(-1,20),F(1,2),F(-1,2),F(2),F(-2)]
    model=json.loads((STAGE/'results/generalisation_protocol.json').read_text())['models'][-1];w,d,s,mg,rg,lam=parameters(model,False)
    for j in range(model['outer_dimension']):
        rates.extend([s[0][j],sum(si[j] for si in s[:2]),2*s[-1][j]])
    rows=[]
    for rate in sorted(set(rates)):
        for kind in (1,2):
            actual=transport.rho1(rate) if kind==1 else transport.second(rate)
            true_measure=direct_mgf(cells,mass,rate,kind)
            expected=(sum((ball(abs_upper(c['_first']))*ball(rate*F(c['actual_node'])).exp() for c in transport.cells),arb(0)) if kind==1 else true_measure)
            al,au=endpoints(actual);el,eu=endpoints(expected)
            if kind==1 and endpoints(true_measure)[1]>au:raise AssertionError('Positive weights failed domination')
            if not al<=el<=eu<=au:raise AssertionError(('Transform missed literal768bit value',rate,kind,float(al),float(au),float(el),float(eu)))
            rows.append({'rate':str(rate),'measure':kind,'production':serial(actual),'literal_768bit_reference':serial(expected),'true_measure_reference':serial(true_measure),'passed':True})
    return {'schema':'qbp-independent-production-transport-audit-v1','status':'passed','transform_checks':len(rows),'precision_bits':768,'records':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'target_source_sha256':hashlib.sha256((STAGE/'tensor_transfer.py').read_bytes()).hexdigest(),'dependency_sha256':{'higherdim_preparation.py':hashlib.sha256((STAGE/'higherdim_preparation.py').read_bytes()).hexdigest()},'scope':'Every checked productionρ1/ρ2transform contains a literal768bit Gaussian-cell sum, including extremely small nonzero rates requiring the positive support fallback. No financial price is evaluated.'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=STAGE);ap.add_argument('--out',type=Path,default=STAGE/'results/tensor_transport_transform_audit.json');args=ap.parse_args()
    r=audit();args.out.write_text(json.dumps(r,indent=2)+'\n');print({k:r[k] for k in ['status','transform_checks']})
