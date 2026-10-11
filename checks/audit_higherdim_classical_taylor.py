"""Independent series and coefficient checks for ordinary-call Taylor prices."""
import sys,json,hashlib,argparse
from pathlib import Path
from math import factorial
from fractions import Fraction as F
from flint import arb,arb_series,ctx
sys.dont_write_bytecode=True
RESEARCH_ROOT=Path(sys.argv[sys.argv.index('--root')+1]) if '--root' in sys.argv else Path(__file__).resolve().parents[1]
if not (RESEARCH_ROOT/'results/generalisation_protocol.json').exists():RESEARCH_ROOT=Path('/private/tmp/qbp-higherdim-curvature-step')
sys.path.insert(0,str(RESEARCH_ROOT))
import higherdim_classical_taylor as impl
from boundary_grid_certificate import ball,endpoints,serial,coefficient_certificate
from taylor_pricing import model_parameters,dot
REPO=impl.REPO

def overlap(a,b):
    al,au=endpoints(a);bl,bu=endpoints(b)
    return al<=bu and bl<=au

def audit():
    ctx.prec=768;identities=[]
    for lam in [F(1,20),F(1,5),F(2,5)]:
        sigma=ball(lam);t=arb_series([arb(1),arb(1)],prec=8)
        for m in [F(-2,5),F(-1,10),F(0),F(1,20),F(1,5),F(2,5)]:
            x=arb_series([ball(m)],prec=8)+t.log();cdf=lambda z:(1+(z/arb(2).sqrt()).erf())/2
            P=100*((x+sigma*sigma/2).exp()*cdf(x/sigma+sigma)-cdf(x/sigma))
            for order in range(2,7):
                direct=impl.C(order,ball(m),sigma);derived=P[order]*factorial(order)
                if not overlap(direct,derived):raise AssertionError((lam,m,order))
                identities.append({'lambda':str(lam),'m':str(m),'order':order,'gaussian_boundary':serial(direct),'direct_series_derivative':serial(derived),'passed':True})
    models=json.loads((REPO/'results/generalisation_protocol.json').read_text())['models'];coefficients=[]
    production=json.loads((RESEARCH_ROOT/'results/higherdim_classical_taylor.json').read_text())
    for model,row in zip(models,production['records']):
        ref=coefficient_certificate(model,None,'vanilla',1)
        actual=row['orders'][0]['partial_sum'];alo,ahi=F(actual['lower_exact_fraction']),F(actual['upper_exact_fraction']);refc=ref['coefficient'];blo,bhi=F(refc['lower_exact_fraction']),F(refc['upper_exact_fraction'])
        if not alo<=bhi or not blo<=ahi:raise AssertionError('Coefficient disagrees')
        coefficients.append({'model_id':model['model_id'],'original_coefficient':actual,'independent_coefficient':refc,'passed':True})
    result={'schema':'qbp-ordinary-call-taylor-independent-audit-v1','status':'passed','kernel_series_checks':len(identities),'geometric_coefficient_checks':len(coefficients),'kernel_series_precision_bits':768,'coefficient_precision_bits':512,'kernel_records':identities,'coefficient_records':coefficients,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'target_sha256':hashlib.sha256((RESEARCH_ROOT/'higherdim_classical_taylor.py').read_bytes()).hexdigest(),'scope':'Boundary recurrence throughorder6 checked against direct Arb power-series differentiation of P(m+log(1+x)); all9 c0+c1 values independently checked against covariance-preserving Gaussian tilt coefficients. No basket-price answer is read.'}
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=RESEARCH_ROOT);ap.add_argument('--out',type=Path,default=RESEARCH_ROOT/'results/higherdim_classical_taylor_audit.json');args=ap.parse_args()
    r=audit();args.out.write_text(json.dumps(r,indent=2)+'\n');print({k:r[k] for k in ['status','kernel_series_checks','geometric_coefficient_checks']})
