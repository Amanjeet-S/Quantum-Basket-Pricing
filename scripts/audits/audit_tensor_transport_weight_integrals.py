"""Independent direct integration checks of the Gaussian transport MGF.
The transform quotient is checked against positive weighted interval sums,
not against another implementation of the same quotient.
"""
import sys,json,hashlib,argparse,time
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
sys.dont_write_bytecode=True
REPO=Path(sys.argv[sys.argv.index('--root')+1]) if '--root' in sys.argv else Path(__file__).resolve().parents[1]
if not (REPO/'results/generalisation_protocol.json').exists():REPO=Path('/private/tmp/qbp-higherdim-curvature-step')
sys.path.insert(0,str(REPO))
from boundary_grid_certificate import ball,endpoints,serial,interval,gaussian_cells,Phi,transport_weight_mass

def continuous_mgf(v,L=6):
    v=ball(v);return (v*v/2).exp()*(Phi(ball(L)-v)-Phi(-ball(L)-v))

def unambiguous_mgf(cells,v,mass):
    v=F(v);N=len(cells)
    rounded=sum((mass*ball(F(1,N))*ball(v*F(c['actual_node'])).exp() for c in cells),arb(0))
    linear=sum((c['_first']*ball(v*F(c['actual_node'])).exp() for c in cells),arb(0))
    cont=continuous_mgf(v)
    direct=sum((c['_second']/2 for c in cells),arb(0)) if v==0 else (cont-rounded-ball(v)*linear)/ball(v*v)
    return direct,cont,rounded,linear

def positive_weight_integral(cells,v,segments=16):
    total=arb(0)
    for c in cells:
        node=F(c['actual_node']);left=F(c['left']['lower_exact_fraction']);right=F(c['right']['upper_exact_fraction'])
        for outer,e,ecdf,side in [(left,c['_left'],c['_left_cdf'],'left'),(right,c['_right'],c['_right_cdf'],'right')]:
            lo,hi=(outer,node) if side=='left' else (node,outer)
            for j in range(segments):
                a=lo+(hi-lo)*F(j,segments);b=lo+(hi-lo)*F(j+1,segments)
                weight=transport_weight_mass(a,b,e,ecdf)
                total+=weight*(ball(v)*interval(a,b)).exp()
    return total

def audit():
    ctx.prec=768;start=time.perf_counter();rows=[];checks=0
    for bits in [2,3,5,10]:
        mass,cells=gaussian_cells(F(6),bits,28,768)
        for v in [F(0),F(1,1000),F(-1,1000),F(1,20),F(-1,20),F(1,2),F(-1,2),F(2),F(-2)]:
            direct,cont,rounded,linear=unambiguous_mgf(cells,v,mass)
            positive=positive_weight_integral(cells,v,16 if bits<10 else 4)
            dl,du=endpoints(direct);pl,pu=endpoints(positive)
            if dl<0 or not pl<=dl<=du<=pu:
                raise AssertionError((bits,v,float(dl),float(du),float(pl),float(pu)))
            opposite,_,_,_=unambiguous_mgf(cells,-v,mass)
            ol,ou=endpoints(opposite)
            if du<ol or ou<dl:raise AssertionError('Even MGF symmetry failed')
            rows.append({'b_y':bits,'v':str(v),'rho2_mgf':serial(direct),'positive_weight_integral_enclosure':serial(positive),'continuous_mgf':serial(cont),'rounded_subprobability_mgf':serial(rounded),'signed_first_moment_transform':serial(linear),'passed':True});checks+=1
    # A tensor product exponential has exact separable means. The telescope
    # must use the raw mass-a transforms, including each untouched coordinate.
    mass,cells=gaussian_cells(F(6),3,28,768);tensor=[]
    for vec in [(F(1,20),F(1,2)),(F(2),F(-1,2),F(1,1000)),(F(0),F(2),F(-2),F(1,20))]:
        cs=[];ds=[];errors=[]
        for v in vec:
            rho,c,d,lin=unambiguous_mgf(cells,v,mass);cs.append(c);ds.append(d);errors.append(ball(v)*lin+ball(v*v)*rho)
        original=arb(1);rounded=arb(1)
        for c,d in zip(cs,ds):original*=c;rounded*=d
        telescope=arb(0)
        for j,e in enumerate(errors):
            product=e
            for k in range(len(vec)):
                if k<j:product*=ds[k]
                elif k>j:product*=cs[k]
            telescope+=product
        diff=original-rounded;err=telescope-diff
        if not err.contains(0):raise AssertionError('Tensor telescope failed')
        tensor.append({'rates':list(map(str,vec)),'continuous_product':serial(original),'rounded_product':serial(rounded),'telescope':serial(telescope),'difference':serial(diff),'roundoff_residual':serial(err),'passed':True})
    return {'schema':'qbp-independent-tensor-transport-mgf-review-v1','status':'passed','precision_bits':768,'scalar_checks':checks,'tensor_checks':len(tensor),'records':rows,'tensor_records':tensor,'elapsed_seconds_diagnostic':time.perf_counter()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'identity':'M_cont(v) - M_rounded(v) = v sum_c mu1,c exp(vc) + v^2 M_rho2(v). Both endpoint laws are subprobabilities of mass a. At v=0, M_rho2(0)=.5 sum_c mu2,c.','independence':'The quotient transform is enclosed by direct positive weighted integration over the Gaussian transport density, with no basket price input.'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=REPO);ap.add_argument('--out',type=Path,default=REPO/'results/tensor_transport_weight_audit.json');args=ap.parse_args()
    result=audit();args.out.write_text(json.dumps(result,indent=2)+'\n');print({k:result[k] for k in ['status','scalar_checks','tensor_checks','elapsed_seconds_diagnostic']})
