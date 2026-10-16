"""Direct power-series audit of tensor first/second curvature envelopes."""
import sys,json,hashlib,argparse
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_series,ctx
sys.dont_write_bytecode=True
STAGE=Path(sys.argv[sys.argv.index('--root')+1]) if '--root' in sys.argv else Path(__file__).resolve().parents[1]
if not (STAGE/'results/generalisation_protocol.json').exists():STAGE=Path('/private/tmp/qbp-higherdim-curvature-step')
sys.path.insert(0,str(STAGE))
from boundary_grid_certificate import parameters,ball,endpoints,serial,abs_upper
from boundary_derivative_certificate import payoff_norms
from tensor_transfer import geometry_polynomials

def cdf(s):return (1+(s/arb(2).sqrt()).erf())/2

def kernel(s,lam,cap):
    l=ball(lam);g=100*(s+l*l/2).exp();x=s/l
    if cap is None:return g*cdf(x+l)-100*cdf(x),g*cdf(x+l)
    H=ball(cap);a=(1+H/100).log();xh=(s-a)/l
    return g*(cdf(x+l)-cdf(xh+l))-100*(cdf(x)-cdf(xh))+H*cdf(xh),g*(cdf(x+l)-cdf(xh+l))

def geometry(model,y,j):
    w,d,s,mg,rg,lam=parameters(model,False)
    T=arb(0);J=arb(0);JJ=arb(0);m=ball(mg+sum(v*z for v,z in zip(rg,y)))
    for wi,di,si in zip(w,d,s):
        e=ball(wi)*(ball(di+sum(v*z for v,z in zip(si,y)))).exp();T+=e;J+=ball(si[j])*e;JJ+=ball(si[j]**2)*e
    return T,T-1,J,JJ,m,rg[j],lam

def direct_series(model,y,j,cap):
    w,d,s,mg,rg,lam=parameters(model,False);x=arb_series([0,1],prec=4)
    T=arb_series([0],prec=4);m=arb_series([ball(mg+sum(v*z for v,z in zip(rg,y)))],prec=4)+ball(rg[j])*x
    for wi,di,si in zip(w,d,s):
        T+=ball(wi)*(arb_series([ball(di+sum(v*z for v,z in zip(si,y)))],prec=4)+ball(si[j])*x).exp()
    pA,_=kernel(m+T.log(),lam,cap);pG,p1=kernel(m,lam,cap);R1=pA-pG;R2=R1-(T-1)*p1
    return R2[1],2*R2[2],R1[1],2*R1[2]

def polynomial_value(poly,model,y):
    w,d,s,mg,rg,lam=parameters(model,False)
    result=arb(0)
    for alpha,c in poly.items():
        exponent=sum(F(a)*(di+sum(v*z for v,z in zip(si,y))) for a,di,si in zip(alpha,d,s))
        result+=ball(c)*ball(exponent).exp()
    return result

def audit():
    ctx.prec=768;models=json.loads((STAGE/'results/generalisation_protocol.json').read_text())['models'];transfer=json.loads((STAGE/'results/tensor_transfer.json').read_text());rows=[];polys=0
    for model in models:
        p=model['outer_dimension'];points=[tuple(F(0) for _ in range(p)),tuple(F(j+1,10) for j in range(p)),tuple(F((-1)**j,2) for j in range(p)),tuple(F(2*(-1)**j) for j in range(p))]
        for y in points:
            for j in range(p):
                T,U,J,JJ,m,r,lam=geometry(model,y,j);directpol={'U2':U*U,'J2':J*J,'U_JJ':U*JJ,'U_J2':U*J*J,'U2_JJ':U*U*JJ,'U2_J2':U*U*J*J,'U4':U**4}
                for key,poly in geometry_polynomials(model,j).items():
                    a=polynomial_value(poly,model,y);b=directpol[key];al,au=endpoints(a);bl,bu=endpoints(b)
                    if au<bl or bu<al:raise AssertionError('Geometry polynomial disagrees')
                    polys+=1
                for cap in [F(1),F(10),F(25),None]:
                    rec=next(r0 for r0 in transfer['records'] if r0['model_id']==model['model_id'] and r0['cap']==(None if cap is None else str(cap)));n=rec['norms'];M2,M3,M4,D1,D2,D31=[ball(n[k]) for k in ['M2','M3','M4','C2_prime','C2_second','C3_prime']];rr=ball(abs(r));aj=abs(J)
                    f1=M2*U*aj+rr*D1*U*U/2+M3*U*U*aj/6
                    f2=M2*(J*J+U*JJ)+2*rr*D1*U*aj+arb(2)/3*M3*U*J*J+rr*rr*D2*U*U/2+rr*D31*U*U*aj/3+M3*U*U*JJ/6+M4*U*U*J*J/12
                    a,b,c,d=direct_series(model,y,j,cap)
                    if abs_upper(a)>endpoints(f1)[1] or abs_upper(b)>endpoints(f2)[1]:raise AssertionError('R2 derivative exceeds envelope')
                    norms=payoff_norms(lam,cap)
                    if cap is not None:
                        M1,P2,P3=[ball(norms[k]) for k in ['M1_upper','scalar_cap_prime_first_upper','scalar_cap_prime_second_upper']]
                        extra1=M1*aj+rr*P2*U;extra2=M1*JJ+2*rr*P2*aj+rr*rr*P3*U
                    else:
                        forward=100*(m+ball(lam*lam/2)).exp()
                        extra1=forward*aj+rr*(forward+M2)*U;extra2=forward*JJ+2*rr*(forward+M2)*aj+rr*rr*(forward+3*M2+M3)*U
                    if abs_upper(c)>endpoints(f1+extra1)[1] or abs_upper(d)>endpoints(f2+extra2)[1]:raise AssertionError('R1 derivative exceeds envelope')
                    rows.append({'model_id':model['model_id'],'cap':None if cap is None else str(cap),'coordinate':j,'point':list(map(str,y)),'R2_first_derivative':serial(a),'R2_second_derivative':serial(b),'R2_first_envelope':serial(f1),'R2_second_envelope':serial(f2),'R1_first_derivative':serial(c),'R1_second_derivative':serial(d),'passed':True})
    return {'schema':'qbp-independent-pointwise-curvature-envelope-audit-v1','status':'passed','precision_bits':768,'derivative_records':len(rows),'derivative_enclosure_checks':len(rows)*4,'geometry_polynomial_checks':polys,'records':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Analytic norm inequalities independently derived; direct Arb power-series differentiation verifies first and second residual envelopes at768bit points in every frozen model/payoff. This numerical check supplements the proof, not a global proof by sampling.'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=STAGE);ap.add_argument('--out',type=Path,default=STAGE/'results/tensor_curvature_envelope_audit.json');args=ap.parse_args()
    r=audit();args.out.write_text(json.dumps(r,indent=2)+'\n');print({k:r[k] for k in ['status','derivative_records','derivative_enclosure_checks','geometry_polynomial_checks']})
