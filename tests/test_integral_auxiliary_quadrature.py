"""Thirty-one exact polynomial/law/norm checks for auxiliary error caps."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,sys,unittest
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import integral_auxiliary_quadrature as q


def beta_case(degree,width):
    N=1<<width;h=F(1,N)
    mass=[F((N-j)**2-(N-j-1)**2,N*N) for j in range(N)]
    assert mass==[2*h*(1-F(2*j+1,2*N)) for j in range(N)]
    exact=F(2,(degree+1)*(degree+2))
    midpoint=sum(p*F(2*j+1,2*N)**degree for j,p in enumerate(mass))
    L1=degree;L2=degree*(degree-1) if degree>=2 else 0
    bound=h*h*(F(L1,6)+F(L2,24));error=abs(exact-midpoint)
    assert error<=bound
    return {'kind':'Beta2 polynomial','degree':degree,'bits':width,'exact_error':str(error),'upper':str(bound),'passed':True}


def strike_case(degree):
    N=8;h=F(1,N);exact=F(1,degree+1)
    midpoint=sum(h*F(2*j+1,2*N)**degree for j in range(N))
    L2=degree*(degree-1) if degree>=2 else 0
    bound=h*h*F(L2,24);error=abs(exact-midpoint)
    assert error<=bound
    return {'kind':'Uniform strike polynomial','degree':degree,'bits':3,'exact_error':str(error),'upper':str(bound),'passed':True}


def norm_case(lam):
    # Coefficients of (A2'-z A2)=lambda(A3+2A2), all exact.
    A2=[-1/lam,-1/lam**2,F(0)]
    A3=[2/lam-1/lam**3,3/lam**2,1/lam**3]
    derivative=[-1/lam**2,1/lam,1/lam**2]
    assert derivative==[lam*(a3+2*a2) for a2,a3 in zip(A2,A3)]
    assert q.joint_norm(4,lam)>0 and q.hermite_norm(3)>0
    result=q.auxiliary_bias({'lambda':str(lam)},25,8,8,F(1,2),F(1,10),F(1,100))
    assert F(result['total_auxiliary_bias_upper'])<=F(result['linear_total_auxiliary_bias_upper'])
    return {'kind':'Exact norm recurrence','lambda':str(lam),'passed':True}


def support_case():
    u=[F(0),F(1,4),F(1,2)];p=[F(1,4),F(1,2),F(1,4)]
    m3=sum(a*x**3 for a,x in zip(p,u));m4=sum(a*x**4 for a,x in zip(p,u))
    assert m4<=max(u)*m3
    b=q.auxiliary_bias({'lambda':'1/5'},10,8,8,max(u),sum(a*x*x for a,x in zip(p,u)),m3,m4)
    assert F(b['fourth_moment_upper_used'])==m4
    return {'kind':'Finite-support fourth moment','moment3':str(m3),'moment4':str(m4),'passed':True}


checks=[(f'test_beta_degree{d}_bits{b}',lambda d=d,b=b:beta_case(d,b)) for d in range(9) for b in [1,2]]
checks += [(f'test_uniform_degree{d}',lambda d=d:strike_case(d)) for d in range(9)]
checks += [(f'test_norm_scale{i}',lambda l=l:norm_case(l)) for i,l in enumerate([F(1,20),F(1,5),F(2,5)])]
checks += [('test_finite_support_moment',support_case)]
assert len(checks)==31
class AuxiliaryQuadratureTests(unittest.TestCase):pass
for name,function in checks:
    def method(self,function=function):self.assertTrue(function()['passed'])
    setattr(AuxiliaryQuadratureTests,name,method)


if __name__=='__main__':
    result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(AuxiliaryQuadratureTests))
    if not result.wasSuccessful():raise SystemExit(1)
    rows=[function() for name,function in checks]
    output={'status':'passed','exact_checks':31,'checks':rows,
            'helper_sha256':hashlib.sha256((ROOT/'integral_auxiliary_quadrature.py').read_bytes()).hexdigest(),
            'test_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'norm_input_sha256':hashlib.sha256((ROOT/'results/gaussian_derivative_norms.json').read_bytes()).hexdigest()}
    (ROOT/'results/integral_auxiliary_quadrature_tests.json').write_text(json.dumps(output,indent=2)+'\n')
