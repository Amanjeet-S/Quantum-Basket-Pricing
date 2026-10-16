"""Exact rational/quadratic-field quadrature and finite-normalisation checks."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import boundary_auxiliary_certificate as aux
from boundary_value_certificate import gauss_constants
from boundary_value_oracle import build_graph
from generalisation_models import models

class Q6:
    def __init__(self,a=0,b=0):self.a,self.b=F(a),F(b)
    def __add__(self,x):
        x=x if isinstance(x,Q6) else Q6(x);return Q6(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return Q6(-self.a,-self.b)
    def __sub__(self,x):return self+-Q6(x) if not isinstance(x,Q6) else self+-x
    def __mul__(self,x):
        x=x if isinstance(x,Q6) else Q6(x);return Q6(self.a*x.a+6*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __pow__(self,n):
        p=Q6(1)
        for _ in range(n):p=p*self
        return p
    def __eq__(self,x):
        x=x if isinstance(x,Q6) else Q6(x);return self.a==x.a and self.b==x.b


def main():
    passed=0
    zm,zp=Q6(F(2,5),F(-1,10)),Q6(F(2,5),F(1,10))
    wm,wp=Q6(F(1,2),F(1,18)),Q6(F(1,2),F(-1,18))
    moment=lambda k:F(2,(k+1)*(k+2))
    for k in range(4):assert wm*zm**k+wp*zp**k==moment(k);passed+=1
    # Orthogonality and squared norm, evaluated by distinct exact coefficients.
    pi=(F(1,10),F(-4,5),F(1))
    for k in (0,1):assert sum(a*moment(j+k) for j,a in enumerate(pi))==0;passed+=1
    norm=sum(a*b*moment(i+j) for i,a in enumerate(pi) for j,b in enumerate(pi))
    assert norm==F(1,300);passed+=1
    assert Q6(moment(4))-(wm*zm**4+wp*zp**4)==Q6(F(1,300));passed+=1
    for b in range(1,9):
        N=1<<b;h=F(1,N)
        cell=lambda j:F(2*N-2*j-1,N*N)
        assert sum(cell(j) for j in range(N))==1
        assert sum(cell(j)*(F(2*j+1,2*N)) for j in range(N))-moment(1)==h*h/6
        assert sum(cell(j)*(F(2*j+1,2*N))**2 for j in range(N))-moment(2)==h*h/12
        passed+=3
    for q in (3,5,8,28):
        g=gauss_constants(q)
        assert F(g['w_minus']['actual'])+F(g['w_plus']['actual'])==1
        assert 0<=F(g['z_minus']['actual'])<F(g['z_plus']['actual'])<=1
        passed+=2
    m=models()[0]
    for payoff in ('capped','vanilla'):
        for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
            for enc in ('centred','positive_part','negative_part') if payoff=='capped' and kind!='positive_r1' else ('centred',):
                g=build_graph(m,25 if payoff=='capped' else None,kind=kind,payoff=payoff,U_star=0,output_encoding=enc)
                expected=F(1,2) if payoff=='capped' and kind!='positive_r1' and enc=='centred' else F(0)
                assert g.nodes[g.output].meta==expected;passed+=1
    result={'status':'passed','checks':passed,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'coverage':'Exact independent Q(sqrt6) degree3 rule, degree4 sharp error, polynomial orthogonality/norm, finite Beta cell mass and linear/quadratic examples, dyadic complement weights and zero-support graph conventions.'}
    target=Path(__file__).resolve().parents[1]/'results/boundary_auxiliary_tests.json';target.write_text(json.dumps(result,indent=2)+'\n')
    print(passed,'exact checks passed')

if __name__=='__main__':main()
