import unittest,json,sys
from pathlib import Path
from fractions import Fraction
from math import comb
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from taylor_pricing import signed_stirling,compositions,model_parameters,terms_for_geometry,kernel_derivatives
from continuous_quadrature import frac_arb,exact_endpoint
from flint import arb,arb_series,ctx
ROOT=Path(__file__).resolve().parents[1]

class TaylorContracts(unittest.TestCase):
    def test_signed_stirling_operator(self):
        st=signed_stirling(9)
        self.assertEqual(st[3],[0,2,-3,1])
        for k in range(10):
            for x in [-5,-1,0,1,3,11]:
                falling=1
                for r in range(k):falling*=x-r
                self.assertEqual(falling,sum(s*x**j for j,s in enumerate(st[k])))
    def test_composition_counts(self):
        for n in [2,4,8]:
            for j in range(10):
                values=list(compositions(j,n))
                self.assertEqual(len(values),comb(n+j-1,j))
                self.assertEqual(len(values),len(set(values)))
                self.assertTrue(all(sum(a)==j and min(a)>=0 for a in values))
    def test_exact_model_centring(self):
        for model in json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']:
            w,d,s,mg,rg,lam=model_parameters(model)
            self.assertEqual(sum(w[i]*d[i] for i in range(len(w))),0)
            self.assertTrue(all(sum(w[i]*s[i][j] for i in range(len(w)))==0 for j in range(len(rg))))
    def test_probability_polynomial_mass(self):
        for n in [2,4,8]:
            model=next(m for m in json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models'] if m['assets']==n)
            for row in terms_for_geometry(model,4):self.assertEqual(sum(t[1] for t in row),1)
    def test_kernel_formal_derivative(self):
        ctx.prec=192;ctx.cap=9;m=arb_series([arb('-.2'),1]);s=arb('0.05');h=arb('0.1');l=(1+h).log()
        phi=lambda z:(1+(z/arb(2).sqrt()).erf())/2
        base=100*((m+s*s/2).exp()*(phi(m/s+s)-phi((m-l)/s+s))-(phi(m/s)-phi((m-l)/s))+h*phi((m-l)/s))
        from math import factorial
        values=kernel_derivatives(arb('-.2'),s,10,8)
        for j in range(9):self.assertTrue(values[j].overlaps(base[j]*factorial(j)))
    def test_archived_certificates_and_nulls(self):
        archive=json.loads((ROOT/'results/taylor_pricing.json').read_text())
        self.assertEqual(len(archive['configurations']),27)
        self.assertEqual(archive['successful_prices'],25)
        self.assertEqual(archive['unresolved_prices'],2)
        for row in archive['configurations']:
            self.assertEqual([o['order'] for o in row['orders']],list(range(9)))
            successful=[]
            for o in row['orders']:
                interval=o['undiscounted_price_interval'];lo=Fraction(interval['lower_exact_fraction']);hi=Fraction(interval['upper_exact_fraction'])
                self.assertLessEqual(lo,hi);self.assertGreaterEqual(lo,0);self.assertLessEqual(hi,row['cap'])
                self.assertEqual(o['meets_undiscounted_radius_001'],(hi-lo)/2<=Fraction(1,100))
                if o['meets_undiscounted_radius_001']:successful.append(o['order'])
            self.assertEqual(row['first_certifying_order'],successful[0] if successful else None)

if __name__=='__main__':unittest.main()
