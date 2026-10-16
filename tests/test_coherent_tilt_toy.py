"""Meaningful finite LCU identities and no-arbitrage model checks."""
from pathlib import Path
from itertools import product
from collections import Counter
import json,math,sys,unittest
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import coherent_tilt_toy as toy
mp=toy.mp


class CoherentTiltTests(unittest.TestCase):
    def test_general_signed_multinomial_identity(self):
        # Nonzero centred log offsets, arbitrary s, and every finite point.
        for s in range(1,5):
            for y1,y2 in product(toy.Y,repeat=2):
                total=mp.mpf(0)
                for j in range(s+1):
                    for alpha1 in range(j+1):
                        alpha2=j-alpha1
                        coefficient=(-1)**(s-j)*math.comb(s,j)*math.comb(j,alpha1)/mp.mpf(2)**j
                        v=(alpha1-alpha2)*mp.mpf('.125')
                        offset=(alpha1-alpha2)*mp.mpf('.03')
                        total+=coefficient*mp.exp(offset+v*(y1+y2))
                u=mp.cosh(mp.mpf('.03')+mp.mpf('.125')*(y1+y2))-1
                self.assertLess(abs(total-u**s),mp.mpf('1e-90'))

    def test_literal_financial_identity_and_forward(self):
        result=toy.financial_identity('1/4')
        forward=mp.mpf(result['forward'])
        self.assertLess(abs(mp.mpf(result['model_asset_forward_check'])-forward),mp.mpf('1e-90'))
        self.assertLess(abs(mp.mpf(result['undiscounted_price'])-mp.mpf(result['classical_order1_price'])-mp.mpf(result['correction_mean'])),mp.mpf('1e-90'))
        self.assertTrue(0<mp.mpf(result['undiscounted_price'])<10)

    def test_retained_actual_rotations_and_amplification_counts(self):
        compiled=json.loads((ROOT/'results/coherent_tilt_toy_compiled.json').read_text())
        counts=Counter(g[0] for g in compiled['chronological_gates'])
        self.assertEqual(counts['T']+counts['Tdg'],4416)
        self.assertEqual(compiled['aggregate_phase_pi'],'0')
        certificate=json.loads((ROOT/'results/coherent_tilt_toy_rotation_certificates.json').read_text())
        self.assertEqual(len(certificate['rotation_certificates']),11)
        self.assertEqual(certificate['actual_rotation_uses'],42)
        audit=json.loads((ROOT/'results/coherent_tilt_toy_audit.json').read_text())
        self.assertEqual(audit['amplification_actual_T_count'],39*4416+19*56)
        self.assertEqual(audit['direct_base_preparation']['actual_T_count'],736)
        self.assertGreater(audit['ideal_A_query_ratio_importance_to_direct_range_diagnostic'],3)

    def test_vanilla_price_curve_from_same_positive_measure(self):
        # Calls, rather than call spreads, obey convexity in strike.
        delta=mp.mpf('.25');lam=mp.mpf('.2');mu=mp.mpf(toy.financial_identity(delta)['mu']);spot=mp.mpf(100)
        forward=mp.mpf(toy.financial_identity(delta)['forward'])
        Phi=lambda x:(1+mp.erf(x/mp.sqrt(2)))/2
        values=[]
        for strike in map(mp.mpf,[80,100,120]):
            price=mp.mpf(0)
            for j1,j2 in product(range(4),repeat=2):
                m=mu+mp.log(mp.cosh(delta*(toy.Y[j1]+toy.Y[j2])))
                a=(m-mp.log(strike/spot))/lam
                call=spot*mp.exp(m+lam*lam/2)*Phi(a+lam)-strike*Phi(a)
                price+=toy.weights[j1]*toy.weights[j2]*call
            self.assertTrue(max(forward-strike,0)<=price<=forward)
            values.append(price)
        self.assertTrue(values[0]>=values[1]>=values[2])
        self.assertGreaterEqual(values[0]-2*values[1]+values[2],0)


if __name__=='__main__':unittest.main()
