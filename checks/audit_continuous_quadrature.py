"""Independent positive strike-integral checks of scalar cap formulae."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys
import time
from flint import acb,arb,ctx

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import continuous_quadrature as cq


def main():
    start=time.perf_counter();ctx.prec=192
    path=ROOT/'results/continuous_quadrature.json';result=json.loads(path.read_text())
    assert result['source_sha256']==hashlib.sha256((ROOT/'continuous_quadrature.py').read_bytes()).hexdigest()
    models={m['model_id']:m for m in result['models']};checks=[];calls=0
    for row in result['configurations']:
        model=models[row['model_id']];cap=row['cap']
        w,r,mu,lam,_=cq.model_checks(model)
        la=cq.frac_arb(lam);h=arb(cap)/100;sqrt2=arb(2).sqrt()
        for y in (-8,-3,-1,0,1,3,8):
            logs=[cq.frac_arb(mu[i]+r[i][0]*y) for i in range(2)]
            ma=sum(cq.frac_arb(w[i])*logs[i].exp() for i in range(2)).log()
            mg=sum(cq.frac_arb(w[i])*logs[i] for i in range(2))
            for kind,m in [('arithmetic',ma),('geometric',mg)]:
                count=[0]
                def positive(u,analytic):
                    count[0]+=1
                    strike=(1+h*u).log(analytic=analytic)
                    return cap*(1+((m-strike)/(la*sqrt2)).erf())/2
                positive_ball=acb.integral(positive,0,1,rel_tol=arb('1e-28'),abs_tol=arb('1e-28'),
                    eval_limit=100000,depth_limit=40)
                assert positive_ball.is_finite() and positive_ball.imag.contains(0)
                closed=cq.capped_conditional(acb(m),la,cap)
                assert closed.real.overlaps(positive_ball.real)
                # The nonnegative strike representation is checked without
                # clipping either interval to manufacture agreement.
                assert closed.real.lower()<=cap and closed.real.upper()>=0
                calls+=count[0]
                checks.append({'case_id':row['case_id'],'y':y,'basket':kind,
                    'closed_kernel_ball':closed.real.str(45,more=True),
                    'positive_strike_integral_ball':positive_ball.real.str(45,more=True),
                    'enclosures_overlap':True,'positive_integrand_calls':count[0]})
        interval=row['undiscounted_price_interval']
        lo,hi=Fraction(interval['lower_exact_fraction']),Fraction(interval['upper_exact_fraction'])
        assert hi-lo<Fraction(1,1000)
        assert Fraction(interval['lower_decimal'])<=lo and Fraction(interval['upper_decimal'])>=hi
        # Independently check the normal-tail erfc bound against Mills' bound.
        L=row['normal_truncation']
        erfc_tail=cap*(arb(L)/arb(2).sqrt()).erfc()
        mills=2*cap*(-arb(L)**2/2).exp()/((2*arb.pi()).sqrt()*L)
        assert erfc_tail.upper()<mills.lower()
        assert Fraction(row['normal_tail_price_bound_upper_exact_fraction'])>=cq.exact_endpoint(erfc_tail,True)
        for check in result['independent_representation_checks']:
            if check['case_id']==row['case_id']:
                ref=check['reference']['undiscounted_price_interval']
                assert max(lo,Fraction(ref['lower_exact_fraction']))<=min(hi,Fraction(ref['upper_exact_fraction']))
                break
        else:raise AssertionError('Missing direct representation check')
    out={'status':'passed','scalar_cap_kernel_checks':len(checks),'positive_strike_integrand_calls':calls,
        'continuous_price_certificates':len(result['configurations']),'checked_real_y_values':[-8,-3,-1,0,1,3,8],
        'precision_bits':192,'checks':checks,'result_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'audit_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'claim_limit':'Independent positive integral representation checks at 126 scalar states, with full retained price certificates from validated outer quadrature. Both use the same Arb arithmetic library.'}
    (ROOT/'results/continuous_quadrature_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
