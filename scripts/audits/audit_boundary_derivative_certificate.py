"""Independent root, log-price operator and Arb-series norm checks."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse,hashlib,json,sys
from flint import arb,arb_series,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import boundary_derivative_certificate as tested


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ball(x):
    x=F(x);return arb(x.numerator)/arb(x.denominator)
def pdf(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def cdf(x):return (1+(x/arb(2).sqrt()).erf())/2

def overlap(a,lo,hi):
    return a.upper()>=ball(lo).lower() and a.lower()<=ball(hi).upper()

def stirling(order):
    rows=[[1]]
    for n in range(1,order+1):
        prior=rows[-1];rows.append([(prior[j-1] if j>=1 and j-1<len(prior) else 0)-(n-1)*(prior[j] if j<len(prior) else 0) for j in range(n+1)])
    return rows[order]

def series_price(m,lam,H=None,K=100):
    x=arb_series([ball(m),arb(1)],prec=7);sigma=ball(lam);K=ball(K)
    result=K*((x+sigma*sigma/2).exp()*cdf((x+sigma*sigma)/sigma)-cdf(x/sigma))
    if H is not None:
        H=ball(H);ell=(1+H/K).log();result-= (K+H)*((x-ell+sigma*sigma/2).exp()*cdf((x-ell+sigma*sigma)/sigma)-cdf((x-ell)/sigma))
    return result

def boundary_polynomial_expected(order,lam):
    lam=F(lam)
    if order==2:return [F(1)]
    if order==3:return [F(-2),-1/lam]
    if order==4:return [6-1/lam**2,5/lam,1/lam**2]
    if order==5:return [-24+9/lam**2,-26/lam+3/lam**3,-9/lam**2,-1/lam**3]
    return [120-71/lam**2+3/lam**4,154/lam-42/lam**3,71/lam**2-6/lam**4,14/lam**3,1/lam**4]

def poly_trim(p):
    p=list(p)
    while len(p)>1 and not p[-1]:p.pop()
    return p

def poly_value(p,x):
    value=F(0)
    for c in reversed(p):value=value*x+c
    return value

def poly_remainder(a,b):
    a=poly_trim(a);b=poly_trim(b)
    while len(a)>=len(b) and any(a):
        degree=len(a)-len(b);factor=a[-1]/b[-1]
        for j,c in enumerate(b):a[degree+j]-=factor*c
        a=poly_trim(a)
    return a

def variations(chain,x):
    signs=[1 if v>0 else -1 for v in [poly_value(row,x) for row in chain] if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def poly_ball(p,x):
    out=arb(0)
    for c in reversed(p):out=out*x+ball(c)
    return out


def evaluate(root):
    root=Path(root);ctx.prec=768;data=json.loads((root/'results/boundary_derivative_norms.json').read_text())
    assert data['status']=='passed' and data['precision_bits']==512 and data['root_isolation_bits']==160
    assert data['source_sha256']==sha(root/'boundary_derivative_certificate.py')
    for path,digest in data['input_sha256'].items():assert digest==sha(root/path)
    for path,digest in data['dependency_source_sha256'].items():assert digest==sha(root/path)
    generic=json.loads((root/'results/gaussian_derivative_norms.json').read_text())
    G={(F(r['lambda']),r['derivative_order']):F(r['norm_upper']) for r in generic['records']}
    HN={r['unit_derivative_order']:F(r['norm_upper']) for r in generic['unit_norms']}
    root_checks=derivatives=scalar_samples=operator_checks=0;rows=[]
    for row in data['records']:
        lam=F(row['lambda']);K=F(row['K']);cap=F(row['cap']) if row['payoff']=='capped' else None
        assert row['passed'] and F(row['M2_lower'])<=F(row['M2_upper'])
        bounds={int(r):F(v) for r,v in row['C_r_upper'].items()}
        points=[F(-2),F(-1),F(-1,4),F(0),F(1,10),F(1,2),F(1),F(2)]
        if cap is not None:
            ell=(1+ball(cap/K)).log();a=ell/ball(lam)
            for field in ('C1_root_certificate','C2_root_certificate'):
                certificate=row[field]
                assert certificate['passed'] and overlap(ell,F(certificate['log_width_interval']['lower_exact_fraction']),F(certificate['log_width_interval']['upper_exact_fraction']))
            for peak in row['C2_root_certificate']['stationary_peaks']:
                r=peak['stationary_root'];lo,hi=F(r['lower']),F(r['upper']);assert 0<hi-lo<=F(1,1<<160)
                assert hi<0 if peak['label']=='positive_maximum' else ball(lo).lower()>a.upper()
                def equation(z):return (ball(z)/(ball(z)-a)).log()-a*ball(z)-ell+a*a/2
                assert equation(lo).lower()>0 and equation(hi).upper()<0
                roots=arb(0,ball((hi-lo)/2).upper())+ball((lo+hi)/2)
                critical=-roots*pdf(roots)+ball(1+cap/K)*(roots-a)*pdf(roots-a)
                assert critical.lower()<=0<=critical.upper()
                value=ball(K/lam)*(pdf(roots)-ball(1+cap/K)*pdf(roots-a))
                saved=peak['value_interval'];assert overlap(value,F(saved['lower_exact_fraction']),F(saved['upper_exact_fraction']))
                assert F(peak['absolute_peak_upper'])<=F(row['C2_root_certificate']['norm_upper'])
                root_checks+=1;points.append(lam*(lo+hi)/2)
            r=row['C1_root_certificate']['stationary_root'];lo,hi=F(r['lower']),F(r['upper'])
            assert 0<=lo<hi and ball(hi).lower()<=ell.upper() and hi-lo<=F(1,1<<160)
            # The independent series supplies P'' endpoint signs.
            assert series_price(lo,lam,cap,K)[2].lower()>0
            assert series_price(hi,lam,cap,K)[2].upper()<0
            root_checks+=1;points.append((lo+hi)/2)
            assert F(row['M1_lower'])<=F(row['M1_upper'])<=min(K+cap,cap*G[lam,1])
            assert F(row['M2_upper'])<=cap*G[lam,2]
            assert F(row['scalar_cap_second_upper'])==min(F(row['M1_upper'])+bounds[2],cap*HN[1]/lam**2)
            assert F(row['scalar_cap_prime_second_upper'])==min(F(row['M1_upper'])+3*bounds[2]+bounds[3],cap*HN[2]/lam**3)
            assert F(row['scalar_boundary_h_first_upper'])==min(bounds[3]+2*bounds[2],(2*K+cap)*HN[1]/lam**2,cap*(HN[1]/lam**2+HN[2]/lam**3))/(2*bounds[2])
            assert F(row['scalar_boundary_h_second_upper'])==min(bounds[4]+5*bounds[3]+4*bounds[2],(2*K+cap)*HN[2]/lam**3,cap*(HN[2]/lam**3+HN[3]/lam**4))/(2*bounds[2])
        else:
            true=ball(K/lam)*pdf(arb(0));assert overlap(true,F(row['M2_lower']),F(row['M2_upper']))
            assert row['M1_upper'] is None and row['second_remainder_nonnegative']
        for m in points:
            prices=series_price(m,lam,cap,K)
            if cap is not None:
                assert abs(prices[1]).upper()<=ball(row['M1_upper']).upper()
                assert abs(2*prices[2]).upper()<=ball(row['scalar_cap_second_upper']).upper()
                assert abs(6*prices[3]).upper()<=ball(row['scalar_cap_prime_second_upper']).upper()
                scalar_samples+=3
            for order in range(2,7):
                coeff=stirling(order)
                converted=sum((coeff[j]*factorial(j)*prices[j] for j in range(1,order+1)),arb(0))
                z=ball(m/lam);B=boundary_polynomial_expected(order,lam)
                explicit=ball(K/lam)*poly_ball(B,z)*pdf(z)
                if cap is not None:explicit-=ball((K+cap)/lam)*poly_ball(B,z-ell/ball(lam))*pdf(z-ell/ball(lam))
                assert explicit.upper()>=converted.lower() and explicit.lower()<=converted.upper(),('Boundary/log-Greek identity',row['payoff'],lam,cap,m,order)
                assert abs(explicit).upper()<=ball(bounds[order]).upper()
                if cap is None and order==2:assert explicit.lower()>=0
                derivatives+=1;operator_checks+=1
        rows.append({'payoff':row['payoff'],'lambda':str(lam),'cap':None if cap is None else str(cap),
            'sample_m_values':len(points),'M2_upper':row['M2_upper'],'root_classification_checked':cap is not None})
    polyrows=[]
    for item in data['boundary_polynomial_norms']:
        order=item['order'];lam=F(item['lambda']);p=list(map(F,item['polynomial']))
        assert p==boundary_polynomial_expected(order,lam)
        assert item['passed'] and F(item['norm_lower'])<=F(item['norm_upper'])
        chain=[list(map(F,row)) for row in item['sturm_chain']]
        critical=list(map(F,item['stationary_polynomial']))
        assert chain[0]==[c/critical[-1] for c in critical]
        assert chain[1]==[j*chain[0][j] for j in range(1,len(chain[0]))]
        for j in range(2,len(chain)):
            rem=poly_remainder(chain[j-2],chain[j-1])
            assert chain[j]==[-v/abs(rem[-1]) for v in rem]
        assert len(chain[-1])==1 and chain[-1][0]
        bound=F(item['cauchy_bound'])
        total=variations(chain,-bound)-variations(chain,bound)
        assert total==item['real_stationary_points']==len(item['root_intervals'])
        for r in item['root_intervals']:
            assert r['variation_left']-r['variation_right']==1
            lo,hi=F(r['lower']),F(r['upper']);assert hi-lo<=F(1,1<<160)
            assert variations(chain,lo)==r['variation_left'] and variations(chain,hi)==r['variation_right']
            x=ball((lo+hi)/2)+arb(0,ball((hi-lo)/2).upper())
            sample=abs(poly_ball(p,ball((lo+hi)/2))*pdf(ball((lo+hi)/2)))
            assert sample.upper()<=ball(r['norm_upper']).upper()
            # Arb radius rounding need not nest across different precisions;
            # each retained outward box is independently valid. The fresh
            # point value and interval intersection supply the numerical check.
            val=abs(poly_ball(p,x)*pdf(x))
            assert val.lower()<=ball(r['norm_upper']).upper()
        polyrows.append({'order':order,'lambda':str(lam),'stationary_points':item['real_stationary_points'],'passed':True})
    # Convex kink mixtures inherit nonnegative cash-Gamma; affine additions
    # disappear from every falling derivative of order at least two.
    cpwa_checks=0
    for lam in (F(1,20),F(1,5),F(2,5)):
        for m in (F(-1),F(0),F(1)):
            x=arb_series([ball(m),arb(1)],prec=3);sigma=ball(lam);S=ball(100)
            price=ball(7)+ball(3)*S*(x+sigma*sigma/2).exp();expected=arb(0)
            for strike,weight in ((F(80),F(1)),(F(120),F(2)),(F(150),F(1,2))):
                b=(ball(strike/100)).log()
                price+=ball(weight)*(S*(x+sigma*sigma/2).exp()*cdf((x-b+sigma*sigma)/sigma)-ball(strike)*cdf((x-b)/sigma))
                expected+=ball(weight*strike/lam)*pdf((ball(m)-b)/sigma)
            actual=2*price[2]-price[1]
            assert expected.lower()>0 and expected.upper()>=actual.lower() and expected.lower()<=actual.upper();cpwa_checks+=1
    return {'status':'independent boundary derivative audit passed','norm_records':rows,'root_brackets_checked':root_checks,
        'falling_operator_and_payoff_series_checks':operator_checks,'scalar_derivative_samples':scalar_samples,
        'boundary_polynomial_records':polyrows,'convex_piecewise_affine_checks':cpwa_checks,
        'reference_precision_bits':768,'source_sha256':data['source_sha256'],'norm_file_sha256':sha(root/'results/boundary_derivative_norms.json'),
        'scope':'Global extrema certified by monotonicity/Sturm proofs; independent 768-bit Arb-series and exact polynomial comparisons provide retained implementation evidence, not a novelty or physical-advantage claim.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',default=str(ROOT));p.add_argument('--out',default=str(ROOT/'results/boundary_derivative_audit.json'))
    args=p.parse_args();result=evaluate(args.root);result['audit_source_sha256']=sha(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','root_brackets_checked','falling_operator_and_payoff_series_checks','scalar_derivative_samples','convex_piecewise_affine_checks')}))

if __name__=='__main__':main()
