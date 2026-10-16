"""Independent exact polynomial, Sturm completeness and Arb norm audit.

Polynomial arithmetic is delegated to Flint and coefficients are rebuilt
from closed Hermite coefficients and falling-factorial products. Neither
the producer's Sturm chain nor its isolator is used as the source of truth.
Run in the existing qf-certification runtime. All output is staged locally.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
from flint import fmpq, fmpq_poly, arb, ctx


def fq(x):
    x=F(x)
    return fmpq(x.numerator,x.denominator)


def fraction(x):
    return F(int(x.numerator),int(x.denominator))


def poly(coefficients):
    return fmpq_poly([fq(x) for x in coefficients])


def coefficients(p):
    return tuple(fraction(x) for x in p.coeffs()) or (F(0),)


def independent_hermite(n):
    c=[F(0)]*(n+1)
    for k in range(n//2+1):
        c[n-2*k]=F((-1)**k*math.factorial(n),
                    2**k*math.factorial(k)*math.factorial(n-2*k))
    return poly(c)


def falling_factorial(n):
    answer=poly([1]);z=poly([0,1])
    for j in range(n):answer*=z-j
    return answer


def independent_norm_polynomial(r,scale):
    lam=F(scale);s=falling_factorial(r);p=poly([0])
    for l in range(1,r+1):
        p+=independent_hermite(l-1)*s[l]*fq(F((-1)**(l-1),1)/lam**l)
    return p


def squarefree(p):
    q,rem=divmod(p,p.gcd(p.derivative()))
    assert rem.is_zero()
    return q/q.leading_coefficient()


def independent_sturm(p):
    p=squarefree(p);q=p.derivative();chain=[p,q]
    while not chain[-1].is_constant():
        _,remainder=divmod(chain[-2],chain[-1])
        assert not remainder.is_zero()
        chain.append(-remainder)
    return chain


def sign(v):return (v>0)-(v<0)


def variation_signs(values):
    values=[x for x in values if x]
    return sum(x!=y for x,y in zip(values,values[1:]))


def variation(chain,x):
    return variation_signs([sign(p(fq(x))) for p in chain])


def infinity_variation(chain,positive=True):
    return variation_signs([sign(p.leading_coefficient())*
                           (1 if positive or p.degree()%2==0 else -1)
                           for p in chain])


def ball_fraction(x):return arb(fq(x))


def exact_endpoint(x,upper):
    a,b=(x.upper() if upper else x.lower()).man_exp()
    a,b=int(a),int(b)
    return F(a*(1<<b)) if b>=0 else F(a,1<<(-b))


def interval_box(lo,hi):
    # Union of separately enclosed exact endpoints, no midpoint-radius reuse.
    return ball_fraction(lo).union(ball_fraction(hi))


def value(p,z):
    # Independent ascending monomial evaluation.
    answer=arb(0);power=arb(1)
    for k in range(p.degree()+1):
        answer+=ball_fraction(fraction(p[k]))*power
        power*=z
    return answer


def function_value(p,z):
    return value(p,z)*(-z*z/2).exp()/(2*arb.pi()).sqrt()


def own_old_triangle(r,lam):
    b=[1/(2*arb.pi()).sqrt(),1/arb.pi()]
    for j in range(2,r):b.append((j-1)*b[j-2])
    s=falling_factorial(r)
    return sum((abs(int(s[l]))*b[l-1]/ball_fraction(lam)**l
                for l in range(1,r+1)),arb(0))


def check_chain(record,chain):
    saved=[poly(row) for row in record['sturm_chain']]
    assert len(saved)==len(chain)
    for actual,expected in zip(saved,chain):
        quotient,remainder=divmod(actual,expected)
        assert remainder.is_zero() and quotient.is_constant()
        assert quotient[0]>0


def audit_profile(record,p,kind='coupled'):
    q=p.derivative()-p.left_shift(1)
    assert coefficients(poly(record['polynomial']))==coefficients(p)
    assert coefficients(poly(record['stationary_polynomial']))==coefficients(q)
    assert coefficients(poly(record['squarefree_stationary']))==coefficients(squarefree(q))
    chain=independent_sturm(q);check_chain(record,chain)
    _,factors=q.factor()
    rational_roots=[fraction(-factor[0]/factor[1])
                    for factor,multiplicity in factors if factor.degree()==1]
    complete=infinity_variation(chain,False)-infinity_variation(chain,True)
    assert complete==record['real_stationary_points']==len(record['root_intervals'])
    B=F(record['cauchy_bound'])
    assert variation(chain,-B)-variation(chain,B)==complete
    assert q(fq(-B)) and q(fq(B))
    lows=[];ups=[];boxes=[];prior_hi=None
    ctx.prec=640
    for root in record['root_intervals']:
        lo,hi=F(root['lower']),F(root['upper'])
        assert lo<hi and hi-lo==F(root['width'])
        assert hi-lo<=F(1,1<<record['root_isolation_bits'])
        if prior_hi is not None:assert prior_hi<=lo
        prior_hi=hi
        assert q(fq(lo)) and q(fq(hi))
        assert variation(chain,lo)-variation(chain,hi)==1
        assert variation(chain,lo)==root['variation_left']
        assert variation(chain,hi)==root['variation_right']
        # Reconstruct the producer's actual precision box and check its
        # endpoints against exact rationals. This checks representation,
        # not the root location (already established independently above).
        ctx.prec=record['precision_bits']
        mid=(lo+hi)/2
        producer_box=ball_fraction(mid)+arb(0,ball_fraction((hi-lo)/2).upper())
        prod_lo=exact_endpoint(producer_box,False)
        prod_hi=exact_endpoint(producer_box,True)
        assert prod_lo<=lo<=hi<=prod_hi
        assert producer_box.contains(fq(lo)) and producer_box.contains(fq(hi))
        ctx.prec=640
        box=interval_box(lo,hi)
        assert exact_endpoint(box,False)<=lo and exact_endpoint(box,True)>=hi
        sample=abs(function_value(p,ball_fraction(mid)))
        # At every root Q=0; bound variation via the derivative Q*phi.
        derivative_up=exact_endpoint(abs(function_value(q,box)),True)
        midpoint_lo=max(F(0),exact_endpoint(sample,False))
        critical_up=exact_endpoint(sample,True)+derivative_up*(hi-lo)/2
        direct_up=exact_endpoint(abs(function_value(p,box)),True)
        critical_up=min(critical_up,direct_up)
        if p.degree()==0:
            constant_global=ball_fraction(abs(fraction(p[0])))/(2*arb.pi()).sqrt()
            critical_up=min(critical_up,exact_endpoint(constant_global,True))
        exact_roots=[z for z in rational_roots if lo<=z<=hi]
        assert len(exact_roots)<=1
        if exact_roots:
            critical_up=min(critical_up,exact_endpoint(abs(function_value(p,ball_fraction(exact_roots[0]))),True))
        assert midpoint_lo>=F(root['midpoint_value_lower'])
        assert critical_up<=F(root['absolute_function_upper'])
        lows.append(midpoint_lo);ups.append(critical_up)
        boxes.append({'lower':str(lo),'upper':str(hi),
                      'producer_Arb_lower':str(prod_lo),'producer_Arb_upper':str(prod_hi),
                      'complete_root_count':'1','independent_value_lower':str(midpoint_lo),
                      'independent_critical_norm_upper':str(critical_up)})
    assert complete>0
    norm_lo=max(lows);norm_up=max(ups)
    assert norm_lo>=F(record['norm_lower'])
    assert norm_up<=F(record['stationary_norm_upper'])
    if kind=='coupled':
        old=exact_endpoint(own_old_triangle(record['derivative_order'],F(record['lambda'])),True)
        assert old<=F(record['old_triangle_constant_upper'])
        own_selected=min(old,norm_up)
        assert own_selected<=F(record['norm_upper'])
    else:
        own_selected=norm_up
        # Unit r=0 is known exactly up to the Arb enclosure of phi(0).
        if p.degree()==0:own_selected=min(own_selected,exact_endpoint(1/(2*arb.pi()).sqrt(),True))
        assert own_selected<=F(record['norm_upper'])
    return {'kind':kind,'derivative_order':record.get('derivative_order'),
            'hermite_order':record.get('hermite_order',record.get('unit_derivative_order')), 'lambda':record.get('lambda'),
            'real_stationary_points':complete,'all_real_roots_accounted_for':True,
            'all_Arb_root_boxes_contain_exact_rational_endpoints':True,
            'independent_norm_lower':str(norm_lo),'independent_norm_upper':str(own_selected),
            'root_audit':boxes,'passed':True}


def custom_root_cases(producer_path):
    spec=importlib.util.spec_from_file_location('producer_under_test',producer_path)
    producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    cases=[]
    # Rational stationary roots and a repeated stationary root.
    for name,p,expected,exact_roots in [
        ('linear_signal_rational_roots',poly([0,1]),2,[F(-1),F(1)]),
        ('quadratic_flat_stationary_root',poly([2,0,1]),1,[F(0)]),
        ('constant_gaussian',poly([1]),1,[F(0)])]:
        q=p.derivative()-p.left_shift(1);chain=independent_sturm(q)
        assert infinity_variation(chain,False)-infinity_variation(chain,True)==expected
        intervals,*_=producer.isolate(coefficients(q),bits=96)
        assert len(intervals)==expected
        for lo,hi in intervals:
            assert variation(chain,lo)-variation(chain,hi)==1
            assert sum(lo<=root<=hi for root in exact_roots)==1
        content,factors=q.factor_squarefree()
        multiplicities=[int(power) for factor,power in factors]
        if name=='quadratic_flat_stationary_root':assert multiplicities==[3]
        cases.append({'case':name,'stationary_polynomial':list(map(str,coefficients(q))),
                      'distinct_real_roots':expected,'squarefree_multiplicities':multiplicities,
                      'passed':True})
    return cases


def main():
    project=Path(__file__).resolve().parent.parent
    if not (project/'gaussian_derivative_norm.py').exists():
        project=Path('/private/tmp/qbp-gaussian-norm-step')
    parser=argparse.ArgumentParser()
    parser.add_argument('--producer',type=Path,default=project/'gaussian_derivative_norm.py')
    parser.add_argument('--input',type=Path,default=project/'results/gaussian_derivative_norms.json')
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('gaussian_derivative_norm_independent_audit.json'))
    args=parser.parse_args();data=json.loads(args.input.read_text())
    records=[]
    for row in data['records']:
        p=independent_norm_polynomial(row['derivative_order'],F(row['lambda']))
        records.append(audit_profile(row,p))
    units=data.get('unit_norms',data.get('unit_hermite_records',data.get('hermite_records',[])))
    for row in units:
        order=row.get('unit_derivative_order',row.get('hermite_order',row.get('degree')))
        records.append(audit_profile(row,independent_hermite(order)*((-1)**order),'unit_Hermite'))
    threeway=[]
    if units:
        assert len(units)==9
        for row in data['records']:
            r=row['derivative_order'];lam=F(row['lambda']);s=falling_factorial(r)
            component=sum(abs(fraction(s[l]))*F(units[l-1]['norm_upper'])/lam**l
                          for l in range(1,r+1))
            assert component==F(row['exact_component_triangle_upper'])
            assert F(row['norm_upper'])/component==F(row['joint_to_exact_component_triangle_upper_ratio'])
            threeway.append({'lambda':str(lam),'derivative_order':r,
                             'exact_component_triangle_upper':str(component),
                             'joint_norm_upper':row['norm_upper'],'passed':True})
    assert data['source_sha256']==hashlib.sha256(args.producer.read_bytes()).hexdigest()
    result={'status':'passed','input_sha256':hashlib.sha256(args.input.read_bytes()).hexdigest(),
            'producer_source_sha256':hashlib.sha256(args.producer.read_bytes()).hexdigest(),
            'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'independent_implementation':'Flint fmpq_poly, falling factorial Stirling coefficients, closed monomial Hermite coefficients; own Euclidean Sturm chain and infinity signs; Arb640 midpoint plus critical-derivative upper envelopes',
            'coupled_records':len(data['records']),'unit_Hermite_records':len(units),
            'profiles':records,'threeway_opponents':threeway,'custom_cases':custom_root_cases(args.producer),
            'scope':'Sharpness over positive threshold measures is proved separately; no specific-cap sharpness or theoretical priority is asserted.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','coupled_records':len(data['records']),
                      'unit_records':len(units),'stationary_roots':sum(x['real_stationary_points'] for x in records),
                      'output':str(args.output)}))


if __name__=='__main__':main()
