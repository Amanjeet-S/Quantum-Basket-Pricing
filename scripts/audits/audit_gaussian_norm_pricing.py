"""Independent retightening, moments, partial sums and actual stopping counters."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import factorial,comb,prod
from collections import defaultdict
import hashlib,json,sys
from flint import arb,fmpq,fmpq_poly,ctx
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT


def ball(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))


def endpoint(x,upper):
    a,b=(x.upper() if upper else x.lower()).man_exp();a,b=int(a),int(b)
    return F(a*(1<<b)) if b>=0 else F(a,1<<(-b))


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def parameters(model):
    w=tuple(map(F,model['weights']));mu=tuple(map(F,model['log_drifts']));R=tuple(tuple(map(F,row)) for row in model['loadings']);lam=F(model['lambda'])
    mg=dot(w,mu);rg=tuple(sum(wi*r[j] for wi,r in zip(w,R)) for j in range(len(R[0])))
    d=tuple(m-mg for m in mu);v=tuple(tuple(x-y for x,y in zip(r,rg)) for r in R)
    assert sum(w)==1 and min(w)>0 and dot(w,d)==0
    assert all(sum(wi*row[j] for wi,row in zip(w,v))==0 for j in range(len(rg)))
    rate=F(model['risk_free_rate'])-F(model['dividend_yield'])
    assert all(m+(dot(row,row)+lam*lam)/2==rate for m,row in zip(mu,R))
    return w,d,v,mg,rg,lam


def own_compositions(j,n):
    for bars in combinations(range(j+n-1),n-1):
        boundaries=(-1,)+bars+(j+n-1,)
        yield tuple(boundaries[k+1]-boundaries[k]-1 for k in range(n))


def own_geometry(params,maximum):
    w,d,v,mg,rg,lam=params;rows=[];raw=[];moments=[]
    for j in range(maximum+1):
        terms=[]
        for alpha in own_compositions(j,len(w)):
            coefficient=F(factorial(j),prod(factorial(a) for a in alpha))
            for wi,a in zip(w,alpha):coefficient*=wi**a
            exposure=tuple(sum(a*row[k] for a,row in zip(alpha,v)) for k in range(len(rg)))
            exponent=dot(alpha,d)+dot(exposure,exposure)/2
            shift=dot(rg,exposure)
            terms.append((ball(coefficient)*ball(exponent).exp(),shift))
        assert len(terms)==comb(len(w)+j-1,j)
        rows.append(terms);raw.append(sum((x for x,shift in terms),arb(0)))
        moments.append(sum(((-1)**(j-l))*comb(j,l)*raw[l] for l in range(j+1)))
    return {'rows':rows,'raw':raw,'moments':moments,'exponentials':sum(map(len,rows))}


def hermite_value(order,z):
    # Closed monomial coefficients, rather than the producer's recurrence.
    coefficients=[F(0)]*(order+1)
    for k in range(order//2+1):coefficients[order-2*k]=F((-1)**k*factorial(order),2**k*factorial(k)*factorial(order-2*k))
    value=arb(0)
    for coefficient in reversed(coefficients):value=value*z+ball(coefficient)
    return value


def cdf_derivatives(z,sigma,maximum):
    values=[(1+(z/arb(2).sqrt()).erf())/2]
    density=(-z*z/2).exp()/(2*arb.pi()).sqrt()
    for k in range(1,maximum+1):values.append((-1)**(k-1)*hermite_value(k-1,z)*density/sigma**k)
    return values


def own_kernel(m,sigma,cap,maximum):
    # Different representation: derivatives of two separate call prices.
    result=[arb(0)]*(maximum+1);exp_price=100*(m+sigma*sigma/2).exp()
    for strike,sign in [(F(100),1),(F(100+cap),-1)]:
        z=(m-ball(strike/100).log())/sigma
        a=cdf_derivatives(z+sigma,sigma,maximum);b=cdf_derivatives(z,sigma,maximum)
        for k in range(maximum+1):
            first=exp_price*sum((comb(k,l)*a[l] for l in range(k+1)),arb(0))
            result[k]+=sign*(first-ball(strike)*b[k])
    return result


def own_coefficients(params,cap,k,geometry):
    w,d,v,mg,rg,lam=params;sigma=ball(lam*lam+dot(rg,rg)).sqrt();tilted=[]
    for row in geometry['rows'][:k+1]:
        values=[arb(0)]*(k+1)
        for factor,shift in row:
            derivatives=own_kernel(ball(mg+shift),sigma,cap,k)
            for l in range(k+1):values[l]+=factor*derivatives[l]
        tilted.append(values)
    answer=tilted[0][0]
    for degree in range(1,k+1):
        polynomial=fmpq_poly([1]);x=fmpq_poly([0,1])
        for a in range(degree):polynomial*=x-a
        coefficient=arb(0)
        for l in range(1,degree+1):
            total=sum(((-1)**(degree-j))*comb(degree,j)*tilted[j][l] for j in range(degree+1))
            coefficient+=int(polynomial[l])*total
        answer+=coefficient/factorial(degree)
    return answer


def fractions_interval(record):return F(record['lower_exact_fraction']),F(record['upper_exact_fraction'])


def main():
    ctx.prec=512
    pricing_path=SOURCE/'results/gaussian_norm_pricing.json';data=json.loads(pricing_path.read_text())
    for name,expected in data['source_hashes'].items():assert hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==expected
    for name,expected in data['input_hashes'].items():assert hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==expected
    norm=json.loads((SOURCE/'results/gaussian_derivative_norms.json').read_text())
    norm_audit=json.loads((SOURCE/'results/gaussian_derivative_norm_independent_audit.json').read_text())
    assert norm_audit['input_sha256']==hashlib.sha256((SOURCE/'results/gaussian_derivative_norms.json').read_bytes()).hexdigest()
    D={(F(row['lambda']),row['derivative_order']):F(row['norm_upper']) for row in norm['records']}
    protocol=json.loads((SOURCE/'results/generalisation_protocol.json').read_text());models={m['model_id']:m for m in protocol['models']}
    old=json.loads((SOURCE/'results/taylor_pricing.json').read_text());old_rows={row['case_id']:row for row in old['configurations']}
    params={key:parameters(model) for key,model in models.items()}
    geometries={}
    for param in params.values():
        w,d,v,mg,rg,lam=param;key=(w,d,v,rg)
        if key not in geometries:geometries[key]=own_geometry(param,9)
    assert sum(row['exponentials'] for row in geometries.values())==25080
    retightened=[];degrees=0
    for row in data['retightened_paths']['configurations']:
        baseline=old_rows[row['case_id']];param=params[baseline['model_id']];w,d,v,mg,rg,lam=param;G=geometries[(w,d,v,rg)]
        lower=upper=F(0);good=[]
        for order,prior in zip(row['orders'],baseline['orders']):
            k=order['order'];a,b=fractions_interval(prior['coefficient_interval']);lower+=a;upper+=b
            assert F(order['coefficient_interval_sum']['lower'])==lower and F(order['coefficient_interval_sum']['upper'])==upper
            r=k+1;independent_rem=max(F(0),endpoint(row['cap']*ball(D[(lam,r)])*G['moments'][r]/factorial(r),True))
            saved_rem=F(order['remainder_upper']);assert independent_rem<=saved_rem
            expected=(max(F(0),lower-saved_rem),min(F(row['cap']),upper+saved_rem))
            assert fractions_interval(order['interval'])==expected
            certified=(expected[1]-expected[0])/2<=F(1,100)
            assert order['certified']==certified
            if certified:good.append(k)
            degrees+=1
        assert row['certified']==bool(good)
        assert row['first_certifying_order']==(good[0] if good else None)
        retightened.append({'case_id':row['case_id'],'degrees_checked':len(row['orders']),
                            'coefficient_interval_sums_exact':True,'moments_independently_recomputed':True,
                            'first_certifying_order':row['first_certifying_order'],'passed':True})

    hierarchy=[];maxima=defaultdict(lambda:-1);total_kernel=total_geometry=0
    for row in data['executed_hierarchy']['configurations']:
        model_id=row['case_id'].rsplit('_cap',1)[0];param=params[model_id];w,d,v,mg,rg,lam=param;key=(w,d,v,rg);G=geometries[key]
        before=maxima[key];maximum=1;expected_kernel=0
        if row['selected'].get('one_sided_monotonicity',False):
            expected_kernel=1
            sigma=ball(lam*lam+dot(rg,rg)).sqrt();partial=own_kernel(ball(mg),sigma,row['cap'],0)[0]
            phi_slope=1/(2*arb.pi()).sqrt()/ball(lam)
            slope=min(endpoint(phi_slope,True),F(100+row['cap'],row['cap']))
            independent_gap=endpoint(row['cap']*ball(slope)*G['raw'][1].log(),True)
            assert independent_gap<=F(row['selected']['remainder_upper'])
            lower,upper=fractions_interval(row['selected']['interval'])
            assert lower<=endpoint(partial,False) and endpoint(partial,True)+independent_gap<=upper
        else:
            for screen in row['screens']:
                k=screen['order'];maximum=max(maximum,k+1)
                own_rem=max(F(0),endpoint(row['cap']*ball(D[(lam,k+1)])*G['moments'][k+1]/factorial(k+1),True))
                assert own_rem<=F(screen['remainder_upper'])
                assert screen['remainder_qualified']==(F(screen['remainder_upper'])<=F(1,100))
                assert screen['coefficient_evaluated']==(screen['remainder_qualified'] or k==8)
            for evaluated in row['evaluated_orders']:
                k=evaluated['order'];expected_kernel+=comb(len(w)+k,k)
                partial=own_coefficients(param,row['cap'],k,G);rem=F(evaluated['remainder_upper'])
                lo=max(F(0),endpoint(partial,False)-rem);hi=min(F(row['cap']),endpoint(partial,True)+rem)
                saved_lo,saved_hi=fractions_interval(evaluated['interval'])
                assert saved_lo<=lo<=hi<=saved_hi
        maxima[key]=max(before,maximum)
        expected_geo=comb(len(w)+maxima[key],maxima[key])-(comb(len(w)+before,before) if before>=0 else 0)
        assert row['kernel_calls']==expected_kernel and row['new_geometry_exponentials']==expected_geo
        total_kernel+=expected_kernel;total_geometry+=expected_geo
        hierarchy.append({'case_id':row['case_id'],'actual_coefficient_partial_sum_independently_enclosed':True,
                          'kernel_calls':expected_kernel,'new_geometry_exponentials':expected_geo,'passed':True})
    assert total_kernel==793==data['executed_hierarchy']['kernel_calls']
    assert total_geometry==760==data['executed_hierarchy']['geometry_exponentials']

    # Instrument actual coefficient calls in an independent fresh stopping run.
    sys.path.insert(0,str(SOURCE));import taylor_order_benchmark as policy;import gaussian_norm_pricing as producer
    ctx.prec=224;original_radius=policy.radius;original_kernel=policy.kernel_derivatives;actual_calls=0;banks={}
    def counted(*args,**kwargs):
        nonlocal actual_calls
        actual_calls+=1;return original_kernel(*args,**kwargs)
    policy.radius=producer.joint_radius;policy.kernel_derivatives=counted
    try:
        for model in protocol['models']:
            param=policy.model_parameters(model);w,d,v,mg,rg,lam=param;bank=banks.setdefault((w,d,v,rg),policy.Geometry((w,d,v,rg)))
            for cap in protocol['caps']:
                before=actual_calls;actual=policy.evaluate(model,cap,'monotone_then_screen',bank)
                saved=next(row for row in data['executed_hierarchy']['configurations'] if row['case_id']==actual['case_id'])
                assert actual_calls-before==actual['kernel_calls']==saved['kernel_calls']
                assert actual['selected']['interval']==saved['selected']['interval']
    finally:policy.radius=original_radius;policy.kernel_derivatives=original_kernel
    assert actual_calls==793 and sum(bank.exponentials for bank in banks.values())==760
    result={'status':'passed','retightened_cases':len(retightened),'retightened_degree_intervals':degrees,
            'independent_geometry_exponential_terms':25080,'retightening_new_kernel_calls':0,
            'hierarchy_actual_kernel_calls':actual_calls,'hierarchy_geometry_exponential_terms':total_geometry,
            'precision_bits_independent':512,'retightened_profiles':retightened,'hierarchy_profiles':hierarchy,
            'input_sha256':hashlib.sha256(pricing_path.read_bytes()).hexdigest(),'source_hashes_verified':data['source_hashes'],
            'input_hashes_verified':data['input_hashes'],'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'independence':'Exact model reconstruction, stars-and-bars compositions, independent high-precision moments, closed monomial Hermites, falling-factorial Stirling coefficients, and separate-call-price derivative formula. Instrumented producer replay checks actual stopping counts after the independent mathematical checks.'}
    (ROOT/'results/gaussian_norm_pricing_independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','retightened_intervals':degrees,'hierarchy_calls':actual_calls,'hierarchy_geometry':total_geometry}))


if __name__=='__main__':main()
