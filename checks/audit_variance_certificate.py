"""Independent exact-map, positive-series moment and confidence audit."""
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from statistics import NormalDist
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time

D = Decimal
ROOT = Path(__file__).resolve().parents[1]
QREF = 48
REF_ERROR = D('1e-130')
SOURCE_ROOT = ROOT


def dec(f):
    f = Fraction(f)
    return D(f.numerator)/D(f.denominator)


@lru_cache(None)
def independent_map(bits,q):
    norm = NormalDist();n = 2**bits
    result = []
    for j in range(n):
        k = j if j>=n//2 else n-j-1
        v = round((norm.inv_cdf((k+.5)/n)+5)*2**q)-5*2**q
        result.append(v if j>=n//2 else -v)
    return result


@lru_cache(None)
def independent_moments(bits):
    # Different traversal and 66 terms rather than the generator's64.
    values = independent_map(bits,QREF)
    sums = [0]*67
    for a in reversed(values[len(values)//2:]):
        square,power = a*a,1
        for k in range(67):
            sums[k] += power
            power = power*square
    return [Fraction(2*v,2**bits*2**(2*k*QREF)) for k,v in enumerate(sums)],max(values)


@lru_cache(None)
def reference_mgf(bits,q,t):
    t = Fraction(t)
    if t==0:return D(1),D(1)
    moments,maximum = independent_moments(bits)
    # Reverse Horner evaluation differs from the generator's forward recurrence.
    square = dec(t)**2
    value = dec(moments[-1]/math.factorial(132))
    for k in range(65,-1,-1):
        value = value*square+dec(moments[k]/math.factorial(2*k))
    a = abs(dec(t))*D(maximum)/D(2**QREF)
    tail = a**134/D(math.factorial(134))/(1-a*a/D(135*136))
    assert tail < REF_ERROR
    delta = D(1)/D(2**(q+1))+D(1)/D(2**(QREF+1))
    return (value-REF_ERROR)*(-abs(dec(t))*delta).exp()-REF_ERROR,\
           (value+REF_ERROR)*(abs(dec(t))*delta).exp()+REF_ERROR


def expmoment(c,v,b,q):
    lo,hi = D(1),D(1)
    for x in v:
        a,z = reference_mgf(b,q,x)
        lo *= a;hi *= z
    return dec(c).exp()*lo-REF_ERROR,dec(c).exp()*hi+REF_ERROR


def geometry(model):
    w = [Fraction(v) for v in model['weights']]
    r = [[Fraction(v) for v in row] for row in model['loadings']]
    mu = [Fraction(v) for v in model['log_drifts']]
    rg = [sum(w[i]*r[i][j] for i in range(len(w))) for j in range(len(r[0]))]
    mg = sum(x*y for x,y in zip(w,mu))
    return w,r,mu,rg,mg


def independent_square(model,b,q,kind,zbits=None):
    w,r,mu,rg,mg = geometry(model)
    if kind=='relative':
        r = [[a-z for a,z in zip(row,rg)] for row in r]
        mu = [a-mg for a in mu]
        rg = [Fraction(0)]*len(rg);mg = Fraction(0)
    aa,ag,gg = D(0),D(0),D(0)
    for i in range(len(w)):
        for j in range(len(w)):
            _,hi = expmoment(mu[i]+mu[j],[a+z for a,z in zip(r[i],r[j])],b,q)
            aa += dec(w[i]*w[j])*hi
        lo,_ = expmoment(mu[i]+mg,[a+z for a,z in zip(r[i],rg)],b,q)
        ag += dec(w[i])*lo
    _,gg = expmoment(2*mg,[2*x for x in rg],b,q)
    square = max(D(0),aa-2*ag+gg)+REF_ERROR
    if kind=='relative':return square
    lam = Fraction(model['lambda'])
    if zbits is None:factor = dec(lam*lam).exp()
    else:_,factor = reference_mgf(zbits,q,2*lam)
    return D(model['spot'])**2*square*factor+REF_ERROR


def density():
    # Positive arcsine series independent of the generator's Machin series.
    term,total,k = D('.5'),D('.5'),0
    while term>D('1e-260'):
        term *= D((2*k+1)**2)/D(8*(k+1)*(2*k+3))
        total += term;k += 1
    return 1/(12*total).sqrt()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config',type=Path,default=SOURCE_ROOT/'results/architecture_resource_candidates.json')
    p.add_argument('--certificates',type=Path,default=ROOT/'results/classical_variance_certificates.json')
    p.add_argument('--source-root',type=Path,default=SOURCE_ROOT)
    p.add_argument('--output',type=Path,default=ROOT/'results/classical_variance_audit.json')
    a = p.parse_args();data = json.loads(a.config.read_text());cert = json.loads(a.certificates.read_text())
    assert cert['configuration_sha256']==hashlib.sha256(a.config.read_bytes()).hexdigest()
    for name,digest in cert['source_hashes'].items():assert digest==hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    sys.path.insert(0,str(a.source_root))
    import encoding_resources as installed
    models = {m['model_id']:m for m in data['models']}
    key = lambda r:(r['case_id'],r['method'],r['q'],r['exp_address_bits'],r['log_address_bits'],r['probability_address_bits'])
    sources = {key(r):r for r in data['configurations']}
    started = time.perf_counter();maps = []
    for bits,q in sorted({(r['b_y'],r['q']) for r in cert['configurations']}|{(r['b_z'],r['q']) for r in cert['configurations'] if r['b_z']}):
        ours = independent_map(bits,q)
        actual,_ = installed.quantile_table(bits,q)
        assert ours==[v-5*2**q for v in actual]
        ref = independent_map(bits,QREF)
        displacement = Fraction(1,2**(q+1))+Fraction(1,2**(QREF+1))
        assert all(abs(Fraction(x,2**q)-Fraction(y,2**QREF))<=displacement for x,y in zip(ours,ref))
        maps.append({'bits':bits,'q':q,'entries':len(ours),'exact_installed_map':True,'nearest_rounding_displacement':True})
    results = []
    with localcontext() as ctx:
        ctx.prec = 300
        phi0 = density()
        for row in cert['configurations']:
            src = sources[key(row)];model = models[row['model_id']]
            cap = D(row['cap']);h = cap/100
            slope = phi0/D(model['lambda'])
            if row['method'] in ('full_global','full_cap'):slope=min(slope,(1+h)/h)
            relative = independent_square(model,row['b_y'],row['q'],'relative')
            candidates = {'range':D(1)}
            if row['method']!='direct':candidates['log_gap']=slope*slope*relative
            if row['method'] in ('direct','full_global','full_cap'):
                candidates['positive_AG_gap']=independent_square(model,row['b_y'],row['q'],'absolute',row['b_z'] if row['method']=='direct' else None)/cap**2
            for name,value in candidates.items():
                assert value<=D(row['exact_grid_second_moment_candidates'][name]),(row['case_id'],name)
            error = D.from_float(src['arithmetic_price_error_upper'])/cap
            assert abs(error-D(row['body_pointwise_normalised_error_upper_decimal']))<D('1e-200')
            raw = min(D(1),(min(candidates.values()).sqrt()+error)**2)
            assert raw<=D(row['finite_output_second_moment_upper_decimal'])
            expected = min(D('.25'),D(row['finite_output_second_moment_upper_decimal']),D(row['legacy_mean_variance_upper_decimal']))
            variance = D(row['finite_output_variance_upper_decimal'])
            assert variance>=expected
            for label,budget in row['confidence_budgets'].items():
                alpha = D('.01') if label=='single_price' else D(1)/2700
                target = D('.005')/cap;n = budget['certified_samples']
                fail = 2*(-D(n)*target*target/(2*variance+2*target/3)).exp() if budget['selected_bound']=='Bernstein' else 2*(-2*D(n)*target*target).exp()
                assert fail+REF_ERROR<=alpha
            assert row['confidence_budgets']['single_price']['certified_samples']<=row['single_price_legacy_samples']
            results.append({'candidate_id':row['candidate_id'],'case_id':row['case_id'],'method':row['method'],'finite_moment_enclosure_conservative':True,'pointwise_body_transfer_checked':True,'confidence_budgets_checked':True})
            if len(results)%25==0:print('Independent finite moment rows',len(results),flush=True)
    out = {'status':'independent exact installed maps,66-term rational coefficients,300-digit moment and fixed-stage confidence audit passed',
        'configuration_sha256':hashlib.sha256(a.config.read_bytes()).hexdigest(),
        'certificate_sha256':hashlib.sha256(a.certificates.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reference_enclosure':str(REF_ERROR),'production_table_checks':maps,
        'configurations':results,'elapsed_seconds':time.perf_counter()-started,
        'qualification':'the arithmetic pointwise premise is the retained independent graph certificate; numerical audit is conditional on correctly rounded Decimal exp/ln/sqrt'}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print('Independent finite variance audit passed:',len(results))


if __name__=='__main__':main()
