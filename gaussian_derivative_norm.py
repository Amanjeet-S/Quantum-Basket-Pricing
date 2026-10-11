"""Certified cancellation-preserving Gaussian derivative norms.

Exact rational polynomial arithmetic and a proved Sturm chain isolate every
real stationary point. Complex-root routines and floating-point roots are not
used as inclusion certificates. Generation uses the existing Arb runtime;
reading retained norm bounds requires only standard Python.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,math,time
ROOT=Path(__file__).resolve().parent

def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return tuple(p or [F(0)])
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([x*c for x in a])
def shift(a):return trim([0]+list(a))
def derivative(a):return trim([i*a[i] for i in range(1,len(a))])
def evaluate(a,x):
    v=0
    for c in reversed(a):v=v*x+c
    return v
def divrem(a,b):
    a=list(trim(a));b=trim(b)
    if b==(F(0),):raise ZeroDivisionError('Zero polynomial')
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and a!=[F(0)]:
        k=len(a)-len(b);c=a[-1]/b[-1];q[k]+=c
        for i,x in enumerate(b):a[i+k]-=c*x
        a=list(trim(a))
    return trim(q),trim(a)
def gcd(a,b):
    while b!=(F(0),):a,b=b,divrem(a,b)[1]
    return scale(a,1/a[-1])
def squarefree(a):
    d=gcd(a,derivative(a));q,rem=divrem(a,d);assert rem==(F(0),)
    return scale(q,1/q[-1])
def sturm(a):
    a=squarefree(trim(a));seq=[a,derivative(a)]
    while seq[-1]!=(F(0),):
        rem=divrem(seq[-2],seq[-1])[1]
        if rem==(F(0),):break
        # Positive rescaling controls integer growth and preserves variations.
        seq.append(scale(rem,-1/abs(rem[-1])))
    assert len(seq[-1])==1 and seq[-1][0]
    return tuple(seq)
def variations(seq,x):
    values=[evaluate(p,x) for p in seq];signs=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def isolate(a,bits=120):
    a=squarefree(a);seq=sturm(a);d=len(a)-1
    M=max([abs(x/a[-1]) for x in a[:-1]] or [F(0)])
    bound=F(2+math.ceil(M));left,right=-bound,bound
    assert evaluate(a,left) and evaluate(a,right)
    total=variations(seq,left)-variations(seq,right)
    pending=[(left,right,total)];out=[];splits=0;target=F(1,1<<bits)
    while pending:
        lo,hi,count=pending.pop()
        if not count:continue
        if count==1 and hi-lo<=target:
            out.append((lo,hi));continue
        mid=(lo+hi)/2
        if evaluate(a,mid)==0:
            # At most d roots; d+1 interior alternatives guarantee a non-root.
            mid=next(lo+(hi-lo)*F(j,d+2) for j in range(1,d+2) if evaluate(a,lo+(hi-lo)*F(j,d+2)))
        nl=variations(seq,lo)-variations(seq,mid)
        nr=variations(seq,mid)-variations(seq,hi)
        assert 0<=nl<=count and 0<=nr<=count and nl+nr==count
        pending.extend([(mid,hi,nr),(lo,mid,nl)]);splits+=1
    out.sort();assert len(out)==total
    assert all(variations(seq,lo)-variations(seq,hi)==1 for lo,hi in out)
    assert all(out[j][1]<=out[j+1][0] for j in range(len(out)-1))
    return out,seq,bound,splits

def hermites(n):
    h=[(F(1),)]
    if n:h.append((F(0),F(1)))
    for j in range(1,n):h.append(add(shift(h[-1]),scale(h[-2],-j)))
    return h
def signed_stirling(n):
    table=[[1]]
    for r in range(1,n+1):
        prior=table[-1];row=[0]*(r+1)
        for l in range(1,r+1):row[l]=(prior[l-1] if l-1<len(prior) else 0)-(r-1)*(prior[l] if l<len(prior) else 0)
        table.append(row)
    return table
def norm_polynomial(r,lam):
    lam=F(lam)
    if not isinstance(r,int) or isinstance(r,bool) or r<1 or lam<=0:raise ValueError('Positive integer derivative and scale required')
    h=hermites(r-1);s=signed_stirling(r)[r];p=(F(0),)
    for l in range(1,r+1):p=add(p,scale(h[l-1],F(s[l]*(-1)**(l-1))/lam**l))
    return p

def certify(r,lam,bits=120,precision=384,unit_index=None):
    from flint import arb,ctx
    from continuous_quadrature import frac_arb,exact_endpoint
    ctx.prec=precision;lam=F(lam)
    p=norm_polynomial(r,lam) if unit_index is None else scale(hermites(unit_index)[unit_index],(-1)**unit_index)
    q=add(derivative(p),scale(shift(p),-1))
    roots,seq,B,splits=isolate(q,bits)
    records=[];lower=F(0);upper=F(0);normaliser=1/(2*arb.pi()).sqrt()
    def f(x):
        v=arb(0)
        for c in reversed(p):v=v*x+frac_arb(c)
        return v*(-x*x/2).exp()*normaliser
    for lo,hi in roots:
        box=frac_arb((lo+hi)/2)+arb(0,frac_arb((hi-lo)/2).upper())
        z=abs(f(box));sample=abs(f(frac_arb((lo+hi)/2)))
        up=exact_endpoint(z,True);low=max(F(0),exact_endpoint(sample,False))
        upper=max(upper,up);lower=max(lower,low)
        records.append({'lower':str(lo),'upper':str(hi),'width':str(hi-lo),
                        'variation_left':variations(seq,lo),'variation_right':variations(seq,hi),
                        'absolute_function_upper':str(up),'midpoint_value_lower':str(low)})
    assert 0<lower<=upper
    # The old analytic triangle bound is independently valid. Taking its upper
    # enclosure as a second cap prevents harmless rounding from worsening r=1.
    s=signed_stirling(r)[r];b=[normaliser,1/arb.pi()]
    for j in range(2,r):b.append((j-1)*b[j-2])
    old=sum(abs(s[l])*b[l-1]/frac_arb(lam)**l for l in range(1,r+1)) if unit_index is None else b[unit_index]
    old_up=exact_endpoint(old,True);selected=min(upper,old_up)
    assert lower<=selected
    return {'kind':'joint' if unit_index is None else 'unit_hermite_density',
      'unit_derivative_order':unit_index,'derivative_order':r,'lambda':str(lam),'polynomial':[str(x) for x in p],
      'stationary_polynomial':[str(x) for x in q],'squarefree_stationary':[str(x) for x in squarefree(q)],
      'sturm_chain':[[str(x) for x in row] for row in seq],'cauchy_bound':str(B),
      'root_isolation_bits':bits,'root_intervals':records,'real_stationary_points':len(roots),
      'subdivision_steps':splits,'precision_bits':precision,'norm_lower':str(lower),
      'stationary_norm_upper':str(upper),'old_triangle_constant_upper':str(old_up),
      'norm_upper':str(selected),'norm_ratio_upper_to_old_upper':str(selected/old_up),
      'relative_norm_enclosure_width':str((upper-lower)/lower),
      'passed':True}

@lru_cache(None)
def profiles():
    data=json.loads((ROOT/'results/gaussian_derivative_norms.json').read_text())
    return {(F(r['lambda']),r['derivative_order']):r for r in data['records']}
def upper_norm(r,lam):return F(profiles()[(F(lam),r)]['norm_upper'])

def main():
    start=time.perf_counter();records=[]
    units=[certify(j+1,F(1),unit_index=j) for j in range(9)]
    for lam in [F(1,20),F(1,5),F(2,5)]:
        for r in range(1,10):
            record=certify(r,lam)
            s=signed_stirling(r)[r]
            exact_triangle=sum(F(abs(s[l]))*F(units[l-1]['norm_upper'])/lam**l for l in range(1,r+1))
            record['exact_component_triangle_upper']=str(exact_triangle)
            record['norm_upper']=str(min(F(record['norm_upper']),exact_triangle))
            record['norm_ratio_upper_to_old_upper']=str(F(record['norm_upper'])/F(record['old_triangle_constant_upper']))
            record['joint_to_exact_component_triangle_upper_ratio']=str(F(record['norm_upper'])/exact_triangle)
            assert F(record['norm_lower'])<=F(record['norm_upper'])<=exact_triangle
            records.append(record)
            print('lambda',lam,'r',r,'ratio',float(F(record['norm_ratio_upper_to_old_upper'])),flush=True)
    out={'schema':'qbp-coupled-gaussian-derivative-norm-v1','status':'passed',
         'records':records,'unit_norms':units,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'All real stationary points isolated with exact rational Sturm chains; Arb384-bit outward Hermite-Gaussian evaluations. Sharp over positive threshold measures of fixed mass, not necessarily over the specific cap measure.',
         'generation_seconds_diagnostic':time.perf_counter()-start}
    (ROOT/'results/gaussian_derivative_norms.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
