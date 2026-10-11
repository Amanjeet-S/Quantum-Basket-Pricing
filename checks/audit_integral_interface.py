"""Exact signed Born counts, coupling endpoints and integral identities."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import comb,factorial
import hashlib,json,sys
from flint import fmpq_poly,fmpq
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT


def poly(c):return fmpq_poly([fmpq(F(x).numerator,F(x).denominator) for x in c])


def hermite(n):
    c=[F(0)]*(n+1)
    for k in range(n//2+1):c[n-2*k]=F((-1)**k*factorial(n),2**k*factorial(k)*factorial(n-2*k))
    return poly(c)


def A(r,lam):
    s=poly([1]);x=poly([0,1])
    for j in range(r):s*=x-j
    result=poly([0])
    for l in range(1,r+1):
        scale=F((-1)**(l-1),1)/lam**l
        result+=hermite(l-1)*s[l]*fmpq(scale.numerator,scale.denominator)
    return result


def no_arbitrage_parameter_audit():
    models=json.loads((SOURCE/'results/generalisation_protocol.json').read_text())['models']
    rows=[]
    for model in models:
        original=[[F(x) for x in row]+[F(model['lambda'])] for row in model['loadings']]
        matrix=[row[:] for row in original];rank=0
        for column in range(len(matrix[0])):
            pivot=next((j for j in range(rank,len(matrix)) if matrix[j][column]),None)
            if pivot is None:continue
            matrix[rank],matrix[pivot]=matrix[pivot],matrix[rank]
            scale=matrix[rank][column];matrix[rank]=[x/scale for x in matrix[rank]]
            for j in range(len(matrix)):
                if j!=rank:
                    factor=matrix[j][column]
                    matrix[j]=[x-factor*y for x,y in zip(matrix[j],matrix[rank])]
            rank+=1
            if rank==len(matrix):break
        assert rank==len(original)
        rate=F(model['risk_free_rate'])-F(model['dividend_yield'])
        assert all(F(mu)+sum(x*x for x in row)/2==rate for mu,row in zip(model['log_drifts'],original))
        rows.append({'model_id':model['model_id'],'stock_count':len(original),'Brownian_dimension':len(original[0]),'exact_volatility_rank':rank,'risk_neutral_log_drift_exact':True,'passed':True})
    return rows


def main():
    born=[];threshold_cases=rounding_cases=0
    for q in range(1,5):
        N=1<<(q+1)
        for a in range(N+1):
            successes=sum(J<a for J in range(N));p=F(successes,N)
            assert p==F(a,N) and 2*p-1==F(a,1<<q)-1
            assert a.bit_length()<=q+2
            threshold_cases+=1
        eta=F(1,10)
        for true,offset in product(map(F,['-1','-3/4','-1/8','0','1/8','3/4','1']),[-eta,F(0),eta]):
            approximation=max(F(-1),min(F(1),true+offset))
            threshold=(1<<q)*(approximation+1);a=threshold.numerator//threshold.denominator
            decoded=F(a,1<<q)-1
            assert F(0)<=approximation-decoded<F(1,1<<q)
            assert abs(decoded-true)<=eta+F(1,1<<q)
            rounding_cases+=1
        born.append({'q':q,'uniform_bits':q+1,'threshold_output_bits':q+2,
                     'zero_and_saturated_endpoints_checked':True,'passed':True})
    coupling=0
    for order,width in [(1,1),(2,1),(2,2),(3,2)]:
        N=1<<width;h=F(1,2*N)
        for cells in product(range(N),repeat=order):
            midpoint=min(F(2*j+1,2*N) for j in cells)
            for ends in product((0,1),repeat=order):
                continuous=min(F(j+e,N) for j,e in zip(cells,ends))
                assert abs(midpoint-continuous)<=h
                coupling+=1
    remainders=[]
    for r in range(1,7):
        for u in [F(0),F(1,3),F(1)]:
            exact=(1+u)**6-sum(F(comb(6,j))*u**j for j in range(r))
            derivative_expectation=F(factorial(6),factorial(6-r))*sum(
                F(comb(6-r,l))*u**l*F(factorial(l)*factorial(r),factorial(r+l))
                for l in range(7-r))
            integral=u**r*derivative_expectation/factorial(r)
            assert integral==exact
            remainders.append({'order':r,'U':str(u),'exact_remainder':str(exact),'passed':True})
    data=json.loads((SOURCE/'results/gaussian_derivative_norms.json').read_text())
    norm={(F(row['lambda']),row['derivative_order']):F(row['norm_upper']) for row in data['records']}
    M=[F(row['norm_upper']) for row in data['unit_norms']];strike=[]
    for lam in [F(1,20),F(1,5),F(2,5)]:
        for r in range(1,9):
            current=A(r,lam);Q=current.derivative()-current.left_shift(1)
            assert Q==(A(r+1,lam)+r*current)*fmpq(lam.numerator,lam.denominator)
            s=poly([1]);x=poly([0,1])
            for j in range(r):s*=x-j
            tri=sum(F(abs(int(s[l])))*M[l]/lam**l for l in range(1,r+1))
            joint=lam*(norm[lam,r+1]+r*norm[lam,r])
            strike.append({'lambda':str(lam),'order':r,'Hermite_component_upper':str(tri),
                           'joint_recurrence_upper':str(joint),'selected_upper':str(min(tri,joint)),
                           'recurrence_exact':True,'passed':True})
    result={'status':'passed','no_arbitrage_parameters':no_arbitrage_parameter_audit(),'threshold_cases':threshold_cases,'signed_rounding_cases':rounding_cases,
            'signed_Born_profiles':born,'midpoint_coupling_endpoint_cases':coupling,
            'exact_Beta_integral_remainders':remainders,'strike_norm_recurrences':strike,
            'norm_input_sha256':hashlib.sha256((SOURCE/'results/gaussian_derivative_norms.json').read_bytes()).hexdigest(),
            'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Encoding interface and mathematical identities. The financial value oracle is not emitted by this audit.'}
    (ROOT/'results/integral_remainder_interface_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','thresholds':threshold_cases,'rounding':rounding_cases,'coupling':coupling,'strike_recurrences':len(strike)}))


if __name__=='__main__':main()
