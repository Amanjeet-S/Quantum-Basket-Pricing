"""Independent exact-cell and financial scaling audit of auxiliary bounds."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def polynomial_integral(coefficients,lo=F(0),hi=F(1)):
    return sum((F(c)*(hi**(j+1)-lo**(j+1))/F(j+1) for j,c in enumerate(coefficients)),F(0))
def value(coefficients,x):return sum((F(c)*x**j for j,c in enumerate(coefficients)),F(0))
def derivative_envelope(coefficients,order):
    result=F(0)
    for j,c in enumerate(coefficients):
        if j<order:continue
        factor=1
        for k in range(order):factor*=j-k
        result+=abs(F(c))*factor
    return result


def evaluate(root):
    root=Path(root);sys.path.insert(0,str(root));import integral_auxiliary_quadrature as helper
    norms=json.loads((root/'results/gaussian_derivative_norms.json').read_text())
    assert norms['status']=='passed'
    joint={(r['derivative_order'],F(r['lambda'])):F(r['norm_upper']) for r in norms['records'] if r['passed']}
    hermite={r['unit_derivative_order']:F(r['norm_upper']) for r in norms['unit_norms'] if r['passed']}
    laws=[];polynomials=[];moment_records=[];financial=[];rejections=0
    coefficients_list=[(0,),(1,),(0,1),(0,-1),(0,0,1),(1,-2,3,-4),
        (-2,1,0,5,-3,2),(1,-1,1,-1,1,-1,1)]
    for N in (1,2,3,4,7,8,16,32):
        h=F(1,N);counts=Counter(min(i,j) for i in range(N) for j in range(N))
        weights=[]
        for cell in range(N):
            lo,hi=cell*h,(cell+1)*h;mid=(lo+hi)/2
            exact_mass=polynomial_integral((2,-2),lo,hi)
            assert exact_mass==F(counts[cell],N*N)
            first=polynomial_integral((-2*mid,2+2*mid,-2),lo,hi)
            # Expand (z-mid)^2 times 2(1-z) independently.
            second=polynomial_integral((2*mid*mid,-4*mid-2*mid*mid,2+4*mid,-2),lo,hi)
            assert first==-h**3/6
            assert second==(1-mid)*h**3/6
            weights.append(exact_mass)
        assert sum(weights)==1
        assert sum(w*F(2*j+1,2*N) for j,w in enumerate(weights))-F(1,3)==h*h/6
        assert F(1,N)*sum(F(2*j+1,2*N)**2 for j in range(N))-F(1,3)==-h*h/12
        laws.append({'cells':N,'minimum_pair_outcomes_enumerated':N*N,
            'exact_cell_mass_and_first_second_moments_checked':True,
            'affine_Beta_error_fraction':str(h*h/6),'quadratic_uniform_error_fraction':str(h*h/12)})
        for coefficients in coefficients_list:
            # Integral of 2(1-z)f(z), independently expanded.
            weighted=[F(0)]*(len(coefficients)+1)
            for j,c in enumerate(coefficients):weighted[j]+=2*c;weighted[j+1]-=2*c
            exact_beta=polynomial_integral(weighted)
            beta_mid=sum(w*value(coefficients,F(2*j+1,2*N)) for j,w in enumerate(weights))
            L1=derivative_envelope(coefficients,1);L2=derivative_envelope(coefficients,2)
            beta_bound=h*h*(L1/6+L2/24);beta_error=abs(exact_beta-beta_mid)
            assert beta_error<=beta_bound
            exact_uniform=polynomial_integral(coefficients)
            uniform_mid=h*sum(value(coefficients,F(2*j+1,2*N)) for j in range(N))
            uniform_bound=h*h*L2/24;uniform_error=abs(exact_uniform-uniform_mid)
            assert uniform_error<=uniform_bound
            polynomials.append({'cells':N,'coefficients':list(coefficients),
                'Beta_error_fraction':str(beta_error),'Beta_bound_fraction':str(beta_bound),
                'uniform_error_fraction':str(uniform_error),'uniform_bound_fraction':str(uniform_bound)})
    support=[F(0),F(1,8),F(1,4),F(3,8)]
    for weights in ((F(1,4),)*4,(F(1,2),F(1,4),F(1,8),F(1,8))):
        for relative in ((F(1),)*4,(F(1,2),F(1,3),F(1,4),F(1,5))):
            moments={j:sum(w*r*u**j for w,r,u in zip(weights,relative,support)) for j in (2,3,4)}
            assert moments[4]<=max(support)*moments[3]
            moment_records.append({'moment2_fraction':str(moments[2]),'moment3_fraction':str(moments[3]),
                'moment4_fraction':str(moments[4]),'weighted_likelihood':relative!=(F(1),)*4})
            for lam in (F(1,20),F(1,5),F(2,5)):
                for H in (F(1),F(10),F(25)):
                    for bits in (1,4,8):
                        D3,D4=joint[3,lam],joint[4,lam]
                        E=hermite[1]/lam+hermite[2]/lam**2
                        F2=hermite[2]/lam+hermite[3]/lam**2
                        hz=hs=F(1,1<<bits);m2,m3,m4=moments[2],moments[3],moments[4]
                        actual=helper.auxiliary_bias({'lambda':str(lam)},H,bits,bits,max(support),m2,m3,m4)
                        bz=H*hz*hz*(D3*m3/12+D4*m4/48)
                        bs=H**3*m2*hs*hs*(F2/lam**2+E/lam)/(48*100**2)
                        zl=H*D3*m3*hz/4;sl=H**2*E*m2*hs/(4*100*lam)
                        assert F(actual['E2_upper'])==E and F(actual['F2_upper'])==F2
                        assert F(actual['D3_upper'])==D3 and F(actual['D4_upper'])==D4
                        assert F(actual['quadratic_beta_bias_upper'])==bz
                        assert F(actual['quadratic_strike_bias_upper'])==bs
                        assert F(actual['total_auxiliary_bias_upper'])==min(bz,zl)+min(bs,sl)
                        assert F(actual['fourth_moment_upper_used'])==m4
                        default=helper.auxiliary_bias({'lambda':str(lam)},H,bits,bits,max(support),m2,m3)
                        assert F(default['fourth_moment_upper_used'])==max(support)*m3
                        assert F(default['quadratic_beta_bias_upper'])>=bz
                        financial.append({'lambda_fraction':str(lam),'cap_fraction':str(H),'bits':bits,
                            'quadratic_Beta_bias_fraction':str(bz),'quadratic_strike_bias_fraction':str(bs),
                            'selected_auxiliary_bias_fraction':actual['total_auxiliary_bias_upper']})
    zero=helper.auxiliary_bias({'lambda':'1/5'},25,8,8,0,0,0)
    assert F(zero['total_auxiliary_bias_upper'])==0
    for kwargs in ({'cap':0},{'U_star':-1},{'moment2_upper':-1},{'moment3_upper':-1},
                   {'moment4_upper':-1},{'beta_bits':True},{'beta_bits':0},{'strike_bits':F(3,2)},{'K':0}):
        arguments={'model':{'lambda':'1/5'},'cap':25,'beta_bits':8,'strike_bits':8,'U_star':F(3,8),
            'moment2_upper':F(1,16),'moment3_upper':F(1,128)};arguments.update(kwargs)
        try:helper.auxiliary_bias(**arguments)
        except ValueError:rejections+=1
        else:raise AssertionError(('Invalid auxiliary inputs accepted',kwargs))
    return {'status':'independent auxiliary quadrature audit passed','exact_cell_laws':laws,
        'signed_polynomial_comparisons':polynomials,'finite_support_moment_comparisons':moment_records,
        'financial_norm_and_scaling_comparisons':financial,'input_rejections':rejections,
        'source_sha256':sha(root/'integral_auxiliary_quadrature.py'),
        'norm_input_sha256':sha(root/'results/gaussian_derivative_norms.json'),
        'proof_review':'Beta weighted first cell moment is nonzero; affine sharpness confirms its h^2/6 term. Uniform quadratic sharpness confirms h^2/24. Financial factors H/2, U powers and strike H/K chain rule independently checked.',
        'scope':'Established midpoint error mechanism specialised to actual emitted laws and separately certified finite financial moments; no generic theorem priority claim'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',default='/private/tmp/qbp-integral-value-step')
    parser.add_argument('--out',default=str(Path(__file__).resolve().parents[1]/'results/integral_auxiliary_quadrature_audit.json'))
    args=parser.parse_args();result=evaluate(args.root);result['audit_source_sha256']=sha(__file__)
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'cell_laws':len(result['exact_cell_laws']),
        'signed_polynomials':len(result['signed_polynomial_comparisons']),
        'financial_scalings':len(result['financial_norm_and_scaling_comparisons']),'input_rejections':result['input_rejections']}))

if __name__=='__main__':main()
