"""Exact Beta/Gauss rule, fixed-parameter and scalar-kernel audit.

Published Gaussian quadrature mathematics is checked through exact quadratic
field identities and its Hermite-interpolation remainder. No pricing answer
enters the auxiliary or scalar verification.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import argparse,gzip,hashlib,json,sys,time
from flint import arb,ctx


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator

def ends(x):
    def endpoint(v):
        m,e=map(int,v.man_exp());return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
    return endpoint(x.lower()),endpoint(x.upper())
def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def Phi(x):return (1+(x/arb(2).sqrt()).erf())/2

def qadd(x,y):return x[0]+y[0],x[1]+y[1]
def qmul(x,y):return x[0]*y[0]+6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def qpow(x,n):
    out=(F(1),F(0))
    for _ in range(n):out=qmul(out,x)
    return out

def independent_scalar(c):
    p={k:F(v) for k,v in c['params'].items()};name=c['name'];lam=ball(p['lambda']) if 'lambda' in p else None
    K=ball(p.get('K',100));offset=ball(c['output_offset']);records=[]
    for n in c['nodes']:
        x=ball(n['input'])
        if name=='exp':v=x.exp()
        elif name=='log':v=x.log()
        elif name=='inv_square':v=1/(x*x)
        elif name=='square':v=x*x
        elif name in ('boundary_h','boundary_density'):
            value=K*phi(x/lam)/lam
            if name=='boundary_h':
                H=ball(p['cap']);ell=(1+H/K).log();value-=(K+H)*phi((x-ell)/lam)/lam
                v=(1+value/ball(p['M2']))/2
            else:v=value/ball(p['M2'])
        elif name in ('cap','cap_prime','vanilla_call','vanilla_call_prime'):
            growth=K*(x+lam*lam/2).exp();lower=x/lam
            if name.startswith('vanilla'):
                v=growth*Phi(lower+lam)
                if name=='vanilla_call':v-=K*Phi(lower)
            else:
                H=ball(p['cap']);ell=(1+H/K).log();high=(x-ell)/lam
                v=growth*(Phi(lower+lam)-Phi(high+lam))
                if name=='cap':v-=K*(Phi(lower)-Phi(high));v+=H*Phi(high)
                if p.get('normalised',F(0)):v/=H
        else:raise ValueError(('Unexpected scalar',name))
        v-=offset;rl,rh=ends(v);lo,hi=F(n['lower']),F(n['upper'])
        assert lo<=rh and rl<=hi,('Scalar enclosure discrepancy',name,n['input'])
        error=max(abs(F(n['integer'],1<<c['q'])-lo),abs(F(n['integer'],1<<c['q'])-hi))
        assert error==F(n['error_upper']);records.append(error)
    largest=max(records);assert largest==F(c['nodal_error_upper'])
    h=F(c['span'])/(1<<c['table_bits']);rem=F(c['second_derivative_upper'])*h*h/8;ulp=F(1,1<<c['q'])
    assert rem==F(c['linear_interpolation_error_upper']) and ulp==F(c['interpolation_floor_error_upper'])
    assert largest+rem+ulp==F(c['output_error_upper'])
    span=F(c['span']);assert not span.numerator&(span.numerator-1) and not span.denominator&(span.denominator-1)
    assert c['fractional_remainder_bits']==c['q']+span.numerator.bit_length()-span.denominator.bit_length()-c['table_bits']
    assert len(c['nodes'])==(1<<c['table_bits'])+1
    values=[n['integer'] for n in c['nodes']];assert min(values)==c['value_min_integer'] and max(values)==c['value_max_integer']
    assert all(type(v) is int and v>=0 for v in values)
    return len(values)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',default='/private/tmp/qbp-boundary-correction-step');p.add_argument('--out',default='/private/tmp/qbp-boundary-grid-audit/results/boundary_auxiliary_audit.json');args=p.parse_args()
    root=Path(args.root);start=time.perf_counter();sys.path.insert(0,str(root));ctx.prec=768
    import boundary_auxiliary_certificate as module
    import boundary_value_certificate as scalars
    before={f:sha(root/f) for f in ('boundary_auxiliary_certificate.py','boundary_value_certificate.py','boundary_value_oracle.py','boundary_derivative_certificate.py','results/boundary_derivative_norms.json')}
    # Work in Q(sqrt(6)), retaining exact first four moments and the positive
    # squared monic norm. The rational remainder is independent of rounding.
    zm,zp=(F(2,5),F(-1,10)),(F(2,5),F(1,10))
    wm,wp=(F(1,2),F(1,18)),(F(1,2),F(-1,18))
    moments=[]
    for k in range(5):
        moment=qadd(qmul(wm,qpow(zm,k)),qmul(wp,qpow(zp,k)));assert moment[1]==0
        truth=F(2,(k+1)*(k+2));error=truth-moment[0]
        if k<=3:assert error==0
        else:assert error==F(1,300)
        moments.append({'degree':k,'exact_Beta_moment':str(truth),'Gauss_moment':str(moment[0]),'difference':str(error)})
    polynomial=[F(1,10),F(-4,5),F(1)];square=[F(0)]*5
    for i,a in enumerate(polynomial):
        for j,b in enumerate(polynomial):square[i+j]+=a*b
    norm=sum(v*F(2,(j+1)*(j+2)) for j,v in enumerate(square));assert norm==F(1,300)
    # Check midpoint-min distribution and the exact within-cell transport
    # moments. Their Taylor bound is independent of a financial price.
    beta_records=[]
    for bits in (1,3,8):
        N=1<<bits;h=F(1,N);sumfirst=F(0);sumsecond=F(0)
        for j in range(N):
            a=j*h;b=(j+1)*h;c=(a+b)/2;mass=2*((b-a)-(b*b-a*a)/2)
            assert mass==F(2*N-2*j-1,N*N)
            first=2*((b*b-a*a)/2-(b**3-a**3)/3)-c*mass
            second=2*((b**3-a**3)/3-(b**4-a**4)/4)-2*c*2*((b*b-a*a)/2-(b**3-a**3)/3)+c*c*mass
            assert first==-h**3/6 and second>=0;sumfirst+=first;sumsecond+=second
        assert sumfirst==-h*h/6 and sumsecond==h*h/12
        for k in range(7):
            grid=sum(F(2*N-2*j-1,N*N)*(F(2*j+1,2*N)**k) for j in range(N));truth=F(2,(k+1)*(k+2))
            derivative=F(k);second=F(k*(k-1));bound=h*h*(derivative/6+second/24)
            assert abs(grid-truth)<=bound
        beta_records.append({'bits':bits,'signed_first_moment_sum':str(sumfirst),'second_moment_sum':str(sumsecond),'polynomial_degrees_checked':list(range(7))})
    constants=[]
    for q in (8,16,28,40):
        c=scalars.gauss_constants(q);ctx.prec=768;root6=arb(6).sqrt();truths={'z_minus':(4-root6)/10,'z_plus':(4+root6)/10,'w_minus':arb(1)/2+root6/18,'w_plus':arb(1)/2-root6/18}
        for key,true in truths.items():
            row=c[key];lo,hi=F(row['ideal_lower']),F(row['ideal_upper']);tl,th=ends(true);assert lo<=tl<=th<=hi
            actual=F(row['actual']);assert actual.denominator&(actual.denominator-1)==0 and actual.denominator<=1<<q
            assert 0<=actual<=1 and F(row['error_upper'])==max(abs(actual-lo),abs(actual-hi))
        assert F(c['w_minus']['actual'])+F(c['w_plus']['actual'])==1
        constants.append(c)
    norms=json.loads((root/'results/boundary_derivative_norms.json').read_text());domainrows=json.loads((root/'results/integral_finite_input_domains.json').read_text())['domains'];domains={x['model_id']:x for x in domainrows}
    momentdata=json.loads((root/'results/boundary_finite_moments.json').read_text());moments_by_mid={x['model_id']:x for x in momentdata['records']}
    modeldata=json.loads((root/'results/generalisation_protocol.json').read_text())['models'];models={m['model_id']:m for m in modeldata}
    original=json.loads((root/'results/original_q28_quantile_nodes.json').read_text())['certificate'];finitechecks=0;calls=[]
    # Direct nonnegative U powers are independently evaluated on actual Y.
    for mid,row in moments_by_mid.items():
        ctx.prec=768;m=models[mid];w=list(map(F,m['weights']));mu=list(map(F,m['log_drifts']));load=[F(r[0]) for r in m['loadings']];mg=sum(a*b for a,b in zip(w,mu));rg=sum(a*b for a,b in zip(w,load));d=[v-mg for v in mu];s=[v-rg for v in load]
        sums={r:arb(0) for r in (2,3,4,6)}
        for node in original['nodes']:
            y=F(node['actual_coordinate']);U=sum((ball(wi)*ball(di+si*y).exp() for wi,di,si in zip(w,d,s)),arb(0))-1
            for r in sums:sums[r]+=U**r
        for r,total in sums.items():
            ref=total/1024;rl,rh=ends(ref);saved=row['centred_U_moments'][str(r)];assert F(saved['lower_exact_fraction'])<=rl<=rh<=F(saved['upper_exact_fraction']);finitechecks+=1
        assert row['classical_conditional_exponential_evaluations']==2048
    for normrow in norms['records']:
        lam=F(normrow['lambda']);cap=None if normrow['payoff']=='vanilla' else F(normrow['cap']);mid=next(mid for mid,m in models.items() if m['assets']==2 and F(m['lambda'])==lam)
        mu={r:F(moments_by_mid[mid]['centred_U_moments'][str(r)]['upper_exact_fraction']) for r in (2,3,4,6)};U=F(domains[mid]['U_star'])
        M2=F(normrow['M2_upper']);M3=F(normrow['C_r_upper']['3']);M4=F(normrow['C_r_upper']['4']);M6=F(normrow['C_r_upper']['6'])
        for kind in ('boundary_r2','stable_gauss_r2','raw_direct_r2','positive_r1'):
            b=module.auxiliary_bias(lam,cap,U,mu[2],mu[3],kind=kind,MU4=mu[4],MU6=mu[6]);ctx.prec=768
            if kind=='boundary_r2':expected=(M3*mu[3]/12+M4*mu[4]/48)/(1<<16)
            elif kind=='stable_gauss_r2':
                c=b['constants'];node=sum(F(c['w_'+s]['ideal_upper'])*F(c['z_'+s]['error_upper']) for s in ('minus','plus'));weight=sum(F(c['w_'+s]['error_upper']) for s in ('minus','plus'))
                trunc=M6*mu[6]/14400;encoded=M3*mu[3]*node/2+M2*mu[2]*weight/2;expected=trunc+encoded
                assert F(b['Gauss_remainder_upper'])==trunc and F(b['fixed_parameter_bias_upper'])==encoded
                assert not b['M6_support_bound_used']
            else:expected=F(0)
            assert F(b['bias_upper'])==expected
            calls.append({'lambda':str(lam),'cap':None if cap is None else str(cap),'kind':kind,'bias_upper':str(expected),'decimal_diagnostic':float(expected)})
    # Each new scalar certificate is independently evaluated at higher
    # precision. The complete graph/circuit expansion is audited separately.
    seedpath=root/'results/boundary_production_profiles.json.gz';family=json.load(gzip.open(seedpath,'rt'));scalar_nodes=0;scalar_certificates=0;quantiles=0
    for digest,c in family['scalar_certificate_archive'].items():
        assert digest==hashlib.sha256(json.dumps(c,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        if c['schema']=='qbp-integral-quantile-v1':
            assert c==original;quantiles+=1;continue
        assert c['schema']=='qbp-boundary-scalar-v1' and c['q']==28 and c['table_bits']==12
        assert c['source_sha256']==sha(root/'boundary_value_certificate.py') and c['norm_source_sha256']==sha(root/'results/boundary_derivative_norms.json')
        ctx.prec=768;scalar_nodes+=independent_scalar(c);scalar_certificates+=1
    after={name:sha(root/name) for name in before};assert before==after
    result={'status':'passed','audit_source_sha256':sha(Path(__file__)),'source_and_norm_sha256':before,'production_archive_sha256':sha(seedpath),'finite_moment_archive_sha256':sha(root/'results/boundary_finite_moments.json'),
        'reference_precision_bits':768,'exact_Gauss_moments':moments,'squared_monic_Gauss_norm':'1/300','probability_quadrature_remainder_factor':'1/7200','monetary_quadrature_remainder_factor':'1/14400',
        'Beta_exact_cell_transport_checks':beta_records,'dyadic_Gauss_precisions_checked':[8,16,28,40],'actual_node_weight_records':constants,
        'finite_U_moment_intervals_checked':finitechecks,'auxiliary_budget_calls_checked':len(calls),'auxiliary_budgets':calls,
        'scalar_certificates_checked':scalar_certificates,'scalar_ordinates_checked':scalar_nodes,'quantile_certificates_preserved':quantiles,'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Exact Gauss-Jacobi/Beta moments and remainder factors, explicitly paid dyadic constant errors, direct nonnegative finite dispersion moments and higher-precision scalar ordinate/formula checks. Whole Clifford/T composition and complete monetary confidence plans are separately required.'}
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('auxiliary_budgets','actual_node_weight_records')}))
if __name__=='__main__':main()
