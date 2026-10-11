"""Independent Gaussian cell and automatic-differentiation transfer audit."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import argparse,gzip,hashlib,json,sys,time
from flint import arb,arb_series,ctx

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/boundary_grid_independent_audit.json'
PREC=768

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(x):
    x=F(x);return arb(x.numerator)/x.denominator

def ends(x):
    def ep(z):
        m,e=map(int,z.man_exp());return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
    return ep(x.lower()),ep(x.upper())
def saved_ball(d):
    lo,hi=F(d['lower_exact_fraction']),F(d['upper_exact_fraction'])
    return ball((lo+hi)/2)+arb(0,ball((hi-lo)/2).upper())
def contains(d,ref,label):
    lo,hi=F(d['lower_exact_fraction']),F(d['upper_exact_fraction']);rl,rh=ends(ref)
    assert lo<=rl<=rh<=hi,(label,lo,rl,rh,hi)
def overlap(d,ref,label):
    lo,hi=F(d['lower_exact_fraction']),F(d['upper_exact_fraction']);rl,rh=ends(ref)
    assert lo<=rh and rl<=hi,(label,lo,rl,rh,hi)
def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def Phi(x):return (1+(x/arb(2).sqrt()).erf())/2

def model_parameters(m):
    w=tuple(map(F,m['weights']));loads=tuple(F(x[0]) for x in m['loadings']);mu=tuple(map(F,m['log_drifts']));lam=F(m['lambda'])
    rg=sum(a*b for a,b in zip(w,loads));mg=sum(a*b for a,b in zip(w,mu));d=tuple(a-mg for a in mu);s=tuple(a-rg for a in loads)
    assert len(w)==2 and sum(w)==1 and min(w)>0 and lam>0
    assert sum(w[i]*d[i] for i in range(2))==0 and sum(w[i]*s[i] for i in range(2))==0
    rate=F(m['risk_free_rate'])-F(m['dividend_yield']);assert all(mu[i]+(loads[i]**2+lam**2)/2==rate for i in range(2))
    return w,d,s,mg,rg,lam


def kernel(m,lam,cap,payoff):
    K=100;lb=ball(lam);x=m/lb;e=(m+lb*lb/2).exp()
    if payoff=='vanilla':return K*(e*Phi(x+lb)-Phi(x))
    H=ball(cap);a=(1+H/K).log();xh=(m-a)/lb
    return K*e*(Phi(x+lb)-Phi(xh+lb))-K*(Phi(x)-Phi(xh))+H*Phi(xh)
def prime(m,lam,cap,payoff):
    lb=ball(lam);e=100*(m+lb*lb/2).exp();x=m/lb
    if payoff=='vanilla':return e*Phi(x+lb)
    a=(1+ball(cap)/100).log();return e*(Phi(x+lb)-Phi((m-a)/lb+lb))


def jets(y,m,cap,payoff):
    w,d,s,mg,rg,lam=model_parameters(m);ys=arb_series([ball(y),arb(1)],prec=3)
    T=sum((ball(wi)*(ball(di)+ball(si)*ys).exp() for wi,di,si in zip(w,d,s)),arb_series([arb(0)],prec=3))
    g=ball(mg)+ball(rg)*ys;a=g+T.log();R1=kernel(a,lam,cap,payoff)-kernel(g,lam,cap,payoff)
    R2=R1-(T-1)*prime(g,lam,cap,payoff)
    return {1:(R1[1],2*R1[2]),2:(R2[1],2*R2[2])}


def boundary_data(L,N):
    low=Phi(-ball(L));M=Phi(ball(L))-low;boundaries=[-ball(L)]
    for j in range(1,N):boundaries.append(arb(2).sqrt()*(2*(low+M*j/N)-1).erfinv())
    boundaries.append(ball(L));return low,M,boundaries

def weight_primitive(x):return ((x*x+1)*Phi(x)+x*phi(x))/2

def tail_reference(model,power,L,extra=F(0)):
    w,d,s,mg,rg,lam=model_parameters(model);v=arb(0)
    # Expand each power directly in the two frozen exponentials. Gaussian
    # exponential tilting computes the outside interval, without a price.
    for j in range(power+1):
        for k in range(j+1):
            coeff=F((-1)**(power-j)*comb(power,j)*comb(j,k))*w[0]**k*w[1]**(j-k)
            a=k*s[0]+(j-k)*s[1]+extra;b=k*d[0]+(j-k)*d[1]
            tail=Phi(ball(a-L))+Phi(ball(-a-L))
            v+=ball(coeff)*ball(b+a*a/2).exp()*tail
    return v


def main():
    global ROOT,OUT
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(ROOT));parser.add_argument('--out',default=str(OUT));parser.add_argument('--original-quantile',default=None)
    args=parser.parse_args();ROOT=Path(args.root);OUT=Path(args.out)
    start=time.perf_counter();ctx.prec=PREC;ctx.cap=3
    path=ROOT/'results/boundary_financial_grid_transfer.json.gz';data=json.load(gzip.open(path,'rt'))
    assert data['source_sha256']==sha(ROOT/'boundary_grid_certificate.py')
    assert set(data['input_sha256'])=={'generalisation_protocol.json'}
    for fn,dig in data['input_sha256'].items():assert dig==sha(ROOT/'results'/fn)
    models={m['model_id']:m for m in json.loads((ROOT/'results/generalisation_protocol.json').read_text())['models']}
    old_path=Path(args.original_quantile) if args.original_quantile else ROOT/'results/original_q28_quantile_nodes.json'
    old_record=json.loads(old_path.read_text());oldq=old_record['certificate']
    assert old_record['certificate_sha256']==hashlib.sha256(json.dumps(oldq,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert oldq['schema']=='qbp-integral-quantile-v1' and oldq['q']==28 and oldq['b_y']==10 and F(oldq['L'])==6
    expected={(mid,cap,payoff,order) for mid,m in models.items() if m['outer_dimension']==1 for payoff,caps in [('capped',(1,10,25)),('vanilla',(None,))] for cap in caps for order in (1,2)}
    seen=set();basis=0;segments=0;derivative_nodes=0;cdf_checks=0;memo={};reference_cells={};checked=[]
    for record in data['records']:
        key=record['model_id'],record['cap'],record['payoff'],record['order'];assert key in expected and key not in seen;seen.add(key)
        model=models[key[0]];cap,payoff,order=key[1:];L=F(record['L']);N=1<<record['b_y'];q=record['q']
        assert L==6 and N==1024 and q==28 and record['cell_count']==1024 and record['precision_bits']==512
        if (L,N) not in reference_cells:reference_cells[(L,N)]=boundary_data(L,N)
        low,M,bounds=reference_cells[(L,N)];contains(record['mass_interval'],M,'retained mass')
        linear=arb(0);quadratic=F(0);subcount=0
        for j,cell in enumerate(record['cells']):
            assert cell['cell']==j;a,b=bounds[j],bounds[j+1];c=F(cell['actual_node']);assert cell['actual_node_integer']==int(c*(1<<q))
            assert c==F(oldq['nodes'][j]['offset_integer'],1<<q)-L
            if j==0:assert F(cell['left']['lower_exact_fraction'])==F(cell['left']['upper_exact_fraction'])==-L
            else:
                pl=Phi(ball(cell['left']['lower_exact_fraction']));ph=Phi(ball(cell['left']['upper_exact_fraction']));target=low+M*j/N
                assert ends(pl)[1]<=ends(target)[0] and ends(ph)[0]>=ends(target)[1];cdf_checks+=1
            if j==N-1:assert F(cell['right']['lower_exact_fraction'])==F(cell['right']['upper_exact_fraction'])==L
            else:
                pl=Phi(ball(cell['right']['lower_exact_fraction']));ph=Phi(ball(cell['right']['upper_exact_fraction']));target=low+M*(j+1)/N
                assert ends(pl)[1]<=ends(target)[0] and ends(ph)[0]>=ends(target)[1];cdf_checks+=1
            target=low+M*ball(F(2*j+1,2*N));ulp=F(1,1<<q)
            assert ends(Phi(ball(c-ulp/2)))[1]<=ends(target)[0] and ends(Phi(ball(c+ulp/2)))[0]>=ends(target)[1]
            assert F(cell['left']['upper_exact_fraction'])<=c<=F(cell['right']['lower_exact_fraction']);basis+=1
            mass=M/N;first=phi(a)-phi(b)-ball(c)*mass
            second=(1+ball(c*c))*mass+(a-2*ball(c))*phi(a)-(b-2*ball(c))*phi(b)
            contains(cell['mass'],mass,'cell mass');contains(cell['first_moment'],first,'cell first moment');contains(cell['second_moment'],second,'cell second moment')
            jetkey=(key[0],cap,payoff,c)
            if jetkey not in memo:memo[jetkey]=jets(c,model,cap,payoff)
            firstd,secondd=memo[jetkey][order];contains(cell['first_derivative_at_actual_node'],firstd,'AD nodal first derivative');derivative_nodes+=1
            contribution=firstd*first;contains(cell['signed_linear_contribution'],contribution,'signed linear contribution');linear+=contribution
            boundsmax=F(0);weighted=F(0);sideparts={side:[] for side in ('left','right')}
            for sub in cell['derivative_subcells']:
                side=sub['side'];u,v=F(sub['left']),F(sub['right']);assert u<v and v-u<=F(record['maximum_derivative_subinterval_width']);sideparts[side].append((u,v))
                mid=(u+v)/2;jetkey=(key[0],cap,payoff,mid)
                if jetkey not in memo:memo[jetkey]=jets(mid,model,cap,payoff)
                contains(sub['second_derivative_interval'],memo[jetkey][order][1],'AD interval midpoint second derivative')
                saved=sub['second_derivative_interval'];maximum=max(abs(F(saved['lower_exact_fraction'])),abs(F(saved['upper_exact_fraction'])))
                assert maximum==F(sub['absolute_second_derivative_upper']);boundsmax=max(boundsmax,maximum)
                e=a if side=='left' else b;ep=low+M*(j if side=='left' else j+1)/N
                weight=weight_primitive(ball(v))-weight_primitive(ball(u))-phi(e)*ball(v-u)-ep*ball((v*v-u*u)/2)
                overlap(sub['gaussian_transport_weight_mass'],weight,'transport weight mass')
                weight_upper=F(sub['gaussian_transport_weight_mass']['upper_exact_fraction']);assert weight_upper>=0
                rem=maximum*weight_upper;assert rem==F(sub['remainder_upper']);weighted+=rem;subcount+=1;segments+=1
            for side,parts in sideparts.items():
                assert parts and all(x[1]==y[0] for x,y in zip(parts,parts[1:]))
                if side=='left':assert parts[0][0]==F(cell['left']['lower_exact_fraction']) and parts[-1][1]==c
                else:assert parts[0][0]==c and parts[-1][1]==F(cell['right']['upper_exact_fraction'])
            assert boundsmax==F(cell['second_derivative_cell_upper'])
            whole=boundsmax*F(cell['second_moment']['upper_exact_fraction'])/2
            assert whole==F(cell['whole_cell_quadratic_upper']) and weighted==F(cell['weighted_transport_remainder_upper'])
            rem=min(whole,weighted);assert rem==F(cell['quadratic_remainder_upper']);quadratic+=rem
        contains(record['first_moment_sum_interval'],linear,'signed first moment sum')
        llo,lhi=F(record['first_moment_sum_interval']['lower_exact_fraction']),F(record['first_moment_sum_interval']['upper_exact_fraction'])
        linear_upper=max(abs(llo),abs(lhi));assert linear_upper==F(record['signed_first_moment_error_upper'])
        assert quadratic==F(record['quadratic_remainder_upper']) and subcount==record['derivative_subinterval_count']
        tail=record['tail_certificate'];w,d,s,mg,rg,lam=model_parameters(model)
        if order==2:
            moment=tail_reference(model,2,L);contains(tail['U_moment_outside'],moment,'outside U squared moment')
            refM=ball(100 if payoff=='vanilla' else 200+F(cap))/(ball(lam)*(2*arb.pi()).sqrt())
            norm=F(tail['boundary_M2_upper']);assert norm>=ends(refM)[1]
            tailupper=norm*F(tail['U_moment_outside']['upper_exact_fraction'])/2
        else:
            moment=tail_reference(model,1,L,rg);contains(tail['weighted_U_moment_outside'],moment,'outside tilted U moment')
            unbounded=100*ball(mg+lam*lam/2).exp()*moment;assert F(tail['vanilla_delta_tail_upper'])>=ends(unbounded)[1]
            mass_tail=2*Phi(-ball(L));Htail=F(cap)*ends(mass_tail)[1] if payoff=='capped' else None
            # The primary tail uses its own outward mass bound, whose exact
            # dyadic value is reproducible from its recorded endpoint rule.
            scale=1<<768;value=ends(2*Phi(-ball(L)))[1]*scale;k=value.numerator//value.denominator
            rounded=F(k+int(value.denominator!=1),scale)
            if payoff=='capped':assert F(tail['tail_error_upper'])>=min(Htail,F(tail['vanilla_delta_tail_upper']))
            tailupper=F(tail['tail_error_upper'])
        assert F(record['tail_upper'])==F(tail['tail_error_upper'])==tailupper
        total=linear_upper+quadratic+tailupper;assert total==F(record['financial_transfer_upper'])
        assert record['status']==('passed' if total<F(1,200) else 'retained allocation failure')
        checked.append({'model_id':key[0],'cap':cap,'payoff':payoff,'order':order,'transfer_upper':str(total),'transfer_decimal_diagnostic':float(total),'cells':N,'segments':subcount})
        print('Independently checked',key,float(total),flush=True)
    assert seen==expected
    coeffchecks=0
    for coefficient in data['coefficients']:
        m=models[coefficient['model_id']];w,d,s,mg,rg,lam=model_parameters(m);sig=(ball(lam*lam+rg*rg)).sqrt();payoff=coefficient['payoff'];cap=coefficient['cap']
        # Explicit Gaussian convolution, independent of the primary helper.
        def p0(x):
            xx=x/sig;e=100*(x+sig*sig/2).exp()
            if payoff=='vanilla':return e*Phi(xx+sig)-100*Phi(xx)
            a=(1+ball(cap)/100).log();xh=(x-a)/sig;return e*(Phi(xx+sig)-Phi(xh+sig))-100*(Phi(xx)-Phi(xh))+ball(cap)*Phi(xh)
        def p1(x):
            xx=x/sig;e=100*(x+sig*sig/2).exp()
            if payoff=='vanilla':return e*Phi(xx+sig)
            a=(1+ball(cap)/100).log();return e*(Phi(xx+sig)-Phi((x-a)/sig+sig))
        c0=p0(ball(mg));c1=sum((ball(wi)*ball(di+si*si/2).exp()*p1(ball(mg+rg*si)) for wi,di,si in zip(w,d,s)),arb(0))-p1(ball(mg))
        contains(coefficient['c0'],c0,'geometric c0');contains(coefficient['c1'],c1,'geometric c1');contains(coefficient['coefficient'],c0+c1,'geometric coefficient sum');coeffchecks+=1
    assert len(data['unsupported_scope'])==18
    result={'status':'passed','certificate_sha256':sha(path),'construction_source_sha256':sha(ROOT/'boundary_grid_certificate.py'),'audit_source_sha256':sha(Path(__file__)),
        'original_quantile_record_sha256':sha(old_path),'original_quantile_archive_sha256':old_record['original_archive_sha256'],'reference_precision_bits':PREC,'records_checked':len(checked),'CDF_brackets_checked':cdf_checks,'actual_rounded_nodes_checked':basis,
        'nodal_first_derivatives_checked_by_Arb_series_AD':derivative_nodes,'segment_second_derivatives_checked_by_Arb_series_AD':segments,
        'unique_automatic_differentiation_arguments':len(memo),'geometric_coefficient_certificates_checked':coeffchecks,'unsupported_higher_dimensional_cases_retained':18,
        'maximum_transfer_upper':str(max(F(r['transfer_upper']) for r in checked)),'records':checked,'elapsed_seconds_diagnostic':time.perf_counter()-start,
        'scope':'Independent higher-precision CDF bracketing, moment identities, actual q28 quantile rounding, closed-kernel series automatic differentiation, exact segment partition/resource-free transport accounting, tilted-Gaussian tails and geometric coefficients. No continuous arithmetic basket price is an input.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}))
if __name__=='__main__':main()
