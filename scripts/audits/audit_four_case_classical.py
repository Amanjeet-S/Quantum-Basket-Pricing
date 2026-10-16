"""Independent full-score and Gaussian-transfer replay.

No producer or producer numerical helper is imported.  Exact arithmetic,
fresh 1024-bit cells and the 256-bit difference of conditional vanilla
calls check the complete retained evidence, including every sampled score.
"""
from array import array
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys,time
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(x):
    x=F(x);return arb(x.numerator)/x.denominator
def endpoint(x,up):
    m,e=map(int,(x.upper() if up else x.lower()).man_exp())
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def bounds(x):
    if not x.is_finite():raise ArithmeticError('Nonfinite independent ball')
    return endpoint(x,False),endpoint(x,True)
def Phi(x):return (1+(x/arb(2).sqrt()).erf())/2
def phi(x):return (-x*x/2).exp()/(2*arb.pi()).sqrt()
def check_ball(record,x,name):
    lo,hi=bounds(x)
    if not F(record['lower_exact_fraction'])<=lo<=hi<=F(record['upper_exact_fraction']):
        raise AssertionError('Producer interval does not enclose higher-precision '+name)
def fraction_record(a,b):return {'lower_exact_fraction':str(a),'upper_exact_fraction':str(b)}

def fresh_cells(stored):
    ctx.prec=1024;L=A(6);low=Phi(-L);a=Phi(L)-low;h=a/1024
    boundaries=[-L]
    for k in range(1,1024):boundaries.append(arb(2).sqrt()*(2*(low+a*k/1024)-1).erfinv())
    boundaries.append(L)
    b1=F(0);b2=arb(0);nodes=[];checks=0
    for k,old in enumerate(stored['actual_cells']):
        midpoint=arb(2).sqrt()*(2*(low+a*A(F(2*k+1,2048)))-1).erfinv()
        ml,mh=bounds(midpoint);a0=ml*(1<<28)+F(1,2);a1=mh*(1<<28)+F(1,2)
        il=a0.numerator//a0.denominator;ih=a1.numerator//a1.denominator
        if il!=ih:raise AssertionError('Ambiguous independent midpoint rounding')
        c=F(il,1<<28);left,right=boundaries[k:k+2]
        if old['cell']!=k or F(old['actual_node'])!=c or old['actual_node_integer']!=il:
            raise AssertionError('Actual grid node changed')
        if not bounds(left)[1]<=c<=bounds(right)[0]:raise AssertionError('Node outside independently reconstructed cell')
        pa,pb=phi(left),phi(right)
        first=pa-pb-A(c)*h
        second=(1+A(c*c))*h+(left-2*A(c))*pa-(right-2*A(c))*pb
        if bounds(second)[0]<0:raise AssertionError('Independent second cell moment negative')
        expected={'left':left,'right':right,'mass':h,'quantile_midpoint':midpoint,
            'first_moment':first,'second_moment':second,
            'left_cdf':low+a*k/1024,'right_cdf':low+a*(k+1)/1024}
        for key,value in expected.items():check_ball(old[key],value,str(k)+' '+key);checks+=1
        b1+=max(map(abs,bounds(first)));b2+=second/2;nodes.append(c)
    check_ball(stored['mass_interval'],a,'retained Gaussian mass')
    check_ball(stored['rho2_mass_interval'],b2,'rho2 total mass')
    if b1>F(stored['rho1_mass_upper']):raise AssertionError('rho1 upper too small')
    return a,b1,b2,nodes,checks+3

class AlternateScores:
    def __init__(self,model,nodes):
        ctx.prec=256;self.lam=A(F(model['lambda']));self.ell=A(F(101,100)).log()
        self.constants=[A(F(w))*(A(F(mu))).exp() for w,mu in zip(model['weights'],model['log_drifts'])]
        self.tables=[[[A(F(r)*c).exp() for c in nodes] for r in row] for row in model['loadings']]
    def value(self,labels):
        ctx.prec=256;sumexp=arb(0)
        for const,row in zip(self.constants,self.tables):
            v=const
            for j,k in enumerate(labels):v=v*row[j][k]
            sumexp=sumexp+v
        m=sumexp.log();growth=100*(m+self.lam*self.lam/2).exp()
        x=m/self.lam;xh=(m-self.ell)/self.lam
        ordinary=growth*Phi(x+self.lam)-100*Phi(x)
        upper_strike=growth*Phi(xh+self.lam)-101*Phi(xh)
        lo,hi=bounds(ordinary-upper_strike)
        lo=max(F(0),lo);hi=min(F(1),hi)
        if lo>hi:raise AssertionError('Alternate formula violates bounded-payoff premise')
        return lo,hi

def binary(record,typecode):
    path=RESULTS/record['file']
    if sha(path)!=record['sha256'] or path.stat().st_size!=record['bytes']:raise AssertionError('Changed binary trace '+str(path))
    a=array(typecode);a.frombytes(path.read_bytes())
    if sys.byteorder!='little':a.byteswap()
    expected=record.get('count',record.get('coordinate_count'))
    if len(a)!=expected:raise AssertionError('Binary record count differs')
    return a

def verify_complete(record,meanlo,meanhi,t):
    Fl=F(t['retained_mass_interval']['lower_exact_fraction']);Fh=F(t['retained_mass_interval']['upper_exact_fraction']);beta=F(t['financial_grid_transfer_upper'])
    lo=max(F(0),Fl*meanlo-beta);hi=min(F(1),Fh*meanhi+beta+1-Fl)
    if record['undiscounted_continuous_price_interval']!=fraction_record(lo,hi):raise AssertionError('Complete continuous endpoints differ')
    if F(record['undiscounted_radius_upper'])!=(hi-lo)/2 or record['certifies_undiscounted_cent']!=((hi-lo)/2<=F(1,100)):
        raise AssertionError('Radius or tolerance admission differs')
    ctx.prec=1024;discount=A(F(-3,100)).exp();check_ball(record['discount_interval'],discount,'discount')
    dl=F(record['discount_interval']['lower_exact_fraction']);dh=F(record['discount_interval']['upper_exact_fraction'])
    if record['discounted_continuous_price_interval']!=fraction_record(dl*lo,dh*hi):raise AssertionError('Discounted endpoints differ')
    if F(record['discounted_radius_upper'])!=(dh*hi-dl*lo)/2:raise AssertionError('Discounted radius differs')

def main():
    start=time.perf_counter();full=RESULTS/'four_case_classical_pricing.json';data=json.loads(full.read_text())
    protocol_path=RESULTS/'four_case_classical_protocol.json';protocol=json.loads(protocol_path.read_text())
    if sha(protocol_path)!=data['protocol_sha256'] or sha(ROOT/'four_case_classical_pricing.py')!=protocol['producer_sha256']:
        raise AssertionError('Frozen producer/protocol hash changed')
    if data['input_sha256']!=protocol['input_sha256']:raise AssertionError('Provenance copies differ')
    for name,value in protocol['input_sha256'].items():
        if sha(ROOT/name)!=value:raise AssertionError('Input hash differs '+name)
    source=json.loads((RESULTS/'dispersion_separation_protocol.json').read_text())
    norms=json.loads((RESULTS/'boundary_derivative_norms.json').read_text())
    expected=list(protocol['case_ids'])
    if [r['model_id'] for r in data['records']]!=expected or len(expected)!=4:raise AssertionError('Case selection differs')
    a,b1,b2,nodes,cell_checks=fresh_cells(data['gaussian_grid'])
    output=[];sample_scores=0;deterministic_scores=0
    for row in data['records']:
        m=next(m for m in source['models'] if m['model_id']==row['model_id']);p=m['outer_dimension'];lam=F(m['lambda'])
        if sum(map(F,m['weights']))!=1 or min(map(F,m['weights']))<=0:raise AssertionError('Invalid basket weights')
        rate=F(m['risk_free_rate'])-F(m['dividend_yield'])
        for drift,loading in zip(m['log_drifts'],m['loadings']):
            if F(drift)+(lam*lam+sum(F(v)**2 for v in loading))/2!=rate:raise AssertionError('No-arbitrage drift equality fails')
        nr=next(n for n in norms['records'] if n['payoff']=='capped' and F(n['cap'])==1 and F(n['lambda'])==lam and F(n['K'])==100)
        M1=F(nr['M1_upper']);M2=F(nr['M2_upper']);t=row['transfer_certificate']
        if F(t['global_M1_upper'])!=M1 or F(t['global_M2_upper'])!=M2:raise AssertionError('Shared payoff norm changed')
        ctx.prec=1024;beta=arb(0)
        for j,dr in enumerate(t['coordinate_derivative_bounds']):
            Aj=max(abs(F(r[j])) for r in m['loadings']);d1=M1*Aj;d2=(M1+M2)*Aj*Aj
            if dr!= {'coordinate':j,'maximum_full_asset_loading':str(Aj),'first_derivative_upper':str(d1),'second_derivative_upper':str(d2)}:
                raise AssertionError('Full-loading derivative envelope differs')
            beta+=A(d1*b1)+A(d2)*b2
        beta*=a**(p-1)
        if bounds(beta)[1]>F(t['financial_grid_transfer_upper']):raise AssertionError('Transfer upper fails 1024-bit audit')
        check_ball(t['retained_mass_interval'],a**p,'full retained mass')
        if F(t['bounded_tail_upper'])!=1-F(t['retained_mass_interval']['lower_exact_fraction']):raise AssertionError('Full bounded tail differs')
        score=AlternateScores(m,nodes);det=row['complete_exactsum'];det_exact_mean=None
        if det is not None:
            if p!=1:raise AssertionError('Unexpected deterministic tensor sum')
            lower=binary(det['scores_lower'],'Q');upper=binary(det['scores_upper'],'Q')
            if len(lower)!=1024 or len(upper)!=1024:raise AssertionError('Incomplete exact finite sum')
            for k,(l,u) in enumerate(zip(lower,upper)):
                lo,hi=score.value((k,))
                if not F(l,1<<40)<=lo<=hi<=F(u,1<<40):raise AssertionError('Exact-sum score enclosure fails '+str(k))
            lo=F(sum(lower),1024*(1<<40));hi=F(sum(upper),1024*(1<<40))
            if det['finite_grid_mean_interval']!=fraction_record(lo,hi):raise AssertionError('Exact-sum rational mean differs')
            verify_complete(det,lo,hi,t);deterministic_scores+=1024;det_exact_mean=(lo,hi)
        mc=row['complete_monte_carlo'];words=binary(mc['random_words'],'H')
        lower=binary(mc['scores_lower'],'Q');upper=binary(mc['scores_upper'],'Q')
        n=protocol['draws_per_case']
        if n!=100000 or len(words)!=n*p or len(lower)!=n or len(upper)!=n:raise AssertionError('Sampling rule/count differs')
        for i,(l,u) in enumerate(zip(lower,upper)):
            labels=tuple(words[i*p+j]&1023 for j in range(p));lo,hi=score.value(labels)
            if not F(l,1<<40)<=lo<=hi<=F(u,1<<40):raise AssertionError('Sample score enclosure fails '+str(i))
        sample_scores+=n;ml=F(sum(lower),n*(1<<40));mh=F(sum(upper),n*(1<<40));eps=F(protocol['statistical_epsilon'])
        if mc['sample_mean_interval']!=fraction_record(ml,mh):raise AssertionError('Exact sample mean differs')
        lo=max(F(0),ml-eps);hi=min(F(1),mh+eps)
        if mc['finite_grid_mean_confidence_interval']!=fraction_record(lo,hi):raise AssertionError('Finite-grid confidence endpoints differ')
        ctx.prec=1024;failure=2*A(-2*n*eps*eps).exp()
        if bounds(failure)[1]>F(mc['hoeffding_failure_upper']) or F(mc['hoeffding_failure_upper'])>F(1,3600):raise AssertionError('Hoeffding confidence allocation fails')
        if F(mc['maximum_outward_score_width'])!=F(max(u-l for u,l in zip(upper,lower)),1<<40):raise AssertionError('Score width differs')
        verify_complete(mc,lo,hi,t)
        primary=det or mc
        output.append({'model_id':row['model_id'],'sample_score_replays':n,
          'exactsum_score_replays':0 if det is None else 1024,
          'primary_undiscounted_price_interval':primary['undiscounted_continuous_price_interval'],
          'primary_undiscounted_radius_upper':primary['undiscounted_radius_upper'],
          'primary_certifies_cent':primary['certifies_undiscounted_cent'],
          'n2_realised_monte_carlo_interval_contains_deterministic_finite_grid_interval':None if det_exact_mean is None else lo<=det_exact_mean[0]<=det_exact_mean[1]<=hi})
        print('audited every score',row['model_id'],flush=True)
    if sample_scores!=400000 or deterministic_scores!=1024:raise AssertionError('Complete replay count differs')
    result={'schema':'qbp-four-case-complete-classical-independent-audit-v1','status':'passed',
      'auditor_sha256':sha(Path(__file__)),'full_result_sha256':sha(full),'protocol_sha256':sha(protocol_path),
      'producer_sha256':sha(ROOT/'four_case_classical_pricing.py'),'input_sha256':protocol['input_sha256'],
      'independent_transfer_precision_bits':1024,'independent_score_precision_bits':256,
      'producer_imports':False,'alternate_score_formula':'Difference of separately evaluated conditional vanilla calls at strikes 100 and 101',
      'gaussian_cell_and_mass_interval_checks':cell_checks,
      'sample_score_replays':sample_scores,'deterministic_exactsum_score_replays':deterministic_scores,
      'all_four_primary_prices_certify_cent':all(r['primary_certifies_cent'] for r in output),
      'records':output,'confidence_limit':'Arithmetic replay does not establish physical independence or uniformity of operating-system random words. The simultaneous sampling guarantee is conditional on the explicit ideal-word premise.',
      'data':'All financial parameters and traces are synthetic; no LSEG observations or derived market data.',
      'audit_seconds_diagnostic':time.perf_counter()-start}
    target=RESULTS/'four_case_classical_independent_audit.json'
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('independent audit passed',flush=True)

if __name__=='__main__':main()
