"""Certified geometric exercise-set lower and arithmetic-geometric upper.

Published bound specialised to the stipulated GBM model, with an original
project implementation of certified scalar optimisation. Mixed-sign Gaussian
covariances are handled by strict convexity of S(z); no new bound is claimed.
The capped contract is a valid independently bounded call spread, rather than
an assertion that the source paper's native capped bound is optimised.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time,platform
import flint
from flint import arb,ctx
from boundary_grid_certificate import parameters,ball,interval,serial,endpoints,Phi,phi
ROOT=Path(__file__).resolve().parent
PRECISION=512
STEPS=80
MAX_DOUBLINGS=80

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sign(x):
    lo,hi=endpoints(x)
    return 1 if lo>0 else -1 if hi<0 else 0

def bracket_monotone(function,initial_lo=F(-1),initial_hi=F(1)):
    lo,hi=F(initial_lo),F(initial_hi);tests=[]
    for k in range(MAX_DOUBLINGS):
        left,right=function(ball(lo)),function(ball(hi));tests.append({'lo':str(lo),'hi':str(hi),'left':serial(left),'right':serial(right)})
        if sign(left)<0 and sign(right)>0:return lo,hi,tests
        if sign(left)>=0:lo-=max(F(1),hi-lo)
        if sign(right)<=0:hi+=max(F(1),hi-lo)
    return None,None,tests

def bisect_monotone(function,lo,hi):
    history=[]
    for k in range(STEPS):
        mid=(lo+hi)/2;value=function(ball(mid));sgn=sign(value)
        item={'mid':str(mid),'score_interval':serial(value),'certified_sign':sgn}
        if sgn<0:lo=mid
        elif sgn>0:hi=mid
        else:
            # Preserve a bracket when the exact zero is near a dyadic midpoint.
            l=(3*lo+hi)/4;r=(lo+3*hi)/4;lv,rv=function(ball(l)),function(ball(r))
            item['quarter_left']={'x':str(l),'score':serial(lv)};item['quarter_right']={'x':str(r),'score':serial(rv)}
            if sign(lv)<0 and sign(rv)>0:lo,hi=l,r
            else:history.append(item);return lo,hi,history,'retained unresolved midpoint sign'
        history.append(item)
    return lo,hi,history,'proved monotone root bracket'

class Family:
    def __init__(self,model):
        self.model=model;self.w,self.d,self.s,self.mg,self.rg,self.lam=parameters(model,False)
        self.sigma=ball(self.lam*self.lam+sum(v*v for v in self.rg)).sqrt()
        rows=model['loadings'];mu=tuple(map(F,model['log_drifts']))
        self.cov=tuple(self.lam*self.lam+sum(F(x)*y for x,y in zip(row,self.rg)) for row in rows)
        self.a=tuple(ball(x)/self.sigma for x in self.cov)
        self.forward=tuple(100*ball(wi)*ball(mi+(sum(F(x)*F(x) for x in row)+self.lam*self.lam)/2).exp() for wi,mi,row in zip(self.w,mu,rows))
        self.EA=sum(self.forward,arb(0));self.EG=100*(ball(self.mg)+self.sigma*self.sigma/2).exp()
        self.counters={'S_queries':0,'Sprime_queries':0,'score_exponential_evaluations':0,'exercise_CDF_evaluations':0,'geometric_call_CDF_evaluations':0}
        signs={sign(ball(q)) for q in self.cov}
        self.classification='mixed' if -1 in signs and 1 in signs else 'nonnegative' if 1 in signs else 'nonpositive' if -1 in signs else 'zero'
        self.minimum=None
        if self.classification=='mixed':
            lo,hi,search=bracket_monotone(self.Sprime)
            if lo is None:self.minimum={'status':'minimum bracket unresolved','bracket_search':search}
            else:
                lo,hi,history,status=bisect_monotone(self.Sprime,lo,hi)
                self.minimum={'status':status,'lo':str(lo),'hi':str(hi),'S_minimum_interval':serial(self.S(interval(lo,hi))),
                  'left_derivative_interval':serial(self.Sprime(ball(lo))),'right_derivative_interval':serial(self.Sprime(ball(hi))),
                  'bracket_search':search,'bisection':history}

    def S(self,z):
        self.counters['S_queries']+=1;self.counters['score_exponential_evaluations']+=len(self.a)
        return sum((A*(a*z-a*a/2).exp() for A,a in zip(self.forward,self.a)),arb(0))
    def Sprime(self,z):
        self.counters['Sprime_queries']+=1;self.counters['score_exponential_evaluations']+=len(self.a)
        return sum((A*a*(a*z-a*a/2).exp() for A,a in zip(self.forward,self.a)),arb(0))
    def L(self,z,k):
        self.counters['exercise_CDF_evaluations']+=len(self.a)+1
        return sum((A*Phi(a-z) for A,a in zip(self.forward,self.a)),arb(0))-ball(k)*Phi(-z)
    def geometric(self,k):
        x=(ball(self.mg)-(ball(k)/100).log())/self.sigma
        self.counters['geometric_call_CDF_evaluations']+=2
        return self.EG*Phi(x+self.sigma)-ball(k)*Phi(x)
    def call(self,k):
        k=F(k);start=time.perf_counter();endpoint=self.EA-ball(k)
        # Every exact fixed threshold yields a valid lower bound irrespective
        # of whether the global optimisation branch is resolved.
        fixed=self.L(ball(0),k);candidates=[{'kind':'positive_part_of_forward_endpoint','value':serial(endpoint)},
          {'kind':'zero_endpoint','value':serial(arb(0))},{'kind':'fixed_z_zero','value':serial(fixed)}]
        root=None;status='endpoint optimum';complete=True
        score=lambda z:self.S(z)-ball(k)
        if self.classification=='nonnegative':
            zero_mass=sum((A for A,q in zip(self.forward,self.cov) if q==0),arb(0))
            margin=ball(k)-zero_mass
            if sign(margin)>0:
                lo,hi,search=bracket_monotone(score)
                if lo is not None:
                    lo,hi,history,bst=bisect_monotone(score,lo,hi)
                    root={'lo':str(lo),'hi':str(hi),'left_score':serial(score(ball(lo))),'right_score':serial(score(ball(hi))),
                      'bracket_search':search,'bisection':history,'bisection_status':bst,'monotonicity':'Every a_i>=0, at least one positive; S strictly increasing.'}
                else:status='root bracket unresolved; fixed-threshold lower';complete=False
            elif sign(margin)==0 and any(q==0 for q in self.cov):
                # The interval cannot decide a zero-covariance asymptotic mass.
                status='asymptotic zero-covariance comparison unresolved';complete=False
        elif self.classification=='mixed':
            minimum=self.minimum
            if 'S_minimum_interval' not in minimum:status='minimum unresolved; fixed-threshold lower';complete=False
            else:
                sm=interval(F(minimum['S_minimum_interval']['lower_exact_fraction']),F(minimum['S_minimum_interval']['upper_exact_fraction']))
                relation=sign(sm-ball(k))
                if relation<0:
                    lo=F(minimum['hi']);right=lo+1;search=[]
                    if sign(score(ball(lo)))>=0:status='right-of-minimum target sign unresolved';complete=False
                    else:
                        for j in range(MAX_DOUBLINGS):
                            value=score(ball(right));search.append({'x':str(right),'score':serial(value)})
                            if sign(value)>0:break
                            right=lo+2*(right-lo)
                        else:status='upper-root bracket unresolved';complete=False
                        if complete:
                            hi=right;lo,hi,history,bst=bisect_monotone(score,lo,hi)
                            root={'lo':str(lo),'hi':str(hi),'left_score':serial(score(ball(lo))),'right_score':serial(score(ball(hi))),
                              'bracket_search':search,'bisection':history,'bisection_status':bst,
                              'monotonicity':'S strictly convex; bracket lies to the right of its independently certified minimum. Upper root is finite local maximum of L; lower root is a local minimum.'}
                elif relation==0:status='tangent/minimum target comparison unresolved; fixed-threshold lower';complete=False
                # If relation>0, S everywhere exceeds target and L decreases.
        if root:
            critical=self.L(interval(F(root['lo']),F(root['hi'])),k);mid=(F(root['lo'])+F(root['hi']))/2
            candidates.extend([{'kind':'upper_critical_root','value':serial(critical)},{'kind':'root_midpoint_fixed_threshold','threshold':str(mid),'value':serial(self.L(ball(mid),k))}])
            status='complete exercise-family optimisation with endpoint comparison'
        lower=max(F(0),*(F(c['value']['lower_exact_fraction']) for c in candidates))
        family_upper=max(F(0),*(F(c['value']['upper_exact_fraction']) for c in candidates)) if complete else None
        geometric=self.geometric(k);upper_bound=geometric+self.EA-self.EG;ulo,uhi=endpoints(upper_bound)
        if lower>uhi:raise ArithmeticError('Optimised lower and proven AG upper contradict')
        return {'model_id':self.model['model_id'],'strike':str(k),'covariance_sign_classification':self.classification,
          'critical_root':root,'optimisation_status':status,'optimisation_complete':complete,
          'candidate_values':candidates,'exercise_family_supremum_interval':None if family_upper is None else serial(interval(lower,family_upper)),
          'call_price_interval':serial(interval(lower,uhi)),'geometric_call_interval':serial(geometric),
          'AG_gap_interval':serial(self.EA-self.EG),'AG_upper_interval':serial(upper_bound),
          'computation_seconds_diagnostic':time.perf_counter()-start}


def main():
    ctx.prec=PRECISION;start=time.perf_counter();protocol=ROOT/'results/generalisation_protocol.json';models=json.loads(protocol.read_text())['models']
    calls=[];records=[];families=[]
    for model in models:
        family=Family(model);cases={}
        for strike in (100,101,110,125):cases[strike]=family.call(strike);calls.append(cases[strike])
        base=cases[100]['call_price_interval'];blo,bhi=map(F,(base['lower_exact_fraction'],base['upper_exact_fraction']))
        for cap in (F(1),F(10),F(25),None):
            if cap is None:lo,hi=blo,bhi
            else:
                uppercall=cases[int(100+cap)]['call_price_interval'];alo,ahi=map(F,(uppercall['lower_exact_fraction'],uppercall['upper_exact_fraction']))
                lo=max(F(0),blo-ahi);hi=min(cap,bhi-alo)
            if lo>hi:raise ArithmeticError('Call-spread enclosure empty')
            radius=(hi-lo)/2;price=interval(lo,hi);discount=ball(-F(model['risk_free_rate'])*F(model['maturity'])).exp()
            records.append({'model_id':model['model_id'],'payoff':'vanilla' if cap is None else 'capped','cap':None if cap is None else str(cap),
              'continuous_price_interval':serial(price),'undiscounted_radius_upper':str(radius),'certifies_undiscounted_cent':radius<=F(1,100),
              'discounted_price_interval':serial(price*discount),'bound_scope':'Optimised geometric exercise-set lower and arithmetic-geometric call upper. Cap uses independent valid call-spread bounds, not a claim of optimised native cap bounds.'})
        families.append({'model_id':model['model_id'],'sigma_G_interval':serial(family.sigma),'geometric_covariances_exact':list(map(str,family.cov)),
          'standardised_covariance_intervals':list(map(serial,family.a)),'weighted_forward_intervals':list(map(serial,family.forward)),
          'EA_interval':serial(family.EA),'EG_interval':serial(family.EG),'classification':family.classification,'minimum_certificate':family.minimum,'operation_counts':family.counters})
    result={'schema':'qbp-caldana-geometric-conditioned-benchmark-v1','status':'completed with every bound and unresolved branch retained','precision_bits':PRECISION,
      'bisection_steps':STEPS,'families':families,'call_bounds':calls,'records':records,'generation_seconds_diagnostic':time.perf_counter()-start,
      'source_sha256':sha(__file__),'protocol_sha256':sha(protocol),'dependencies':{n:sha(ROOT/n) for n in ('boundary_grid_certificate.py',)},
      'source_attribution':'Caldana et al. (2016), accepted manuscript supplied by root, Propositions1.1/2.1, Corollary2.2 and Eq27. Formula/optimisation specialisation and exact covariance classification are project implementation; no new generic bound claimed.',
      'python_version':platform.python_version(),'python_flint_version':flint.__version__,
      'scope':'Four strikes for each of nine established frozen GBM models. Geometric conditional event optimisation, including mixed-sign covariances, and valid AG upper. No basket-price answer, financial grid, price simulation or new dataset. Cases failing cent accuracy are retained.'}
    out=ROOT/'results/caldana_conditioning_benchmark.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print('Caldana/AG specialisationcertified',sum(r['certifies_undiscounted_cent'] for r in records),'of',len(records),'seconds',result['generation_seconds_diagnostic'],flush=True)
if __name__=='__main__':main()
