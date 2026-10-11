"""Certified full-batch comparator based on published modified IQAE.

The confidence allocation and greedy search follow Fukuzawa et al. (2023).
A fixed rational radius, outward inverses, explicit finite ceilings and a
financial-prior bootstrap are checked consequences, not an unchanged use of
its published query coefficient. No physical quantum device is called here.
"""
from dataclasses import dataclass
from decimal import Decimal,localcontext
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import math
from certified_scalar import (ceil_fraction,floor_fraction,pi_rational_interval,
    probability_interval,inverse_probability_interval,amplitude_interval,
    sqrt_upper_fraction,INVERSE_TOLERANCE)

RADIUS=F(69,1000)
ENDPOINT_ALLOWANCE=F(1,10**60)


def quadrant(m,lo,hi):
    j=floor_fraction(m*lo)
    return j if m*hi<=j+1 else None


@lru_cache(None)
def cover_certificate():
    critical={s:[probability_interval(F(j,s)) for j in range(1,s)] for s in (3,5,7)}
    spans=[]
    for a in critical[3]:
        for b in critical[5]:
            for c in critical[7]:
                spans.append(max(a[0],b[0],c[0])-min(a[1],b[1],c[1]))
    margin=min(spans)-2*RADIUS-8*INVERSE_TOLERANCE
    assert margin>0
    return {'triple_span_lower':min(spans),'expanded_probability_span_upper':2*RADIUS+8*INVERSE_TOLERANCE,
            'strict_margin_lower':margin,'critical_enclosures':critical,'checked_triples':len(spans)}


def log_upper(value):
    with localcontext() as ctx:
        ctx.prec=80
        number=Decimal(value.numerator)/Decimal(value.denominator)
        upper=number.ln()+Decimal('1e-60')
    if not 0<upper<100:raise ValueError('Declared logarithm domain requires 0<ln(value)<100')
    return F(upper)


@dataclass(frozen=True)
class Plan:
    epsilon: F
    alpha: F
    prior_upper: F
    guarded_prior_upper: F
    gamma_upper: F
    angle_tolerance: F
    normalised_angle_width_target: F
    multiplier_upper: F
    maximum_odd_multiplier: int
    initial_multiplier: int
    maximum_rounds: int
    maximum_shots: int
    maximum_Q_calls: int
    shots_by_multiplier: tuple

    @property
    def total_shots_upper(self):return self.maximum_shots
    @property
    def total_grover_calls_upper(self):return self.maximum_Q_calls
    @property
    def grover_power_upper(self):return max(0,(self.maximum_odd_multiplier-1)//2)
    def shot_table(self):return dict(self.shots_by_multiplier)

    def cost_upper(self,a_cost,q_cost):
        """Maximum linear cost over a conservative set of admissible chains.

        These are mathematical scheduling bounds using the supplied exact
        emitted A/Q costs, not another gate compilation or asymptotic count.
        """
        if self.maximum_rounds==0:return {'cost_upper':0,'witness_multipliers':[],'witness_shots':0,'witness_Q_calls':0}
        table=self.shot_table();last=self.maximum_odd_multiplier
        prefix=[None]*(last+1);parents={};costs={}
        best=None;best_m=None
        for m in range(1,last+1):
            if m%2 and m in table:
                n=table[m];local=n*(a_cost+q_cost*((m-1)//2))
                if m==self.initial_multiplier:
                    current=local;parent=None
                elif m>=3*self.initial_multiplier and prefix[m//3] is not None:
                    prior_cost,prior_m=prefix[m//3];current=prior_cost+local;parent=prior_m
                else:current=None
                if current is not None:
                    costs[m]=current;parents[m]=parent
                    if best is None or current>best:best=current;best_m=m
            prefix[m]=None if best is None else (best,best_m)
        chain=[];cursor=best_m
        while cursor is not None:chain.append(cursor);cursor=parents[cursor]
        chain.reverse()
        return {'cost_upper':best,'witness_multipliers':chain,
                'witness_shots':sum(table[m] for m in chain),
                'witness_Q_calls':sum(table[m]*((m-1)//2) for m in chain)}


@lru_cache(maxsize=500)
def plan(epsilon,alpha=F(1,100),prior_upper=F(1)):
    epsilon,alpha,prior=F(epsilon),F(alpha),F(prior_upper)
    if not F(1,1000000)<=epsilon<F(1,2):raise ValueError('Declared epsilon range is1e-6<=epsilon<.5')
    if not 0<alpha<1 or not 0<=prior<=1:raise ValueError('Invalid alpha or proved prior')
    cover_certificate()
    guard=min(F(1),prior+F(1,10**18))
    variance=guard*(1-guard) if guard<=F(1,2) else F(1,4)
    gamma=2*sqrt_upper_fraction(variance)
    if prior<=2*epsilon:
        return Plan(epsilon,alpha,prior,guard,gamma,F(0),F(0),F(0),0,1,0,0,0,())
    theta_epsilon=(epsilon-ENDPOINT_ALLOWANCE)/gamma
    pi_upper=pi_rational_interval()[1]
    target=4*theta_epsilon/pi_upper;bound=1/target
    upper=floor_fraction(bound)
    if upper%2==0:upper-=1
    if upper>1000000:raise ValueError('Declared finite scheduling table limit exceeded')
    hi=inverse_probability_interval(prior)[1]
    initial=floor_fraction(1/hi)
    if initial%2==0:initial-=1
    initial=max(1,initial)
    assert initial<=upper
    rounds=1;m=initial
    while m*3<=bound:m*=3;rounds+=1
    table=[]
    for m in range(initial,upper+1,2):
        if m!=initial and m<3*initial:continue
        failure=(2*alpha/3)*m/bound
        shots=ceil_fraction(log_upper(2/failure)/(2*RADIUS*RADIUS))
        table.append((m,shots))
    interim=Plan(epsilon,alpha,prior,guard,gamma,theta_epsilon,target,bound,upper,initial,rounds,0,0,tuple(table))
    shot_cap=interim.cost_upper(1,0)['cost_upper']
    query_cap=interim.cost_upper(0,1)['cost_upper']
    return Plan(epsilon,alpha,prior,guard,gamma,theta_epsilon,target,bound,upper,initial,rounds,shot_cap,query_cap,tuple(table))


def find_next_multiplier(current,lo,hi):
    maximum=floor_fraction(1/(hi-lo))
    if maximum%2==0:maximum-=1
    for candidate in range(maximum,3*current-1,-2):
        if quadrant(candidate,lo,hi) is not None:return candidate
    raise AssertionError('The certified three/five/seven cover failed')


def estimate(count_ones,epsilon,alpha=F(1,100),prior_upper=F(1)):
    setup=plan(F(epsilon),F(alpha),F(prior_upper));prior=setup.prior_upper
    lo=F(0);hi=inverse_probability_interval(prior)[1];m=setup.initial_multiplier;history=[]
    def result(status,interval):
        a,b=interval
        return {'status':status,'amplitude_interval':interval,'estimate':(a+b)/2,
                'plan':setup,'rounds':history,'total_shots':sum(x['shots'] for x in history),
                'total_Q_calls':sum(x['shots']*x['k'] for x in history)}
    def current_interval():
        a,b=amplitude_interval(lo,hi);return a,min(b,prior)
    initial=current_interval()
    if initial[1]-initial[0]<=2*setup.epsilon:return result('certified',initial)
    table=setup.shot_table()
    for index in range(setup.maximum_rounds):
        j=quadrant(m,lo,hi)
        if j is None:raise AssertionError('Selected power crossed a quadrant')
        n=table[m];k=(m-1)//2;count=count_ones(k,n)
        if isinstance(count,bool) or count!=int(count):raise ValueError('Backend count must be an integer')
        count=int(count)
        if not 0<=count<=n:raise ValueError('Backend count outside shot range')
        observed=F(count,n);a=max(F(0),observed-RADIUS);b=min(F(1),observed+RADIUS)
        ql,qu=(a,b) if j%2==0 else (1-b,1-a)
        left=inverse_probability_interval(ql)[0];right=inverse_probability_interval(qu)[1]
        lo=max(lo,(j+left)/m);hi=min(hi,(j+right)/m)
        history.append({'m':m,'k':k,'shots':n,'ones':count,'quadrant':j,
                        'conditional_failure_upper':(2*setup.alpha/3)*m/setup.multiplier_upper,
                        'normalised_angle_interval':(lo,hi)})
        if lo>hi:return result('confidence_failure_empty_intersection',(F(0),F(0)))
        amplitude=current_interval()
        if amplitude[0]>amplitude[1]:return result('confidence_failure_prior_intersection',(F(0),F(0)))
        if amplitude[1]-amplitude[0]<=2*setup.epsilon:return result('certified',amplitude)
        if hi-lo<=setup.normalised_angle_width_target:
            raise AssertionError('The proved monetary stopping implication failed')
        m=find_next_multiplier(m,lo,hi)
        if m>setup.maximum_odd_multiplier:raise AssertionError('The proved multiplier bound failed')
    raise AssertionError('The proved maximum round count failed')


def rational_record(value):
    if isinstance(value,F):return {'numerator':str(value.numerator),'denominator':str(value.denominator),'float':float(value)}
    if isinstance(value,Plan):
        out={k:rational_record(v) for k,v in value.__dict__.items() if k!='shots_by_multiplier'}
        out['shot_table_entries']=len(value.shots_by_multiplier)
        out['shot_table_sha256']=sha256(repr(value.shots_by_multiplier).encode()).hexdigest()
        return out
    if isinstance(value,dict):return {k:rational_record(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [rational_record(v) for v in value]
    return value
