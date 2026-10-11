"""Certified dyadic finite-value quantum mean planning.

The finite output's second moment supplies bucket *mean* priors. It is never
assigned to the variance of the original Bernoulli comparison flag. The
published dyadic reduction is implemented with the project's separately
proved scalar controllers. Resource caps use supplied emitted A/Q profiles.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import math
import full_batch_miqae
import iterative_estimator
from certified_scalar import sqrt_upper_fraction

MIN_EPSILON = F(1, 10**6)


@dataclass(frozen=True)
class Bucket:
    index: int
    weight: F
    prior: F
    top_atom: bool


def sharp_grid_mean_upper(q, second_moment_upper, lower_integer=1, upper_integer=None):
    """Sharp mean bound on {0,l,l+2**-q,...,h} from a second moment.

    This is the inverse lower convex hull of the exact finite squared-value
    graph. It is attainable by a zero/l or adjacent-grid-point mixture. The
    support class is all allowed register integers, not a claim that every
    such value is attained by the particular financial circuit.
    """
    if isinstance(q,bool) or not isinstance(q,int) or not 1<=q<=64:
        raise ValueError('The declared register range is 1 <= q <= 64')
    v=F(second_moment_upper);upper_integer=(1 << q) if upper_integer is None else upper_integer
    if any(isinstance(a,bool) or not isinstance(a,int) for a in (lower_integer,upper_integer)):
        raise ValueError('Support endpoints must be register integers')
    if not 0<=v<=1 or not 1<=lower_integer<=upper_integer<=1 << q:
        raise ValueError('Invalid finite support or second moment')
    delta=F(1,1 << q);lo=lower_integer*delta;hi=upper_integer*delta
    if v<lo*lo:return v/lo
    if v>=hi*hi:return hi
    scaled=v*(1 << (2*q))
    j=math.isqrt(scaled.numerator//scaled.denominator)
    assert lower_integer<=j<upper_integer
    x,y=j*delta,(j+1)*delta
    return (v+x*y)/(x+y)


def validate(q, mean_upper, second_moment_upper, tolerance, alpha):
    if isinstance(q, bool) or not isinstance(q, int) or not 1 <= q <= 64:
        raise ValueError('The declared register range is 1 <= q <= 64')
    u, v, t, a = map(F, (mean_upper, second_moment_upper, tolerance, alpha))
    if not 0 <= u <= 1 or not 0 <= v <= 1:
        raise ValueError('Moment upper bounds must lie in [0,1]')
    if not 0 < t < F(1,2) or not 0 < a < 1:
        raise ValueError('Invalid tolerance or confidence')
    return u, v, t, a


def partition(q, omitted_prefix_bits, mean_upper, second_moment_upper, sharp_priors=False,
              max_integer=None):
    """Exact weights, bucket priors and omitted-prefix mean enclosure.

    Omitted values satisfy 0 <= a < 2**s. Their represented midpoint has
    worst-case mean error tail_upper/2. Retained buckets are s,...,q-1,
    followed by the separate atom a=2**q.
    """
    s = omitted_prefix_bits
    if isinstance(q,bool) or not isinstance(q,int) or not 1<=q<=64:
        raise ValueError('The declared register range is 1 <= q <= 64')
    maximum=(1 << q) if max_integer is None else max_integer
    if isinstance(maximum,bool) or not isinstance(maximum,int) or not 0<=maximum<=1 << q:
        raise ValueError('The support maximum must be a register integer in [0,2**q]')
    if not isinstance(s, int) or isinstance(s, bool) or not 0 <= s <= q:
        raise ValueError('The omitted prefix must lie between zero and q')
    u, v = F(mean_upper), F(second_moment_upper)
    if not 0 <= u <= 1 or not 0 <= v <= 1:
        raise ValueError('Invalid proved moments')
    tail = min(F(min((1 << s)-1,maximum), 1 << q), u, sqrt_upper_fraction(v))
    buckets = []
    for b in range(s, q):
        weight = F(1 << (b+1), 1 << q)
        # Every value in the bucket is at most 1-2**(-b-1) after scaling.
        highest=min((1 << (b+1))-1,maximum)
        prior=(F(0) if highest<1 << b else
               min(F(highest,1 << (b+1)),2*v/(weight*weight),u/weight))
        if sharp_priors and prior>0:
            prior=min(prior,sharp_grid_mean_upper(q,v,1 << b,highest)/weight)
        buckets.append(Bucket(b, weight, prior, False))
    buckets.append(Bucket(q,F(1),min(v,u) if maximum==1 << q else F(0),True))
    return tuple(buckets), tail


def default_prefix(q, mean_upper, second_moment_upper, tolerance, tail_fraction=F(1,4),
                   max_integer=None):
    """Largest omitted prefix whose midpoint error fits the tail allowance."""
    t, fraction = F(tolerance), F(tail_fraction)
    if not 0 <= fraction < 1:
        raise ValueError('The tail budget fraction must lie in [0,1)')
    eligible = [s for s in range(q+1)
                if partition(q,s,mean_upper,second_moment_upper,max_integer=max_integer)[1]/2 <= fraction*t]
    return max(eligible)


def capped_allocation(buckets, budget, weights, enforce_probability_floor=True):
    """Exact finite water filling with each zero-shot midpoint as a cap.

    Each active probability tolerance is at least MIN_EPSILON. Those lower
    allowances are charged before distributing the remaining budget. This
    respects the declared scalar implementation range without a hidden
    extrapolation. No global gate-cost optimality is asserted.
    """
    budget = F(budget)
    if len(buckets) != len(weights) or any(F(w) <= 0 for w in weights):
        raise ValueError('Each bucket requires a positive rational weight')
    caps = [b.weight*b.prior/2 for b in buckets]
    floors = ([min(cap,b.weight*MIN_EPSILON) for b,cap in zip(buckets,caps)]
              if enforce_probability_floor else [F(0) for cap in caps])
    if sum(floors,F(0)) > budget:
        raise ValueError('The bucket tolerance falls outside the scalar implementation range')
    if sum(caps,F(0)) <= budget:
        return tuple(caps)
    values = list(floors)
    remaining = budget-sum(values,F(0))
    free = {i for i in range(len(values)) if values[i] < caps[i]}
    while free:
        denominator = sum((F(weights[i]) for i in free),F(0))
        increments = {i:remaining*F(weights[i])/denominator for i in free}
        saturated = {i for i in free if increments[i] >= caps[i]-values[i]}
        if not saturated:
            for i in free: values[i] += increments[i]
            remaining = F(0)
            break
        for i in saturated:
            remaining -= caps[i]-values[i]
            values[i] = caps[i]
        free -= saturated
    assert sum(values,F(0)) <= budget
    assert all(floor <= value <= cap for floor,value,cap in zip(floors,values,caps))
    if any(value < cap and value/b.weight < MIN_EPSILON
           for b,value,cap in zip(buckets,values,caps)):
        raise ValueError('The bucket tolerance falls outside the scalar implementation range')
    return tuple(values)


def allocation_weights(buckets, resources, style='equal', metric='logical_T_count'):
    """Positive exact weights for a reproducible finite candidate search.

    The cost/moment choice follows the continuous surrogate sum d_b/e_b:
    d_b=(c_A+c_Q)*r_b*sqrt(prior_b). Its minimiser without caps is
    e_b proportional to sqrt(d_b). Rounded rational weights simply specify
    a candidate allocation; the actual discrete scalar caps choose winners.
    """
    if style == 'equal': return tuple(F(1) for b in buckets)
    if style == 'dyadic': return tuple(sqrt_upper_fraction(b.weight) for b in buckets)
    if style != 'cost_moment': raise ValueError('Unknown allocation style')
    result=[]
    for b in buckets:
        if b.prior==0:
            result.append(F(1,1 << 100));continue
        profile=resources[b.index]
        c=max(1,F(profile['A'].get(metric,0))+F(profile['Q'].get(metric,0)))
        d=c*b.weight*sqrt_upper_fraction(b.prior)
        # A zero prior is already a zero-shot atom; retain a positive weight.
        result.append(max(F(1,1 << 100),sqrt_upper_fraction(d)))
    return tuple(result)


def _controller(name):
    if name == 'miqae': return full_batch_miqae
    if name == 'three_five': return iterative_estimator
    raise ValueError('Unknown scalar controller')


def scalar_cost(setup, controller, a_cost, q_cost):
    if controller == 'miqae':
        return setup.cost_upper(a_cost,q_cost)['cost_upper']
    return setup.total_shots_upper*a_cost+setup.total_grover_calls_upper*q_cost


def build_plan(q, mean_upper, second_moment_upper, tolerance, alpha,
               resources, controller='miqae', allocation='equal',
               omitted_prefix_bits=None, tail_fraction=F(1,4),
               metric='logical_T_count', sharp_priors=False,max_integer=None):
    """Plan independent bucket means with one rigorous aggregate guarantee.

    resources[b] contains {'A':{metric:cost}, 'Q':{metric:cost}, 'qubits':n}.
    All per-operation metric costs must be nonnegative exact numbers.
    A/Q denote emitted circuits for this bucket, including numerical-value
    computation and uncomputation, fresh comparison bits, bucket selection,
    reflection work and any measurement/feedforward instrument resources.
    """
    u,v,t,alpha=validate(q,mean_upper,second_moment_upper,tolerance,alpha)
    maximum=(1 << q) if max_integer is None else max_integer
    if isinstance(maximum,bool) or not isinstance(maximum,int) or not 0<=maximum<=1 << q:
        raise ValueError('The support maximum must be a register integer in [0,2**q]')
    # The entire finite expectation has the stronger joint prior. This is
    # strictly separate from adding the potentially loose bucket priors.
    whole_upper=min(u,sqrt_upper_fraction(v),F(maximum,1 << q))
    if sharp_priors and maximum:whole_upper=min(whole_upper,sharp_grid_mean_upper(q,v,upper_integer=maximum))
    if whole_upper <= 2*t:
        return {'q':q,'controller':controller,'allocation':allocation,
                'status':'whole_prior_zero_shot','mean_upper':u,'max_integer':maximum,
                'second_moment_upper':v,'tolerance':t,'alpha':alpha,
                'estimate_without_shots':whole_upper/2,
                'total_error_upper':whole_upper/2,'buckets':(),
                'aggregate_metrics_upper':{},'total_shots_upper':0,
                'total_Q_calls_upper':0,'peak_qubits_upper':0}
    u=whole_upper
    s=default_prefix(q,u,v,t,tail_fraction,maximum) if omitted_prefix_bits is None else omitted_prefix_bits
    buckets,tail=partition(q,s,u,v,sharp_priors,maximum)
    tail_error=tail/2
    if tail_error >= t: raise ValueError('Omitted-prefix error exhausts the target')
    weights=allocation_weights(buckets,resources,allocation,metric)
    errors=capped_allocation(buckets,t-tail_error,weights,
                             enforce_probability_floor=(allocation!='equal'))
    active=[i for i,(b,e) in enumerate(zip(buckets,errors)) if e < b.weight*b.prior/2]
    per_alpha=alpha/len(active) if active else F(0)
    module=_controller(controller)
    records=[];aggregate={};shots=queries=peak=0
    for i,(b,e) in enumerate(zip(buckets,errors)):
        record={'index':b.index,'weight':b.weight,'prior_upper':b.prior,
                'top_atom':b.top_atom,'mean_error_allowance':e,
                'probability_tolerance':e/b.weight,'alpha':per_alpha if i in active else F(0)}
        if i not in active:
            record.update(status='prior_zero_shot',estimate_without_shots=b.prior/2,
                          maximum_shots=0,maximum_Q_calls=0,metric_caps={})
        else:
            profile=resources[b.index]
            setup=module.plan(e/b.weight,per_alpha,b.prior)
            metrics=set(profile['A'])|set(profile['Q'])
            caps={}
            for name in metrics:
                ca,cq=profile['A'].get(name,0),profile['Q'].get(name,0)
                if F(ca)<0 or F(cq)<0: raise ValueError('Negative emitted resource cost')
                caps[name]=scalar_cost(setup,controller,ca,cq)
                aggregate[name]=aggregate.get(name,0)+caps[name]
            n=setup.total_shots_upper
            k=setup.total_grover_calls_upper
            shots+=n;queries+=k
            peak=max(peak,int(profile.get('qubits',0)))
            record.update(status='estimated',scalar_plan=setup,maximum_shots=n,
                          maximum_Q_calls=k,metric_caps=caps)
        records.append(record)
    total_error=tail_error+sum((e for e in errors),F(0))
    assert total_error<=t
    assert sum((r['alpha'] for r in records),F(0))<=alpha
    return {'q':q,'controller':controller,'allocation':allocation,
            'status':'bucket_estimation','mean_upper':u,'second_moment_upper':v,
            'max_integer':maximum,
            'tolerance':t,'alpha':alpha,'omitted_prefix_bits':s,
            'omitted_mean_interval':(F(0),tail),'omitted_midpoint':tail/2,
            'omitted_error_upper':tail_error,'total_error_upper':total_error,
            'buckets':tuple(records),'aggregate_metrics_upper':aggregate,
            'total_shots_upper':shots,'total_Q_calls_upper':queries,
            'peak_qubits_upper':peak,'estimate_clipping_interval':(F(0),whole_upper)}


def choose_plan(q, mean_upper, second_moment_upper, tolerance, alpha, resources,
                controller='miqae', metric='logical_T_count',
                allocations=('equal','dyadic','cost_moment'), prefixes=None,
                max_integer=None,sharp_priors=False):
    """Return the cheapest completely evaluated finite candidate and all peers."""
    if prefixes is None:
        default=default_prefix(q,mean_upper,second_moment_upper,tolerance,max_integer=max_integer)
        prefixes=tuple(sorted({0,default,max(0,default-1)}))
    candidates=[];unsupported=[]
    for s in prefixes:
        for allocation in allocations:
            try:
                candidate=build_plan(q,mean_upper,second_moment_upper,tolerance,alpha,
                    resources,controller,allocation,s,metric=metric,
                    max_integer=max_integer,sharp_priors=sharp_priors)
            except ValueError as err:
                unsupported.append({'prefix':s,'allocation':allocation,'reason':str(err)})
                continue
            candidates.append(candidate)
    if not candidates: raise ValueError('No declared-range candidate is admissible')
    winner=min(candidates,key=lambda p:(p['aggregate_metrics_upper'].get(metric,0),
                p['total_Q_calls_upper'],p['total_shots_upper'],
                p.get('omitted_prefix_bits',0),p['allocation']))
    return {'winner':winner,'candidates':tuple(candidates),'unsupported':tuple(unsupported),
            'selection_metric':metric,'selection_scope':'minimum over the retained finite candidate set only'}


def plan(q,v2,u,e,alpha,bucketcosts,controller='miqae',**kwargs):
    """Requested moment-first interface; return one fully certified allocation."""
    aliases={'three_five':'three_five','3/5':'three_five','miqae':'miqae','full_batch_miqae':'miqae'}
    if controller not in aliases: raise ValueError('Unknown scalar controller')
    return build_plan(q,u,v2,e,alpha,bucketcosts,aliases[controller],**kwargs)


def execute(setup, count_ones, resources):
    """Run the certified bucket controller against an ideal backend contract.

    count_ones(index,k,n) must return n fresh independent measurements of
    Q_index**k A_index|0>. An exact emulator is useful for checks, but it is
    not a physical quantum experiment or independent market observation.
    """
    if setup['status']=='whole_prior_zero_shot':
        upper=2*setup['estimate_without_shots']
        return {'status':'certified','mean_interval':(F(0),upper),
                'estimate':upper/2,'bucket_runs':(), 'total_shots':0,
                'total_Q_calls':0,'actual_metrics':{}}
    module=_controller(setup['controller'])
    lo,hi=setup['omitted_mean_interval']
    runs=[];shots=queries=0;metrics={}
    for record in setup['buckets']:
        index=record['index'];weight=record['weight']
        if record['status']=='prior_zero_shot':
            a,b=F(0),record['prior_upper']
            run={'index':index,'status':'prior_zero_shot','amplitude_interval':(a,b)}
        else:
            run=module.estimate(lambda k,n:count_ones(index,k,n),
                    record['probability_tolerance'],record['alpha'],record['prior_upper'])
            run['index']=index
            if run['status']!='certified':
                return {'status':'confidence_failure_bucket','failed_bucket':index,
                        'mean_interval':None,'estimate':None,'bucket_runs':tuple(runs+[run])}
            a,b=run['amplitude_interval'];n=run['total_shots']
            k=run.get('total_Q_calls',run.get('total_grover_calls'))
            shots+=n;queries+=k
            profile=resources[index]
            for name in set(profile['A'])|set(profile['Q']):
                metrics[name]=metrics.get(name,0)+n*profile['A'].get(name,0)+k*profile['Q'].get(name,0)
        lo+=weight*a;hi+=weight*b;runs.append(run)
    global_lo,global_hi=setup['estimate_clipping_interval']
    lo=max(lo,global_lo);hi=min(hi,global_hi)
    if lo>hi:
        return {'status':'confidence_failure_prior_intersection','mean_interval':None,
                'estimate':None,'bucket_runs':tuple(runs)}
    assert hi-lo<=2*setup['tolerance']
    assert shots<=setup['total_shots_upper'] and queries<=setup['total_Q_calls_upper']
    assert all(value<=setup['aggregate_metrics_upper'][name] for name,value in metrics.items())
    return {'status':'certified','mean_interval':(lo,hi),'estimate':(lo+hi)/2,
            'bucket_runs':tuple(runs),'total_shots':shots,'total_Q_calls':queries,
            'actual_metrics':metrics}


def exact_bucket_means(integer_values, q):
    """Independent-law diagnostic, never used to supply a production prior."""
    values=list(integer_values)
    if not values or any(isinstance(a,bool) or not isinstance(a,int) or not 0<=a<=1 << q for a in values):
        raise ValueError('Invalid finite-law integer outputs')
    means=[F(0) for _ in range(q+1)]
    for a in values:
        if a==1 << q:means[q]+=F(1,len(values))
        elif a:
            b=a.bit_length()-1
            means[b]+=F(a,(1 << (b+1))*len(values))
    return tuple(means)


def rational_record(value):
    if isinstance(value,F):
        return {'numerator':str(value.numerator),'denominator':str(value.denominator),'float':float(value)}
    if isinstance(value,full_batch_miqae.Plan): return full_batch_miqae.rational_record(value)
    if isinstance(value,iterative_estimator.Plan): return iterative_estimator.decimal_record(value)
    if isinstance(value,dict): return {str(k):rational_record(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [rational_record(v) for v in value]
    return value
