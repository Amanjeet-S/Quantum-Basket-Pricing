"""A conservative, confidence-certified 3/5 amplitude-estimation controller.

The backend contract is count_ones(k, shots) for fresh independent measurements
of Q**k A|0>. No quantum device is invoked by the included simulation adapter.
Classical branch decisions use rational intervals, not floating-point phases.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
import math

F = Fraction
RADIUS = F(1, 25)
INVERSE_TOLERANCE = F(1, 1 << 70)


def ceil_fraction(value):
    return -((-value.numerator)//value.denominator)


def floor_fraction(value):
    return value.numerator//value.denominator


@lru_cache(None)
def pi_rational_interval(terms=80):
    def atan_interval(den):
        x=F(1,den);total=F(0)
        for j in range(terms):
            total+=(-1 if j%2 else 1)*x**(2*j+1)/(2*j+1)
        remainder=x**(2*terms+1)/(2*terms+1)
        return total-remainder,total+remainder
    al,au=atan_interval(5);bl,bu=atan_interval(239)
    return 16*al-4*bu,16*au-4*bl


def decimal_fraction(value, rounding):
    with localcontext() as ctx:
        ctx.prec=80;ctx.rounding=rounding
        return Decimal(value.numerator)/Decimal(value.denominator)


@lru_cache(maxsize=20000)
def probability_interval(u):
    """Enclose f(u)=sin²(pi*u/2) for rational u in [0,1].

    Directed Decimal interval arithmetic encloses the cosine Taylor polynomial.
    Its decreasing alternating tail is enclosed by the next positive term.
    Pi itself has an exact rational Machin-series enclosure.
    """
    u=F(u)
    if not 0<=u<=1:raise ValueError('Probability input must lie in [0,1]')
    if u==0:return F(0),F(0)
    if u==1:return F(1),F(1)
    if u==F(1,2):return F(1,2),F(1,2)
    pl,pu=pi_rational_interval()
    with localcontext() as ctx:
        ctx.prec=80;ctx.rounding=ROUND_FLOOR
        xl=decimal_fraction(pl,ROUND_FLOOR)*decimal_fraction(u,ROUND_FLOOR)
        x2l=xl*xl
        ctx.rounding=ROUND_CEILING
        xu=decimal_fraction(pu,ROUND_CEILING)*decimal_fraction(u,ROUND_CEILING)
        x2u=xu*xu
        tl=tu=Decimal(1);sl=su=Decimal(1)
        for j in range(1,45):
            den=Decimal((2*j-1)*(2*j))
            ctx.rounding=ROUND_FLOOR;tl=tl*x2l/den
            ctx.rounding=ROUND_CEILING;tu=tu*x2u/den
            if j%2:
                ctx.rounding=ROUND_FLOOR;sl=sl-tu
                ctx.rounding=ROUND_CEILING;su=su-tl
            else:
                ctx.rounding=ROUND_FLOOR;sl=sl+tl
                ctx.rounding=ROUND_CEILING;su=su+tu
        # At j=45 the remaining terms decrease because |x|<pi<4.
        ctx.rounding=ROUND_CEILING
        tail=tu*x2u/Decimal(89*90)
        ctx.rounding=ROUND_FLOOR;cl=sl-tail
        ctx.rounding=ROUND_CEILING;cu=su+tail
        ctx.rounding=ROUND_FLOOR;lo=(Decimal(1)-cu)/2
        ctx.rounding=ROUND_CEILING;hi=(Decimal(1)-cl)/2
    return max(F(0),F(lo)),min(F(1),F(hi))


@lru_cache(maxsize=20000)
def inverse_probability_interval(p):
    """Enclose the inverse of increasing sin²(pi*u/2) on [0,1]."""
    p=F(p)
    if not 0<=p<=1:raise ValueError('Probability must lie in [0,1]')
    if p in (F(0),F(1),F(1,2)):return p,p
    lo,hi=F(0),F(1);eta=INVERSE_TOLERANCE
    while hi-lo>eta:
        mid=(lo+hi)/2;fl,fu=probability_interval(mid)
        if fu<p:lo=mid
        elif fl>p:hi=mid
        else:
            # |f(x)-f(y)| >= sin²(pi*|x-y|/2) >= |x-y|².
            # Therefore a probability enclosure of width <=eta² implies
            # inverse coordinate error <=eta, without dividing by a slope.
            if fu-fl>eta*eta:raise ArithmeticError('Cosine enclosure too wide')
            return max(lo,mid-eta),min(hi,mid+eta)
    return lo,hi


def quadrant(m, lo, hi):
    """Return j if [m*lo,m*hi] lies in one closed [j,j+1] cell."""
    a,b=m*lo,m*hi;j=floor_fraction(a)
    if b<=j+1:return j
    return None


def amplitude_interval(lo,hi):
    return probability_interval(lo)[0],probability_interval(hi)[1]


@dataclass(frozen=True)
class Plan:
    epsilon: Fraction
    alpha: Fraction
    maximum_rounds: int
    shots_per_round: int
    target_multiplier: int
    multiplier_upper: int
    grover_power_upper: int
    total_shots_upper: int
    total_grover_calls_upper: int
    prior_upper: Fraction
    guarded_prior_upper: Fraction
    amplitude_angle_lipschitz_upper: Fraction
    initial_multiplier: int
    inverse_coordinate_tolerance: Fraction=INVERSE_TOLERANCE


def sqrt_upper_fraction(value,bits=100):
    """An exact rational upper square-root, without a libm assumption."""
    value=F(value)
    if value<0:raise ValueError('Negative square-root input')
    scale=1<<bits;integer=ceil_fraction(value*scale*scale)
    root=math.isqrt(integer)
    if root*root<integer:root+=1
    return F(root,scale)


def plan(epsilon,alpha=F(1,100),prior_upper=F(1)):
    epsilon,alpha=F(epsilon),F(alpha)
    if not F(1,1000000)<=epsilon<F(1,2):
        raise ValueError('Declared implementation range: 1e-6 <= epsilon < .5')
    if not 0<alpha<1:raise ValueError('Failure probability must lie in (0,1)')
    prior=F(prior_upper)
    if not 0<=prior<=1:raise ValueError('Prior must lie in [0,1]')
    guard=min(F(1),prior+F(1,10**18))
    variance=guard*(1-guard) if guard<=F(1,2) else F(1,4)
    gamma=2*sqrt_upper_fraction(variance)
    if prior<=2*epsilon:
        return Plan(epsilon,alpha,0,0,0,0,0,0,0,prior,guard,gamma,1)
    upper_t=inverse_probability_interval(prior)[1]
    initial=floor_fraction(1/upper_t)
    if initial%2==0:initial-=1
    initial=max(1,initial)
    target=max(1,ceil_fraction(gamma/(5*epsilon)))
    rounds=1;power=initial
    while power<target:rounds+=1;power*=3
    with localcontext() as ctx:
        ctx.prec=80
        arg=Decimal(2*rounds*alpha.denominator)/Decimal(alpha.numerator)
        # Decimal ln is correctly rounded, not directed. Add an explicit
        # envelope before the exact integer ceiling; the argument here is
        # finite and its logarithm is below 100 for declared confidence runs.
        logarithm=arg.ln()+Decimal('1e-60')
        if not 0<logarithm<100:raise ValueError('Declared logarithm envelope requires 0<log(2R/alpha)<100')
        shots=ceil_fraction(F(logarithm)/(2*RADIUS*RADIUS))
    maximum_m=initial if initial>=target else min(5*target,initial*5**(rounds-1))
    maximum_queries=((initial-1)//2 if rounds==1 else min((3*maximum_m)//4,
        (initial*((5**rounds-1)//4)-rounds)//2))
    return Plan(epsilon,alpha,rounds,shots,target,maximum_m,(maximum_m-1)//2,
                rounds*shots,shots*maximum_queries,prior,guard,gamma,initial)


def estimate(count_ones,epsilon,alpha=F(1,100),prior_upper=F(1)):
    """Return a rigorous finite-amplitude confidence interval and query record.

    Coverage is at least 1-alpha under ideal independent quantum measurements.
    Empty intersections are a flagged failure outcome, whose probability is
    already bounded by alpha, rather than an unjustified confidence repair.
    """
    prior=F(prior_upper);setup=plan(epsilon,alpha,prior)
    if not 0<=prior<=1:raise ValueError('Prior must be a proved probability bound')
    lo=F(0);hi=inverse_probability_interval(prior)[1];m=setup.initial_multiplier;history=[]
    def result(status, interval):
        a,b=interval
        return {'status':status,'amplitude_interval':(a,b),'estimate':(a+b)/2,
                'plan':setup,'rounds':history,'total_shots':sum(x['shots'] for x in history),
                'total_grover_calls':sum(x['shots']*x['k'] for x in history)}
    initial=amplitude_interval(lo,hi)
    initial=(initial[0],min(initial[1],prior))
    if initial[1]-initial[0]<=2*setup.epsilon:return result('certified',initial)
    for round_id in range(setup.maximum_rounds):
        j=quadrant(m,lo,hi)
        if j is None:raise AssertionError('Amplification interval crossed a quadrant')
        k=(m-1)//2;n=setup.shots_per_round
        supplied=count_ones(k,n)
        if isinstance(supplied,bool) or supplied!=int(supplied):raise ValueError('Backend count must be an integer')
        successes=int(supplied)
        if not 0<=successes<=n:raise ValueError('Invalid backend count')
        phat=F(successes,n);pl=max(F(0),phat-RADIUS);pu=min(F(1),phat+RADIUS)
        ql,qu=(pl,pu) if j%2==0 else (1-pu,1-pl)
        left=inverse_probability_interval(ql)[0]
        right=inverse_probability_interval(qu)[1]
        newlo=max(lo,(j+left)/m);newhi=min(hi,(j+right)/m)
        history.append({'round':round_id+1,'m':m,'k':k,'shots':n,'ones':successes,
                        'quadrant':j,'bernoulli_interval':(pl,pu),
                        'theta_normalised_interval':(newlo,newhi)})
        if newlo>newhi:return result('confidence_failure_empty_intersection',(F(0),F(0)))
        lo,hi=newlo,newhi
        ai=amplitude_interval(lo,hi);ai=(ai[0],min(ai[1],prior))
        if ai[0]>ai[1]:return result('confidence_failure_prior_intersection',(F(0),F(0)))
        if ai[1]-ai[0]<=2*setup.epsilon:return result('certified',ai)
        factors=[s for s in (3,5) if quadrant(s*m,lo,hi) is not None]
        if not factors:raise AssertionError('The proved 3/5 quadrant cover failed')
        m*=factors[0]
    raise AssertionError('The proved deterministic round bound failed')


def decimal_record(value):
    """Serialisation only; control decisions never use these floats."""
    if isinstance(value,Fraction):return {'numerator':value.numerator,'denominator':value.denominator,'float':float(value)}
    if isinstance(value,Plan):return {k:decimal_record(v) for k,v in value.__dict__.items()}
    if isinstance(value,dict):return {k:decimal_record(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [decimal_record(v) for v in value]
    return value


def original_constructor_resources(row,alpha=F(1,100),prior_upper=F(1)):
    """Complete worst-case costs using the original emitted seven-T profile.

    An improved architecture can substitute its checked A and Q profiles;
    the shot/power controller and proof do not assume their particular cost.
    """
    cap=int(row['cap']);setup=plan(F(1,200*cap),alpha,prior_upper)
    rr=row['resources'];shots=setup.total_shots_upper;qs=setup.total_grover_calls_upper
    aT=rr['logical_T_count_exact_7T_Toffoli'];qT=rr['grover_iterate_T_count_exact_7T_Toffoli']
    ac=rr['expanded_gate_counts'];rc=rr['grover_zero_reflection']['gate_counts'];calls=shots+2*qs
    logical={'X':calls*ac['X']+qs*rc.get('X',0),
             'CNOT':calls*ac['CNOT'],
             'CCX':calls*ac['Toffoli']+qs*rc.get('CCX',0),
             'CZ':qs*rc.get('CZ',0),'Z':qs,
             'H':calls*rr['includes_finite_uniform_preparation_H_gates']}
    # A standard exact seven-T CCX has six CNOT and two H gates. The assigned
    # three-T-layer schedule comes from the constructor's declared profile;
    # the counts here do not claim an optimal Clifford/depth compilation.
    expanded=dict(logical);ccx=expanded.pop('CCX')
    expanded['CNOT']+=6*ccx;expanded['H']+=2*ccx
    expanded['T']=4*ccx;expanded['Tdg']=3*ccx
    reflection_gate_upper=sum(rc.values())
    a_depth=rr['state_preparation_depth_including_uniform_H_layer_upper']
    q_depth=2*a_depth+reflection_gate_upper+1
    return {'plan':setup,'confidence_scope':'one price at 1-alpha; simultaneous 27 prices requires alpha_total/27',
            'statistical_currency_tolerance':F(1,200),'finite_amplitude_tolerance':setup.epsilon,
            'clean_A_invocations_upper':shots+2*qs,'grover_invocations_upper':qs,
            'flag_measurements_upper':shots,'fresh_shot_resets_upper':shots,
            'T_count_upper':shots*aT+qs*qT,
            'logical_gate_counts_upper':logical,'Clifford_T_gate_counts_upper':expanded,
            'serial_all_logical_gate_depth_upper':shots*a_depth+qs*q_depth,
            'maximum_shot_all_logical_gate_depth_upper':a_depth+setup.grover_power_upper*q_depth,
            'logical_qubits_upper':rr['grover_iterate_logical_qubits_upper'],
            'one_shot_T_depth_upper':rr['logical_T_depth_serial_3_layers_per_Toffoli']+setup.grover_power_upper*rr['grover_iterate_T_depth_serial_upper'],
            'serial_T_depth_all_shots_upper':shots*rr['logical_T_depth_serial_3_layers_per_Toffoli']+qs*rr['grover_iterate_T_depth_serial_upper'],
            'explicit_7T_sequence_T_depth_serial_upper':shots*rr['logical_T_depth_serial_3_layers_per_Toffoli']+qs*rr['grover_iterate_T_depth_serial_upper']+2*logical['CCX'],
            'explicit_7T_sequence_T_depth_maximum_shot_upper':rr['logical_T_depth_serial_3_layers_per_Toffoli']+setup.grover_power_upper*rr['grover_iterate_T_depth_serial_upper']+2*(ac['Toffoli']+setup.grover_power_upper*(2*ac['Toffoli']+rc.get('CCX',0))),
            'explicit_7T_sequence_all_gate_dependency_depth_serial_upper':shots*a_depth+qs*q_depth+14*logical['CCX'],
            'explicit_7T_sequence_all_gate_dependency_depth_maximum_shot_upper':a_depth+setup.grover_power_upper*q_depth+14*(ac['Toffoli']+setup.grover_power_upper*(2*ac['Toffoli']+rc.get('CCX',0))),
            'explicit_7T_sequence_all_gate_depth_serial_upper':shots*rr['includes_finite_uniform_preparation_H_gates']+logical['X']+logical['CNOT']+15*logical['CCX']+logical['CZ']+logical['Z']+2*qs*rr['includes_finite_uniform_preparation_H_gates'],
            'H_gates_uniform_preparation_upper':(shots+2*qs)*rr['includes_finite_uniform_preparation_H_gates'],
            'rotation_synthesis_error':0,'precision_scope':'exact logical Clifford+T gates; physical hardware error is not assumed zero experimentally'}
