"""Complete finite-confidence canonical QPE amplitude-estimation plans.

The amplitude theorem is Brassard et al. (2002), Theorem 12. Financial
priors, exact binary ceilings, finite arithmetic allowances and emitted
controlled-Grover/QFT charges are project-specific specialisations. The
sampling adapter is an ideal diagnostic, not physical quantum execution.
"""
from dataclasses import dataclass, asdict
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import cmath
import math

from certified_scalar import (ceil_fraction, pi_rational_interval,
                              probability_interval, sqrt_upper_fraction)

PI_UPPER = F(355, 113)
POSTPROCESS_ALLOWANCE = F(1, 10**60)
QFT_NORM_ERROR = F(1, 64)
CONFIDENCE_NORM_GRID_BITS = 80
assert pi_rational_interval()[1] < PI_UPPER


@lru_cache(None)
def log_upper(value, terms=80):
    """Exact rational enclosure of ln(value) for value >= 1.

    Binary reduction leaves a mantissa in [1,2); the positive atanh
    series has z <= 1/3 and its omitted tail is geometrically bounded.
    """
    value = F(value)
    if value < 1:
        raise ValueError('Logarithm input must be at least one')
    exponent = 0
    while value >= 2:
        value /= 2
        exponent += 1
    def series(x):
        z = (x-1)/(x+1)
        total = sum((2*z**(2*j+1)/F(2*j+1)
                     for j in range(terms)), F(0))
        tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return total+tail
    return exponent*series(F(2))+series(value)


def binomial_failure_upper(repetitions, good_probability_lower):
    """Exact probability of at least half bad outcomes for odd repeats."""
    n, s = repetitions, F(good_probability_lower)
    if n < 1 or n % 2 == 0 or not F(1, 2) < s <= 1:
        raise ValueError('Odd positive repetitions and majority success required')
    bad = 1-s
    return sum((F(comb(n,j))*bad**j*s**(n-j)
                for j in range((n+1)//2,n+1)), F(0))


@lru_cache(None)
def confidence_repetitions(alpha, norm_error=QFT_NORM_ERROR,
                           confidence='exact_binomial'):
    alpha, norm_error = F(alpha), F(norm_error)
    if not 0 < alpha < 1 or not 0 <= norm_error <= QFT_NORM_ERROR:
        raise ValueError('Invalid failure probability or QFT norm allowance')
    success = 8/PI_UPPER**2-norm_error
    gap = success-F(1, 2)
    assert gap > 0
    hoeffding = max(1, ceil_fraction(log_upper(1/alpha)/(2*gap*gap)))
    if hoeffding % 2 == 0:
        hoeffding += 1
    if confidence == 'hoeffding':
        n = hoeffding
    elif confidence == 'exact_binomial':
        n = 1
        while binomial_failure_upper(n,success) > alpha:
            n += 2
    else:
        raise ValueError('Unsupported finite confidence rule')
    return n, success, binomial_failure_upper(n,success), hoeffding


@dataclass(frozen=True)
class Plan:
    epsilon: F
    alpha: F
    prior_upper: F
    phase_register_bits: int
    M: int
    repetitions: int
    variance_upper: F
    amplitude_error_upper: F
    postprocess_allowance: F
    qft_norm_error_allowance: F
    good_probability_lower: F
    binomial_failure_upper: F
    hoeffding_repetitions_upper: int
    confidence_rule: str

    @property
    def zero_shots(self):
        return self.repetitions == 0

    @property
    def total_shots_upper(self):
        return self.repetitions

    @property
    def total_controlled_grover_calls(self):
        return self.repetitions*(self.M-1)

    @property
    def total_grover_calls_upper(self):
        return self.total_controlled_grover_calls

    def record(self):
        return {k:(str(v) if isinstance(v,F) else v)
                for k,v in asdict(self).items()}


@lru_cache(maxsize=4096)
def plan(epsilon, alpha=F(1, 100), prior_upper=F(1),
         confidence='exact_binomial', norm_error=QFT_NORM_ERROR):
    """Minimal supported power-of-two M, followed by fresh odd repeats.

    The zero-shot midpoint has deterministic error at most prior_upper/2.
    All active plans reserve a certified 1e-60 output arithmetic allowance.
    """
    epsilon,alpha,prior,norm_error=map(F,(epsilon,alpha,prior_upper,norm_error))
    if not F(1,10**6) <= epsilon < F(1,2):
        raise ValueError('Declared probability tolerance is 1e-6 <= epsilon < .5')
    if not 0 < alpha < 1 or not 0 <= prior <= 1:
        raise ValueError('Invalid failure probability or certified prior')
    if not 0 <= norm_error <= QFT_NORM_ERROR:
        raise ValueError('Inverse QFT error exceeds its declared allowance')
    if confidence not in ('exact_binomial','hoeffding'):
        raise ValueError('Unsupported finite confidence rule')
    variance = prior*(1-prior) if prior <= F(1,2) else F(1,4)
    if prior <= 2*epsilon:
        return Plan(epsilon,alpha,prior,0,1,0,variance,prior/2,F(0),
                    F(0),F(1),F(0),0,confidence)
    root = sqrt_upper_fraction(variance)
    bits,M = 1,2
    while True:
        error = 2*PI_UPPER*root/M+PI_UPPER**2/(M*M)
        if error+POSTPROCESS_ALLOWANCE <= epsilon:
            break
        bits += 1
        M *= 2
    if bits > 32:
        raise ValueError('Phase-register certificate is limited to 32 bits')
    repeats,success,fail,hoeffding=confidence_repetitions(alpha,norm_error,confidence)
    assert fail <= alpha
    return Plan(epsilon,alpha,prior,bits,M,repeats,variance,error,
                POSTPROCESS_ALLOWANCE,norm_error,success,fail,hoeffding,confidence)


def compiled_plan(epsilon,alpha=F(1,100),prior_upper=F(1),
                  confidence='exact_binomial',profile_getter=None):
    """Use the actual retained inverse-QFT error in the confidence constant.

    The binary M depends on amplitude tolerance and prior, not Fourier error.
    Its actual width certificate is therefore available before confidence
    amplification, without a circular choice of M or gate synthesis.
    """
    preliminary=plan(epsilon,alpha,prior_upper,confidence=confidence)
    if preliminary.zero_shots:
        return preliminary,None
    if profile_getter is None:
        from qft_synthesis import inverse_qft_profile
        profile_getter=inverse_qft_profile
    profile=profile_getter(preliminary.phase_register_bits)
    if profile['phase_register_bits'] != preliminary.phase_register_bits:
        raise ValueError('Compiled Fourier profile width is inconsistent')
    actual_error=F(profile['norm_error_upper'])
    if not 0<=actual_error<=QFT_NORM_ERROR:
        raise ValueError('Compiled Fourier profile has an unsupported norm error')
    scale=1<<CONFIDENCE_NORM_GRID_BITS
    error=F(ceil_fraction(actual_error*scale),scale)
    # Outward rounding keeps the guarantee while avoiding unnecessarily huge
    # denominator powers in exact binomial confidence records. The grid also
    # represents the maximum admissible 1/64 bound exactly.
    assert actual_error<=error<=QFT_NORM_ERROR
    setup=plan(epsilon,alpha,prior_upper,confidence=confidence,norm_error=error)
    assert setup.M==preliminary.M
    return setup,profile


def output_interval(y,M):
    """Outward exact-rational enclosure of sin^2(pi*y/M)."""
    if M < 2 or M & (M-1) or not 0 <= y < M:
        raise ValueError('Invalid binary phase outcome')
    folded=min(y,M-y)
    lo,hi=probability_interval(F(2*folded,M))
    assert hi-lo <= 2*POSTPROCESS_ALLOWANCE
    return lo,hi


def certified_median(setup, outcomes):
    """Median of outward-enclosed output midpoints, projected onto the prior."""
    if setup.zero_shots:
        if list(outcomes):
            raise ValueError('Zero-shot plans have no measurement outcomes')
        return {'estimate':setup.prior_upper/2,
                'interval':(F(0),setup.prior_upper),'outcomes':()}
    outcomes=tuple(outcomes)
    if len(outcomes) != setup.repetitions:
        raise ValueError('Complete fresh repetition transcript required')
    estimates=[]
    for y in outcomes:
        lo,hi=output_interval(y,setup.M)
        estimates.append(min(setup.prior_upper,max(F(0),(lo+hi)/2)))
    estimates.sort()
    estimate=estimates[len(estimates)//2]
    return {'estimate':estimate,
            'interval':(max(F(0),estimate-setup.epsilon),
                        min(setup.prior_upper,estimate+setup.epsilon)),
            'outcomes':outcomes}


def ideal_phase_probability(y,M,omega):
    """Dirichlet-kernel probability; floating-point diagnostic only."""
    delta=math.remainder(omega-y/M,1.0)
    if abs(delta) < 1e-15:
        return 1.0
    return (math.sin(math.pi*M*delta)/(M*math.sin(math.pi*delta)))**2


def ideal_amplitude_distribution(p,M):
    """Reduced phase-register law for the two Grover eigenphases."""
    if not 0 <= p <= 1:
        raise ValueError('Probability outside [0,1]')
    omega=math.asin(math.sqrt(p))/math.pi
    if p == 0:
        return [1.0]+[0.0]*(M-1)
    if p == 1:
        return [float(y==M//2) for y in range(M)]
    return [(ideal_phase_probability(y,M,omega)
             +ideal_phase_probability(y,M,1-omega))/2 for y in range(M)]


def execute_ideal(setup,p,rng,maximum_distribution_size=1<<20):
    """Seeded ideal simulation of independent QPE outcomes, never a QPU call."""
    p=float(p)
    if not 0 <= p <= float(setup.prior_upper):
        raise ValueError('Simulation probability violates the certified prior')
    if setup.zero_shots:
        result=certified_median(setup,())
    else:
        if setup.M > maximum_distribution_size:
            raise ValueError('Diagnostic categorical law exceeds its explicit size cap')
        probabilities=ideal_amplitude_distribution(p,setup.M)
        total=sum(probabilities)
        if abs(total-1) > 1e-9:
            raise ArithmeticError('Floating-point Dirichlet diagnostic did not normalise')
        # The ideal random source contract is independent categorical outcomes.
        # The seed itself does not prove the physical independence premise.
        cumulative=[];s=0.0
        for value in probabilities:
            s+=value/total
            cumulative.append(s)
        import bisect
        outcomes=[min(setup.M-1,bisect.bisect_left(cumulative,rng.random()))
                  for _ in range(setup.repetitions)]
        result=certified_median(setup,outcomes)
    result['execution_scope']='seeded ideal Dirichlet categorical simulation; no physical quantum device'
    return result


def ideal_inverse_qft_gates(bits,wires=None):
    """Little-endian inverse DFT, including final output-order swaps."""
    wires=tuple(range(bits)) if wires is None else tuple(wires)
    if bits < 1 or len(wires) != bits or len(set(wires)) != bits:
        raise ValueError('Invalid Fourier-register wires')
    forward=[]
    for j in range(bits-1,-1,-1):
        forward.append(('H',wires[j]))
        for k in range(j-1,-1,-1):
            forward.append(('CP',F(1,1<<(j-k)),wires[k],wires[j]))
    for j in range(bits//2):
        forward.append(('SWAP',wires[j],wires[bits-1-j]))
    return [(('CP',-g[1],*g[2:]) if g[0]=='CP' else g)
            for g in reversed(forward)]


def inverse_qft_symbolic_counts(bits):
    if not 1 <= bits <= 32:
        raise ValueError('Supported phase register has 1 through 32 bits')
    return {'H':bits,'controlled_phase':bits*(bits-1)//2,
            'SWAP':bits//2,'CX_for_SWAP':3*(bits//2),
            'CX_for_controlled_phase':bits*(bits-1),
            'single_Rz_rotations_before_synthesis':3*bits*(bits-1)//2}


def charge(setup,a_profile,controlled_q_profile,inverse_qft_profile):
    """Charge emitted complete profiles, with no financial table recomputation.

    Each supplied profile has metrics in `gate_counts`, or direct numeric
    metrics in `metrics`; the dictionaries are additive counts/serial caps.
    The controlled-Q profile already includes its global-minus Z, controlled
    good predicate and controlled zero reflection. It is never inferred by
    treating Q as an uncontrolled invocation.
    """
    if setup.zero_shots:
        return {'repetitions':0,'controlled_Q_calls':0,'metrics':{},
                'peak_logical_qubits':0,'phase_register_measurements':0,
                'complete_reinitialisations':0,'source':'zero-shot prior midpoint'}
    if F(inverse_qft_profile['norm_error_upper']) > setup.qft_norm_error_allowance:
        raise ValueError('Supplied emitted inverse QFT exceeds the plan norm allowance')
    if inverse_qft_profile['phase_register_bits'] != setup.phase_register_bits:
        raise ValueError('Fourier profile width does not match the plan')
    def metrics(profile):
        return profile.get('metrics',profile.get('gate_counts',{}))
    ap,qp,fp=map(metrics,(a_profile,controlled_q_profile,inverse_qft_profile))
    keys=set(ap)|set(qp)|set(fp)
    per={k:ap.get(k,0)+(setup.M-1)*qp.get(k,0)+fp.get(k,0)
         for k in keys}
    h_key='gate_H' if 'gate_H' in keys else 'H'
    per[h_key]=per.get(h_key,0)+setup.phase_register_bits
    if 'expanded_primitive_operations' in keys:
        per['expanded_primitive_operations']+=setup.phase_register_bits
    if 'uniform_H' in keys:
        per['uniform_H']+=setup.phase_register_bits
    # Include terminal phase observations in the total observation metric;
    # their separate field is an explanation of that included component.
    per['measurements']=per.get('measurements',0)+setup.phase_register_bits
    totals={k:setup.repetitions*v for k,v in per.items()}
    # Financial scratch, reflection scratch and synthesis ancillas may all
    # be reused between invocations. Adding supplied maxima is conservative.
    financial=max(a_profile.get('logical_qubits',a_profile.get('qubits',0)),
                  controlled_q_profile.get('logical_qubits',controlled_q_profile.get('qubits',0)))
    extra=inverse_qft_profile.get('extra_ancillas',0)
    return {'repetitions':setup.repetitions,
            'controlled_Q_calls':setup.total_controlled_grover_calls,
            'per_repetition_metrics':per,'metrics':totals,
            'peak_logical_qubits':financial+setup.phase_register_bits+extra,
            'phase_register_measurements':setup.repetitions*setup.phase_register_bits,
            'complete_reinitialisations':setup.repetitions,
            'source':'one emitted A, M-1 controlled Q and one emitted inverse QFT per fresh repetition',
            'depth_scope':'supplied depth metrics are serial additive upper bounds; no parallel scheduling claim'}
