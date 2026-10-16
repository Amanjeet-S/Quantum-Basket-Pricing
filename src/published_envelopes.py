"""Published IQAE/MIQAE theorem envelopes, separate from my finite controller.

These are costs of the source algorithms under their stated ideal premises.
They are not a transfer of a published coefficient to full_batch_miqae.py.
The MIQAE full-batch integer rounding corollary is labelled explicitly.
"""
from fractions import Fraction as F
from decimal import Decimal,localcontext
from certified_scalar import (pi_rational_interval,probability_interval,
    sqrt_upper_fraction,ceil_fraction,floor_fraction)
from full_batch_miqae import log_upper


def log_lower(x):
    with localcontext() as ctx:
        ctx.prec=80
        value=Decimal(x.numerator)/Decimal(x.denominator)
        lower=value.ln()-Decimal('1e-60')
    return F(lower)


def financial_angle_parameter(epsilon,prior_upper=F(1)):
    epsilon,prior=F(epsilon),F(prior_upper)
    if not 0<epsilon<F(1,2) or not 0<=prior<=1:raise ValueError('Invalid tolerance/prior')
    if prior<=2*epsilon:return {'zero_shots':True,'gamma':F(0),'epsilon_angle':F(0)}
    variance=prior*(1-prior) if prior<=F(1,2) else F(1,4)
    gamma=2*sqrt_upper_fraction(variance)
    return {'zero_shots':False,'gamma':gamma,'epsilon_angle':epsilon/gamma}


def theorem_envelopes(epsilon,alpha=F(1,100),prior_upper=F(1)):
    epsilon,alpha=F(epsilon),F(alpha)
    if not 0<alpha<1:raise ValueError('Invalid failure probability')
    parameters=financial_angle_parameter(epsilon,prior_upper)
    if parameters['zero_shots']:
        empty={'total_shots_upper':0,'total_grover_calls_upper':0,'grover_power_upper':0,'maximum_rounds':0}
        return {'financial_postprocessing':parameters,'iqae_published':dict(empty),'miqae_published_continuous_shot_term':dict(empty),'miqae_full_batch_integer_corollary':dict(empty)}
    eta=parameters['epsilon_angle']
    if eta>=F(1,4):raise ValueError('Declared financial source-envelope range requires physical angle tolerance<.25 radians')
    pi_upper=pi_rational_interval()[1]
    l2_lower=log_lower(F(2));log2_upper=log_upper(pi_upper/(4*eta))/l2_lower
    logarithm=log_upper((2/alpha)*log2_upper)
    rounds_iqae=max(1,ceil_fraction(log_upper(pi_upper/(8*eta))/l2_lower))
    sin_squared_upper=probability_interval(F(1,7))[1]
    sine_upper=sqrt_upper_fraction(sin_squared_upper)
    denominator=1-2*sine_upper
    assert denominator>0
    max_shots_round_upper=32*logarithm/(denominator*denominator)
    structural_rounds=1;source_multiplier=F(1);multiplier_bound=pi_upper/(4*eta)
    while 2*source_multiplier<=multiplier_bound:source_multiplier*=2;structural_rounds+=1
    # At cumulative N>=Nmax, the sufficient source circle-growth condition
    # holds. A pool can cross that threshold by at most one100-shot batch.
    iqae_shots=structural_rounds*(ceil_fraction(max_shots_round_upper)+100)
    iqae_queries=ceil_fraction(50*logarithm/eta)
    # Choose Nshots=100 only after proving this integer is in the source range.
    # The lower bound32*ln(...) suffices because the source multiplier is>32.
    assert 32*log_lower((2/alpha)*(log_lower(pi_rational_interval()[0]/(4*eta))/log_upper(F(2))))>=100
    iqae={'source':'Grinko et al.2021, final journal HTML Theorem1 equations6-8 and Algorithms1-2',
          'algorithm_scope':'published algorithm with base Nshots100 and its no-overshooting rule; theorem-based envelope, not a retained QPU execution',
          'maximum_rounds':structural_rounds,'published_maximum_rounds':rounds_iqae,'base_shots':100,
          'published_reported_Q_term_upper':iqae_queries,
          'independent_structural_Q_cap_upper':ceil_fraction((max_shots_round_upper+100)*multiplier_bound),
          'shot_rounding_scope':'structural doubling-round bound includes initialk0; cumulative Nmax threshold plus one100-shot batch perround; source reported50 query term kept distinct',
          'maximum_shots_per_round_real_upper':max_shots_round_upper,
          'total_shots_upper':iqae_shots,'total_grover_calls_upper':iqae_queries,
          'grover_power_upper':ceil_fraction(pi_upper/(8*eta)),
          'classical_assumption':'exact source-algorithm real postprocessing; coefficients not assigned to my modified controller'}
    bound=pi_upper/(4*eta)
    rounds_miqae=1;m=1
    while m*3<=bound:m*=3;rounds_miqae+=1
    delta_lower=probability_interval(F(3,7))[0]-F(1,4)
    assert delta_lower>F(13,100)
    coefficient=2/(delta_lower*delta_lower)
    initial_real_shots_upper=coefficient*log_upper(3*bound/alpha)
    assert initial_real_shots_upper<F(1000000), 'Declared integer inverse-padding guard domain exceeded'
    published_queries=ceil_fraction(62*log_upper(6/alpha)/eta)
    published={'source':'Fukuzawa et al.2023, published ALENEX chapter metadata verified; actual text consulted arXivv4 sections2-3',
        'algorithm_scope':'source full-batch choice, with its real Nmax expressions; source query term is retained separately from integer implementation',
        'maximum_rounds_including_k0':rounds_miqae,
        'initial_maximum_real_shots_upper':initial_real_shots_upper,
        'total_shots_upper':ceil_fraction(rounds_miqae*initial_real_shots_upper),
        'total_grover_calls_upper':published_queries,'grover_power_upper':ceil_fraction((bound-1)/2),
        'source_version_limit':'final publisher chapter full text was not accessed; no implication that its final text was read'}
    # Ni=ceil(Nmax_i)+one numerical guard shot can cost at mostNmax_i+3
    # in the retained domain. Σki<=3M/4 and at mostR rounds charge this
    # explicitly. My rational-radius controller has its own DP budget instead.
    integer={**published,'algorithm_scope':'documented full-batch source method with upward integer shots and one guard shot; a derived resource corollary, not the unchanged62 bound',
        'total_shots_upper':rounds_miqae*ceil_fraction(initial_real_shots_upper+3),
        'total_grover_calls_upper':published_queries+ceil_fraction(9*bound/4),
        'integer_padding_extra_queries_upper':ceil_fraction(9*bound/4),
        'guard_scope':'source real bounds plus<=3 additional shots perround; Nmax<=1e6 makes one extra guard shot decrease the probability span by>6e-8, exceeding inversepadding8*2^-70 and rational square-root rounding; monetary endpoints still require their explicit allowance'}
    return {'financial_postprocessing':parameters,'iqae_published':iqae,
        'miqae_published_continuous_shot_term':published,'miqae_full_batch_integer_corollary':integer}


def linear_gate_budget(caps,clean_A_T,grover_T):
    shots=caps['total_shots_upper'];queries=caps['total_grover_calls_upper']
    return {'total_shots_upper':shots,'total_grover_calls_upper':queries,
        'A_and_adjoint_calls_upper':shots+2*queries,
        'T_count_upper':shots*clean_A_T+queries*grover_T,
        'flag_measurements_upper':shots,'scope':'logical gates of the same finite oracle; no comparison of physical runtime units'}
