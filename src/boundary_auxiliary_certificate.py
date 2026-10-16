"""Paid deterministic auxiliary allowances for the finite financial law."""
from fractions import Fraction as F
from boundary_derivative_certificate import payoff_norms
from boundary_value_certificate import gauss_constants


def auxiliary_bias(lam,cap,U_star,MU2,MU3,beta_bits=8,q=28,kind='boundary_r2',MU4=None,MU6=None,K=100):
    U,M2u,M3u=map(F,(U_star,MU2,MU3))
    if min(U,M2u,M3u)<0 or beta_bits<1:raise ValueError('Nonnegative finite-support moments required')
    norm=payoff_norms(lam,cap,K);M2=F(norm['M2_upper']);M3=F(norm['C_r_upper']['3'])
    M4=F(norm['C_r_upper']['4']);M6=F(norm['C_r_upper']['6'])
    M4u=U*M3u if MU4 is None else F(MU4)
    M6u=U**3*M3u if MU6 is None else F(MU6)
    if min(M4u,M6u)<0:raise ValueError('Nonnegative moments required')
    if kind=='boundary_r2':
        h=F(1,1<<beta_bits);error=h*h*(M3*M3u/12+M4*M4u/48)
        return {'kind':kind,'beta_bits':beta_bits,'Beta_cell_first_moment':'-h^3/6',
            'Beta_cell_second_moment_sum':'h^2/12','bias_upper':str(error),
            'finite_U_moment_2_upper':str(M2u),'finite_U_moment_3_upper':str(M3u),
            'finite_U_moment_4_upper':str(M4u),'M4_support_bound_used':MU4 is None,
            'premise':'Moments of the ideal financial U on the actual outer dyadic law; scalar and arithmetic discrepancies are paid separately.'}
    if kind=='stable_gauss_r2':
        constants=gauss_constants(q)
        node=F(0);weight=F(0)
        for suffix in ('minus','plus'):
            zr=constants['z_'+suffix];wr=constants['w_'+suffix]
            node+=F(wr['ideal_upper'])*F(zr['error_upper'])
            weight+=F(wr['error_upper'])
        trunc=M6*M6u/14400;fixed=M3*M3u*node/2+M2*M2u*weight/2
        return {'kind':kind,'bias_upper':str(trunc+fixed),'Gauss_degree':'3',
            'monic_orthogonal_polynomial':'z^2-(4/5)z+1/10','squared_norm':'1/300',
            'Gauss_remainder_upper':str(trunc),'fixed_parameter_bias_upper':str(fixed),
            'finite_U_moment_6_upper':str(M6u),'M6_support_bound_used':MU6 is None,
            'constants':constants,
            'premise':'All moments refer to the identical ideal U on the actual outer dyadic law. Actual dyadic weights sum one. Constant encoding error is not charged twice in the arithmetic graph.'}
    if kind in ('raw_direct_r2','positive_r1'):return {'kind':kind,'bias_upper':'0','auxiliary_free':True}
    raise ValueError('Unknown correction kind')
