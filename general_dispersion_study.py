"""Continuous-law diagnostic for the distribution-general dispersion theorem.

The theorem is the proof. Positive two-dimensional quadrature and refinement
are diagnostics, not certified integration or finite-quantum execution.
Risk-neutral drifts are normalised by each law's own exponential moments.
"""
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parent
D=Decimal
PAD=D('1e-100')
EPSILONS=[1,.5,.25,.125,.0625,.03125,.015625]
LAWS=('Gaussian','standardised_uniform','Rademacher','symmetric_variance_gamma')


def mgf_decimal(law,t):
    t=D(t)
    if not t:return D(1)
    if law=='Gaussian':return (t*t/2).exp()
    if law=='Rademacher':return (t.exp()+(-t).exp())/2
    if law=='standardised_uniform':
        a=D(3).sqrt()*t
        return (a.exp()-(-a).exp())/(2*a)
    if law=='symmetric_variance_gamma':
        nu=D('.4');base=1-nu*t*t/2
        assert base>0
        return (-base.ln()/nu).exp()
    raise ValueError(law)


def law_quadrature(law,order):
    if law=='Gaussian':
        x,w=np.polynomial.hermite.hermgauss(order)
        return math.sqrt(2)*x,w/math.sqrt(math.pi)
    if law=='standardised_uniform':
        x,w=np.polynomial.legendre.leggauss(order)
        return math.sqrt(3)*x,w/2
    if law=='Rademacher':return np.array([-1.,1.]),np.array([.5,.5])
    # Golub-Welsch Jacobi matrix for normalised Gamma(5/2,scale.4).
    # The first eigenvector components squared sum to one, so no Gamma
    # normalising constant is introduced as an unchecked float.
    gamma_order=order//2;alpha=1.5
    n=np.arange(gamma_order,dtype=float)
    diagonal=2*n+alpha+1
    k=np.arange(1,gamma_order,dtype=float)
    off=np.sqrt(k*(k+alpha))
    jacobi=np.diag(diagonal)+np.diag(off,1)+np.diag(off,-1)
    x,eigenvectors=np.linalg.eigh(jacobi)
    g=.4*x;wg=eigenvectors[0,:]**2
    z,wz=law_quadrature('Gaussian',order)
    y=(np.sqrt(g)[:,None]*z[None,:]).ravel()
    weight=(wg[:,None]*wz[None,:]).ravel()
    return y,weight


def fractional_residual(y,mu,t,order,along_order):
    """Positive density integration avoids subtraction at small dispersion."""
    k,wk=np.polynomial.legendre.leggauss(order)
    strike_log=np.log1p(.1*(k+1)/2)
    a,wa=np.polynomial.legendre.leggauss(along_order)
    along=(a+1)/2;along_weight=wa/2
    z=t*np.asarray(y)
    gap=np.log1p(2*np.sinh(z/2)**2)  # log cosh(z), stable at zero.
    out=np.empty_like(gap)
    for start in range(0,len(gap),1024):
        j=gap[start:start+1024]
        arguments=(mu+j[:,None,None]*along[None,:,None]-strike_log[None,None,:])/.2
        density=np.exp(-arguments*arguments/2)/math.sqrt(2*math.pi)
        out[start:start+1024]=j/.2*np.einsum('ijk,j,k->i',density,along_weight,wk/2)
    return out


def main():
    b=D('.2');lam=D('.2');m0=D('.01')-lam*lam/2
    x,w=np.polynomial.legendre.leggauss(128)
    derivative=float(np.sum(np.exp(-((float(m0)-np.log1p(.1*(x+1)/2))/.2)**2/2)/math.sqrt(2*math.pi)/.2*w/2))
    fourth={'Gaussian':3.,'standardised_uniform':1.8,'Rademacher':1.,'symmetric_variance_gamma':4.2}
    variance_square={k:v-1 for k,v in fourth.items()}
    records=[];quadrature_checks=[]
    for law in LAWS:
        y,weights=law_quadrature(law,64);fine_y,fine_w=law_quadrature(law,128)
        quadrature_checks.append({'law':law,'coarse_weight_sum':float(weights.sum()),'fine_weight_sum':float(fine_w.sum()),
            'coarse_E_Y2':float(weights@(y*y)),'fine_E_Y2':float(fine_w@(fine_y*fine_y)),
            'coarse_E_Y4':float(weights@(y**4)),'fine_E_Y4':float(fine_w@(fine_y**4)),
            'analytical_E_Y4':fourth[law]})
        for eps in EPSILONS:
            with localcontext() as ctx:
                ctx.prec=220
                ee=D(str(eps));t=b*ee
                moment=mgf_decimal(law,t)
                drift=m0-moment.ln()
                M=(mgf_decimal(law,2*t)/2-2*moment+D('1.5'))
                assert M>=0
                ratio=M/(ee**4)
                relative_limit=b**4*D(str(fourth[law]))/4
                # S_T=100exp(mu+-epsilon*b*Y+lambda*Z).
                risk_neutral_mean=100*drift.exp()*moment*(lam*lam/2).exp()
                assert abs(risk_neutral_mean-100*D('.01').exp())<PAD
            residual=fractional_residual(y,float(drift),float(t),48,12)
            refined=fractional_residual(fine_y,float(drift),float(t),96,24)
            mean=float(weights@residual);fm=float(fine_w@refined)
            variance=float(weights@((residual-mean)**2));fv=float(fine_w@((refined-fm)**2))
            c1=.02*derivative
            cvar=.0004*derivative*derivative*variance_square[law]
            records.append({'law':law,'dispersion_parameter':eps,'risk_neutral_log_drift_decimal':str(drift),
                'risk_neutral_stock_mean_decimal':str(risk_neutral_mean),
                'relative_gap_second_moment_decimal':str(M),
                'relative_gap_second_moment_numerical_enclosure_radius_decimal':str(PAD),
                'relative_gap_second_moment_divided_by_epsilon4_decimal':str(ratio),
                'relative_gap_second_moment_limit_decimal':str(relative_limit),
                'fractional_mean_coarse_diagnostic':mean,'fractional_mean_refined_diagnostic':fm,
                'fractional_variance_coarse_diagnostic':variance,'fractional_variance_refined_diagnostic':fv,
                'fractional_mean_divided_by_epsilon2_diagnostic':fm/(eps**2),'mean_limit_coefficient_diagnostic':c1,
                'fractional_variance_divided_by_epsilon4_diagnostic':fv/(eps**4),'variance_limit_coefficient_diagnostic':cvar,
                'canonical_flag_variance_diagnostic':fm*(1-fm),
                'flag_variance_divided_by_epsilon2_diagnostic':fm*(1-fm)/(eps**2),
                'fractional_mean_absolute_refinement_diagnostic':abs(fm-mean),
                'fractional_variance_absolute_refinement_diagnostic':abs(fv-variance),
                'degenerate_fractional_variance':law=='Rademacher',
                'qualification':'continuous-law positive quadrature diagnostic; no certified integration or current finite-grid asymptotic claim'})
    out={'status':'four-law weak-dispersion diagnostics with correct risk-neutral normalisation',
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'models':{'assets':2,'weights':['.5','.5'],'centred_loadings':['-.2','.2'],'lambda':'.2','strike':100,'cap':10,'rate':'.03','dividend_yield':'.02','variance_gamma_nu':'.4'},
        'common_conditioned_kernel_derivative_diagnostic':derivative,
        'quadrature_checks':quadrature_checks,'configurations':records,
        'refinement_orders':{'latent': [64,128], 'strike': [48,96], 'along_log_gap': [12,24]},
        'domains':{'VG_mgf':'|t|<sqrt(2/.4); largest required |t|=.4',
            'weak_theorem_norm_exponential':'2epsilon0B=.4; exponential moment at .4+delta withdelta=.1 finite, so degree4 integrability follows',
            'financial_scope':'Gaussian and variance-gamma give established financial latent model classes; bounded uniform/two-point laws are theoretical controls; no market calibration'},
        'novelty':'models and controls have published antecedents; the general theorem and encoding consequence require source comparison before a novelty claim'}
    (ROOT/'results/general_dispersion_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Completed',len(records),'continuous-law diagnostics across',len(LAWS),'laws.')


if __name__=='__main__':main()
