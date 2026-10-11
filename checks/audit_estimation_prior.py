"""Independent positive-series density and direct-moment prior check."""
from decimal import Decimal, localcontext
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = Decimal
REF_ERROR = D('1e-220')


def density_reference():
    # pi = 6 arcsin(1/2). Positive binomial-integral series; the successive
    # ratio is below 1/4, so the remaining sum is bounded by4/3 next_term.
    # This differs from the generator's alternating Machin series.
    with localcontext() as ctx:
        ctx.prec = 300
        term, total, n = D('.5'), D('.5'), 0
        while True:
            nxt = term*D((2*n+1)**2)/D(8*(n+1)*(2*n+3))
            if nxt < D('1e-260'):
                assert n < 1000
                pi_lo = 6*total-D('1e-280')
                pi_hi = 6*(total+D(4)*nxt/3)+D('1e-280')
                centre = 1/(2*((pi_lo+pi_hi)/2)).sqrt()
                error = max(centre-1/(2*pi_hi).sqrt(),
                            1/(2*pi_lo).sqrt()-centre)+D('1e-270')
                assert error < REF_ERROR
                return centre
            total += nxt; term = nxt; n += 1


def independent_mean(model, cap, phi0):
    with localcontext() as ctx:
        ctx.prec = 300
        fweights = list(map(Fraction, model['weights']))
        fr = [list(map(Fraction, rr)) for rr in model['loadings']]
        fmu = list(map(Fraction, model['log_drifts']))
        fmg = sum(w*m for w,m in zip(fweights,fmu))
        frg = [sum(fweights[i]*fr[i][j] for i in range(len(fr)))
               for j in range(model['outer_dimension'])]
        lam = Fraction(model['lambda'])
        cvt = lambda f: D(f.numerator)/D(f.denominator)
        ea = sum(cvt(w)*D(model['spot'])*
                 (cvt(mu+(lam*lam+sum(x*x for x in rr))/2)).exp()
                 for w,mu,rr in zip(fweights,fmu,fr))
        eg = D(model['spot'])*(cvt(fmg+(lam*lam+sum(x*x for x in frg))/2)).exp()
        egap = (sum(cvt(w)*(cvt(mu-fmg+
                    sum((a-b)**2 for a,b in zip(rr,frg))/2)).exp()
                    for w,mu,rr in zip(fweights,fmu,fr))).ln()
        h = D(cap)/D(model['strike'])
        slope = min(phi0/cvt(lam),(1+h)/h)
        return {'ag': ea-eg, 'J0': egap, 'slope': slope,
                'slope_mean': D(cap)*slope*egap,
                'mean': min(D(cap),ea-eg,D(cap)*slope*egap)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=ROOT/'results/generalisation_resources.json')
    parser.add_argument('--priors', type=Path, default=ROOT/'results/estimation_priors.json')
    parser.add_argument('--output', type=Path, default=ROOT/'results/estimation_prior_audit.json')
    args = parser.parse_args()
    config = json.loads(args.config.read_text()); priors = json.loads(args.priors.read_text())
    assert priors['configuration_sha256'] == hashlib.sha256(args.config.read_bytes()).hexdigest()
    assert priors['script_sha256'] == hashlib.sha256((ROOT/'estimation_prior.py').read_bytes()).hexdigest()
    models = {m['model_id']:m for m in config['models']}
    paired = {(r['case_id'],r['method']):r for r in config['configurations']}
    phi0 = density_reference(); records = []
    with localcontext() as ctx:
        ctx.prec = 300
        for row in priors['configurations']:
            source = paired[row['case_id'],row['method']]
            model = models[source['model_id']]
            reference = independent_mean(model,row['cap'],phi0)
            envelope = D(row['decimal_numerical_allowance'])
            assert abs(reference['ag']-D(row['arithmetic_geometric_mean_currency_centre_decimal']))+REF_ERROR <= envelope
            assert abs(reference['slope_mean']-D(row['slope_Jensen_mean_currency_centre_decimal']))+REF_ERROR <= envelope
            assert D(row['mean_currency_upper_decimal'])+REF_ERROR >= reference['mean']
            # Required sample budget uses exact promised0.005, not a binary
            # approximation that can exceed the statistical allocation.
            target = D('.005')/D(row['cap'])
            assert D(row['statistical_probability_target_decimal']) >= target-REF_ERROR
            assert D(row['statistical_probability_target_decimal']) <= target+REF_ERROR
            upper = D(row['finite_mean_probability_upper_decimal'])
            residual_bias = D.from_float(source['arithmetic_price_error_upper'])+D.from_float(source['distribution_price_bound'])
            expected_upper = min(D(1),(reference['mean']+residual_bias)/D(row['cap']))
            assert upper+REF_ERROR >= expected_upper
            variance = D(row['finite_variance_upper_decimal'])
            expected_var = upper*(1-upper) if upper <= D('.5') else D('.25')
            assert variance+REF_ERROR >= expected_var
            alpha = D(row['failure_probability']); n = row['classical_certified_samples']
            if row['classical_selected_bound'] == 'Bernstein':
                failure = 2*(-D(n)*target*target/(2*variance+2*target/3)).exp()+REF_ERROR
            else:
                failure = 2*(-2*D(n)*target*target).exp()+REF_ERROR
            assert failure <= alpha
            records.append({'case_id':row['case_id'],'method':row['method'],
                'whole_reference_enclosure_within_generator_allowance':True,
                'finite_amplitude_prior_conservative':True,
                'statistical_allocation_exact':True,
                'selected_classical_samples':n,
                'independent_failure_upper_decimal':str(failure)})
    output = {'status':'independent direct-lognormal moments, positive arcsine density and confidence budgets passed',
        'configuration_sha256':hashlib.sha256(args.config.read_bytes()).hexdigest(),
        'prior_sha256':hashlib.sha256(args.priors.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reference_numerical_enclosure':str(REF_ERROR), 'configurations':records,
        'qualification':'conditional on correctly rounded Decimal operations; analytical confidence guarantee concerns independent exact finite draws'}
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print('Independent financial priors and confidence budgets passed:',len(records))


if __name__ == '__main__':
    main()
