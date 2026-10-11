"""Analytical residual-amplitude prior and certified classical sample budgets.

The bounded-payoff and lognormal-moment argument is established mathematics.
This is a project-specific parameter certificate, not a new quantum theorem.
"""
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path

D = Decimal
ROOT = Path(__file__).resolve().parent
PAD = D('1e-100')


def exact_dispersion(model):
    weights = list(map(Fraction, model['weights']))
    loadings = [list(map(Fraction, row)) for row in model['loadings']]
    assert all(w > 0 for w in weights) and sum(weights) == 1
    variances = [sum(x*x for x in row) for row in loadings]
    geometric = [sum(weights[i]*loadings[i][j] for i in range(len(weights)))
                 for j in range(model['outer_dimension'])]
    value = sum(w*v for w,v in zip(weights, variances))-sum(x*x for x in geometric)
    assert value >= 0
    assert value == Fraction(model['weighted_residual_dispersion'])
    rate = Fraction(model['risk_free_rate'])-Fraction(model['dividend_yield'])
    lam = Fraction(model['lambda'])
    assert all(Fraction(mu) == rate-(v+lam*lam)/2
               for mu,v in zip(model['log_drifts'], variances))
    return value


def normal_density_at_zero():
    with localcontext() as ctx:
        ctx.prec = 220
        def atan(x):
            term, total, n = x, x, 0
            while abs(term) > D('1e-230'):
                n += 1
                term *= -x*x
                total += term/D(2*n+1)
            return total
        pi = 16*atan(D(1)/5)-4*atan(D(1)/239)
        return 1/(2*pi).sqrt()


def mean_bound(model, cap):
    dispersion = exact_dispersion(model)
    with localcontext() as ctx:
        ctx.prec = 220
        dd = D(dispersion.numerator)/D(dispersion.denominator)
        # In this model the spot equals the common normalising strike.
        assert model['spot'] == model['strike'] == 100
        growth = (D(model['risk_free_rate'])-D(model['dividend_yield'])).exp()
        ag_centre = D(model['spot'])*growth*(1-(-dd/2).exp())
        weights = list(map(Fraction, model['weights']))
        loadings = [list(map(Fraction, row)) for row in model['loadings']]
        rg = [sum(weights[i]*loadings[i][j] for i in range(len(weights)))
              for j in range(model['outer_dimension'])]
        def decimal_fraction(x):
            return D(x.numerator)/D(x.denominator)
        exponential_mean = sum(decimal_fraction(w)*
            decimal_fraction(dispersion/2-sum((a-b)*b for a,b in zip(rr,rg))).exp()
            for w,rr in zip(weights,loadings))
        jensen = exponential_mean.ln()
        assert jensen >= 0
        h = D(cap)/D(model['strike'])
        slope = min(normal_density_at_zero()/D(model['lambda']), (1+h)/h)
        slope_jensen_centre = D(cap)*slope*jensen
        upper = min(D(cap), ag_centre+PAD, slope_jensen_centre+PAD)
        return {'dispersion_exact_numerator': str(dispersion.numerator),
                'dispersion_exact_denominator': str(dispersion.denominator),
                'mean_currency_upper_decimal': str(upper),
                'arithmetic_geometric_mean_currency_centre_decimal': str(ag_centre),
                'slope_Jensen_mean_currency_centre_decimal': str(slope_jensen_centre),
                'Jensen_log_gap_mean_upper_decimal': str(jensen+PAD),
                'conditional_kernel_slope_upper_decimal': str(slope+PAD),
                'decimal_numerical_allowance': str(PAD),
                'numerical_assumption': 'correctly rounded Decimal exp/ln/sqrt at220 significant digits; Machin alternating-series pi with 1e-230 stop; no exponent overflow/underflow'}


def prior_for_row(model, row, failure_probability=D('.01')):
    assert 0 < failure_probability < 1
    mean = mean_bound(model, row['cap'])
    with localcontext() as ctx:
        ctx.prec = 220
        cap = D(row['cap'])
        # Recorded component bounds are outward binary floats. Their decimal
        # strings round-trip; use exact float values rather than shortened
        # strings when forming a new certified budget.
        arithmetic = D.from_float(row['arithmetic_price_error_upper'])
        distribution = D.from_float(row['distribution_price_bound'])
        residual_bias = arithmetic+distribution
        upper = min(D(1), (D(mean['mean_currency_upper_decimal'])+residual_bias)/cap+PAD)
        assert upper >= 0
        variance = upper*(1-upper)+PAD if upper <= D('.5') else D('.25')
        target_currency = D(str(row['price_tolerance']))/2
        assert D.from_float(row['reserved_statistical_price_allowance']) >= target_currency
        target_probability = target_currency/cap
        assert target_probability > 0
        log_term = (2/failure_probability).ln()+PAD
        # Classical samples are the finite fractional outputs in [0,1].
        bernstein = ((2*variance+2*target_probability/3)*log_term/
                     (target_probability*target_probability)+PAD)
        hoeffding = log_term/(2*target_probability*target_probability)+PAD
        n_bernstein = int(bernstein.to_integral_value(rounding=ROUND_CEILING))
        n_hoeffding = int(hoeffding.to_integral_value(rounding=ROUND_CEILING))
        bound_at_n = 2*(-D(n_bernstein)*target_probability*target_probability/
                        (2*variance+2*target_probability/3)).exp()+PAD
        assert bound_at_n <= failure_probability
        chosen = min(n_bernstein, n_hoeffding)
        method = 'Bernstein' if n_bernstein <= n_hoeffding else 'Hoeffding'
        chosen_failure = (bound_at_n if method == 'Bernstein' else
                          2*(-2*D(chosen)*target_probability*target_probability).exp()+PAD)
        assert chosen_failure <= failure_probability
        return {'case_id': row['case_id'], 'method': row['method'], 'cap': row['cap'],
                **mean, 'finite_residual_bias_upper_decimal': str(residual_bias+PAD),
                'finite_mean_probability_upper_decimal': str(upper),
                'finite_variance_upper_decimal': str(variance),
                'statistical_currency_target_decimal': str(target_currency),
                'statistical_probability_target_decimal': str(target_probability),
                'failure_probability': str(failure_probability),
                'classical_Bernstein_samples': n_bernstein,
                'classical_Hoeffding_samples': n_hoeffding,
                'classical_certified_samples': chosen,
                'classical_selected_bound': method,
                'classical_selected_failure_bound_decimal': str(chosen_failure),
                'failure_bound_at_Bernstein_budget_decimal': str(bound_at_n),
                'confidence_scope': 'independent exact finite-output draws, mathematical/numerical assumptions and source deterministic bias certificate; no QPU noise model',
                'cost_scope': 'number of classical finite-function evaluations; no host operation or runtime comparison with logical T gates'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, default=ROOT/'results/generalisation_resources.json')
    parser.add_argument('--output', type=Path, default=ROOT/'results/estimation_priors.json')
    parser.add_argument('--failure-probability', type=D, default=D('.01'))
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    models = {m['model_id']: m for m in config['models']}
    rows = [prior_for_row(models[row['model_id']], row, args.failure_probability)
            for row in config['configurations']]
    result = {'status': 'analytical model-specific amplitude/variance certificate and classical finite-output budgets',
              'configuration_sha256': hashlib.sha256(args.config.read_bytes()).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'configurations': rows,
              'contribution_limit': 'known Lipschitz/AM-GM/lognormal moments and Bernstein concentration; useful project parameter, theoretical novelty unverified'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'configurations': len(rows),
         'mean_bound_currency': rows[0]['mean_currency_upper_decimal'],
         'Bernstein_sample_range': [min(r['classical_Bernstein_samples'] for r in rows),
                                    max(r['classical_Bernstein_samples'] for r in rows)]}, indent=2))


if __name__ == '__main__':
    main()
