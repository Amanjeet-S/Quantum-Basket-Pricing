"""Price-independent higher-dimensional finite-law preparation.

All financial parameters and the actual rounded Gaussian quantile law are
unchanged. Separability avoids Cartesian tensor enumeration. Cached empirical
MGFs are enclosed by finite exponential series with an explicit remainder.
The method is established multinomial/Gaussian mathematics, not a novelty claim.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import copy, hashlib, json, time
from flint import arb, ctx
from continuous_quadrature import frac_arb, serialise_interval
from integral_value_oracle import model_parameters
from taylor_pricing import compositions
from boundary_grid_certificate import coefficient_certificate, endpoints, serial
from boundary_derivative_certificate import payoff_norms

ROOT = Path(__file__).resolve().parent
PRECISION = 512
SERIES_ORDER = 120


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def qnodes(L=6, b_y=10, q=28):
    """Return exactly rounded midpoint quantiles, with compressed provenance.

    Symmetry halves the costly inverse-error-function work. Every positive
    coordinate's enclosure lies inside one strict rounding bin, so ties cannot
    alter the mirrored negative coordinate. Existing small tables are audited
    independently below, rather than silently assumed identical.
    """
    L = F(L)
    if L <= 0 or type(b_y) is not int or not 1 <= b_y <= 20 or type(q) is not int or q < 2:
        raise ValueError('Positive L, b_y in 1 through 20 and q at least two required')
    scale = 1 << q
    if (L * scale).denominator != 1:
        raise ValueError('The retained offset must be exactly dyadic')
    ctx.prec = PRECISION
    start = time.perf_counter()
    N = 1 << b_y
    sqrt2 = arb(2).sqrt()
    lower = (1 + (-frac_arb(L) / sqrt2).erf()) / 2
    mass = 1 - 2 * lower
    integers = [0] * N
    maximum_error = F(0)
    last = F(0)
    least_margin = None
    for j in range(N // 2, N):
        probability = lower + frac_arb(F(2 * j + 1, 2 * N)) * mass
        y = (2 * probability - 1).erfinv() * sqrt2
        yl, yu = endpoints(y)
        centre = (yl + yu) / 2
        z = centre * scale + F(1, 2)
        integer = z.numerator // z.denominator
        low_bin, high_bin = F(2 * integer - 1, 2 * scale), F(2 * integer + 1, 2 * scale)
        if not low_bin < yl <= yu < high_bin:
            raise ArithmeticError('Quantile enclosure does not prove its unique rounding bin')
        coordinate = F(integer, scale)
        if not 0 <= coordinate <= L or coordinate < last:
            raise ArithmeticError('Rounded positive quantiles failed retained range or monotonicity')
        last = coordinate
        maximum_error = max(maximum_error, abs(coordinate - yl), abs(coordinate - yu))
        margin = min(yl - low_bin, high_bin - yu)
        least_margin = margin if least_margin is None else min(least_margin, margin)
        integers[j], integers[N - 1 - j] = integer, -integer
    ys = tuple(F(i, scale) for i in integers)
    digest = hashlib.sha256(json.dumps(integers, separators=(',', ':')).encode()).hexdigest()
    return ys, {'schema': 'qbp-compressed-rounded-quantiles-v1', 'L': str(L), 'b_y': b_y, 'q': q,
        'node_count': N, 'inverse_error_function_evaluations': N // 2,
        'maximum_node_error_upper': str(maximum_error), 'minimum_strict_rounding_bin_margin_lower': str(least_margin),
        'coordinate_min': str(ys[0]), 'coordinate_max': str(ys[-1]), 'symmetric_actual_nodes': True,
        'actual_signed_coordinate_integer_sha256': digest, 'mass_interval': serial(mass),
        'precision_bits': PRECISION, 'generation_seconds_diagnostic': time.perf_counter() - start,
        'scope': 'Uniform labels on the exactly rounded truncated-normal midpoint quantiles. Symmetry is a proved construction, not a distribution approximation.'}


class EmpiricalMGF:
    """Enclose sum_i weight_i exp(v*y_i), including signed weights.

    This general helper supports transport-measure MGFs. A caller passing Arb
    weights retains their enclosures. No assertion of positivity or symmetry
    is needed; the absolute remainder is charged using their total variation.
    """
    def __init__(self, coordinates, weights=None, order=SERIES_ORDER):
        ctx.prec = PRECISION
        self.coordinates = tuple(map(F, coordinates))
        if not self.coordinates or type(order) is not int or order < 2:
            raise ValueError('Nonempty coordinates and a positive series order required')
        N = len(self.coordinates)
        weights = tuple(weights) if weights is not None else (F(1, N),) * N
        if len(weights) != N:
            raise ValueError('Coordinate and weight counts differ')
        self.order = order
        self.radius = max(map(abs, self.coordinates))
        self.total_variation = arb(0)
        self.coefficients = [arb(0) for _ in range(order + 1)]
        for y, weight in zip(self.coordinates, weights):
            value = weight if isinstance(weight, arb) else frac_arb(F(weight))
            self.total_variation += abs(value)
            self.coefficients[0] += value
            ya = frac_arb(y)
            for r in range(1, order + 1):
                value *= ya / r
                self.coefficients[r] += value
        self.cache = {}
        self.query_count = 0
        self.maximum_remainder_upper = F(0)

    def __call__(self, rate):
        rate = F(rate)
        self.query_count += 1
        if rate not in self.cache:
            x = frac_arb(rate)
            value = arb(0)
            for coefficient in reversed(self.coefficients):
                value = value * x + coefficient
            argument = frac_arb(abs(rate) * self.radius)
            remainder = self.total_variation * argument.exp() * argument ** (self.order + 1) / factorial(self.order + 1)
            rup = endpoints(remainder)[1]
            value += arb(0, frac_arb(rup).upper())
            self.maximum_remainder_upper = max(self.maximum_remainder_upper, rup)
            self.cache[rate] = value
        return self.cache[rate]


class FiniteMGF:
    """Certified positive symmetric finite-law MGF without tensor enumeration."""
    def __init__(self, L=6, b_y=10, q=28, order=SERIES_ORDER):
        ctx.prec = PRECISION
        if type(order) is not int or order < 2 or order % 2:
            raise ValueError('An even exponential-series order is required')
        start = time.perf_counter()
        self.coordinates, self.node_certificate = qnodes(L, b_y, q)
        self.order = order
        self.radius = max(map(abs, self.coordinates))
        self.coefficients = [arb(0) for _ in range(order // 2 + 1)]
        # Each y and -y pair has identical even powers. The normalisation is
        # therefore 2/N and odd empirical moments vanish exactly.
        positive = self.coordinates[len(self.coordinates) // 2:]
        normalisation = F(1, len(positive))
        for y in positive:
            value = frac_arb(normalisation)
            self.coefficients[0] += value
            square = frac_arb(y * y)
            for r in range(1, order // 2 + 1):
                value *= square / ((2 * r - 1) * (2 * r))
                self.coefficients[r] += value
        self.cache = {F(0): arb(1)}
        self.query_count = 0
        self.maximum_remainder_upper = F(0)
        self.elapsed = time.perf_counter() - start
        self.moment_geometry_cache = {}
        self.domain_geometry_cache = {}

    def __call__(self, rate):
        rate = abs(F(rate))
        self.query_count += 1
        if rate not in self.cache:
            square = frac_arb(rate * rate)
            value = arb(0)
            for coefficient in reversed(self.coefficients):
                value = value * square + coefficient
            argument = frac_arb(rate * self.radius)
            # A complete even polynomial through degree order has the generic
            # exponential remainder bound below. It is deliberately absolute;
            # positivity also proves that the omitted cosh tail is nonnegative.
            remainder = argument.exp() * argument ** (self.order + 1) / factorial(self.order + 1)
            rup = endpoints(remainder)[1]
            lo, hi = endpoints(value)
            self.cache[rate] = frac_arb((lo + hi + rup) / 2) + arb(0, frac_arb((hi + rup - lo) / 2).upper())
            self.maximum_remainder_upper = max(self.maximum_remainder_upper, rup)
        return self.cache[rate]

    def profile(self):
        return {'node_certificate': self.node_certificate, 'series_order': self.order,
            'empirical_even_power_multiplications': len(self.coordinates) // 2 * (self.order // 2),
            'distinct_MGF_rates': len(self.cache), 'MGF_requests': self.query_count,
            'polynomial_coefficient_evaluations': (len(self.cache) - 1) * (self.order // 2 + 1),
            'remainder_exponential_evaluations': len(self.cache) - 1,
            'distinct_moment_geometries': len(self.moment_geometry_cache),
            'distinct_finite_domain_geometries': len(self.domain_geometry_cache),
            'maximum_MGF_remainder_upper': str(self.maximum_remainder_upper),
            'generation_seconds_diagnostic': self.elapsed,
            'scope': 'Paid deterministic shared preparation. Rates are exact rational model combinations; the finite law is unchanged. No financial payoff or basket-price answer is an input.'}


@lru_cache(maxsize=16)
def get_mgf(L=6, b_y=10, q=28, order=SERIES_ORDER):
    return FiniteMGF(F(L), b_y, q, order)


def product_grid_moments(model, L=6, b_y=10, q=28, max_order=6, mgf=None):
    """Multinomial factorisation of exact product-law T and U moments."""
    if type(max_order) is not int or not 1 <= max_order <= 6:
        raise ValueError('Moment orders one through six are supported')
    ctx.prec = PRECISION
    start = time.perf_counter()
    w, d, s, mg, rg, lam = model_parameters(model)
    mgf = mgf or get_mgf(L, b_y, q)
    if (F(mgf.node_certificate['L']), mgf.node_certificate['b_y'], mgf.node_certificate['q']) != (F(L), b_y, q):
        raise ValueError('The finite MGF belongs to a different actual input law')
    geometry_key = (w, d, s, rg, mg + lam * lam / 2, max_order)
    if geometry_key in mgf.moment_geometry_cache:
        record = copy.deepcopy(mgf.moment_geometry_cache[geometry_key])
        record['model_id'] = model['model_id']
        record['reused_exact_geometry'] = True
        record['additional_multinomial_exponential_evaluations'] = 0
        record['generation_seconds_diagnostic'] = time.perf_counter() - start
        return record
    raw, term_counts, exponentials = [], [], 0
    for r in range(max_order + 1):
        total = arb(0)
        terms = 0
        for alpha in compositions(r, len(w)):
            coefficient = F(factorial(r), prod(factorial(a) for a in alpha)) * prod(wi ** a for wi, a in zip(w, alpha))
            vector = tuple(sum(alpha[i] * s[i][j] for i in range(len(w))) for j in range(len(rg)))
            exponent = sum(a * di for a, di in zip(alpha, d))
            total += frac_arb(coefficient) * frac_arb(exponent).exp() * prod(mgf(v) for v in vector)
            terms += 1
        raw.append(total)
        term_counts.append(terms)
        exponentials += terms
    centred = {}
    for r in range(1, max_order + 1):
        value = sum(((-1) ** (r - j) * comb(r, j) * raw[j] for j in range(r + 1)), arb(0))
        lo, hi = endpoints(value)
        if hi < 0:
            raise ArithmeticError('A nonnegative U moment has a negative certified upper bound')
        centred[str(r)] = serialise_interval(max(F(0), lo), hi)
    weighted_T = sum((frac_arb(wi) * frac_arb(di).exp() * prod(mgf(x + y) for x, y in zip(rg, si)) for wi, di, si in zip(w, d, s)), arb(0))
    level = prod(mgf(x) for x in rg)
    gap = 100 * frac_arb(mg + lam * lam / 2).exp() * (weighted_T - level)
    glo, ghi = endpoints(gap)
    if ghi < 0:
        raise ArithmeticError('The nonnegative correlated underlying gap failed certification')
    record = {'schema': 'qbp-higherdim-factorised-finite-moments-v1', 'status': 'passed',
        'model_id': model['model_id'], 'assets': model['assets'], 'outer_dimension': len(rg),
        'L': str(F(L)), 'b_y': b_y, 'q': q, 'distribution': 'quantile',
        'centred_U_moments': centred, 'raw_T_moments': [serial(x) for x in raw],
        'correlated_underlying_gap_interval': serialise_interval(max(F(0), glo), ghi),
        'multinomial_terms_by_order': term_counts,
        'multinomial_exponential_evaluations': exponentials,
        'gap_exponential_evaluations': len(w) + 1,
        'shared_geometry_source_model_id': model['model_id'], 'reused_exact_geometry': False,
        'additional_multinomial_exponential_evaluations': exponentials + len(w) + 1,
        'avoided_tensor_labels': (1 << b_y) ** len(rg),
        'generation_seconds_diagnostic': time.perf_counter() - start,
        'scope': 'Ideal dispersion and correlated underlying gap on the exact independent rounded-node law. No conditional payoff or basket price is evaluated.'}
    mgf.moment_geometry_cache[geometry_key] = copy.deepcopy(record)
    return record


def finite_domain(model, L=6, b_y=10, q=28, mgf=None):
    ctx.prec = PRECISION
    mgf = mgf or get_mgf(L, b_y, q)
    node_error = F(mgf.node_certificate['maximum_node_error_upper'])
    actual_lo, actual_hi = min(mgf.coordinates), max(mgf.coordinates)
    lo, hi = actual_lo - node_error, actual_hi + node_error
    w, d, s, mg, rg, lam = model_parameters(model)
    maximum, vertices = F(1), []
    geometry_key = (w, d, s, lo, hi)
    reused = geometry_key in mgf.domain_geometry_cache
    if reused:
        maximum, vertices, source_model = copy.deepcopy(mgf.domain_geometry_cache[geometry_key])
    else:
        for signs in product((0, 1), repeat=len(rg)):
            y = tuple(hi if sign else lo for sign in signs)
            value = sum((frac_arb(wi) * frac_arb(di + sum(v * z for v, z in zip(si, y))).exp() for wi, di, si in zip(w, d, s)), arb(0))
            bound = endpoints(value)[1]
            maximum = max(maximum, bound)
            vertices.append({'signs': signs, 'T_upper': str(bound)})
        source_model = model['model_id']
        mgf.domain_geometry_cache[geometry_key] = copy.deepcopy((maximum, vertices, source_model))
    mlo = mg + sum(min(v * lo, v * hi) for v in rg)
    mhi = mg + sum(max(v * lo, v * hi) for v in rg)
    return {'model_id': model['model_id'], 'L_transfer': str(F(L)), 'b_y': b_y, 'q': q,
        'coordinate_bounds': [str(lo), str(hi)], 'mG_bounds': [str(mlo), str(mhi)],
        'actual_coordinate_bounds': [str(actual_lo), str(actual_hi)],
        'ideal_quantile_enclosure_padding': str(node_error),
        'T_star': str(maximum), 'U_star': str(maximum - 1), 'vertices_evaluated': len(vertices),
        'vertex_upper_enclosures': vertices, 'covers_actual_dyadic_nodes': True,
        'shared_geometry_source_model_id': source_model, 'reused_exact_geometry': reused,
        'additional_vertex_exponential_evaluations': 0 if reused else len(w) * len(vertices),
        'covers_ideal_quantile_enclosures': True,
        'scope': 'Finite-law oracle normalisation by coordinatewise convexity. Coordinate extrema are padded by the certified maximum rounding displacement, covering both ideal quantile enclosures and actual dyadics. Continuous/grid transfer has its own L box and nodal displacement proof.'}


def preparation(model, L=6, b_y=10, q=28, max_order=6):
    mgf = get_mgf(L, b_y, q)
    moments = product_grid_moments(model, L, b_y, q, max_order, mgf)
    domain = finite_domain(model, L, b_y, q, mgf)
    cases = []
    for cap in (1, 10, 25, None):
        payoff = 'vanilla' if cap is None else 'capped'
        cases.append({'model_id': model['model_id'], 'payoff': payoff, 'cap': cap,
            'coefficients': coefficient_certificate(model, cap, payoff, 1),
            'shared_derivative_norms': payoff_norms(model['lambda'], cap),
            'scope': 'Both quantum and classical competitors receive identical coefficients, finite moments, correlated gap, finite domain and derivative norms.'})
    return {'model_id': model['model_id'], 'moments': moments, 'domain': domain, 'cases': cases}


def run(b_y=10, q=28, output='higherdim_preparation.json'):
    start = time.perf_counter()
    protocol = ROOT / 'results/generalisation_protocol.json'
    dependencies = ['continuous_quadrature.py', 'integral_value_oracle.py', 'taylor_pricing.py',
        'boundary_grid_certificate.py', 'boundary_derivative_certificate.py',
        'results/boundary_derivative_norms.json', 'results/generalisation_protocol.json']
    hashes = {name: sha(ROOT / name) for name in dependencies}
    family = json.loads(protocol.read_text())['models']
    records = []
    for model in family:
        if model['outer_dimension'] not in (3, 7):
            continue
        record = preparation(model, b_y=b_y, q=q)
        records.append(record)
        print(model['model_id'], 'U2', float(F(record['moments']['centred_U_moments']['2']['upper_exact_fraction'])),
            'U6', float(F(record['moments']['centred_U_moments']['6']['upper_exact_fraction'])),
            'Ustar', float(F(record['domain']['U_star'])), flush=True)
    if hashes != {name: sha(ROOT / name) for name in dependencies}:
        raise RuntimeError('A frozen preparation input changed during generation')
    result = {'schema': 'qbp-higherdim-shared-preparation-family-v1', 'status': 'passed',
        'records': records, 'model_count': len(records), 'payoff_case_count': 4 * len(records),
        'L': '6', 'b_y': b_y, 'q': q, 'shared_MGF_profile': get_mgf(6, b_y, q).profile(),
        'input_sha256': hashes, 'source_sha256': sha(__file__),
        'elapsed_seconds_diagnostic': time.perf_counter() - start,
        'scope': 'Twenty-four frozen p3/p7 financial cases. This file supplies paid shared preparation, not a complete monetary-transfer certificate or quantum advantage claim.'}
    (ROOT / 'results' / output).write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--b-y', type=int, default=10)
    parser.add_argument('--q', type=int, default=28)
    parser.add_argument('--output', default='higherdim_preparation.json')
    args = parser.parse_args()
    run(args.b_y, args.q, args.output)
