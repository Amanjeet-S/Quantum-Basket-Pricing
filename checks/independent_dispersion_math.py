"""Independent high-precision algebra for the dispersion certificate audit.

This module imports no production helper. Signed U powers are formed by
polynomial multiplication. The Peano MGF uses a cancellation-safe factorial
moment series, not the production quotient of three transforms.
"""
from fractions import Fraction as F
from math import factorial
from flint import arb, ctx

PRECISION = 1024
SERIES = 224


def ball(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def endpoint(x, top):
    y = x.upper() if top else x.lower()
    m, e = map(int, y.man_exp())
    return F(m * (1 << e)) if e >= 0 else F(m, 1 << -e)


def interval(lo, hi):
    lo, hi = F(lo), F(hi)
    return ball((lo + hi) / 2) + arb(0, ball((hi - lo) / 2).upper())


def serial(x):
    return {"lower_exact_fraction": str(endpoint(x, False)),
            "upper_exact_fraction": str(endpoint(x, True))}


def contains_record(rec, value):
    return (F(rec["lower_exact_fraction"]) <= endpoint(value, False)
            <= endpoint(value, True) <= F(rec["upper_exact_fraction"]))


def positive(x):
    lo, hi = endpoint(x, False), endpoint(x, True)
    if hi < 0:
        raise ArithmeticError("Negative upper endpoint for nonnegative moment")
    return interval(max(F(0), lo), hi)


def phi(x):
    return (-x * x / 2).exp() / (2 * arb.pi()).sqrt()


def Phi(x):
    return (1 + (x / arb(2).sqrt()).erf()) / 2


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def add(p, q):
    out = dict(p)
    for a, c in q.items():
        out[a] = out.get(a, F(0)) + c
    return {a: c for a, c in out.items() if c}


def mul(p, q):
    out = {}
    for a, ca in p.items():
        for b, cb in q.items():
            alpha = tuple(x + y for x, y in zip(a, b))
            out[alpha] = out.get(alpha, F(0)) + ca * cb
    return {a: c for a, c in out.items() if c}


def base_polynomials(w, s, j=None, degree=4):
    zero = (0,) * len(w)
    T = {tuple(int(i == k) for k in range(len(w))): wi for i, wi in enumerate(w)}
    U = add(T, {zero: F(-1)})
    powers = [{zero: F(1)}]
    for _ in range(degree):
        powers.append(mul(powers[-1], U))
    if j is None:
        return powers
    u = {a: c * s[i][j] for i, (a, c) in enumerate(T.items())}
    v = {a: c * s[i][j] ** 2 for i, (a, c) in enumerate(T.items())}
    u2 = mul(u, u)
    return {"U2": powers[2], "J2": u2, "U_JJ": mul(U, v),
            "U_J2": mul(U, u2), "U2_JJ": mul(powers[2], v),
            "U2_J2": mul(powers[2], u2), "U4": powers[4]}


class IndependentTransport:
    def __init__(self, actual):
        ctx.prec = PRECISION
        self.L = F(6)
        self.N = 1024
        assert len(actual) == self.N
        self.mass = Phi(ball(6)) - Phi(ball(-6))
        low = Phi(ball(-6))
        sqrt2 = arb(2).sqrt()
        boundaries = [ball(-6)]
        boundaries.extend(sqrt2 * (2 * (low + self.mass * j / self.N) - 1).erfinv()
                          for j in range(1, self.N))
        boundaries.append(ball(6))
        self.nodes = []
        first = []
        weights = []
        second = []
        self.input_checks = 0
        for j, saved in enumerate(actual):
            a, b = boundaries[j:j + 2]
            midpoint = sqrt2 * (2 * (low + self.mass * ball(F(2 * j + 1, 2 * self.N))) - 1).erfinv()
            scale = 1 << 28
            lo = endpoint(midpoint, False) * scale + F(1, 2)
            hi = endpoint(midpoint, True) * scale + F(1, 2)
            assert lo.numerator // lo.denominator == hi.numerator // hi.denominator
            node = F(lo.numerator // lo.denominator, scale)
            assert node == F(saved["actual_node"])
            assert node * scale == saved["actual_node_integer"]
            assert endpoint(a, True) <= node <= endpoint(b, False)
            y = ball(node)
            cellmass = self.mass / self.N
            m1 = phi(a) - phi(b) - y * cellmass
            m2 = (1 + y * y) * cellmass + (a - 2 * y) * phi(a) - (b - 2 * y) * phi(b)
            atom = F(saved["rho1_positive_atom_weight"])
            assert atom >= max(abs(endpoint(m1, False)), abs(endpoint(m1, True)))
            for key, value in [("left", a), ("right", b), ("mass", cellmass),
                               ("first_moment", m1), ("second_moment", m2),
                               ("left_cdf", low + self.mass * j / self.N),
                               ("right_cdf", low + self.mass * (j + 1) / self.N)]:
                assert contains_record(saved[key], value), (j, key)
                self.input_checks += 1
            self.nodes.append(node)
            first.append(m1)
            weights.append(atom)
            second.append(m2)
        J = SERIES
        rounded = [arb(0)] * (J + 3)
        signed = [arb(0)] * (J + 2)
        rho1 = [arb(0)] * (J + 1)
        for node, mu, weight in zip(self.nodes, first, weights):
            power = arb(1)
            for r in range(J + 3):
                rounded[r] += self.mass / self.N * power
                if r <= J + 1:
                    signed[r] += mu * power
                if r <= J:
                    rho1[r] += ball(weight) * power
                power *= ball(node) / (r + 1)
        normal = [self.mass, arb(0)]
        for r in range(2, J + 3):
            normal.append(arb(0) if r % 2 else (r - 1) * normal[r - 2] - 2 * ball(6) ** (r - 1) * phi(ball(6)))
        normal = [x / factorial(j) for j, x in enumerate(normal)]
        self.rho2 = [normal[r + 2] - rounded[r + 2] - signed[r + 1] for r in range(J + 1)]
        self.rho2mass = positive(sum(second, arb(0)) / 2)
        assert (self.rho2[0] - self.rho2mass).contains(0)
        self.coefficients = {"rounded": rounded[:J + 1], "rho1": rho1, "rho2": self.rho2}
        self.variation = {"rounded": self.mass, "rho1": sum((ball(w) for w in weights), arb(0)), "rho2": self.rho2mass}
        self.cache = {}

    def mgf(self, kind, rate):
        rate = F(rate)
        key = (kind, rate)
        if key not in self.cache:
            if kind == "continuous":
                v = ball(rate)
                self.cache[key] = (v * v / 2).exp() * (Phi(ball(6) - v) - Phi(ball(-6) - v))
            else:
                v = ball(rate)
                val = arb(0)
                for a in reversed(self.coefficients[kind]):
                    val = val * v + a
                x = ball(abs(rate) * self.L)
                error = self.variation[kind] * x.exp() * x ** (SERIES + 1) / factorial(SERIES + 1)
                radius = endpoint(error, True)
                self.cache[key] = positive(interval(endpoint(val, False) - radius,
                                                    endpoint(val, True) + radius))
        return self.cache[key]


class IndependentGeometry:
    def __init__(self, w, d, s, rg, degree=4):
        self.w, self.d, self.s, self.rg = w, d, s, rg
        self.powers = base_polynomials(w, s, degree=degree)
        self.basis = {}
        self.hybrid_cache = {}
        self.full_moments = [self.full(poly) for poly in self.powers]

    def basis_value(self, alpha):
        if alpha not in self.basis:
            v = tuple(sum(alpha[i] * self.s[i][j] for i in range(len(alpha)))
                      for j in range(len(self.rg)))
            b = dot(alpha, self.d)
            self.basis[alpha] = (v, b, ball(b + dot(v, v) / 2).exp(), dot(self.rg, v))
        return self.basis[alpha]

    def full(self, poly):
        return positive(sum((ball(c) * self.basis_value(a)[2] for a, c in poly.items()), arb(0)))

    def outside(self, poly):
        result = arb(0)
        for a, coefficient in poly.items():
            v, b, full, shift = self.basis_value(a)
            retained = full
            for x in v:
                retained *= Phi(ball(6 - x)) - Phi(ball(-6 - x))
            result += ball(coefficient) * (full - retained)
        return positive(result)

    def hybrid(self, poly, j, kind, transport):
        total = arb(0)
        for a, coefficient in poly.items():
            key = (a, j, kind)
            if key not in self.hybrid_cache:
                v, b, full, shift = self.basis_value(a)
                value = ball(b).exp()
                for l, rate in enumerate(v):
                    scalar = ("rho1" if kind == 1 else "rho2") if l == j else "continuous" if l < j else "rounded"
                    value *= transport.mgf(scalar, rate)
                self.hybrid_cache[key] = value
            total += ball(coefficient) * self.hybrid_cache[key]
        return positive(total)


def call_falling(m, sigma, strike, order):
    strike = ball(strike)
    z = (m - (strike / 100).log()) / sigma
    growth = 100 * (m + sigma * sigma / 2).exp()
    if order == 0:
        return growth * Phi(z + sigma) - strike * Phi(z)
    if order == 1:
        return growth * Phi(z + sigma)
    # Independent derivative-operator expansion, as opposed to the
    # producer's recurrence for a polynomial in z at each smoothing scale.
    polynomial = [1]
    for r in range(2, order):
        new = [0] * (len(polynomial) + 1)
        for j, a in enumerate(polynomial):
            new[j] -= r * a
            new[j + 1] += a
        polynomial = new
    hermite = [arb(1), z]
    for j in range(2, len(polynomial)):
        hermite.append(z * hermite[-1] - (j - 1) * hermite[-2])
    value = sum((a * (-1) ** j * hermite[j] / sigma ** j
                 for j, a in enumerate(polynomial)), arb(0))
    return strike * phi(z) / sigma * value


def coefficients(geometry, mg, lam, cap, maximum):
    sigma = ball(lam * lam + dot(geometry.rg, geometry.rg)).sqrt()
    def kernel(m, order):
        val = call_falling(m, sigma, 100, order)
        if cap is not None:
            val -= call_falling(m, sigma, 100 + cap, order)
        return val
    out = [kernel(ball(mg), 0)]
    for k in range(1, maximum + 1):
        val = arb(0)
        for a, coefficient in geometry.powers[k].items():
            v, b, factor, shift = geometry.basis_value(a)
            val += ball(coefficient) * factor * kernel(ball(mg + shift), k)
        out.append(val / factorial(k))
    return out


def envelope_components(norm, r, moments, kind):
    M2, M3, M4 = (F(norm["C_r_upper"][str(k)]) for k in (2, 3, 4))
    if norm["payoff"] == "vanilla":
        D1 = M2 * F(norm["scalar_boundary_density_first_upper"])
        D2 = M2 * F(norm["scalar_boundary_density_second_upper"])
    else:
        D1 = 2 * M2 * F(norm["scalar_boundary_h_first_upper"])
        D2 = 2 * M2 * F(norm["scalar_boundary_h_second_upper"])
    D31 = min(M4 + 3 * M3, D2 + 2 * D1)
    r = abs(r)
    uj = (moments["U2"] * moments["J2"]).sqrt()
    u2j = min(endpoint((moments["U2"] * moments["U2_J2"]).sqrt(), True),
              endpoint((moments["U4"] * moments["J2"]).sqrt(), True))
    if kind == 1:
        values = {"U_Tj": ball(M2) * uj, "r_U2": ball(r * D1 / 2) * moments["U2"],
                  "U2_Tj": ball(M3 / 6) * ball(u2j)}
    else:
        values = {"Tj2_U_Tjj": ball(M2) * (moments["J2"] + moments["U_JJ"]),
                  "r_U_Tj": ball(2 * r * D1) * uj,
                  "U_Tj2": ball(F(2, 3) * M3) * moments["U_J2"],
                  "r2_U2": ball(r * r * D2 / 2) * moments["U2"],
                  "r_U2_Tj": ball(r * D31 / 3) * ball(u2j),
                  "U2_Tjj": ball(M3 / 6) * moments["U2_JJ"],
                  "U2_Tj2": ball(M4 / 12) * moments["U2_J2"]}
    return values, {"M2": M2, "M3": M3, "M4": M4, "C2_prime": D1,
                    "C2_second": D2, "C3_prime": D31}
