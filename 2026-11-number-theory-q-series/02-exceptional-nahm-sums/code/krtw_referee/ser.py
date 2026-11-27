# Independent referee series library: elements of Z[[q]][[x]] truncated at x^X, q^N.
# Represented as list (index = x-degree) of flint.fmpz_poly in q, truncated at degree N.
from flint import fmpz_poly
from fractions import Fraction
import itertools

class S:
    def __init__(self, X, N, c=None):
        self.X, self.N = X, N
        self.c = c if c is not None else [fmpz_poly(0) for _ in range(X + 1)]
    def __add__(self, o):
        return S(self.X, self.N, [a + b for a, b in zip(self.c, o.c)])
    def __sub__(self, o):
        return S(self.X, self.N, [a - b for a, b in zip(self.c, o.c)])
    def __mul__(self, o):
        X, N = self.X, self.N
        r = [fmpz_poly(0) for _ in range(X + 1)]
        for i, a in enumerate(self.c):
            if a.is_zero(): continue
            for j in range(X + 1 - i):
                b = o.c[j]
                if b.is_zero(): continue
                r[i + j] += a.mul_low(b, N)
        return S(X, N, r)
    def __eq__(self, o):
        return all(a == b for a, b in zip(self.c, o.c))
    def iszero(self):
        return all(a.is_zero() for a in self.c)
    def inv(self):
        # inverse, requires constant term (x^0) invertible power series in q with c0=+-1
        X, N = self.X, self.N
        a0 = self.c[0]
        assert abs(int(a0[0])) == 1
        # invert a0 mod q^N by Newton/naive
        inv0 = series_inv(a0, N)
        r = [fmpz_poly(0) for _ in range(X + 1)]
        r[0] = inv0
        for d in range(1, X + 1):
            acc = fmpz_poly(0)
            for i in range(1, d + 1):
                acc += self.c[i].mul_low(r[d - i], N)
            r[d] = -(acc.mul_low(inv0, N))
        return S(X, N, r)
    def maxcoef(self):
        return max((max((abs(int(v)) for v in p.coeffs()), default=0) for p in self.c), default=0)

def series_inv(p, N):
    c = [int(v) for v in p.coeffs()] + [0] * N
    c = c[:N]
    a0 = c[0]
    r = [0] * N
    r[0] = a0  # a0 = +-1, 1/a0 = a0
    for n in range(1, N):
        s = 0
        for i in range(1, n + 1):
            s += c[i] * r[n - i]
        r[n] = -s * a0
    return fmpz_poly(r)

def qpoch_inv_table(d, nmax, N):
    """1/(q^d;q^d)_n mod q^N for n=0..nmax, as fmpz_poly"""
    out = [fmpz_poly([1])]
    cur = fmpz_poly([1])
    for n in range(1, nmax + 1):
        # multiply by 1/(1-q^{dn}) = sum q^{dnj}
        e = d * n
        geo = fmpz_poly([1 if (i % e == 0) else 0 for i in range(N)])
        cur = cur.mul_low(geo, N)
        out.append(cur)
    return out

def monomial(X, N, xd, qe, coef=1):
    s = S(X, N)
    if xd <= X and 0 <= qe < N:
        s.c[xd] = fmpz_poly([0] * qe + [coef])
    elif qe < 0:
        raise ValueError("negative q exponent")
    return s

def from_qpoly_in_x(X, N, terms):
    """terms: dict (xdeg, qexp)->coef"""
    s = S(X, N)
    for (a, b), v in terms.items():
        if a <= X and b < N:
            assert b >= 0
            s.c[a] += fmpz_poly([0] * b + [v])
    return s

def poch(X, N, factors):
    """product over factors of (1 - c x^a q^b), each factor given as (c,a,b); finite product"""
    r = monomial(X, N, 0, 0)
    for (c, a, b) in factors:
        f = monomial(X, N, 0, 0) - monomial(X, N, a, b, c)
        r = r * f
    return r

def poch_inf(X, N, c, a, b, step):
    """(c x^a q^b ; q^step)_inf truncated (b>0 required unless a>0)"""
    facs = []
    j = 0
    while True:
        e = b + step * j
        if a == 0 and e >= N: break
        if a > 0 and (e >= N and True):
            break
        facs.append((c, a, e))
        j += 1
        if j > 10 * N: break
    return poch(X, N, facs)
