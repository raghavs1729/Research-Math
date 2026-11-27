"""Referee 2 -- own exact truncated-series arithmetic (nothing imported from code/, referee/, rr-identities/).

Univariate q-series: flint.fmpz_poly truncated to length L (i.e. mod q^L).
Bivariate series in Z[[x,q]] truncated mod (x^(D+1), q^L): class BS, a list of D+1 fmpz_poly
(entry d = coefficient of x^d).
"""
from flint import fmpz_poly


# ---------------------------------------------------------------- univariate
def mono(e, L, c=1):
    """c*q^e mod q^L (e >= 0)."""
    assert e >= 0
    if e >= L:
        return fmpz_poly()
    return fmpz_poly([0] * e + [c])


def trunc(p, L):
    return p.truncate(L)


def mul(a, b, L):
    return a.mul_low(b, L)


def inv(p, L):
    """1/p mod q^L, p(0) = +-1."""
    c0 = int(p[0])
    assert c0 in (1, -1)
    out = [0] * L
    pc = [int(p[i]) for i in range(L)]
    for n in range(L):
        s = 1 if n == 0 else 0
        for i in range(1, n + 1):
            if pc[i]:
                s -= pc[i] * out[n - i]
        out[n] = s * c0
    return fmpz_poly(out)


def poch(a, step, n, L, sign=1):
    """prod_{j=0}^{n-1} (1 - sign*q^(a+j*step)) mod q^L; n may be 'inf'."""
    r = fmpz_poly([1])
    j = 0
    while (n == 'inf' or j < n):
        e = a + j * step
        if e >= L:
            if n == 'inf':
                break
            j += 1
            continue
        assert e >= 0
        r = r.mul_low(fmpz_poly([1]) - mono(e, L, sign), L)
        j += 1
    return r


def qpoch(n, L):
    """(q;q)_n, 0 for n<0 is NOT returned: caller must treat 1/(q)_n = 0 for n < 0."""
    return poch(1, 1, n, L)


def Qpoch(n, L):
    """(q^2;q^2)_n."""
    return poch(2, 2, n, L)


class InvCache:
    """cached 1/(q;q)_n and 1/(Q;Q)_n mod q^L (0 for n<0)."""

    def __init__(self, L, nmax):
        self.L = L
        self.iq = [inv(qpoch(n, L), L) for n in range(nmax + 1)]
        self.iQ = [inv(Qpoch(n, L), L) for n in range(nmax + 1)]

    def q(self, n):
        return self.iq[n] if n >= 0 else fmpz_poly()

    def Q(self, n):
        return self.iQ[n] if n >= 0 else fmpz_poly()


# ---------------------------------------------------------------- bivariate
class BS:
    __slots__ = ('D', 'L', 'c')

    def __init__(self, D, L, c=None):
        self.D, self.L = D, L
        self.c = c if c is not None else [fmpz_poly() for _ in range(D + 1)]

    @staticmethod
    def const(D, L, p):
        r = BS(D, L)
        r.c[0] = trunc(fmpz_poly(p) if not isinstance(p, fmpz_poly) else p, L)
        return r

    @staticmethod
    def from_terms(D, L, terms):
        """terms: iterable of (xexp, qexp, coeff)."""
        r = BS(D, L)
        for (i, j, c) in terms:
            if i <= D and j < L:
                assert i >= 0 and j >= 0
                r.c[i] = r.c[i] + mono(j, L, c)
        return r

    @staticmethod
    def xseries(D, L, coeffs):
        """sum_d coeffs(d) x^d (q-free)."""
        r = BS(D, L)
        for d in range(D + 1):
            r.c[d] = fmpz_poly([coeffs(d)])
        return r

    def copy(self):
        return BS(self.D, self.L, list(self.c))

    def __add__(self, o):
        return BS(self.D, self.L, [a + b for a, b in zip(self.c, o.c)])

    def __sub__(self, o):
        return BS(self.D, self.L, [a - b for a, b in zip(self.c, o.c)])

    def __neg__(self):
        return BS(self.D, self.L, [-a for a in self.c])

    def scale(self, k):
        return BS(self.D, self.L, [a * k for a in self.c])

    def __mul__(self, o):
        D, L = self.D, self.L
        out = [fmpz_poly() for _ in range(D + 1)]
        nzs = [i for i in range(D + 1) if not self.c[i].is_zero()]
        nzo = [j for j in range(D + 1) if not o.c[j].is_zero()]
        for i in nzs:
            for j in nzo:
                if i + j > D:
                    break
                out[i + j] = out[i + j] + self.c[i].mul_low(o.c[j], L)
        return BS(D, L, out)

    def mulq(self, p):
        return BS(self.D, self.L, [a.mul_low(p, self.L) if not a.is_zero() else a for a in self.c])

    def shift(self, a, b):
        """multiply by x^a q^b, a, b >= 0."""
        assert a >= 0 and b >= 0
        D, L = self.D, self.L
        out = [fmpz_poly() for _ in range(D + 1)]
        for i in range(D + 1 - a):
            if not self.c[i].is_zero():
                out[i + a] = self.c[i].left_shift(b).truncate(L)
        return BS(D, L, out)

    def shift_down_x(self, a):
        """divide by x^a; asserts divisibility (low coefficients zero)."""
        for i in range(a):
            assert self.c[i].is_zero(), "not divisible by x^%d" % a
        out = self.c[a:] + [fmpz_poly() for _ in range(a)]
        return BS(self.D, self.L, out)

    def inverse(self):
        """1/self, requires c[0] invertible (constant term +-1)."""
        D, L = self.D, self.L
        i0 = inv(self.c[0], L)
        out = [fmpz_poly() for _ in range(D + 1)]
        out[0] = i0
        for d in range(1, D + 1):
            s = fmpz_poly()
            for i in range(1, d + 1):
                if not self.c[i].is_zero() and not out[d - i].is_zero():
                    s = s + self.c[i].mul_low(out[d - i], L)
            out[d] = (-s).mul_low(i0, L)
        return BS(D, L, out)

    def is_zero(self):
        return all(a.is_zero() for a in self.c)

    def __eq__(self, o):
        return (self - o).is_zero()

    def first_diff(self, o):
        d = self - o
        for i in range(self.D + 1):
            if not d.c[i].is_zero():
                cs = d.c[i].coeffs()
                j = next(t for t, v in enumerate(cs) if v != 0)
                return (i, j, int(cs[j]))
        return None

    def at_x1(self):
        """sum of x-coefficients (valid only if truncation in x is irrelevant, i.e. D large)."""
        s = fmpz_poly()
        for a in self.c:
            s = s + a
        return s.truncate(self.L)


def bpoch_x(D, L, xpow, qa, qstep, n, sign=1):
    """prod_{j=0}^{n-1} (1 - sign * x^xpow q^(qa+j*qstep)) as BS; n may be 'inf'."""
    r = BS.const(D, L, fmpz_poly([1]))
    j = 0
    while n == 'inf' or j < n:
        e = qa + j * qstep
        if n == 'inf' and (e >= L or xpow > D):
            break
        f = BS.from_terms(D, L, [(0, 0, 1), (xpow, e, -sign)])
        r = r * f
        j += 1
    return r
