"""Exact truncated bivariate series  sum_{n>=0} x^n c_n(q),  c_n in Z[q]/(q^N)  (python-flint).

Used by the krtw_* scripts.  A series is a dict {x-degree: fmpz_poly}; all arithmetic is exact modulo q^N
(N = global truncation order set by setN).  x-degrees are never truncated (every series built here has
x-degree bounded in terms of the q-order, because each power of x carries a positive power of q).
"""
from flint import fmpz_poly

N = 30
_cache = {}

def setN(n):
    global N, _cache
    N = n
    _cache = {}

def tr(p):
    if p.degree() >= N:
        return fmpz_poly(p.coeffs()[:N])
    return p

def mono(e, c=1):
    """c*q^e as fmpz_poly (0 if e >= N)"""
    if e >= N:
        return fmpz_poly([0])
    return fmpz_poly([0] * e + [c])

class S:
    __slots__ = ("d",)
    def __init__(self, d=None):
        self.d = {} if d is None else {k: v for k, v in d.items() if v != 0}
    @staticmethod
    def one():
        return S({0: fmpz_poly([1])})
    @staticmethod
    def xq(a, e, c=1):
        """c * x^a q^e"""
        return S({a: mono(e, c)})
    def __add__(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d[k] + v if k in d else v
        return S(d)
    def __sub__(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d[k] - v if k in d else -v
        return S(d)
    def __neg__(self):
        return S({k: -v for k, v in self.d.items()})
    def __mul__(self, o):
        if isinstance(o, int):
            return S({k: v * o for k, v in self.d.items()})
        if isinstance(o, fmpz_poly):
            return S({k: tr(v * o) for k, v in self.d.items()})
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                p = tr(v1 * v2)
                if p != 0:
                    k = k1 + k2
                    d[k] = d[k] + p if k in d else p
        return S(d)
    __rmul__ = __mul__
    def shift(self, a, e):
        """multiply by x^a q^e"""
        if e >= N:
            return S()
        m = mono(e)
        return S({k + a: tr(v * m) for k, v in self.d.items()})
    def subx(self, c):
        """x -> x q^c  (c >= 0)"""
        d = {}
        for k, v in self.d.items():
            if k * c < N:
                d[k] = tr(v * mono(k * c))
        return S(d)
    def iszero(self):
        return all(v == 0 for v in self.d.values())
    def __eq__(self, o):
        return (self - o).iszero()
    def at1(self):
        """x = 1"""
        p = fmpz_poly([0])
        for v in self.d.values():
            p = p + v
        return p
    def maxx(self):
        return max(self.d) if self.d else -1
    def __repr__(self):
        return " + ".join("x^%d*(%s)" % (k, self.d[k]) for k in sorted(self.d)) or "0"

def qpoch(n, base=1):
    """(q^base; q^base)_n  as fmpz_poly mod q^N"""
    key = ("p", n, base, N)
    if key not in _cache:
        p = fmpz_poly([1])
        for i in range(1, n + 1):
            if base * i < N:
                p = tr(p * fmpz_poly([1] + [0] * (base * i - 1) + [-1]))
        _cache[key] = p
    return _cache[key]

def inv(p):
    """1/p mod q^N for p with constant term 1"""
    c = p.coeffs() + [0] * N
    assert c[0] == 1
    r = [0] * N
    r[0] = 1
    for n in range(1, N):
        s = 0
        for j in range(1, n + 1):
            if c[j]:
                s += c[j] * r[n - j]
        r[n] = -s
    return fmpz_poly(r)

def iqpoch(n, base=1):
    """1/(q^base;q^base)_n ; 0 for n<0"""
    if n < 0:
        return fmpz_poly([0])
    key = ("ip", n, base, N)
    if key not in _cache:
        _cache[key] = inv(qpoch(n, base))
    return _cache[key]

def qbin(n, k, base=1):
    if k < 0 or k > n:
        return fmpz_poly([0])
    key = ("b", n, k, base, N)
    if key not in _cache:
        _cache[key] = tr(tr(qpoch(n, base) * iqpoch(k, base)) * iqpoch(n - k, base))
    return _cache[key]

def xpoch_inf(sign, a, e, base=1):
    """prod_{i>=0} (1 + sign * x^a q^{e+base*i})  mod q^N   (a>=1, e>=1 or (a>=1,e>=0))"""
    r = S.one()
    i = 0
    while e + base * i < N:
        r = r * (S.one() + S.xq(a, e + base * i, sign))
        i += 1
    return r

def xpoch(sign, a, e, n, base=1):
    """prod_{i=0}^{n-1} (1 + sign * x^a q^{e+base*i})"""
    r = S.one()
    for i in range(n):
        if e + base * i >= N:
            break
        r = r * (S.one() + S.xq(a, e + base * i, sign))
    return r

def xpoch_inv(sign, a, e, n, base=1):
    """1/prod_{i=0}^{n-1} (1 + sign * x^a q^{e+base*i}),  e>=1"""
    r = S.one()
    for i in range(n):
        ee = e + base * i
        if ee >= N:
            break
        # 1/(1+s x^a q^ee) = sum_m (-s)^m x^{am} q^{ee m}
        g = S()
        m = 0
        while ee * m < N:
            g = g + S.xq(a * m, ee * m, (-sign) ** m)
            m += 1
        r = r * g
    return r
