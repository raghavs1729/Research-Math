"""Minimal exact rational functions in Q(x, q) on top of flint.fmpz_mpoly (num/den, gcd-reduced)."""
import flint

CTX = flint.fmpz_mpoly_ctx.get(('x', 'q'), 'lex')
X, Q = CTX.gens()
ONE = CTX.from_dict({(0, 0): 1})
ZERO = CTX.from_dict({})


class RF:
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None, reduce=True):
        if isinstance(n, int):
            n = ONE * n
        if d is None:
            d = ONE
        elif isinstance(d, int):
            d = ONE * d
        if reduce:
            if n.is_zero():
                d = ONE
            else:
                g = n.gcd(d)
                if not g.is_one():
                    n = n / g
                    d = d / g
            # normalise sign: leading coeff of d positive
            lc = d.leading_coefficient()
            if lc < 0:
                n, d = -n, -d
        self.n, self.d = n, d

    def __add__(self, o):
        if isinstance(o, int):
            o = RF(o)
        if self.n.is_zero():
            return o
        if o.n.is_zero():
            return self
        g = self.d.gcd(o.d)
        d1 = self.d / g
        d2 = o.d / g
        return RF(self.n * d2 + o.n * d1, self.d * d2)

    __radd__ = __add__

    def __neg__(self):
        return RF(-self.n, self.d, reduce=False)

    def __sub__(self, o):
        if isinstance(o, int):
            o = RF(o)
        return self + (-o)

    def __rsub__(self, o):
        return RF(o) - self

    def __mul__(self, o):
        if isinstance(o, int):
            return RF(self.n * o, self.d)
        if self.n.is_zero() or o.n.is_zero():
            return RF(0)
        g1 = self.n.gcd(o.d)
        g2 = o.n.gcd(self.d)
        return RF((self.n / g1) * (o.n / g2), (self.d / g2) * (o.d / g1), reduce=True)

    __rmul__ = __mul__

    def inv(self):
        return RF(self.d, self.n)

    def __truediv__(self, o):
        if isinstance(o, int):
            return RF(self.n, self.d * o)
        return self * o.inv()

    def __pow__(self, e):
        if e >= 0:
            return RF(self.n ** e, self.d ** e, reduce=False)
        return RF(self.d ** (-e), self.n ** (-e), reduce=False)

    def iszero(self):
        return self.n.is_zero()

    def __eq__(self, o):
        if isinstance(o, int):
            o = RF(o)
        return (self.n * o.d - o.n * self.d).is_zero()

    def subs_x(self, xnew_poly):
        """substitute x -> xnew (an fmpz_mpoly), q -> q"""
        return RF(self.n.compose(xnew_poly, Q), self.d.compose(xnew_poly, Q))

    def factor_str(self):
        def fs(p):
            c, fl = p.factor()
            parts = []
            for f, e in fl:
                s = str(f)
                parts.append('(' + s + ')' + ('^%d' % e if e > 1 else ''))
            return c, parts
        cn, pn = fs(self.n)
        cd, pd = fs(self.d)
        from fractions import Fraction
        c = Fraction(int(cn), int(cd))
        return '%s * %s / [%s]' % (c, ' '.join(pn) if pn else '1', ' '.join(pd) if pd else '1')

    def __repr__(self):
        return '(%s)/(%s)' % (self.n, self.d)


def rx(): return RF(X)
def rq(): return RF(Q)
def mono(i, j):
    """x^i q^j as RF (j may be negative)"""
    if j >= 0:
        return RF(X ** i * Q ** j, ONE, reduce=False) if i >= 0 else RF(Q ** j, X ** (-i))
    return RF(X ** i if i >= 0 else ONE, Q ** (-j) * (ONE if i >= 0 else X ** (-i)))


def poch(a_poly, n, step=1):
    """(a; q^step)_n for a polynomial a (fmpz_mpoly), n >= 0 -> RF"""
    p = ONE
    for i in range(n):
        p = p * (ONE - a_poly * Q ** (step * i))
    return RF(p, ONE, reduce=False)


def qpoch(n, step=1):
    """(q^step; q^step)_n"""
    return poch(Q ** step, n, step)
