"""numeric helpers: evaluate F, G=R*F, and the four terms of Au's Proposition 2.1"""
import re
import mpmath as mp
from wz import pari


def pari_to_func(expr, varnames):
    """compile a PARI rational function into a python function of mp numbers"""
    num = str(pari.numerator(expr)).replace('^', '**')
    den = str(pari.denominator(expr)).replace('^', '**')
    rx = re.compile(r'(?<![\w.])(\d+)(?![\w.])')
    num = rx.sub(r'mpf(\1)', num)
    den = rx.sub(r'mpf(\1)', den)
    code = "lambda %s: (%s)/(%s)" % (",".join(varnames), num, den)
    return eval(code, {'mpf': mp.mpf, '__builtins__': {}})


class WZNum:
    def __init__(self, F, R, extra=()):
        """F Term in n,k (other symbols listed in extra must be given at call time)"""
        self.F = F
        self.vars = ['n', 'k'] + list(extra)
        self.Rf = pari_to_func(R, self.vars)

    def Fv(self, n, k, **kw):
        vals = {'n': mp.mpmathify(n), 'k': mp.mpmathify(k)}
        vals.update({v: mp.mpmathify(x) for v, x in kw.items()})
        return self.F.evalf(vals, mp)

    def Gv(self, n, k, **kw):
        args = [mp.mpmathify(n), mp.mpmathify(k)] + [mp.mpmathify(kw[v]) for v in self.vars[2:]]
        return self.Rf(*args) * self.Fv(n, k, **kw)

    def sumF(self, n, k0=0, **kw):
        return mp.nsum(lambda k: self.Fv(n, k + k0, **kw), [0, mp.inf])

    def sumG(self, a, k=0, n0=0, **kw):
        return mp.nsum(lambda n: self.Gv(n + a, k, **kw), [n0, mp.inf])


def hsum_k(W, n, kmax_factor=40, tol=None, k_start=0, **kw):
    """robust sum_{k>=k_start} F(n,k): direct summation with ratio recurrence past the bulk;
       if the terms decay only algebraically, the tail is extrapolated with mp.nsum."""
    F = W.F
    if not hasattr(W, '_sig'):
        W._sig = pari_to_func(F.ratio('k'), W.vars)
        W._rho = pari_to_func(F.ratio('n'), W.vars)
    tol = tol or mp.mpf(10) ** (-mp.mp.dps + 5)
    extra = [mp.mpmathify(kw[v]) for v in W.vars[2:]]
    n = mp.mpmathify(n)
    k = mp.mpf(k_start)
    t = W.Fv(n, k, **kw)
    s = t
    small = 0
    it = 0
    kbulk = kmax_factor * (abs(n) + 5)
    while True:
        r = W._sig(n, k, *extra)
        t = t * r
        k += 1
        s += t
        it += 1
        if abs(t) < tol * max(abs(s), mp.mpf(10) ** (-mp.mp.dps)):
            small += 1
            if small > 30 and k > kbulk:
                break
        else:
            small = 0
        if abs(t) > mp.mpf(10) ** 300:
            raise RuntimeError("divergent")
        if k > kbulk + 200:
            # algebraic decay or divergence: check ratio -> 1
            if abs(abs(r) - 1) < 0.05:
                tail = mp.nsum(lambda j: W.Fv(n, j, **kw), [k + 1, mp.inf])
                return s + tail
            if abs(r) > 1:
                raise RuntimeError("divergent")
        if it > 200000:
            raise RuntimeError("no convergence")
    return s


def hsum_n(W, a, k=0, n_terms=None, tol=None, **kw):
    """sum_{n>=0} G(n+a, k) by direct evaluation (ratio recurrence for F in n, R evaluated each step)"""
    if not hasattr(W, '_sig'):
        W._sig = pari_to_func(W.F.ratio('k'), W.vars)
        W._rho = pari_to_func(W.F.ratio('n'), W.vars)
    tol = tol or mp.mpf(10) ** (-mp.mp.dps + 5)
    extra = [mp.mpmathify(kw[v]) for v in W.vars[2:]]
    a = mp.mpmathify(a)
    k = mp.mpmathify(k)
    n = a
    f = W.Fv(n, k, **kw)
    s = mp.mpf(0)
    small = 0
    it = 0
    while True:
        g = W.Rf(n, k, *extra) * f
        s += g
        if abs(g) < tol * max(abs(s), mp.mpf(10) ** (-mp.mp.dps)):
            small += 1
            if small > 30:
                break
        else:
            small = 0
        f = f * W._rho(n, k, *extra)
        n += 1
        it += 1
        if n_terms and it >= n_terms:
            break
        if it > 60000 or abs(g) > mp.mpf(10) ** 200:
            raise RuntimeError("no convergence")
    return s
