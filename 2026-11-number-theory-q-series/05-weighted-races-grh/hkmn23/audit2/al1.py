"""Second audit: L_1(q) by Lemma S, re-implemented with mpmath (Hurwitz zeta), own
characters (achars.py), and a DIFFERENT quadrature (tanh-sinh, mpmath.quad) instead of the
author's 64-point trapezoid in Arb.  Also the closed-form sandwich.  Also the constants
b, b0 / c_{-1} from the functional equation (L'/L at 1 via Stieltjes constants, at 2 via
Hurwitz zeta derivatives), compared with the author's Taylor-coefficient values in
l1const.json.  Also c(chi) of Theorem 5(a).
Numerical (mpmath, 40 digits), not ball arithmetic.
"""
import json, sys
import mpmath as mp
from achars import primitive_characters, parity, characters, conductor

mp.mp.dps = 40


def chival(v):
    return mp.expjpi(2 * mp.mpf(v.numerator) / v.denominator) if v is not None else mp.mpc(0)


class Char:
    def __init__(self, vals, f):
        self.f = f
        self.vals = vals
        self.c = [chival(v) if v is not None else 0 for v in vals]
        self.a = parity(vals, f)

    def L_and_Lp(self, s):
        f = self.f
        if f == 1:
            return mp.zeta(s), mp.zeta(s, 1, 1)  # s != 1 enforced by caller
        S0 = mp.mpc(0); S1 = mp.mpc(0)
        for k in range(1, f + 1):
            ck = self.c[k % f]
            if ck == 0:
                continue
            S0 += ck * mp.zeta(s, mp.mpf(k) / f)
            S1 += ck * mp.zeta(s, mp.mpf(k) / f, 1)
        fs = mp.power(f, -s)
        L = fs * S0
        Lp = -mp.log(f) * L + fs * S1
        return L, Lp

    def logder(self, s):
        if self.f > 1 and s == 1:
            return self.logder1()
        L, Lp = self.L_and_Lp(s)
        return Lp / L

    def logder1(self):
        # L(s) = f^{-s} sum chi(k) [ -psi(k/f) - gamma_1(k/f)(s-1) + ... ]   (sum chi = 0)
        f = self.f
        S0 = mp.mpc(0); S1 = mp.mpc(0)
        for k in range(1, f + 1):
            ck = self.c[k % f]
            if ck == 0:
                continue
            a = mp.mpf(k) / f
            S0 += ck * (-mp.digamma(a))
            S1 += ck * (-mp.stieltjes(1, a))
        L = S0 / f
        Lp = -mp.log(f) * L + S1 / f
        return Lp / L

    def V(self, alpha):
        s = mp.mpf(1) / 2 + alpha
        f = self.f
        if f == 1:
            e = s - 1
            if e == 0:
                ld = mp.euler
            elif abs(e) < mp.mpf(10) ** -8:
                g0, g1, g2 = mp.euler, mp.stieltjes(1), mp.stieltjes(2)
                z = 1 / e + g0 - g1 * e + g2 * e ** 2 / 2
                zp = -1 / e ** 2 - g1 + g2 * e
                ld = zp / z + 1 / e
            else:
                ld = self.logder(s) + 1 / e
            return mp.re(ld) + 1 / s - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2
        return mp.re(self.logder(s)) + mp.log(mp.mpf(f) / mp.pi) / 2 + mp.digamma((s + self.a) / 2) / 2

    def S(self):
        g = lambda th: self.V(mp.sqrt(mp.mpf(5) / 4 - mp.cos(th))) / mp.sqrt(mp.mpf(5) / 4 - mp.cos(th))
        return mp.quad(g, [0, mp.pi / 3, mp.pi], method='gauss-legendre') / mp.pi

    def sandwich(self):
        r5 = mp.sqrt(5) / 2
        lo = self.V(r5) / r5
        hi = (2 * self.V(mp.mpf(1) / 2) + mp.mpf(2) / 3 * self.V(mp.mpf(3) / 2)) / 2
        return lo, hi

    def c_closed(self):
        # Theorem 5(a)
        if self.f == 1:
            return mp.log(2 * mp.pi) + mp.log(2) - mp.mpf(1) / 2 - mp.zeta(-1, 1, 1) / mp.zeta(-1)
        d = mp.log(2) if self.a == 0 else 1 - mp.log(2)
        return mp.re(self.logder(2) - self.logder1()) + d

    def fe_constants(self):
        """b and second constant via the functional equation, using chibar."""
        f = self.f
        bar = Char([(-v) % 1 if v is not None else None for v in self.vals], f)
        l1 = bar.logder1(); l2 = bar.logder(2)
        b = -l1 - mp.log(mp.mpf(f) / mp.pi) + mp.euler + mp.log(2)
        b2 = -l2 - mp.log(mp.mpf(f) / mp.pi) + mp.euler + mp.log(2) - 1
        return b, b2


def L1(q, cache={}):
    from sympy import divisors
    tot = mp.mpf(0)
    for f in divisors(q):
        if f not in cache:
            if f == 1:
                cache[f] = [Char([None], 1)]
                cache[f][0].vals = [0]
            else:
                cache[f] = [Char(v, f) for v in primitive_characters(f)]
        for ch in cache[f]:
            tot += ch.S()
    return tot


if __name__ == '__main__':
    ref = json.load(open('../l1const.json'))
    out = {}
    zeta = Char([0], 1)
    print('S(zeta) =', mp.nstr(zeta.S(), 25), ' sandwich', [mp.nstr(x, 12) for x in zeta.sandwich()])
    print('c(zeta) =', mp.nstr(zeta.c_closed(), 20))
    qs = [int(x) for x in sys.argv[1:]] or [1, 5, 8, 12, 15, 20]
    for q in qs:
        v = L1(q)
        r = mp.mpf(ref['moduli'][str(q)]['L1'][0]) if str(q) in ref['moduli'] else None
        print(f'L_1({q}) = {mp.nstr(v, 25)}   author {mp.nstr(r, 25) if r is not None else "-"}   diff {mp.nstr(v - r, 3) if r is not None else "-"}')
        out[q] = str(v)
    # per-character detail for an even character (5) and the imprimitive case (12, 15)
    for f in [5, 12, 15, 3, 4]:
        for vals in primitive_characters(f):
            ch = Char(vals, f)
            S = ch.S(); lo, hi = ch.sandwich()
            b, b2 = ch.fe_constants()
            print(f'  f={f} parity={ch.a} chi(2)={vals[2]} S={mp.nstr(S, 20)} sandwich=[{mp.nstr(lo, 10)},{mp.nstr(hi, 10)}]'
                  f' b={mp.nstr(b, 15)} second={mp.nstr(b2, 15)} c={mp.nstr(ch.c_closed(), 15)}')
    json.dump(out, open('al1.json', 'w'), indent=1)
