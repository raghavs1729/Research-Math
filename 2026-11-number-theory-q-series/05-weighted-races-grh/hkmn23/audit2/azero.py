"""Second audit: independent re-certification of a zero-side witness (Theorem 6) from the author's cached
zero brackets.  Everything except the bracket list is mine:
  * characters (achars.py), matched to the cache files by testing sign changes;
  * every bracket re-certified: the real function u*Lambda(1/2+it,chi) (completed L-function times a
    constant phase, u = eps^{-1/2}, eps from the Gauss sum) has certified opposite signs at the two ends,
    with L(s,chi) = f^{-s} sum chi(k) zeta(s,k/f) (flint Hurwitz zeta, not flint's Dirichlet code);
  * zeta zeros cross-checked against mpmath.zetazero;
  * S(chi) from al1 (mpmath, Lemma S re-implemented), the remainder r and F_Z(y0) in Arb.
Usage: python3 azero.py q a y0 T
"""
import sys, json, math, glob
import flint
from flint import acb, arb, fmpq
import mpmath as mp
from sympy import divisors
from achars import primitive_characters, parity
from al1 import Char

flint.ctx.prec = 96
mp.mp.dps = 30


def cval(v):
    return acb.exp_pi_i(acb(fmpq(2 * v.numerator, v.denominator)))


class PC:
    def __init__(self, vals, f):
        self.vals, self.f, self.a = vals, f, parity(vals, f)
        self.c = [cval(v) if v is not None else acb(0) for v in vals]
        tau = sum((self.c[k] * acb.exp_pi_i(acb(fmpq(2 * k, f))) for k in range(1, f)), acb(0))
        eps = tau / (acb(1j) ** self.a * arb(f).sqrt())
        self.u = 1 / eps.sqrt()

    def Zt(self, t):
        f = self.f
        s = acb(arb(1) / 2, arb(t))
        L = sum((self.c[k % f] * s.zeta(acb(fmpq(k, f))) for k in range(1, f + 1) if self.vals[k % f] is not None), acb(0))
        L *= acb(f) ** (-s)
        E = (acb(f) / acb(arb.pi())) ** ((s + self.a) / 2) * ((s + self.a) / 2).gamma() * L
        return self.u * E

    def sign(self, t):
        z = self.Zt(t)
        assert z.imag.contains(0) or abs(z.imag) < abs(z.real) * 1e-6, ('not real', z)
        if z.real > 0:
            return 1
        if z.real < 0:
            return -1
        return 0


def bracket_ok(ch, g):
    lo, hi = g.lower(), g.upper()
    return ch.sign(lo) * ch.sign(hi) < 0


def main(q, a, y0, T):
    from math import gcd
    N = sum(1 for b in range(q) if gcd(b, q) == 1 and (b * b - a) % q == 0)
    thr = arb(2 * N) / 3
    y0 = arb(y0)
    total_S = arb(0); found = arb(0)
    FZ = acb(0)
    # zeta
    zz = []
    k = 1
    while True:
        batch = acb.zeta_zeros(k, 50)
        stop = False
        for z in batch:
            if z.imag > T:
                stop = True; break
            zz.append(z.imag)
        if stop:
            break
        k += 50
    for i in [0, 9, len(zz) - 1]:
        assert abs(float(zz[i].mid()) - float(mp.im(mp.zetazero(i + 1)))) < 1e-12
    w = lambda g: 1 / ((arb(1) / 4 + g * g) * (arb(9) / 4 + g * g)).sqrt()
    term = lambda g: acb.exp_pi_i(acb(0, 0) + acb(g * y0 / arb.pi())) / (acb(arb(1) / 2, g) * acb(arb(3) / 2, g))
    Sz = arb(mp.nstr(Char([0], 1).S(), 30)) + arb('1e-25') * arb(0, 1)
    total_S += Sz
    for g in zz:
        found += 2 * w(g)
        FZ -= term(g) + term(-g)
    nb = 2 * len(zz)
    report = []
    for f in divisors(q):
        if f < 3:
            continue
        files = sorted(glob.glob(f'../raw/zeros_cache/z_{f}_*_T{int(T)}.json'))
        prims = [PC(v, f) for v in primitive_characters(f)]
        brackets = {fn: [arb(x) for x in json.load(open(fn))] for fn in files}
        # match files to my characters by the first brackets
        match = {}
        for fn, gs in brackets.items():
            cands = [i for i, ch in enumerate(prims) if all(bracket_ok(ch, g) for g in gs[:6])]
            assert len(cands) == 1, (fn, cands)
            match[cands[0]] = fn
        assert len(match) == len(prims), (f, len(match), len(prims))
        for i, ch in enumerate(prims):
            gs = brackets[match[i]]
            # disjointness and certification of every bracket
            for j in range(1, len(gs)):
                assert gs[j - 1] < gs[j]
            bad = sum(1 for g in gs if not bracket_ok(ch, g))
            assert bad == 0, (f, i, bad)
            S = Char(ch.vals, f).S()
            total_S += arb(mp.nstr(S, 30)) + arb('1e-25') * arb(0, 1)
            # positive zeros of chi from its file; negative zeros of chi = - positive zeros of conj(chi)
            conjvals = [(-v) % 1 if v is not None else None for v in ch.vals]
            jb = [k for k, c2 in enumerate(prims) if c2.vals == conjvals][0]
            neg = brackets[match[jb]]
            cb = ch.c[a % f].conjugate()
            for g in gs:
                found += w(g); FZ -= cb * term(g)
            for g in neg:
                found += w(g); FZ -= cb * term(-g)
            nb += len(gs) + len(neg)
            report.append((f, i, len(gs), match[i].split('/')[-1]))
    # characters mod q are in bijection with primitive chars of conductors f | q: the sums above cover all of them
    r = total_S - found
    lower = FZ.real - r
    print(f'({q},{a}) y0={y0.str(15)}: {nb} zeros used (all brackets re-certified), L_1 = {total_S.str(12)}, r = {r.str(8)}')
    print(f'  F_Z(y0) = {FZ.real.str(10)} (imag {FZ.imag.str(3)}),  F(y0) >= {lower.str(10)}  vs (2/3)N = {thr.str(8)}  -> certified: {bool(lower > thr)}')
    for x in report:
        print('   ', x)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]))
