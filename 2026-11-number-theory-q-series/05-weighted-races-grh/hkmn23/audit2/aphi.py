"""Second audit: Phi(q) = F_{q,1}(0) of Theorem 5 from the closed form
c(chi) = Re[(L'/L)(2,chi) - (L'/L)(1,chi)] + delta(parity), recomputed with mpmath
(L'/L(1) via digamma and Stieltjes constants gamma_1(k/f), L'/L(2) via Hurwitz zeta), own characters.
Also the GRH-free comparison with (2/3)N(q,1) and with the author's raw/t2_closed blocks.
"""
import sys, json, glob
import mpmath as mp
from math import gcd
from sympy import divisors, totient
from achars import primitive_characters, parity

mp.mp.dps = 25
_cache = {}


def blocks(f):
    if f in _cache:
        return _cache[f]
    ks = [k for k in range(1, f + 1) if gcd(k, f) == 1]
    dg = {k: mp.digamma(mp.mpf(k) / f) for k in ks}
    st = {k: mp.stieltjes(1, mp.mpf(k) / f) for k in ks}
    z2 = {k: mp.zeta(2, mp.mpf(k) / f) for k in ks}
    z2d = {k: mp.zeta(2, mp.mpf(k) / f, 1) for k in ks}
    tot = mp.mpf(0)
    for vals in primitive_characters(f):
        c = {k: mp.expjpi(2 * mp.mpf(vals[k].numerator) / vals[k].denominator) for k in ks}
        L1 = sum(c[k] * (-dg[k]) for k in ks) / f
        L1p = -mp.log(f) * L1 + sum(c[k] * (-st[k]) for k in ks) / f
        L2 = sum(c[k] * z2[k] for k in ks) / f ** 2
        L2p = -mp.log(f) * L2 + sum(c[k] * z2d[k] for k in ks) / f ** 2
        a = parity(vals, f)
        d = mp.log(2) if a == 0 else 1 - mp.log(2)
        tot += mp.re(L2p / L2 - L1p / L1) + d
    _cache[f] = tot
    return tot


def Phi(q):
    cz = mp.log(2 * mp.pi) + mp.log(2) - mp.mpf(1) / 2 - mp.zeta(-1, 1, 1) / mp.zeta(-1)
    return cz + sum(blocks(f) for f in divisors(q) if f > 1)


def N1(q):
    return sum(1 for b in range(1, q + 1) if gcd(b, q) == 1 and (b * b - 1) % q == 0) if q > 1 else 1


if __name__ == '__main__':
    # author's blocks
    auth = {}
    for fn in glob.glob('../raw/t2_closed_400_*.jsonl'):
        for line in open(fn):
            r = json.loads(line)
            auth[r.get('f', r.get('q'))] = r
    for q in [int(x) for x in sys.argv[1:]]:
        P = Phi(q); thr = mp.mpf(2) / 3 * N1(q)
        print(f'q={q:4d} phi={int(totient(q)):3d} N={N1(q):2d} Phi={mp.nstr(P, 12)}  (2/3)N={mp.nstr(thr, 8)}  Phi-(2/3)N={mp.nstr(P - thr, 8)}  rel={mp.nstr((P - thr) / thr, 4)}')
