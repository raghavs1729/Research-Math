#!/usr/bin/env python3
"""fq_spot.py q1 q2 ... -- independent F_q (AF (9)) for composite q: own enumeration of all
characters mod q (brute-force character table via CRT generators), conductor and primitive
inducing character found by brute force, F(chi') = log(d/pi) + 2 Re L'/L(1,chi') + digamma((1+kappa)/2),
F(principal) = 2 + gamma - log 4pi; R_{q,1} = #{b : b^2 = 1 mod q}.  mpmath, 30 digits."""
import sys, itertools
from math import gcd
import mpmath as mp
mp.mp.dps = 30


def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def units(n):
    return [m for m in range(1, n + 1) if gcd(m, n) == 1]


def gens(q):
    """list of (generator g mod q, order) such that the unit group is the direct product."""
    out = []
    for p, e in factor(q).items():
        pe = p ** e; rest = q // pe
        def lift(x):  # x mod pe, 1 mod rest
            for t in range(q):
                if t % pe == x % pe and t % rest == 1 % rest:
                    return t
        if p == 2:
            if e >= 2:
                out.append((lift(pe - 1), 2))
            if e >= 3:
                out.append((lift(5), 2 ** (e - 2)))
        else:
            ph = pe - pe // p
            g = next(g for g in range(2, pe) if gcd(g, p) == 1 and
                     len({pow(g, k, pe) for k in range(ph)}) == ph)
            out.append((lift(g), ph))
    return out


def chartable(q):
    G = gens(q)
    # discrete log of every unit as a vector over G
    logs = {}
    for ks in itertools.product(*[range(o) for _, o in G]):
        m = 1
        for (g, o), k in zip(G, ks):
            m = m * pow(g, k, q) % q
        logs[m] = ks
    assert len(logs) == len(units(q)) or q == 1
    chars = []
    for js in itertools.product(*[range(o) for _, o in G]):
        val = {m: sum(mp.mpf(j * k) / o for j, k, (_, o) in zip(js, logs[m], G)) for m in logs}
        chars.append({m: mp.expjpi(2 * v) for m, v in val.items()})
    return chars


def conductor(q, ch):
    for d in sorted(x for x in range(1, q + 1) if q % x == 0):
        if all(abs(ch[m] - 1) < 1e-20 for m in ch if m % d == 1 % d):
            return d


def F_q(q):
    F = 2 + mp.euler - mp.log(4 * mp.pi)
    if q <= 2:
        return F, 1
    nprim = 0
    for ch in chartable(q):
        if all(abs(v - 1) < 1e-20 for v in ch.values()):
            continue
        d = conductor(q, ch)
        cp = {}
        for r in units(d):
            m = next(m for m in ch if m % d == r % d)
            cp[r] = ch[m]
        kappa = 0 if abs(ch[q - 1] - 1) < 1e-20 else 1
        L = lambda s: mp.mpf(d) ** (-s) * sum(cp[r] * mp.zeta(s, mp.mpf(r) / d) for r in cp)
        h = mp.mpf("1e-12")
        Lp = (L(1 + h) - L(1 - h)) / (2 * h); L1 = (L(1 + h) + L(1 - h)) / 2
        F += mp.log(mp.mpf(d) / mp.pi) + 2 * mp.re(Lp / L1) + mp.digamma(mp.mpf(1 + kappa) / 2)
    return F, sum(1 for b in units(q) if b * b % q == 1)


for q in map(int, sys.argv[1:]):
    F, R = F_q(q)
    print("q = %4d  F_q = %s  2R_{q,1} = %d  F_q < 2R: %s" % (q, mp.nstr(F, 12), 2 * R, F < 2 * R))
