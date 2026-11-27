#!/usr/bin/env python3
"""
a4_plain_mc.py -- NON-RIGOROUS sanity checks (audit item 3), own code.
 (P) plain float evaluation: composite 20-point Gauss-Legendre (numpy) on [0, 60]
     of prod J0(b_j t) * exp(-V t^2/2) * sin t / t, all listed zeros, V from the
     closed form with mpmath L'/L (numerical differentiation of mp.dirichlet).
 (MC) Monte Carlo with common random phases for the two M values, unlisted zeros
     replaced by N(0, V); paired difference.
Usage: python3 a4_plain_mc.py q M1 M2 nsamples
"""
import sys, math
import numpy as np
from scipy.special import j0
import mpmath as mp

HERE = "/home/user/weighted-prime-races/monotonicity/nonmonotone"


def kron_list(q):
    # chi(n) = (-q/n) via Euler's criterion for prime q=163; generic via flint-free Jacobi
    def jac(a, n):
        a %= n; r = 1
        while a:
            while a % 2 == 0:
                a //= 2
                if n % 8 in (3, 5): r = -r
            a, n = n, a
            if a % 4 == 3 and n % 4 == 3: r = -r
            a %= n
        return r if n == 1 else 0
    D = -q
    out = [0] * q
    for n in range(1, q):
        m, r = n, 1
        while m % 2 == 0:
            m //= 2
            if D % 2 == 0: r = 0; break
            r *= 1 if D % 8 in (1, 7) else -1
        out[n] = r * jac(D, m) if r else 0
    return out


def main():
    q = int(sys.argv[1]); Ms = [float(eval(x)) for x in sys.argv[2:4]]; ns = int(float(sys.argv[4]))
    g = np.array([float(l.split()[0]) for l in open("%s/zeros_q%d.txt" % (HERE, q)) if l.strip() and not l.startswith("#")])
    chi = kron_list(q)
    mp.mp.dps = 30
    res = {}
    for M in Ms:
        s = mp.mpf(M) + mp.mpf(1) / 2
        L = lambda x: mp.dirichlet(x, chi)
        lp = mp.diff(L, s) / L(s)
        sig2 = 4 * M * (mp.log(q / mp.pi) / 2 + mp.digamma((s + 1) / 2) / 2 + lp)
        b = 4 * M / np.sqrt(M * M + g * g)
        V = float(sig2) - float(np.sum(b * b / 2))
        x, w = np.polynomial.legendre.leggauss(20)
        edges = np.linspace(0, 60, 6001)
        tot = 0.0
        for a, c in zip(edges[:-1], edges[1:]):
            t = (a + c) / 2 + (c - a) / 2 * x
            P = np.prod(j0(np.outer(t, b)), axis=1)
            tot += (c - a) / 2 * np.sum(w * P * np.exp(-V * t * t / 2) * np.sin(t) / t)
        res[M] = (b, V)
        print("q=%d M=%g: sigma2=%.12f V=%.8f  plain delta = %.12f" % (q, M, float(sig2), V, 0.5 + tot / math.pi), flush=True)
    rng = np.random.default_rng(12345)
    B = 50000
    diffs = []; hits = {M: [] for M in Ms}
    for _ in range(ns // B):
        cs = np.cos(rng.uniform(0, 2 * np.pi, size=(B, len(g))))
        z = rng.standard_normal(B)
        h = []
        for M in Ms:
            b, V = res[M]
            X = 1 + cs @ b + math.sqrt(V) * z
            h.append(X > 0)
            hits[M].append(h[-1].mean())
        diffs.append((h[1].astype(float) - h[0]).mean())
    d = np.array(diffs)
    for M in Ms:
        a = np.array(hits[M]); print("  MC M=%g: %.6f +- %.6f" % (M, a.mean(), a.std(ddof=1) / math.sqrt(len(a))))
    print("  MC paired delta(%g)-delta(%g) = %+.6f +- %.6f  (%d samples)" % (Ms[1], Ms[0], d.mean(), d.std(ddof=1) / math.sqrt(len(d)), len(d) * B))


if __name__ == "__main__":
    main()
