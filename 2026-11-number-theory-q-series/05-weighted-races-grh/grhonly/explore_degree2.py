#!/usr/bin/env python3
"""
explore_degree2.py -- NUMERICAL EXPLORATION ONLY (status [num]; no certificate).
What would a degree-two multiplier kappa = kappa_0 (1 - L + s (L^2 - E L^2)),
L = sum b_n cos(theta_n), give in Theorem 6.1 (GRH-ONLY.md, Remark 6.3.1)?
The non-resonance conditions of order 6 (head), 5 (head vs other zeros) and
4 (head vs pairs) that a proof would need are NOT verified here, and
||kappa||_inf is bounded by 1 + sum b + |s| ((sum b)^2 + sum b^2/2) (rigorous)
with positivity of kappa only checked on the grid.
"""
import os, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")

def run(zfile, M, VA, VB, N):
    g = np.array([float(l) for l in open(os.path.join(ZDIR, zfile)) if l.strip()])
    T = g[-1] + 0.25
    amp = 4 * M / np.sqrt(M * M + g * g)
    a = amp[:N]
    Tp = T - 4 * g[N - 1]                       # buffer for a degree-2 multiplier
    K = (g > g[N - 1]) & (g <= Tp)
    Vh, VK = (a ** 2).sum() / 2, (amp[K] ** 2).sum() / 2
    n = 8
    th = 2 * np.pi * np.arange(n) / n
    grids = np.meshgrid(*([th] * N), indexing="ij")
    C = [np.cos(x) for x in grids]
    H = sum(a[i] * C[i] for i in range(N))
    res = {}
    for tier, V in (("A", VA), ("B", VB)):
        Vt = V - Vh - VK
        def J(p):
            k0, s, b = p[0], p[1], np.abs(p[2:])
            L = sum(b[i] * C[i] for i in range(N))
            kap = k0 * (1 - L + s * (L * L - (b ** 2).sum() / 2))
            if kap.min() < 0 or k0 <= 0:
                return 1e9
            sup = k0 * (1 + b.sum() + abs(s) * (b.sum() ** 2 + (b ** 2).sum() / 2))
            return np.mean((1 - kap * (1 + H)) ** 2) + np.mean(kap ** 2) * VK + sup ** 2 * Vt
        p = np.array([1 / (1 + V), 0.0] + list(0.6 * a))
        best = J(p); step = 0.1
        while step > 1e-5:
            imp = False
            for i in range(len(p)):
                for sg in (1, -1):
                    q = p.copy(); q[i] += sg * step
                    v = J(q)
                    if v < best - 1e-13:
                        best, p, imp = v, q, True
            if not imp:
                step /= 2
        res[tier] = (best, p)
    return res

for name, zf, VA, VB in (("chi_4", "chi4_zeros.txt", 0.169630835, 0.155573631),
                         ("chi_3", "chi3_zeros.txt", 0.125451412, 0.113234208)):
    for N in (2, 3, 4, 5):
        r = run(zf, 0.5, VA, VB, N)
        print("%s m=0 N=%d  degree-2 exploration: tier A density >= %.4f (s=%.2f), tier B >= %.4f (s=%.2f)   [num, not certified]"
              % (name, N, 1 - r["A"][0], r["A"][1][1], 1 - r["B"][0], r["B"][1][1]))
