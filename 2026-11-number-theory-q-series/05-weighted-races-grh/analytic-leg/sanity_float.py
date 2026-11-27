#!/usr/bin/env python3
"""sanity_float.py -- NON-RIGOROUS float64 sanity test of the inequalities of
lem:diffint(b) (eq:Abound), lem:negstruct (eq:mainlb pieces) and the zone
bounds of lem:analytic-leg, against the actual zeros, at weights where the
quantities are not astronomically small.  Nothing here is part of the proof;
it only guards against a wrong direction / wrong factor in the chain.

Truncation: only the listed zeros (gamma <= 640 resp. 700) are used, so
Phi and D are the truncated ones; every inequality tested is termwise in the
zeros and holds for the truncated product/sums as well.
"""
import math, os
import numpy as np
from scipy.special import j0, j1, erf

HERE = os.path.dirname(os.path.abspath(__file__))
EXPLORE = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))
bad = 0
for q, fn in ((4, "chi4_zeros.txt"), (3, "chi3_zeros.txt")):
    g = np.loadtxt(os.path.join(EXPLORE, fn))
    for M in (40.0, 100.0, 300.0):
        a = 2 * M / np.sqrt(M * M + g * g)
        da = 2 * g * g / (M * M + g * g) ** 1.5
        D = (4 * a * da).sum()
        s2 = (2 * a * a).sum()
        S4 = ((2 * a) ** 4).sum()
        S = g <= 2 * M
        NM, N2M = int((g <= M).sum()), int(S.sum())
        worstA = worstZ = worstMain = -1e9
        for t in np.concatenate([np.linspace(0.6, 16, 400), np.linspace(16, 200, 200)]):
            J = j0(2 * a * t)
            logabs = np.log(np.abs(J) + 1e-300)
            tot = logabs.sum()
            # exact |d/dM Phi| <= sum_g |J1| 2 a' t prod_{g' != g} |J0|
            terms = np.abs(j1(2 * a * t)) * 2 * da * t * np.exp(tot - logabs)
            lhs = terms.sum()
            # eq:Abound
            A = D * t * t / 2 + 4 * t / M * (a * t > 1).sum()
            Aex = (np.minimum(1, a * t) * 2 * da * t).sum()
            prodS = np.exp(logabs[S].sum() - logabs[S].min())      # >= sup_g prod_{S\g}
            supS = np.exp((logabs[S].sum() - logabs[S])).max()
            rhs = A * supS
            worstA = max(worstA, math.log(lhs + 1e-300) - math.log(rhs + 1e-300), Aex - A)
            # zone bound of lem:analytic-leg
            if t <= 16:
                zone = 0.447 ** NM * 0.733 ** (N2M - NM) / 0.447
            else:
                zone = (1.9041 / math.sqrt(t)) ** 9 * 0.733 ** (N2M - 10)
            worstZ = max(worstZ, math.log(supS + 1e-300) - math.log(zone))
        for t in np.linspace(1e-4, 0.125, 200):
            z = 2 * a * t
            Phi = np.exp(np.log(j0(z)).sum())
            B = (j1(z) / j0(z) * 2 * da * t).sum()
            beta = s2 / 2 + S4 / 1088
            worstMain = max(worstMain, math.exp(-beta * t * t) - Phi, D * t * t / 2 - B)
        ok = worstA <= 1e-9 and worstZ <= 1e-9 and worstMain <= 1e-12
        bad += not ok
        print(f"q={q} M={M:5.0f}: N(M)={NM:3d} N(2M)={N2M:3d}  max log(lhs/rhs) eq:Abound = {worstA:9.2f}  "
              f"max log(sup prod / zone bound) = {worstZ:9.2f}  main-term pointwise slack >= {-worstMain:.1e}  "
              f"{'OK' if ok else 'VIOLATION'}")
print("sanity_float:", "OK" if not bad else f"{bad} VIOLATION(S)")
