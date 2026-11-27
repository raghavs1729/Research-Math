#!/usr/bin/env python3
"""
general_pairs.py -- numerics for Section 8 of EXPLICIT.md (general weight pairs m < m').

  Phi_{m,m'}(t) = sum_{p<=t} s(p) log p * p^m (1 - (p/t)^kappa),  kappa = m'-m,  M = m+1/2, M' = m'+1/2
  heuristic (GRH):  Phi(t) = (kappa/2) t^M ( 1/(M M') + sum_{gamma>0} 4 Re[t^{i gamma}/((M+i gamma)(M'+i gamma))] + o(1) )
  Lambda(M,M') = 1/(M M') - sum_{gamma>0} 4/(|M+i gamma| |M'+i gamma|)      (zeros + Riemann-von Mangoldt tail; NOT certified here)

Phi is continuous and monotone between consecutive primes, so its sign on [p0, X] is decided at primes.
Lemma 6:  if Phi >= 0 on [p0, x] then I_{m,m'}(x) = D_m(x) - x^{-kappa} D_{m'}(x) >= Phi(x)/log x.
"""
import os, sys, math
import numpy as np, mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = 1 << 23

def zeros(q):
    return [float(l.split()[0]) for l in open(os.path.join(HERE, '..', '..', 'prototype', 'explore', 'chi%d_zeros.txt' % q)) if l.strip()]

def Lam(q, M, M2, zs):
    T = zs[-1]
    part = sum(4 / math.sqrt((M * M + g * g) * (M2 * M2 + g * g)) for g in zs)
    tail = float(mp.quad(lambda t: 4 / mp.sqrt((M * M + t * t) * (M2 * M2 + t * t)) * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf]))
    return 1 / (M * M2) - part - tail

sv = np.ones(NMAX + 1, dtype=bool); sv[:2] = False
for i in range(2, int(NMAX ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = False
allp = np.nonzero(sv)[0]
for q in (4, 3):
    if q == 4:
        P = allp[allp > 2]; s = np.where(P % 4 == 3, 1.0, -1.0)
    else:
        P = allp[allp != 3]; s = np.where(P % 3 == 2, 1.0, -1.0)
    x = P.astype(np.longdouble); lx = np.log(x); s = s.astype(np.longdouble)
    zs = zeros(q)
    print('q = %d, primes <= 2^23' % q)
    print('  %-12s %-9s %-10s | %-22s %-28s | %-26s' % ("(m,m')", 'Lambda', '(k/2)Lam', 'primes with Phi<=0: #, last', 'min Phi/(k t^M) on [1e4,2^23] (at)', 'primes with I<=0: #, last'))
    for (m, m2) in ((0, 1), (0, 0.5), (-0.25, 0.75), (0, 1.5), (0, 2), (0, 3), (0.5, 1.5), (1, 2), (1, 3), (1.5, 2.5), (2, 3), (3, 4)):
        k = m2 - m; M = m + 0.5; M2 = m2 + 0.5
        A = np.cumsum(s * lx * x ** m); Bv = np.cumsum(s * lx * x ** m2)
        Phi = A - x ** (-k) * Bv
        Dm = np.cumsum(s * x ** m); Dm2 = np.cumsum(s * x ** m2)
        I = Dm - x ** (-k) * Dm2
        badP = np.nonzero(Phi[1:] <= 0)[0]; badI = np.nonzero(I[1:] <= 0)[0]
        nrm = Phi / (k * x ** M)
        sel = x >= 1e4
        j = np.argmin(nrm[sel])
        lam = Lam(q, M, M2, zs)
        print('  %-12s %-9.4f %-10.4f | %-6d %-15s %-10.4f at %-14d | %-6d %-12s'
              % ('(%g,%g)' % (m, m2), lam, lam / 2, len(badP), (int(P[1:][badP[-1]]) if len(badP) else '-'), float(nrm[sel][j]), int(P[sel][j]),
                 len(badI), (int(P[1:][badI[-1]]) if len(badI) else '-')))
