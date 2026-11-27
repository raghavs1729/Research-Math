#!/usr/bin/env python3
"""
mertens_C.py q [X]  -- independent value of log C(q, 1) (Theorem 1.8 / L-Z constant), prime q.

  phi(q) log C(q,1) = -gamma - log(1 - 1/q)
                      + sum_{chi != chi0} [ -log|L(1,chi)| + Re sum_p sum_{k>=2} (chi(p)^k - chi(p))/(k p^k) ]
(derivation: sum_{p<=x, p=1 (q)} log(1-1/p) = (1/phi) sum_chi sum_{p<=x} chi(p) log(1-1/p);
 principal part = -gamma - log log x - log(1-1/q) + o(1) (Mertens);  for chi != chi0,
 sum_p chi(p) log(1-1/p) = -log L(1,chi) + sum_p sum_{k>=2}(chi(p)^k - chi(p))/(k p^k), and for
 a = 1 only real parts survive, so no branch of log L(1,chi) is needed.)
The p-sum is truncated at X; |tail| <= sum_{p>X} 2/(p(p-1)) per character < 2/X, so the
total error in log C is < 2 (phi-1)/(phi X) < 2/X.
"""
import sys
import numpy as np
import mpmath as mp

mp.mp.dps = 30
q = int(sys.argv[1]); X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 7
phi = q - 1
g = next(g for g in range(2, q) if len({pow(g, k, q) for k in range(phi)}) == phi)
ind = {pow(g, k, q): k for k in range(phi)}

sv = np.ones(X + 1, dtype=bool); sv[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
P = np.nonzero(sv)[0]
P = P[P % q != 0]
idx = np.array([ind[int(r)] for r in range(1, q)] , dtype=np.int64)  # ind of residue r at position r-1
ip = idx[(P % q) - 1]
pf = P.astype(np.float64)

tot = -mp.euler - mp.log(1 - mp.mpf(1) / q)
for j in range(1, phi):
    L1 = -sum(mp.expjpi(mp.mpf(2 * j * ind[m]) / phi) * mp.digamma(mp.mpf(m) / q) for m in range(1, q)) / q
    # Re sum_{k>=2} (chi^k - chi)/(k p^k) with chi = e(j ind/phi); sum k up to 40 (p >= 2: 2^-40 negligible
    # relative to the 2/X budget is false for p=2? 2^-41/41 < 1e-14, fine)
    s = 0.0
    th = 2 * np.pi * j * ip / phi
    for k in range(2, 60):
        s += np.sum((np.cos(k * th) - np.cos(th)) / (k * pf ** k))
    tot += -mp.log(abs(L1)) + s
logC = tot / phi
print("q = %d  X = %d  log C(q,1) = %s  (+- %.1e)  C = %s  1/C = %s" %
      (q, X, mp.nstr(logC, 15), 2.0 / X + 1e-12, mp.nstr(mp.exp(logC), 15), mp.nstr(mp.exp(-logC), 15)))
