#!/usr/bin/env python3
"""Tables for the section on Chebyshev's bias (heuristic comparison; no theorem).
 (1) observed fraction of failing primes of each class per half-decade (exact scan) against the
     LI-model probability P(Y >= W(x) + c B(L)) (Gil-Pelaez inversion of prod J0(r xi) exp(-tv xi^2/2)
     where it exceeds 0.02, saddle-point approximation below; averaged over the half-decade in log x);
 (2) residue / non-residue comparison of the last failures.
Usage: bias_tables.py merged1 [merged2 ...]"""
import sys, math, collections
import numpy as np
from scipy.special import j0
import model
CQ = [1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12]
CA = [0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11]
def make_tail(q):
    r, tv, info = model.amps(q)
    xi = np.arange(1, 24001) * 0.0025                      # (0, 60]
    logphi = np.zeros_like(xi)
    for i in range(0, len(r), 200):
        logphi += np.sum(np.log(np.abs(j0(np.outer(xi, r[i:i + 200]))) + 1e-300), axis=1)
    sign = np.ones_like(xi)
    for i in range(0, len(r), 200):
        sign *= np.prod(np.sign(j0(np.outer(xi, r[i:i + 200]))), axis=1)
    phi = sign * np.exp(logphi - tv * xi * xi / 2)
    def tail(w):
        v = 0.5 - np.sum(phi * np.sin(w * xi) / xi) * 0.0025 / math.pi
        if v > 0.02: return v                      # numerical inversion is reliable here
        return min(10 ** model.tail_prob(r, tv, w)[1], 0.02)     # far tail: saddle-point approximation
    return tail
bins = collections.defaultdict(lambda: [0, 0])
for fn in sys.argv[1:]:
    for l in open(fn):
        w = l.split()
        if w and w[0] == "BIN":
            c, b = int(w[1]), int(w[2])
            bins[(c, b // 50)][0] += int(w[3]); bins[(c, b // 50)][1] += int(w[4])
print("(1) fraction of the primes of the class at which the inequality fails, per half-decade: observed / LI model")
hd = [h for h in range(16, 30) if any(bins.get((c, h), (0, 0))[0] for c in range(17))]
print("| (q,a) | c | " + " | ".join(f"10^{h/2:g}–10^{(h+1)/2:g}" for h in hd) + " |")
print("|---|---|" + "---|" * len(hd))
for q in (1, 3, 4, 5, 8, 12):
    tail = make_tail(q)
    for c in range(17):
        if CQ[c] != q: continue
        cs = model.CSQ[(q, CA[c])]
        row = []
        for h in hd:
            n, f = bins.get((c, h), (0, 0))
            if n == 0: row.append("–"); continue
            us = np.linspace(h / 2, (h + 1) / 2, 11)
            pm = float(np.mean([max(tail(model.W_of(10 ** u) + cs * model.B_of(10 ** u)), 0.0) for u in us]))
            obs = f / n
            fmt = lambda v: ("%.3f" % v) if v >= 0.0995 else (("%.1e" % v) if v > 0 else "0")
            row.append(f"{fmt(obs)} / {fmt(pm) if pm > 1e-9 else '<1e-9'}")
        print(f"| ({q},{CA[c]}) | {cs} | " + " | ".join(row) + " |")
