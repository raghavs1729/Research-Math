#!/usr/bin/env python3
"""compare_G.py sievelog zerofile -- G(x) = 12 sqrt(x) log x log f(x;13,1) from primes (psieve.c)
against the leading term of Theorem 1.18 under GRH: -S_1(log x) - 4 (zeros below 400)."""
import sys, re, numpy as np
rows = [l.split() for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
g = np.array([float(d) for _, d in rows]); w = 2 / (0.25 + g * g)
out = []
for l in open(sys.argv[1]):
    if l.startswith('CP'):
        d = dict(kv.split('=') for kv in l.split()[1:]); u = float(d['u']); x = float(d['x'])
        S = np.sum(w * np.cos(g * np.log(x)))
        out.append((u, float(d['G']), -S - 4))
a = np.array(out)
for lo in (12, 16, 18):
    m = a[:, 0] >= lo - 1e-9; diff = a[m, 1] - a[m, 2]
    print("u>=%d: n=%d G range [%.2f,%.2f]  G-(-S-4): mean %.3f rms %.3f max %.3f  corr %.4f" % (
        lo, m.sum(), a[m, 1].min(), a[m, 1].max(), diff.mean(), np.sqrt((diff ** 2).mean()), np.abs(diff).max(),
        np.corrcoef(a[m, 1], a[m, 2])[0, 1]))
    if lo == 18:
        dflip = a[m, 1] - (a[m, 2] + 2 * (4 + a[m, 2]) * 0 )
        # sign-flip control: G vs +S - 4 (wrong sign of the oscillating term)
        S_ = -(a[m, 2] + 4)
        alt = S_ - 4
        print("   control, wrong sign (+S-4): rms %.3f corr %.3f ; control, threshold 0 (-S): rms %.3f" % (
            np.sqrt(((a[m, 1] - alt) ** 2).mean()), np.corrcoef(a[m, 1], alt)[0, 1],
            np.sqrt(((a[m, 1] - (a[m, 2] + 4)) ** 2).mean())))
for r in a[::10]:
    print("u=%.2f G=%.4f -S-4=%.4f" % tuple(r))
