"""Referee 5: size of the term Pi^-_beta(y) of Lemma P0' (PRIMECORR Sec. 3).

AI-produced adversarial referee check, 2026-11-27.

Pi^-_beta(y) = sum chi(n) Lambda(n) n^{-1/2} b(y + log n), chi = chi_{-4}.
(M1) the plateau-free bump on [9,21] used by p7_checks.py (file: "negligible for y >= 6");
(M2) a P7-type window (plateau [9.8,19.2], steps of width g = 0.8), as in r5_explicit_formula.py.
Only y + log n >= y + log 3 enters (chi(2) = 0), so |b| beyond v ~ y+1 matters.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
P = out.append


def chi_lambda(N):
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    lam = np.zeros(N + 1)
    for p in np.nonzero(s)[0]:
        if p == 2:
            continue
        cp = 1.0 if p % 4 == 1 else -1.0
        pk, c = int(p), cp
        while pk <= N:
            lam[pk] = c * np.log(p)
            pk *= int(p)
            c *= cp
    return lam


N = int(np.exp(13.0))
lam = chi_lambda(N)
idx = np.nonzero(lam)[0]
coef = lam[idx] / np.sqrt(idx)
logn = np.log(idx)


def step(x):
    x = np.asarray(x, dtype=float)
    a = np.where(x > 0, np.exp(-1 / np.maximum(x, 1e-300)), 0.0)
    b = np.where(x < 1, np.exp(-1 / np.maximum(1 - x, 1e-300)), 0.0)
    return a / (a + b)


def bump01(x):
    o = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    o[m] = np.exp(-1.0 / (x[m] * (1 - x[m])))
    return o


wins = {
    "p7_checks bump on [9,21]": (lambda r: bump01((np.abs(r) - 9.0) / 12.0) / np.exp(-4.0), 9.0, 21.0),
    "P7-type, plateau [9.8,19.2], g=0.8": (lambda r: step((np.abs(r) - 9.0) / 0.8) * step((20.0 - np.abs(r)) / 0.8), 9.0, 20.0),
}
for name, (beta, lo, hi) in wins.items():
    rg = np.linspace(lo, hi, 40001)
    dr = rg[1] - rg[0]
    bt = beta(rg)
    vt = np.linspace(0, 40, 40001)
    bv = np.array([np.sum(bt * np.cos(rg * v)) * dr / np.pi for v in vt])
    b = lambda v: np.interp(np.abs(v), vt, bv, right=0.0)
    vals = []
    for y in [6.0, 7.5, 9.0]:
        sel = logn <= 40 - y
        vals.append((y, np.sum(coef[sel] * b(y + logn[sel]))))
    P(f"  {name}: max |b(v)| for v >= 7: {np.abs(bv[vt >= 7]).max():.2e}; Pi^-(y): " +
      ", ".join(f"y={y}: {s:+.2e}" for y, s in vals))

open(os.path.join(HERE, "r5_piminus_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
