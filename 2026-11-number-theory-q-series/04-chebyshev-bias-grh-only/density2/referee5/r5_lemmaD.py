"""Referee 5: independent checks of INVLO.md Lemma D (densities of random cosine sums).

AI-produced adversarial referee check, 2026-11-27.

(D1) |J_0(y)| <= sqrt(2/(pi y)) (Nicholson: y (J_0^2 + Y_0^2) increases to 2/pi), mpmath.
(D2) 0 < J_0(y) <= exp(-y^2/4) on [0, j_{0,1}), mpmath, including the endpoint region.
(D3) sup_{y > 2} |J_0(y)| < pi^{-1/2} (used in the proof of (b)).
(D4) Lemma D(a) with exactly three equal coefficients a = 1 (the tightest case):
     max density of cos t1 + cos t2 + cos t3 versus 6/pi^2, by Monte Carlo.
(D5) Lemma D(b) with few equal coefficients (theta = pi^{-1/2} regime) and with spread
     coefficients: true max density (Monte Carlo, max over x) versus the bound.
"""
import os

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
P = out.append
mp.mp.dps = 30

P("(D1) y*(J0^2+Y0^2) increasing to 2/pi, hence |J0(y)| <= sqrt(2/(pi y))")
ys = [mp.mpf(t) for t in ['0.001', '0.01', '0.1', '0.5', '1', '2', '5', '10', '100', '1000', '1e4', '1e5']]
vals = [y * (mp.besselj(0, y) ** 2 + mp.bessely(0, y) ** 2) for y in ys]
inc = all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
P(f"  samples increasing: {inc}; last value {mp.nstr(vals[-1], 12)} vs 2/pi = {mp.nstr(2 / mp.pi, 12)}")
worst = mp.mpf(0)
for y in np.concatenate([np.linspace(0.01, 50, 5001), np.linspace(50, 5000, 20001)]):
    y = mp.mpf(y)
    r = abs(mp.besselj(0, y)) / mp.sqrt(2 / (mp.pi * y))
    worst = max(worst, r)
P(f"  max |J0(y)|/sqrt(2/(pi y)) on [0.01, 5000] = {mp.nstr(worst, 12)} (must be <= 1)")

P("")
P("(D2) J0(y) <= exp(-y^2/4) on [0, j01)")
j01 = mp.besseljzero(0, 1)
worst = mp.mpf(1)
argworst = None
for t in np.linspace(1e-4, 1 - 1e-9, 4001):
    y = j01 * mp.mpf(t)
    d = mp.exp(-y * y / 4) - mp.besselj(0, y)
    if d < worst:
        worst, argworst = d, y
P(f"  j01 = {mp.nstr(j01, 10)}; min of exp(-y^2/4) - J0(y) = {mp.nstr(worst, 6)} at y = {mp.nstr(argworst, 6)}"
  f" (small-y expansion: difference ~ y^4/64 > 0)")
P(f"  J0(2) = {mp.nstr(mp.besselj(0, 2), 8)} <= exp(-1) = {mp.nstr(mp.exp(-1), 8)}")

P("")
P("(D3) sup_{y>2} |J0(y)| versus pi^{-1/2}")
j11 = mp.besseljzero(1, 1)  # first extremum of J0 after 0
P(f"  |J0(j11)| = {mp.nstr(abs(mp.besselj(0, j11)), 8)} at y = {mp.nstr(j11, 6)};"
  f" J0 is decreasing on [2, j11]; J0(2) = {mp.nstr(mp.besselj(0, 2), 6)}; pi^(-1/2) = {mp.nstr(1 / mp.sqrt(mp.pi), 6)}")

rng = np.random.default_rng(20261127)


def mc_max_density(coeffs, nsamp=200_000_000, chunk=10_000_000, binw=None):
    """max over x of the density of sum_k c_k cos(theta_k), by histogram."""
    coeffs = np.asarray(coeffs, dtype=float)
    span = coeffs.sum()
    if binw is None:
        binw = span / 1500
    edges = np.arange(-span - binw, span + 2 * binw, binw)
    h = np.zeros(len(edges) - 1)
    done = 0
    while done < nsamp:
        m = min(chunk, nsamp - done)
        s = np.zeros(m)
        for c in coeffs:
            s += c * np.cos(rng.uniform(0, 2 * np.pi, m))
        h += np.histogram(s, edges)[0]
        done += m
    dens = h / (nsamp * binw)
    i = dens.argmax()
    # crude relative statistical error at the maximum
    return dens[i], 0.5 * (edges[i] + edges[i + 1]), 1 / np.sqrt(h[i])


P("")
P("(D4) Lemma D(a), three equal coefficients a = 1")
d, x0, rel = mc_max_density([1, 1, 1], nsamp=400_000_000, binw=0.004)
P(f"  max density = {d:.4f} (at x ~ {x0:+.3f}, rel. stat. err ~{rel:.1e}); bound 6/pi^2 = {6 / np.pi ** 2:.4f};"
  f" ratio = {d / (6 / np.pi ** 2):.3f}")

P("")
P("(D5) Lemma D(b): bound (2pi)^-1 [2 sqrt(pi) Sigma^-1/2 + 12 theta^(n-3)/(pi a)] vs max density")


def boundDb(c):
    c = np.asarray(c, dtype=float)
    a, b = c.min(), c.max()
    theta = max(np.exp(-a * a / (b * b)), np.pi ** -0.5)
    Sig = np.sum(c ** 2)
    return (2 * np.sqrt(np.pi) / np.sqrt(Sig) + 12 * theta ** (len(c) - 3) / (np.pi * a)) / (2 * np.pi), theta


cases = [("n=3 equal", [1, 1, 1]), ("n=4 equal", [1] * 4), ("n=6 equal", [1] * 6),
         ("n=10 equal", [1] * 10), ("n=5, c in [0.5,1]", list(np.linspace(0.5, 1, 5))),
         ("n=30, c in [0.2,1]", list(rng.uniform(0.2, 1, 30)))]
for name, c in cases:
    bd, th = boundDb(c)
    d, x0, rel = mc_max_density(c, nsamp=100_000_000)
    P(f"  {name:20s}: theta = {th:.4f}; max density = {d:.4f} (x ~ {x0:+.2f}); bound = {bd:.4f}; ok: {d <= bd}")

open(os.path.join(HERE, "r5_lemmaD_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
