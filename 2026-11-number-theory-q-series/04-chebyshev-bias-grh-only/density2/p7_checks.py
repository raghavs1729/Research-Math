"""Numerical check of Lemma P0' (PRIMECORR.md, band decoupling) for chi_{-4}:
for a window beta in C_c^infty supported on a frequency band [s1,s2] (even extension),

    Pi_beta(y) = lim_X sum_{n<=X} chi(n)Lambda(n) n^{-1/2} b(y - log n)
               = Omega_beta(y) - sum_{gamma in band} 2 beta(gamma) cos(gamma y) - Pi^-_beta(y),

with b(v) = (1/pi) int_0^inf beta(r) cos(r v) dr (Schwartz, not compactly supported).
The prime sum is truncated at log n <= y + V; we report the effect of V.
AI-produced (Anthropic's Claude), 2026-11-27.
"""
import os
import numpy as np
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, "..", "..", "prototype", "explore",
                                                  "chi4_zeros.txt")) if l.strip()])


def bump01(x):
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    out[m] = np.exp(-1.0 / (x[m] * (1 - x[m])))
    return out


def beta_band(r, s1, s2):
    """smooth window: 0 outside [s1,s2], plateau-free bump scaled to max 1 (even in r)."""
    r = np.abs(np.asarray(r, dtype=float))
    x = (r - s1) / (s2 - s1)
    return bump01(x) / np.exp(-4.0)


def b_kernel(v, s1, s2, nr=20001):
    rr = np.linspace(s1, s2, nr)
    bb = beta_band(rr, s1, s2)
    v = np.atleast_1d(v)
    # (1/2pi) int_{-inf}^{inf} beta(r) e^{-irv} dr = (1/pi) int_0^inf beta(r) cos(rv) dr
    return np.array([np.trapezoid(bb * np.cos(rr * x), rr) / np.pi for x in v])


def lambda_chi(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    lam = np.zeros(N + 1)
    for p in np.nonzero(sieve)[0]:
        if p == 2:
            continue
        chi_p = 1.0 if p % 4 == 1 else -1.0
        pk, c = p, chi_p
        while pk <= N:
            lam[pk] = c * np.log(p)
            pk *= p
            c *= chi_p
    return lam


out = []
P = out.append
s1, s2 = 9.0, 21.0
inband = Z[(Z > s1) & (Z < s2)]
P(f"band [{s1},{s2}]: ordinates inside: {np.round(inband, 3).tolist()}")
y_list = [6.0, 7.5, 9.0]
Vmax = 4.0
N = int(np.exp(max(y_list) + Vmax)) + 10
lam = lambda_chi(N)
idx = np.nonzero(lam)[0]
logn = np.log(idx)
coef = lam[idx] * idx ** -0.5


def Omega(r):
    return np.log(4 / np.pi) + np.real(digamma((1.5 + 1j * r) / 2))


for y in y_list:
    zero_side = -np.sum(2 * beta_band(inband, s1, s2) * np.cos(inband * y))
    rr = np.linspace(s1, s2, 20001)
    Om = 2 * np.trapezoid(beta_band(rr, s1, s2) * Omega(rr) * np.cos(rr * y), rr) / (2 * np.pi)
    for V in [2.0, 3.0, 4.0]:
        sel = logn <= y + V
        # kernel on a fine grid in v = y - log n
        vv = y - logn[sel]
        grid = np.linspace(vv.min() - 1e-9, vv.max() + 1e-9, 6001)
        kg = b_kernel(grid, s1, s2)
        kv = np.interp(vv, grid, kg)
        Pi = np.sum(coef[sel] * kv)
        P(f"  y={y:4.1f} V={V}: prime sum = {Pi:+.5f}; Omega - zero sum = {Om + zero_side:+.5f};"
          f" diff = {Pi - (Om + zero_side):+.2e}  (Omega = {Om:+.1e})")

open(os.path.join(HERE, "p7_checks_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
