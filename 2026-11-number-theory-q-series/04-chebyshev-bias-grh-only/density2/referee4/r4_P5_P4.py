"""Referee 4: checks for PRIMECORR.md Proposition P5 and Remark P4.

AI-produced, 2026-11-27, adversarial referee (referee 4).

P5:
 (P5a) spectrum Omega*N (Omega = 5, i.e. M = 0.1 rescaled to weight 1/2), m_n = floor(log(n+2)) + 1,
       random perturbations with |gamma'-gamma| < delta*min(1, gamma/2): counting-law constant
       (half-open unit intervals) before/after, and the sigma^2 ratio versus
       [(1+delta)^-2, (1-delta)^-2].
 (P5b) the factor 2 in the counting law is attained: spectrum {k + 0.6}, one shift of 0.45.
 (P5c) absolute continuity: the characteristic function of an LI model is bounded by three factors
       (numerical integral of the 3-factor bound), so the density is bounded.
P4:
 (P4a) the ratio (X/T)^{k/2} X^{-1/2+eps} (log T)^{-k/2} for sample values.
 (P4b) Sum_rho what(gamma/T)^2 versus T log T for chi_{-4} (the "Sum a_rho^2 ~ T log T" of §3).
 (P4c) the relation lengths in INVLO Prop. E's range: L >= C_2 T / (alpha^2 delta log L).
"""
import os
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0

out = []
P = out.append
rng = np.random.default_rng(11)


def counting_constant(g, m, shift=3.0):
    """max over t >= 0 of (sum_{gamma in (t,t+1]} m) / log(t + shift); sup attained with t just below a point."""
    order = np.argsort(g)
    g = g[order]
    m = m[order]
    best = 0.0
    # for half-open windows (t, t+1], the maximal sums occur for t = gamma_i - 1 (window (gamma_i - 1, gamma_i])
    cs = np.concatenate([[0], np.cumsum(m)])
    for i in range(len(g)):
        t = max(g[i] - 1.0, 0.0)
        lo = np.searchsorted(g, t, side='right')
        hi = np.searchsorted(g, t + 1.0, side='right')
        s = cs[hi] - cs[lo]
        best = max(best, s / math.log(t + shift))
    return best


def sigma2(g, m, M=0.5):
    return float(np.sum(m * 8 * M * M / (M * M + g * g)))


P("(P5a) perturbation of Omega*N, Omega = 5")
n = np.arange(1, 2001)
g = 5.0 * n
m = np.floor(np.log(n + 2)) + 1
A0 = counting_constant(g, m)
s0 = sigma2(g, m)
for delta in [0.1, 0.25, 0.49]:
    worstA, rmin, rmax = 0.0, 10.0, 0.0
    for trial in range(20):
        sh = (rng.uniform(-1, 1, len(g)) * delta * np.minimum(1.0, g / 2))
        gp = g + sh
        worstA = max(worstA, counting_constant(gp, m, shift=3.5))
        r = sigma2(gp, m) / s0
        rmin, rmax = min(rmin, r), max(rmax, r)
    P(f"  delta={delta}: counting constant before {A0:.3f}, after (log(t+3.5)) <= {worstA:.3f} (claim <= 2A = {2*A0:.3f});"
      f" sigma^2 ratio in [{rmin:.4f}, {rmax:.4f}] within [{(1+delta)**-2:.4f}, {(1-delta)**-2:.4f}]: "
      f"{(1+delta)**-2 <= rmin and rmax <= (1-delta)**-2}")
P(f"  sigma^2 of the unperturbed stand-in spectrum at weight 1/2: {s0:.4f}")

P("")
P("(P5b) factor 2 attained")
g = np.arange(0, 40) + 0.6
m = np.ones_like(g)
gp = g.copy()
gp[10] += 0.45  # 10.6 -> 11.05 (allowed: |shift| < delta*min(1, gamma/2) with delta = 0.49 -> < 0.49)
c_before = max(np.sum((g > t) & (g <= t + 1)) for t in np.linspace(0, 39, 3901))
c_after = max(np.sum((gp > t) & (gp <= t + 1)) for t in np.linspace(0, 39, 3901))
P(f"  max count in a unit interval: before {c_before}, after {c_after}: the doubling in P5 is sharp up to the"
  f" additive constant in the logarithm.")

P("")
P("(P5c) three-factor bound for the characteristic function")
for a in [0.05, 0.2, 1.0]:
    xi = np.linspace(0, 2e5 / a, 4000001)
    b3 = np.minimum(1.0, (2 / (np.pi * xi[1:] * a)) ** 1.5)
    I = 2 * np.trapezoid(np.concatenate([[1.0], b3]), xi)
    P(f"  a={a}: 2*int_0^X min(1,(2/(pi xi a))^1.5) = {I:.4f} (-> 12/(pi a) = {12/(np.pi*a):.4f}); "
      f"|J0(xi a)|^3 <= bound on grid: {bool(np.all(np.abs(j0(xi[1:200001]*a))**3 <= b3[:200000] + 1e-15))}")

P("")
P("(P4a) Hardy-Littlewood error/main ratio (X/T)^{k/2} X^{-1/2+eps} (log T)^{-k/2}, eps = 0.01, T = 10")
for k in [2, 3, 4]:
    row = []
    for X in [1e6, 1e9, 1e12]:
        row.append(f"X={X:.0e}: {(X/10)**(k/2) * X**(-0.5+0.01) * math.log(10)**(-k/2):.2e}")
    P(f"  k={k}: " + "; ".join(row))
P("  needed error per shift tuple for the main term to dominate: X^{1-k/2} T^{k/2}(log T)^{k/2}:"
  " O(1) for k=2, O(X^{-1/2}) for k=3, O(X^{-1}) for k=4 (so 'o(X^{-1/2}) for k>=3' should read 'O(X^{1-k/2})').")

P("")
P("(P4b) Sum_rho what(gamma/T)^2 versus T log T (chi_{-4}, window v = exp(-1/(1-4u^2)))")
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")) if l.strip()])
xg, wg = leggauss(3000)
ug, wgu = 0.5 * xg, 0.5 * wg
vg = np.exp(-1.0 / (1.0 - 4.0 * ug ** 2))
VN = np.sum(wgu * vg)


def what(xi):
    return ((np.cos(np.outer(np.atleast_1d(xi), ug)) @ (wgu * vg)) / VN) ** 2


for T in [5, 10, 20, 40, 60]:
    s2 = np.sum(what(Z / T) ** 2)
    P(f"  T={T:3d}: s_2(T) = {s2:8.3f}; T log T = {T*math.log(T):8.2f}; ratio {s2/(T*math.log(T)):.4f};"
      f" what(640/T)^2 = {what(np.array([640.0/T]))[0]**2:.1e}")

P("")
P("(P4c) relation lengths in INVLO Prop. E's range (illustrative C_2 = 1)")
for alpha, delta in [(0.0891, 0.1), (0.0891, 0.01)]:
    for T in [1e2, 1e3, 1e4]:
        # smallest L with T <= alpha^2 delta L log L
        L = 3.0
        while alpha ** 2 * delta * L * math.log(L) < T:
            L *= 1.01
        P(f"  alpha={alpha}, delta={delta}, block height T={T:.0e}: Prop E speaks only for L >= {L:.2e};"
          f" the relations may have any length <= L, i.e. far beyond any fixed moment order")

P("")
P("(P4d) H_corr at l = 2 and a double zero above T: Delta_2 >= 4 what(gamma_0/T)^2 versus c_0^2")
for sharp in [1.0, 0.5]:
    vg2 = np.exp(-sharp / (1.0 - 4.0 * ug ** 2))
    VN2 = np.sum(wgu * vg2)
    wh = lambda xi: ((np.cos(np.outer(np.atleast_1d(xi), ug)) @ (wgu * vg2)) / VN2) ** 2
    xi = np.linspace(0, 1, 2001)
    c0 = wh(xi).min()
    xs = np.linspace(1, np.pi, 4001)
    viol = xs[4 * wh(xs) ** 2 >= c0 ** 2]
    P(f"  bump exp(-{sharp}/(1-4u^2)): c_0 = min_[0,1] what = {c0:.4f}; what(2) = {wh(np.array([2.0]))[0]:.4f};"
      f" a double zero at gamma_0 = xi*T violates H_corr(2,T) for all xi in [1, {viol.max():.3f}]")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
