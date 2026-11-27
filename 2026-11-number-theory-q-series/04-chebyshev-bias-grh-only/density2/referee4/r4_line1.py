"""Referee 4: numerics for LINE1.md §4 (Proposition 4.1, Lemma 4.2).

AI-produced, 2026-11-27, adversarial referee (referee 4).

(L1) Prop 4.1: concentration of a cos(phi) (phi uniform) on intervals of length 4 eps, a >= 8 eps:
     sup = arccos(1 - 4eps/a)/pi <= 1/3.
(L2) Lemma 4.2 on the certified zeros of L(s, chi_{-4}) (all ordinates up to 640), U = 600,
     hhat(r) = Phi(r/U), Phi = 1 on [-1/2,1/2], smooth step to 0 at |s| = 1:
     ratio(y) = |S(y)| / N with S two-sided (2 sum cos, N = 2 sum hhat) and one-sided (sum e^{i gamma y}),
     real and imaginary parts separately, for y in [K/U, 2 log U + 2 log log U].  The explicit formula
     controls only the two-sided (cosine) sum.
(L3) the prime spikes for y <= log U: predicted size of the two-sided ratio at y = log n,
     2 pi Lambda(n) n^{-1/2} phi(0) / (2 N/U) with phi(0) = (1/2pi) int Phi, compared with the numerics.
"""
import os
import math
import numpy as np

out = []
P = out.append

P("(L1) arcsine concentration")
for ratio in [8, 10, 20, 100]:
    a = 1.0
    eps = a / ratio
    l = 4 * eps
    # by symmetry and monotonicity of the arcsine density in |x|, the sup is at an endpoint
    sup = math.acos(1 - l / a) / math.pi
    # brute force over interval positions
    phi = np.linspace(0, 2 * np.pi, 2000001)[:-1]
    xs = a * np.cos(phi)
    starts = np.linspace(-a, a - l, 2001)
    xs_sorted = np.sort(xs)
    cnt = (np.searchsorted(xs_sorted, starts + l, side='right') - np.searchsorted(xs_sorted, starts, side='left')) / len(xs)
    P(f"  a/eps = {ratio}: sup_x P(a cos phi in [x, x+4eps]) = {cnt.max():.5f} (formula {sup:.5f}) <= 1/3: {cnt.max() <= 1/3 + 1e-5}")

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")) if l.strip()])
U = 600.0


def step(t):
    t = np.clip(t, 0, 1)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        a = np.where(t > 0, np.exp(-1 / np.where(t > 0, t, 1)), 0.0)
        b = np.where(t < 1, np.exp(-1 / np.where(t < 1, 1 - t, 1)), 0.0)
    return a / (a + b)


def Phi(s):
    s = np.abs(s)
    return np.where(s <= 0.5, 1.0, np.where(s >= 1, 0.0, step((1 - s) / 0.5)))


hz = Phi(Z / U)
N1 = hz.sum()
P("")
P(f"(L2) Lemma 4.2 numerics, U = {U:.0f}: one-sided N = sum_gamma hhat(gamma) = {N1:.2f}")
ymax = 2 * math.log(U) + 2 * math.log(math.log(U))
P(f"  2 log U + 2 log log U = {ymax:.3f}")
ys = np.linspace(1e-3, ymax + 2, 120001)
S2 = np.empty_like(ys)
ReS1 = np.empty_like(ys)
ImS1 = np.empty_like(ys)
for i in range(0, len(ys), 2000):
    Y = ys[i:i + 2000]
    C = np.cos(np.outer(Y, Z)) @ hz
    S = np.sin(np.outer(Y, Z)) @ hz
    S2[i:i + 2000] = 2 * C
    ReS1[i:i + 2000] = C
    ImS1[i:i + 2000] = S
r2 = np.abs(S2) / (2 * N1)
rre = np.abs(ReS1) / N1
rim = np.abs(ImS1) / N1
for K, Cc in [(10, 4), (30, 6), (100, 8)]:
    m = (ys >= K / U) & (ys <= ymax - Cc)
    P(f"  K={K:3d}, C={Cc}: y in [{K/U:.3f}, {ymax-Cc:.2f}]: max two-sided ratio {r2[m].max():.3f};"
      f" max |Im one-sided|/N {rim[m].max():.3f}; max |one-sided|/N {np.sqrt(rre[m]**2 + rim[m]**2).max():.3f}")
for lo, hi in [(0.02, 0.5), (0.5, math.log(U)), (math.log(U), 2 * math.log(U)), (2 * math.log(U), ymax), (ymax, ymax + 2)]:
    m = (ys >= lo) & (ys <= hi)
    j = np.argmax(r2[m])
    P(f"  y in [{lo:.2f},{hi:.2f}]: max two-sided ratio {r2[m].max():.3f} at y = {ys[m][j]:.3f}"
      f" (e^y = {math.exp(ys[m][j]):.1f}); max |Im S1|/N {rim[m].max():.3f}")

P("")
P("(L3) prime spikes for small y (two-sided ratio at y = log n)")
phi0 = (1 / (2 * np.pi)) * np.trapezoid(Phi(np.linspace(-1, 1, 200001)), np.linspace(-1, 1, 200001))
for nn, lamchi in [(3, -math.log(3)), (5, math.log(5)), (7, -math.log(7)), (9, math.log(3)), (13, math.log(13))]:
    y0 = math.log(nn)
    j = np.argmin(np.abs(ys - y0))
    pred = abs(lamchi / math.sqrt(nn) * U * phi0) / (2 * N1)
    P(f"  n={nn:2d}: ratio at y=log n: {r2[j]:.3f}; explicit-formula prediction |Lambda chi n^-1/2 U phi(0)|/(2N) = {pred:.3f}"
      f"  (O(1/log U): the o_U(1) term of Lemma 4.2 decays slowly)")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
