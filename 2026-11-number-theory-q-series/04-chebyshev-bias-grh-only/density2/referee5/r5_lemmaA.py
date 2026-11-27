"""Referee 5: independent checks of INVLO.md Lemma A (needle polynomial and constants).

AI-produced adversarial referee check, 2026-11-27.

Independent construction (not reusing invlo_checks.py):
  G(phi) = F_N(phi - pi/2) + F_N(phi + pi/2),  F_N = Fejer kernel,
  g(cos phi) = G(phi),  q = g/g(0).
(A1) exact rational value of int_{-1}^{1} g^2 ds (Chebyshev expansion), compared with
     the claimed bound 8 pi (N+1) (i.e. int q^2 <= 8 pi/(N+1));
(A2) g(0) in [N+1, N+2];
(A3) min of q on |s| <= 2/(N+1), claimed >= 8/(3 pi^2) > 0.27;
(A4) |q(s)| <= 2 (2|s|)^N for |s| >= 1, evaluated in multiprecision (no overflow);
(A5) 0 <= q <= 2 on [-1,1];
(A6) Gauss-Hermite sharpness constant: central weight * sqrt(N) -> sqrt(pi/2) = 1.2533;
(A7) Gaussian W at x = mu: the Lemma A bound drops below 1 only for N > ~2.0e7.
"""
import os
from fractions import Fraction

import mpmath as mp
import numpy as np
from numpy.polynomial.hermite_e import hermegauss

HERE = os.path.dirname(os.path.abspath(__file__))
out = []
P = out.append


def cheb_coeffs(N):
    """Exact Chebyshev coefficients of g: g = sum_k a_k T_k.
    G(phi) = sum_{|k|<=N} w_k (e^{ik(phi-pi/2)} + e^{ik(phi+pi/2)}), w_k = 1-|k|/(N+1)
           = 2 w_0 + sum_{k>=1} 4 w_k cos(k pi/2) cos(k phi)."""
    a = []
    for k in range(N + 1):
        w = Fraction(N + 1 - k, N + 1)
        ck = [1, 0, -1, 0][k % 4]  # cos(k pi/2)
        a.append(2 * w if k == 0 else 4 * w * ck)
    return a


def int_Tn(n):
    """int_{-1}^{1} T_n(s) ds."""
    return Fraction(0) if n % 2 else Fraction(2, 1 - n * n)


def int_g2(a):
    """exact int_{-1}^1 (sum a_k T_k)^2 ds, using T_j T_k = (T_{j+k} + T_{|j-k|})/2."""
    n = len(a)
    tot = Fraction(0)
    nz = [k for k in range(n) if a[k] != 0]
    for j in nz:
        for k in nz:
            tot += a[j] * a[k] * (int_Tn(j + k) + int_Tn(abs(j - k))) / 2
    return tot


def g_at(a, s):
    """evaluate sum a_k T_k(s) in mpmath (s may be outside [-1,1])."""
    s = mp.mpf(s)
    # Clenshaw recurrence
    b1 = mp.mpf(0)
    b2 = mp.mpf(0)
    for k in range(len(a) - 1, 0, -1):
        b1, b2 = 2 * s * b1 - b2 + mp.mpf(a[k].numerator) / a[k].denominator, b1
    return s * b1 - b2 + mp.mpf(a[0].numerator) / a[0].denominator


def G_closed(phi, N):
    def F(t):
        t = np.asarray(t, dtype=float)
        sn = np.sin(t / 2)
        with np.errstate(divide='ignore', invalid='ignore'):
            v = (np.sin((N + 1) * t / 2) / sn) ** 2 / (N + 1)
        v = np.where(np.abs(sn) < 1e-14, N + 1.0, v)
        return v
    return F(phi - np.pi / 2) + F(phi + np.pi / 2)


P("(A1)-(A5) needle polynomial, exact / multiprecision")
mp.mp.dps = 60
worst_A4 = 0.0
for N in [1, 2, 3, 4, 5, 7, 10, 20, 50, 100, 200, 300, 400]:
    a = cheb_coeffs(N)
    g0 = sum(a[k] * [1, 0, -1, 0][k % 4] for k in range(len(a)))  # T_k(0) = cos(k pi/2)
    I = int_g2(a)
    ratio_int = float(I / g0 ** 2) / (8 * np.pi / (N + 1))
    # A3: min on |s| <= 2/(N+1): s = cos(phi), phi near pi/2
    s0 = 2.0 / (N + 1)
    phis = np.arccos(np.linspace(-s0, s0, 20001))
    qin = G_closed(phis, N) / float(g0)
    # A5: on the whole interval
    phall = np.linspace(0, np.pi, 200001)
    qall = G_closed(phall, N) / float(g0)
    # A4: outside, |s| in [1, 1e3], log grid, both signs
    worst = 0.0
    for s in list(np.logspace(0, 3, 121)) + [1.0000001, 1.001, 1.01]:
        for sg in (1, -1):
            val = abs(g_at(a, sg * s)) / (mp.mpf(g0.numerator) / g0.denominator)
            bound = 2 * (2 * mp.mpf(s)) ** N
            worst = max(worst, float(val / bound))
    worst_A4 = max(worst_A4, worst)
    P(f"  N={N:4d}: g(0)-(N+1)={float(g0) - (N + 1):.4f} in [0,1]: {0 <= float(g0) - (N + 1) <= 1};"
      f" int q^2 / (8pi/(N+1)) = {ratio_int:.4f};  (N+1)*int q^2 = {float(I / g0 ** 2) * (N + 1):.4f};"
      f" min q on |s|<=2/(N+1) = {qin.min():.4f} (>= {8 / (3 * np.pi ** 2):.4f});"
      f" q range on [-1,1] = [{qall.min():.2e}, {qall.max():.4f}];"
      f" max |q|/(2(2|s|)^N) on 1<=|s|<=1e3 = {worst:.3e}")
P(f"  overall worst ratio in (A4): {worst_A4:.3e} (must be <= 1)")

P("")
P("(A6) Gauss-Hermite sharpness: (N+1)-point rule for N(0,1), central weight")
def gauss_hermite_prob(n):
    """Golub-Welsch for the probabilists' Hermite weight (stable for large n)."""
    off = np.sqrt(np.arange(1, n))
    J = np.diag(off, 1) + np.diag(off, -1)
    x, V = np.linalg.eigh(J)
    return x, V[0, :] ** 2


for N in [10, 40, 160, 640, 1280]:
    x, w = gauss_hermite_prob(N + 1)
    w = w / w.sum()
    P(f"  N={N:4d}: max weight = {w.max():.5f}; max weight*sqrt(N) = {w.max() * np.sqrt(N):.4f};"
      f" *sqrt(N+1) = {w.max() * np.sqrt(N + 1):.4f}  (limit sqrt(pi/2) = {np.sqrt(np.pi / 2):.4f})")

P("")
P("(A7) Lemma A bound for a standard Gaussian W at x = mu (b = 0, sigma = 1):")


def boundA(N):
    R = 16 * 2 * np.sqrt(N)
    return 352 * (1 / np.sqrt(2 * np.pi)) * R / (N + 1) + 56 * 16.0 ** (-N)


lo, hi = 1.0, 1e9
for _ in range(200):
    mid = np.sqrt(lo * hi)
    if boundA(mid) > 1:
        lo = mid
    else:
        hi = mid
P(f"  bound(N) = 1 at N = {hi:.4e}; bound(1e7) = {boundA(1e7):.3f}, bound(2.1e7) = {boundA(2.1e7):.4f}")

open(os.path.join(HERE, "r5_lemmaA_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
