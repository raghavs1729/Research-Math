"""Referee 5: numerical checks of PRIMECORR (1.1) and Lemma P0' for chi_{-4}.

AI-produced adversarial referee check, 2026-11-27.
Written from scratch (own window, own quadrature, own sieve).

(X1) (1.1): Pi_T(y) = Omega_T(y) - sum_gamma 2 what(gamma/T) cos(gamma y), with a window
     v(u) = C exp(-1/(2(1/4 - u^2))) on (-1/2, 1/2) (not the author's bump), w = v*v.
(X2) Lemma P0', regularised form, for a P7-type plateau window beta_I: beta_I = 1 on
     [s1+g, s2-g], smooth steps on [s1, s1+g] and [s2-g, s2], s1 = 9, s2 = 20, g = 0.8.
     Margins (8.2, 9.8) and (19.2, 20.8) are zero-free, so (Z) holds; the band holds the
     four zeros 10.24, 12.99, 16.34, 18.29, all with beta_I = 1.
     For eps > 0 (h_eps = (beta e^{i.y}) * phi_eps; g_eps(u) = b(u-y) exp(-eps^2 u^2/2))
     the explicit formula is absolutely convergent and is checked directly.
(X3) The eps -> 0 limit: partial sums of sum chi(n)Lambda(n) n^{-1/2} b(y - log n) over
     log n <= y + V, against Omega_beta(y) - G_beta(iota(y)) - Pi^-_beta(y).
"""
import os

import numpy as np
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))
ZF = os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")
Z = np.array([float(l) for l in open(ZF) if l.strip()])
out = []
P = out.append
P(f"zeros: {len(Z)}, max ordinate {Z.max():.3f}; min gap {np.diff(Z).min():.4f}")


def chi_lambda(N):
    """array lam[n] = chi_{-4}(n) Lambda(n) for n <= N."""
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
        lp = np.log(p)
        while pk <= N:
            lam[pk] = c * lp
            pk *= int(p)
            c *= cp
    return lam


NMAX = int(np.exp(14.3))
lam = chi_lambda(NMAX)
idx = np.nonzero(lam)[0]
coef = lam[idx] / np.sqrt(idx)
logn = np.log(idx)
P(f"prime powers up to {NMAX}")


def Omega_chi(r):
    return np.log(4 / np.pi) + np.real(digamma((1.5 + 1j * np.asarray(r)) / 2))


# ---------------------------------------------------------------- (X1)
P("")
P("(X1) identity (1.1) with v(u) = C exp(-1/(2(1/4-u^2)))")
U = np.linspace(-0.5, 0.5, 40001)[1:-1]
du = U[1] - U[0]
vr = np.exp(-1.0 / (2 * (0.25 - U ** 2)))
vr /= vr.sum() * du


def vhat(xi):
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    return np.array([np.sum(vr * np.cos(x * U)) * du for x in xi])


# w = v * v on the grid of step du, support [-1, 1]
wgrid = np.convolve(vr, vr) * du
tw = (np.arange(len(wgrid)) - (len(wgrid) - 1) / 2) * du
P(f"  int w = {wgrid.sum() * du:.8f}; what(0) = {vhat(0)[0] ** 2:.8f}")
vz_cache = {}
for T in [3.0, 6.0]:
    wz = vhat(Z / T) ** 2
    P(f"  T={T}: what(gamma_max/T) = {wz[-1]:.2e}; c_0 = min_[0,1] what = {min(vhat(np.linspace(0, 1, 11)) ** 2):.4f}")
    rr = np.linspace(0, 60 * T, 60001)
    wr = vhat(rr / T) ** 2
    Om_r = wr * Omega_chi(rr)
    for y in [8.0, 10.0, 12.0, 14.0]:
        sel = np.abs(y - logn) <= 1 / T
        Pi = np.sum(coef[sel] * T * np.interp(T * (y - logn[sel]), tw, wgrid))
        zs = np.sum(2 * wz * np.cos(Z * y))
        Om = np.trapezoid(Om_r * np.cos(rr * y), rr) / np.pi
        P(f"    y={y:5.1f}: Pi_T = {Pi:+.6f}; Omega_T - zero sum = {Om - zs:+.6f}; diff = {Pi - (Om - zs):+.2e};"
          f" Omega_T = {Om:+.1e}")

# ---------------------------------------------------------------- (X2)
P("")
s1, s2, g = 9.0, 20.0, 0.8
band = Z[(Z > s1) & (Z < s2)]
marg = Z[((Z > s1 - g) & (Z < s1 + g)) | ((Z > s2 - g) & (Z < s2 + g))]
P(f"(X2) P7-type window: s1={s1}, s2={s2}, g={g}; zeros in margins: {marg.tolist()};"
  f" band zeros: {np.round(band, 3).tolist()}")


def step(x):
    """C^infty step: 0 for x<=0, 1 for x>=1."""
    x = np.asarray(x, dtype=float)
    a = np.where(x > 0, np.exp(-1 / np.maximum(x, 1e-300)), 0.0)
    b = np.where(x < 1, np.exp(-1 / np.maximum(1 - x, 1e-300)), 0.0)
    return a / (a + b)


def beta(r):
    r = np.abs(np.asarray(r, dtype=float))
    return step((r - s1) / g) * step((s2 - r) / g)


P(f"  beta at band zeros: {np.round(beta(band), 12).tolist()}")
rgrid = np.linspace(s1, s2, 44001)
dr = rgrid[1] - rgrid[0]
bet = beta(rgrid)


def bker(v):
    """b(v) = (1/pi) int_0^inf beta(r) cos(r v) dr (even, real)."""
    v = np.atleast_1d(np.asarray(v, dtype=float))
    res = np.empty(len(v))
    for i0 in range(0, len(v), 2000):
        vv = v[i0:i0 + 2000]
        res[i0:i0 + 2000] = (np.cos(np.outer(vv, rgrid)) @ bet) * dr / np.pi
    return res


vtab = np.linspace(-40, 40, 160001)
btab = bker(vtab)
P(f"  |b(v)| at v = 0, 2, 5, 10, 20, 30: " +
  ", ".join(f"{abs(np.interp(t, vtab, btab)):.2e}" for t in [0, 2, 5, 10, 20, 30]))


def b(v):
    return np.interp(v, vtab, btab, left=0.0, right=0.0)


for y in [7.0, 8.5]:
    P(f"  y = {y}")
    for eps in [0.5, 0.45, 0.4]:
        # prime side (absolutely convergent): sum chi Lambda n^-1/2 [g(log n) + g(-log n)]
        damp = np.exp(-eps ** 2 * logn ** 2 / 2)
        prime = np.sum(coef * (b(logn - y) + b(logn + y)) * damp)
        tailw = np.exp(logn[-1] / 2 - eps ** 2 * logn[-1] ** 2 / 2)
        # zero side: sum over all zeros of h_eps(gamma) = sum_{gamma>0} 2 Re h_eps(gamma)
        ss = np.linspace(-s2, s2, 88001)
        dss = ss[1] - ss[0]
        bs = beta(ss) * np.exp(1j * ss * y)

        def h_eps(r):
            return np.array([np.sum(bs * np.exp(-(x - ss) ** 2 / (2 * eps ** 2))) * dss / np.sqrt(2 * np.pi) / eps
                             for x in np.atleast_1d(r)])
        Zr = Z[Z < s2 + 12 * eps]
        zero = np.sum(2 * np.real(h_eps(Zr)))
        rr = np.linspace(-s2 - 12 * eps, s2 + 12 * eps, 8001)
        Om = np.trapezoid(np.real(h_eps(rr)) * Omega_chi(rr), rr) / (2 * np.pi)
        P(f"    eps={eps}: prime side = {prime:+.7f}; Omega side - zero side = {Om - zero:+.7f};"
          f" diff = {prime - (Om - zero):+.2e} (prime truncation weight e^(u/2 - eps^2 u^2/2) at u_max: {tailw:.1e})")

# ---------------------------------------------------------------- (X3)
P("")
P("(X3) eps -> 0: partial sums over log n <= y + V vs Omega_beta - G_beta - Pi^-_beta")
for y in [6.0, 7.0, 8.0]:
    Gb = np.sum(2 * beta(band) * np.cos(band * y))
    rr = np.linspace(s1, s2, 44001)
    Om = np.trapezoid(beta(rr) * Omega_chi(rr) * np.cos(rr * y), rr) / np.pi
    Pim = np.sum(coef * b(y + logn))
    target = Om - Gb - Pim
    row = []
    for V in [2.0, 4.0, 5.0, 6.0, 14.3 - y]:
        sel = logn <= y + V
        row.append((V, np.sum(coef[sel] * b(y - logn[sel]))))
    P(f"  y={y}: target = {target:+.6f} (G_beta = {Gb:+.5f}, Omega = {Om:+.1e}, Pi^- = {Pim:+.1e}); " +
      "; ".join(f"V={V:.1f}: {s:+.6f} (diff {s - target:+.1e})" for V, s in row))

open(os.path.join(HERE, "r5_explicit_formula_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
