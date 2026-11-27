"""Referee 4b: numerical checks of PRIMECORR.md Lemma P0' and of the empty-band case of Theorem P7.

AI-produced, 2026-11-27, adversarial referee (referee 4, second read).

chi = chi_{-4}; certified zeros up to 640.  Windows of the Theorem P7 type (plateau = 1, smooth
e^{-1/t} transitions of width g), even in r; kernel b(v) = (1/pi) int_0^inf beta(r) cos(rv) dr.
 (R1) the Gaussian-regularised explicit formula used in the proof, for eps = 0.5, 0.7, 1.0
      (absolutely convergent):
        sum_n chi(n)Lambda(n)n^{-1/2} [g_eps(log n) + g_eps(-log n)]
          = (1/2pi) int h_eps Omega_chi - sum_rho h_eps(gamma_rho),
      g_eps(u) = b(u - y) e^{-eps^2 u^2/2},  h_eps = (beta e^{i . y}) * phi_eps.
 (R2) the limit identity of Lemma P0' (both prime terms, b(y - log n) and b(y + log n)) with
      a smooth cutoff in u = log n at several heights:
        sum_n chi Lambda n^{-1/2} [b(y - log n) + b(y + log n)] = Omega_beta(y) - G_beta(iota(y)).
 (R3) the one-level count n_I at y = 0: n_I = (Omega_beta(0) - sum_n 2 chi Lambda n^{-1/2} b(log n))/2.
 (R4) the EMPTY band (s1, s2, g) = (1, 2.5, 0.4): (Z) holds, no zero in [0.6, 2.9]; the band sum
      Pi_I(y) is ~ Omega_beta(y) -> 0, so M_{a,b} = 0 for b >= 1 and (D) holds for every L.
"""
import os
import math
import numpy as np
from scipy.special import digamma

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")) if l.strip()])
out = []
P = out.append


def S(t):
    t = np.clip(t, 0, 1)
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        a = np.where(t > 0, np.exp(-1 / np.where(t > 0, t, 1)), 0.0)
        b = np.where(t < 1, np.exp(-1 / np.where(t < 1, 1 - t, 1)), 0.0)
    return a / (a + b)


def beta(r, s1, s2, g):
    x = np.abs(np.asarray(r, dtype=float))
    return np.where(x <= s1, 0.0, np.where(x < s1 + g, S((x - s1) / g),
                    np.where(x <= s2 - g, 1.0, np.where(x < s2, S((s2 - x) / g), 0.0))))


def Omega_chi(r):
    return np.log(4 / np.pi) + np.real(digamma((1.5 + 1j * r) / 2))


def lam_chi(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    lam = np.zeros(N + 1)
    for p in np.nonzero(sieve)[0]:
        if p == 2:
            continue
        cp = 1.0 if p % 4 == 1 else -1.0
        pk, c = int(p), cp
        while pk <= N:
            lam[pk] = c * math.log(p)
            pk *= int(p)
            c *= cp
    return lam


UMAX = 16.0
NMAX = int(math.exp(UMAX)) + 10
lam = lam_chi(NMAX)
nz = np.nonzero(lam)[0]
un = np.log(nz)
cn = lam[nz] / np.sqrt(nz)
P(f"prime powers with chi != 0 up to e^{UMAX:.0f} = {NMAX}: {len(nz)}")


def kernel_grid(s1, s2, g, vmin, vmax, step=2.5e-4, nr=6001):
    r = np.linspace(s1, s2, nr)
    br = beta(r, s1, s2, g)
    wr = np.full(nr, r[1] - r[0])
    wr[0] *= 0.5
    wr[-1] *= 0.5
    v = np.arange(vmin, vmax + step, step)
    bv = np.empty_like(v)
    for i in range(0, len(v), 4000):
        bv[i:i + 4000] = (np.cos(np.outer(v[i:i + 4000], r)) @ (br * wr)) / np.pi
    return v, bv


def smooth_cut(u, U0, W):
    return 1.0 - S((u - U0) / W)


def check_band(name, s1, s2, g, ylist, cutoffs, eps_list):
    marg = ((Z >= s1 - g) & (Z <= s1 + g)) | ((Z >= s2 - g) & (Z <= s2 + g))
    band = Z[(Z > s1) & (Z < s2)]
    P(f"{name}: band ({s1}, {s2}), margin g = {g}; zeros in margins: {Z[marg].tolist()} -> (Z) holds: {not marg.any()};"
      f" band zeros: {np.round(band, 4).tolist()} (n_I = {len(band)})")
    vmax = UMAX + max(ylist) + 1
    v, bv = kernel_grid(s1, s2, g, -vmax, vmax)
    bfun = lambda x: np.interp(x, v, bv)
    rr = np.linspace(0, s2, 40001)
    br = beta(rr, s1, s2, g)
    for y in ylist:
        Om = np.trapezoid(br * np.cos(rr * y) * Omega_chi(rr), rr) / np.pi
        Gb = np.sum(2 * beta(Z, s1, s2, g) * np.cos(Z * y))
        target = Om - Gb
        # (R2) limit identity, smooth cutoffs
        vals = []
        for U0, W in cutoffs:
            w = smooth_cut(un, U0, W)
            s = np.sum(cn * (bfun(y - un) + bfun(y + un)) * w)
            vals.append((U0, W, s))
        P(f"  y = {y}: Omega_beta = {Om:+.3e}; target Omega - G_beta = {target:+.6f}; smooth-cutoff prime sums: "
          + "; ".join(f"cut [{U0:.0f},{U0+W:.0f}]: {s:+.6f} (diff {s-target:+.1e})" for U0, W, s in vals))
        # (R1) regularised identity
        for eps in eps_list:
            ge = lambda u: bfun(u - y) * np.exp(-eps ** 2 * u ** 2 / 2)
            prime = np.sum(cn * (ge(un) + ge(-un)))
            tail_bound = math.exp(UMAX / 2 - eps ** 2 * UMAX ** 2 / 2)
            # h_eps(gamma) = int beta(s) e^{isy} phi_eps(gamma - s) ds  (complex)
            sgrid = np.linspace(-s2, s2, 40001)
            bs = beta(sgrid, s1, s2, g) * np.exp(1j * sgrid * y)
            ds = sgrid[1] - sgrid[0]

            def h_eps(x):
                x = np.atleast_1d(x)
                res = np.empty(len(x), dtype=complex)
                for i in range(0, len(x), 200):
                    X = x[i:i + 200]
                    phi = np.exp(-(X[:, None] - sgrid[None, :]) ** 2 / (2 * eps ** 2)) / math.sqrt(2 * math.pi * eps ** 2)
                    res[i:i + 200] = (phi * bs[None, :]).sum(axis=1) * ds
                return res
            hz = h_eps(Z)
            zero_side = 2 * np.real(hz.sum())
            rgrid = np.linspace(-(s2 + 12 * eps), s2 + 12 * eps, 8001)
            hr = h_eps(rgrid)
            om_eps = np.real(np.trapezoid(hr * Omega_chi(rgrid), rgrid)) / (2 * math.pi)
            P(f"    (R1) eps = {eps}: prime side {prime:+.8f}; Omega - zeros {om_eps - zero_side:+.8f};"
              f" diff {prime - (om_eps - zero_side):+.1e} (prime tail factor at u={UMAX:.0f}: {tail_bound:.0e})")
    # (R3) one-level count at y = 0
    Om0 = np.trapezoid(br * Omega_chi(rr), rr) / np.pi
    for U0, W in cutoffs:
        w = smooth_cut(un, U0, W)
        s0 = np.sum(cn * 2 * bfun(un) * w)
        P(f"  (R3) y = 0, cut [{U0:.0f},{U0+W:.0f}]: (Omega(0) - prime sum)/2 = {(Om0 - s0)/2:.5f} (n_I = {len(band)})")


check_band("P7-type band", 8.13, 19.87, 1.4, [5.0, 6.0], [(12.0, 3.0), (13.0, 3.0)], [0.5, 0.7, 1.0])
P("")
check_band("(R4) EMPTY band", 1.0, 2.5, 0.4, [5.0, 6.0], [(12.0, 3.0), (13.0, 3.0)], [0.7, 1.0])
P("  -> no zeros in the band: Pi_I(y) -> 0 (Lemma P0'), so M_{a,b} = 0 for b >= 1, mu_b(0) = 0, and (D)")
P("     holds for every L. With 2*1 <= 2.5 <= 3*1 and s1/log(s1+3) = 1/log 4 = o(L), all hypotheses of")
P("     Theorem P7 as stated hold for chi_{-4} under GRH: the theorem would give 'no atoms' unconditionally.")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
