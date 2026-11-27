"""Referee 4: independent numerical checks of PRIMECORR.md (1.1) and P1(e) for chi_{-4}.

AI-produced, 2026-11-27, adversarial referee (referee 4).

Window: v(u) = exp(-0.5/(1-4u^2)) on (-1/2,1/2), normalised (different from the author's bump);
w = v*v (FFT convolution on a fine grid), what = vhat^2 (Gauss-Legendre quadrature).
(E1) (1.1):  Pi_T(y) = Omega_T(y) - sum_{gamma>0} 2 m_gamma what(gamma/T) cos(gamma y),
     with Omega_T computed two ways: (i) (1/2pi) int e^{iry} what(r/T) Omega_chi(r) dr,
     Omega_chi(r) = log(4/pi) + Re digamma((3/2 + i r)/2); (ii) the closed form of DENSITY (I5):
     int T w(T(u-y)) Omega(u) du, Omega(u) = -e^{-3u/2}/(1-e^{-2u}).  (i) = (ii) checks the
     normalisation of the Gamma term.
(E2) P1(e) at y = 0:  2 s_{2j}(T) = (1/2pi) int what(r/T)^{2j} Omega_chi(r) dr
                                    - 2 sum_{n <= e^{2j/T}} Lambda(n) chi(n) n^{-1/2} T w^{*2j}(T log n).
(E3) the zero-counting facts used in INVLO Thm C / Prop E for chi_{-4} up to 640: distinct
     ordinates in (t,2t], max gap (H_chi), all multiplicities 1 in the certified range.
"""
import os
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.signal import fftconvolve
from scipy.special import digamma
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ZF = os.path.join(HERE, "..", "..", "..", "prototype", "explore", "chi4_zeros.txt")
Z = np.array([float(l) for l in open(ZF) if l.strip()])
out = []
P = out.append
AB = 0.5  # bump sharpness


def vraw(u):
    u = np.asarray(u, dtype=float)
    r = np.zeros_like(u)
    m = np.abs(u) < 0.5
    r[m] = np.exp(-AB / (1.0 - 4.0 * u[m] ** 2))
    return r


# normalisation and vhat by Gauss-Legendre on [-1/2,1/2] (integrand smooth, flat at the ends)
xg, wg = leggauss(4000)
ug = 0.5 * xg
wgu = 0.5 * wg
vg = vraw(ug)
VN = np.sum(wgu * vg)


def vhat(xi):
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    res = np.empty_like(xi)
    for i in range(0, len(xi), 2000):
        X = xi[i:i + 2000]
        res[i:i + 2000] = (np.cos(np.outer(X, ug)) @ (wgu * vg)) / VN
    return res


def what(xi):
    return vhat(xi) ** 2


# w = v*v and w^{*2} = v^{*4} on fine grids
h = 2e-6
ugrid = np.arange(-0.5, 0.5 + h / 2, h)
vgrid = vraw(ugrid) / VN
w_grid = fftconvolve(vgrid, vgrid) * h
w_u = np.linspace(-1.0, 1.0, len(w_grid))
w4_grid = fftconvolve(w_grid, w_grid) * h
w4_u = np.linspace(-2.0, 2.0, len(w4_grid))
P(f"bump v = exp(-{AB}/(1-4u^2)); int v = {np.trapezoid(vgrid, ugrid):.12f}; int w = {np.trapezoid(w_grid, w_u):.12f};"
  f" int w*w = {np.trapezoid(w4_grid, w4_u):.10f}")


def w(u):
    return np.interp(u, w_u, w_grid, left=0.0, right=0.0)


def w4(u):
    return np.interp(u, w4_u, w4_grid, left=0.0, right=0.0)


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
            lam[pk] = c * np.log(p)
            pk *= int(p)
            c *= cp
    return lam


P(f"zeros: {len(Z)}, max {Z.max():.3f}; all listed ordinates distinct: {len(np.unique(np.round(Z, 9))) == len(Z)}")

P("")
P("(E1) identity (1.1)")
ymax = 12.0
Tlist = [3.0, 5.0]
NN = int(np.exp(ymax + 1 / min(Tlist))) + 5
lam = lam_chi(NN)
for T in Tlist:
    tail = what(np.array([Z.max() / T]))[0]
    wz = what(Z / T)
    rr = np.linspace(0, 400 * T, 400001)
    wr = what(rr / T)
    Om_r = Omega_chi(rr)
    for y in [6.0, 8.5, 10.0, 12.0]:
        lo = max(2, int(np.exp(y - 1 / T)) - 1)
        hi = min(NN, int(np.exp(y + 1 / T)) + 2)
        n = np.arange(lo, hi + 1)
        Pi = np.sum(lam[n] * n ** -0.5 * T * w(T * (y - np.log(n))))
        zs = np.sum(2 * wz * np.cos(Z * y))
        Om1 = 2 * np.trapezoid(wr * Om_r * np.cos(rr * y), rr) / (2 * np.pi)
        Om2 = quad(lambda u: T * w(T * (u - y)) * (-np.exp(-1.5 * u) / (1 - np.exp(-2 * u))),
                   y - 1 / T, y + 1 / T, limit=400, points=[y])[0]
        P(f"  T={T}, y={y:5.2f}: Pi_T = {Pi:+.8f}; Omega_T - zero sum = {Om1 - zs:+.8f}; diff = {Pi-(Om1-zs):+.2e};"
          f" Omega_T via (i) {Om1:+.3e}, via (ii) {Om2:+.3e}")
    P(f"  T={T}: what(gamma_max/T) = {tail:.1e} (zero-sum truncation scale)")

P("")
P("(E2) P1(e) at y = 0")
for T, j in [(3.0, 1), (1.0, 1), (0.8, 1), (2.0, 2), (1.5, 2)]:
    lhs = 2 * np.sum(what(Z / T) ** (2 * j))
    rr = np.linspace(0, 400 * T, 800001)
    dens = 2 * np.trapezoid(what(rr / T) ** (2 * j) * Omega_chi(rr), rr) / (2 * np.pi)
    nmax = int(np.exp(2 * j / T)) + 1
    lam2 = lam_chi(max(nmax, 3))
    if j == 1:
        kern = lambda u: T * w4(T * u)
    else:
        # w^{*4} = w4 * w4 on [-4, 4]
        w8_grid = fftconvolve(w4_grid, w4_grid) * h
        w8_u = np.linspace(-4.0, 4.0, len(w8_grid))
        kern = lambda u: T * np.interp(T * u, w8_u, w8_grid, left=0.0, right=0.0)
    ps = sum(lam2[n] * n ** -0.5 * kern(np.log(n)) for n in range(2, nmax + 1))
    rhs = dens - 2 * ps
    primes_used = [n for n in range(2, nmax + 1) if lam2[n] != 0 and kern(np.log(n)) != 0]
    P(f"  T={T}, j={j}: 2 s_2j (zeros) = {lhs:.8f}; Gamma term = {dens:.8f}; prime sum term = {-2*ps:+.8f}"
      f" (n used: {primes_used}); RHS = {rhs:.8f}; diff = {lhs-rhs:+.2e}")

P("")
P("(E3) zero-counting facts (chi_{-4}, certified range)")
for t in [20, 40, 80, 160, 320]:
    cnt = np.sum((Z > t) & (Z <= 2 * t))
    P(f"  #(distinct ordinates in ({t},{2*t}]) = {cnt}; (t/2pi) log t = {t/(2*np.pi)*np.log(t):.1f}; ratio {cnt/(t*np.log(t)):.3f} per t log t")
gaps = np.diff(np.concatenate([[0.0], Z]))
P(f"  max gap between consecutive ordinates (incl. 0 -> gamma_1) = {gaps.max():.3f}; so every interval of length"
  f" {3*gaps.max():.2f} in (0,640] contains >= 3 ordinates (H_chi of INVLO Thm C; numerics, not a proof)")

txt = "\n".join(out) + "\n"
open(__file__.replace('.py', '_results.txt'), 'w').write(txt)
print(txt)
