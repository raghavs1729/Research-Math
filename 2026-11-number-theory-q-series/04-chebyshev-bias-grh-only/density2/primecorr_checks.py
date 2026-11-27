"""Numerical check of identity (1.1) of PRIMECORR.md for chi_{-4}:

    Pi_T(y) = Omega_T(y) - sum_{gamma>0} 2 m_gamma what(gamma/T) cos(gamma y),

Pi_T(y) = sum_n chi(n) Lambda(n) n^{-1/2} T w(T(y - log n)),  w = v*v,  v a C^infty bump on
(-1/2,1/2) normalised to integral 1, what = vhat^2.  Zeros: prototype/explore/chi4_zeros.txt
(first 511 ordinates, all simple).  Omega_T(y) = (1/2pi) int e^{iry} what(r/T) Omega(r) dr with
Omega(r) = log(4/pi) + Re digamma((1/2 + 1 + ir)/2) (chi_{-4} is odd, a = 1).

AI-produced (Anthropic's Claude), 2026-11-27.
"""
import os
import numpy as np
from scipy.special import digamma
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.array([float(l) for l in open(os.path.join(HERE, "..", "..", "prototype", "explore",
                                                  "chi4_zeros.txt")) if l.strip()])

# bump v(u) = exp(-1/(1-4u^2)) on (-1/2,1/2), normalised
U = np.linspace(-0.5, 0.5, 20001)[1:-1]
vraw = np.exp(-1.0 / (1.0 - 4.0 * U ** 2))
vnorm = np.trapezoid(vraw, U)
du = U[1] - U[0]


def vhat(xi):
    xi = np.atleast_1d(xi)
    return np.array([np.trapezoid(vraw * np.cos(x * U), U) / vnorm for x in xi])


def what(xi):
    return vhat(xi) ** 2


# w = v*v on a grid (support [-1,1])
vgrid = vraw / vnorm
wgrid = np.convolve(vgrid, vgrid) * du
Ugrid_w = np.linspace(-1.0, 1.0, len(wgrid))


def w(u):
    return np.interp(u, Ugrid_w, wgrid, left=0.0, right=0.0)


# chi_{-4} Lambda(n) up to N
def lambda_chi(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    primes = np.nonzero(sieve)[0]
    lam = np.zeros(N + 1)
    for p in primes:
        if p == 2:
            continue  # chi(2^k) = 0
        chi_p = 1.0 if p % 4 == 1 else -1.0
        pk, c = p, chi_p
        while pk <= N:
            lam[pk] = c * np.log(p)
            pk *= p
            c *= chi_p
    return lam


out = []
P = out.append
ymax = 13.0
Tlist = [4.0, 8.0]
N = int(np.exp(ymax + 1.0 / min(Tlist))) + 10
lam = lambda_chi(N)
ns = np.arange(N + 1)
logn = np.log(np.maximum(ns, 1))
P(f"zeros used: {len(Z)} (max {Z.max():.2f}); primes up to {N}")
for T in Tlist:
    wz = what(Z / T)
    tail = what(np.array([Z.max() / T]))[0]
    P(f"T={T}: what(gamma_max/T) = {tail:.2e} (truncation of the zero sum)")

    def Omega(r):
        return np.log(4 / np.pi) + np.real(digamma((1.5 + 1j * r) / 2))

    for y in [7.0, 9.0, 11.0, 13.0]:
        lo, hi = int(np.exp(y - 1 / T)) - 1, int(np.exp(y + 1 / T)) + 2
        idx = np.arange(max(lo, 2), min(hi, N) + 1)
        Pi = np.sum(lam[idx] * idx ** -0.5 * T * w(T * (y - logn[idx])))
        zsum = np.sum(2 * wz * np.cos(Z * y))
        # Omega_T(y): (1/2pi) int what(r/T) Omega(r) cos(ry) dr (even integrand), r up to 40T
        rr = np.linspace(0, 40 * T, 40001)
        integrand = what(rr / T) * Omega(rr) * np.cos(rr * y)
        Om = 2 * np.trapezoid(integrand, rr) / (2 * np.pi)
        P(f"  y={y:5.1f}: Pi_T = {Pi:+.6f}; Omega_T - zero sum = {Om - zsum:+.6f}; "
          f"difference = {Pi - (Om - zsum):+.2e}; (Omega_T = {Om:+.2e})")

open(os.path.join(HERE, "primecorr_checks_results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
