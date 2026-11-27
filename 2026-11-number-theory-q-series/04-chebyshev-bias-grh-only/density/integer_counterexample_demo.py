#!/usr/bin/env python3
"""
integer_counterexample_demo.py -- [num] demonstration of the mechanics of DENSITY.md
Theorem 9.6 (integer multiplicities do NOT exclude an atom at 0): two steps of the
construction "round a band of multiplicities to integers, then cancel the rounding error
exactly by a correction supported on many tiny equispaced gaps".

Base: Example 9.2 (real nested weights, F_0 = 0 exactly on J = (4pi/3, 2pi)), M = 1/2.
Step j: band B_j = [N_j, N_{j+1}); round m_n (n in B_j) to integers; the change of
Psi (F = (1 - D/M) Psi, Psi even) is the trig polynomial eP; correction
phi_j = (1 - D/M)(eP * beta_j / mean(beta_j)), beta_j = smooth even bumps of period
2pi/N'_j (N'_j > 2 N_{j+1}) and total measure 2 pi rho_j.  Then F_{j+1} = F_j + phi_j:
 - equals F_j outside the gaps, so still 0 on J minus the gaps;
 - its multiplicities m_n = F^(n)(M+in)/(2M) are integers for all n < N_{j+1}.
Everything is computed on a grid of 2^21 points by FFT.  Output: integer_counterexample_demo_results.txt
"""
import math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, "integer_counterexample_demo_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n")

M = 0.5
L = 2 ** 21
t = 2 * math.pi * np.arange(L) / L
def rho_saw(t, N):
    h = 2 * math.pi / N; s = np.mod(t, h)
    return h * np.exp(-M * s) / (1 - math.exp(-M * h)) - 1 / M
A0 = 2 * math.pi / (1 - math.exp(-2 * math.pi * M))
A1 = (2 * math.pi / 3) * math.exp(4 * math.pi * M / 3) / (1 - math.exp(-2 * math.pi * M / 3))
b0 = A1 / (2 * (A1 - A0)); b1 = 0.5 - b0
F = 1 + 2 * M * (b0 * rho_saw(t, 1) + b1 * rho_saw(t, 3))
J = (t > 4 * math.pi / 3 + 1e-3) & (t < 2 * math.pi - 1e-3)

def coeffs(f):                      # c_n = (1/2pi) int f e^{-int}  (used only for the smooth corrections)
    return np.fft.fft(f) / L
def mult(c, n):                     # m_n = c_n (M + i n) / (2M)
    return c[n] * (M + 1j * n) / (2 * M)
NMAX = 5000
m_base = np.array([0.0] + [b0 + (b1 if n % 3 == 0 else 0.0) for n in range(1, NMAX)])   # exact base multiplicities
corr = np.zeros(NMAX, dtype=complex)                                                     # accumulated corrections
def m_total(n):
    return m_base[n] + corr[n]

def bump_periodic(Np, rho):
    """smooth even bumps of width w = 2 pi rho / Np centred at 2 pi k / Np, and their exact derivative"""
    h = 2 * math.pi / Np; w = h * rho
    s = np.mod(t + h / 2, h) - h / 2            # signed distance to nearest grid point
    x = 2 * s / w                                # bump support |x| < 1
    b = np.zeros_like(t); db = np.zeros_like(t); inside = np.abs(x) < 1
    xi = x[inside]
    b[inside] = np.exp(-1 / (1 - xi ** 2))
    db[inside] = b[inside] * (-2 * xi / (1 - xi ** 2) ** 2) * (2 / w)
    return b, db

log("base: Example 9.2, M = 1/2, max|F_0| on J = %.2e" % np.abs(F[J]).max())
Ns = [1, 8, 64]
gaps = np.zeros_like(t, dtype=bool)
n_all = np.arange(1, 200)
for j in range(len(Ns) - 1):
    Nj, Nj1 = Ns[j], Ns[j + 1]
    band = np.arange(Nj, Nj1)
    m = np.array([m_total(n) for n in band])
    assert np.abs(m.imag).max() < 1e-8, "model class violated"
    dm = np.round(m.real) - m.real
    # change of Psi: Psi^(n) = 2 M^2 m_n / (M^2 + n^2) (two-sided, even)
    dPsi = 2 * M * M * dm / (M * M + band ** 2)
    eP = np.zeros_like(t); deP = np.zeros_like(t)
    for n, d in zip(band, dPsi):
        eP += 2 * d * np.cos(n * t); deP -= 2 * d * n * np.sin(n * t)
    Np = 3 * Nj1; rho = 0.1 / (j + 1) ** 2
    beta, dbeta = bump_periodic(Np, rho)
    bm = beta.mean()
    phi_e = eP * beta / bm
    dphi_e = (deP * beta + eP * dbeta) / bm          # exact derivative
    phi = phi_e - dphi_e / M
    F = F + phi
    gaps |= beta > 0
    cphi = coeffs(phi)
    corr += np.array([mult(cphi, n) if n > 0 else 0 for n in range(NMAX)])
    mm = np.array([m_total(n) for n in range(1, Nj1)])
    log("step %d: band [%d,%d), gap period 2pi/%d, gap measure %.4f, ||phi||_2 = %.3e" %
        (j, Nj, Nj1, Np, (beta > 0).mean(), math.sqrt((phi ** 2).mean())))
    log("    multiplicities n < %d: max |m - round(m)| = %.2e, max |Im m| = %.2e" %
        (Nj1, np.abs(mm.real - np.round(mm.real)).max(), np.abs(mm.imag).max()))
    log("    integer values: %s" % np.round(mm.real).astype(int).tolist()[:20])
    S = J & ~gaps
    log("    max|F| on J minus gaps = %.2e   (measure of J minus gaps / 2pi = %.4f)" %
        (np.abs(F[S]).max(), S.mean()))
mhi = np.array([m_total(n).real for n in range(Ns[-1], Ns[-1] + 600)])
log("multiplicities just above the last band (still real, to be rounded at the next step): max |m| = %.2f" % np.abs(mhi).max())
out.close()
