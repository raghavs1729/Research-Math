#!/usr/bin/env python3
"""
integer_counterexample_positive_demo.py -- [num] the construction of DENSITY.md Theorem 9.6
exactly as stated: rounding UP, bump b = c * c~ (nonnegative Fourier transform), K = 4,
N_j = 4^j, N'_j = 3 N_{j+1}, rho_j = 0.05 (j+1)^(-5/4).  Three steps (bands [1,4), [4,16),
[16,64)).  Checks after each step: F = 0 on J minus the gaps (exactly, up to rounding),
multiplicities below N_{j+1} are integers >= 1, all multiplicity changes above N_{j+1} are >= 0.
Output: integer_counterexample_positive_demo_results.txt
"""
import math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, "integer_counterexample_positive_demo_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n")
M = 0.5; L = 2 ** 22
t = 2 * math.pi * np.arange(L) / L
def rho_saw(t, N):
    h = 2 * math.pi / N; s = np.mod(t, h)
    return h * np.exp(-M * s) / (1 - math.exp(-M * h)) - 1 / M
A0 = 2 * math.pi / (1 - math.exp(-2 * math.pi * M))
A1 = (2 * math.pi / 3) * math.exp(4 * math.pi * M / 3) / (1 - math.exp(-2 * math.pi * M / 3))
b0 = A1 / (2 * (A1 - A0)); b1 = 0.5 - b0
F = 1 + 2 * M * (b0 * rho_saw(t, 1) + b1 * rho_saw(t, 3))
J = (t > 4 * math.pi / 3 + 1e-3) & (t < 2 * math.pi - 1e-3)
NMAX = 20000
m_base = np.array([0.0] + [b0 + (b1 if n % 3 == 0 else 0.0) for n in range(1, NMAX)])
corr = np.zeros(NMAX)

# profile b = c * c~ on [-1,1], c(x) = exp(-1/(1-4x^2)) on (-1/2,1/2); tabulate b and b'
xs = np.linspace(-1, 1, 40001); dx = xs[1] - xs[0]
cx = np.zeros_like(xs); ins = np.abs(xs) < 0.5
cx[ins] = np.exp(-1 / (1 - 4 * xs[ins] ** 2))
bt = np.convolve(cx, cx[::-1], mode="same") * dx
dbt = np.gradient(bt, dx)
def profile(x):
    return np.interp(x, xs, bt, left=0, right=0), np.interp(x, xs, dbt, left=0, right=0)

def bumps(Np, rho):
    h = 2 * math.pi / Np; w = h * rho
    s = np.mod(t + h / 2, h) - h / 2
    x = 2 * s / w
    b, db = profile(x)
    return b, db * (2 / w)

log("base: Example 9.2, M = 1/2; base multiplicities in {%.6f, 0.5}; max|F_0| on J = %.1e" % (b0, np.abs(F[J]).max()))
Ns = [1, 4, 16, 64]
gaps = np.zeros(L, dtype=bool)
k = np.fft.fftfreq(L, d=1.0 / L).astype(int)
for j in range(3):
    Nj, Nj1 = Ns[j], Ns[j + 1]
    band = np.arange(Nj, Nj1)
    mcur = m_base[band] + corr[band]
    dm = np.ceil(mcur - 1e-12) - mcur
    dm[dm < 0] = 0
    dPsi = 2 * M * M * dm / (M * M + band ** 2)
    eP = np.zeros(L); deP = np.zeros(L)
    for n, d in zip(band, dPsi):
        eP += 2 * d * np.cos(n * t); deP -= 2 * d * n * np.sin(n * t)
    Np = 3 * Nj1; rho = 0.05 * (j + 1) ** (-1.25)
    beta, dbeta = bumps(Np, rho)
    bm = beta.mean()
    phi_e = eP * beta / bm; dphi_e = (deP * beta + eP * dbeta) / bm
    phi = phi_e - dphi_e / M
    F = F + phi
    gaps |= beta > 0
    c = np.fft.fft(phi_e) / L                          # Psi-coefficients of the correction
    n_idx = np.arange(1, NMAX)
    m_change = (M * M + n_idx ** 2) * np.real(c[n_idx]) / (2 * M * M)
    m_change[: Nj1 - 1] = 0                           # below N_{j+1}: handled by rounding
    low = (M * M + band ** 2) * np.real(c[band]) / (2 * M * M)
    corr[band] += low                                  # = dm (rounding)
    corr[Nj1:NMAX] += m_change[Nj1 - 1:NMAX - 1]
    mm = m_base[1:Nj1] + corr[1:Nj1]
    S = J & ~gaps
    log("step %d: band [%d,%d), gap period 2pi/%d, gap measure/2pi %.4f, ||phi||_2 = %.3e" %
        (j, Nj, Nj1, Np, (beta > 0).mean(), math.sqrt((phi ** 2).mean())))
    log("   multiplicities n < %d: %s ; max dist to Z = %.1e ; min = %.3f" %
        (Nj1, np.round(mm, 6).tolist()[:16], np.abs(mm - np.round(mm)).max(), mm.min()))
    log("   multiplicity changes above N_{j+1}: min = %.2e (should be >= 0 up to FFT error), max = %.2f" %
        (m_change[Nj1 - 1:].min(), m_change[Nj1 - 1:].max()))
    log("   max|F| on J minus gaps = %.1e ; measure(J minus gaps)/2pi = %.4f" % (np.abs(F[S]).max(), S.mean()))
out.close()
