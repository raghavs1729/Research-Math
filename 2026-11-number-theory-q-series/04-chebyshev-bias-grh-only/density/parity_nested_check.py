#!/usr/bin/env python3
"""
parity_nested_check.py -- [num] sanity checks for DENSITY.md section 9 (nested comb classes).

 (1) comb <-> Fourier consistency: F = 1 + 2M sum_k b_k rho_k  (rho_k = exponential sawtooth of
     period 2pi/N_k) agrees with the Fourier series 1 + Re sum m_n 4M/(M+in) e^{int},
     m_n = sum_{k: N_k | n} b_k;
 (2) Example 9.2: nonnegative REAL weights in the nested class N = (1, 3) with F = 0 exactly
     on the arc (4pi/3, 2pi) (atom >= 1/3), for several M;
 (3) an INTEGER nested family (N_k = 3^k, b_k in {-1,0,1}, bounded): empirical measure of
     {|F| < eps} for decreasing eps -> 0 (Theorem 9.1 predicts -> 0).
Output: parity_nested_check_results.txt
"""
import math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, "parity_nested_check_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n")

def rho(t, N, M):
    h = 2 * math.pi / N
    s = np.mod(t, h)
    return h * np.exp(-M * s) / (1 - math.exp(-M * h)) - 1 / M

def F_comb(t, Ns, bs, M):
    return 1 + 2 * M * sum(b * rho(t, N, M) for N, b in zip(Ns, bs))

def F_fourier(t, Ns, bs, M, nmax):
    n = np.arange(1, nmax + 1)
    m = np.zeros(nmax)
    for N, b in zip(Ns, bs):
        m[(n % N) == 0] += b
    coef = m * 4 * M / (M + 1j * n)
    return 1 + np.real(np.exp(1j * np.outer(t, n)) @ coef)

# (1) consistency
M = 0.5; Ns = [1, 3, 9]; bs = [1, -1, 1]
t = np.linspace(0.05, 2 * math.pi - 0.05, 400)
t = t[np.min(np.abs(np.mod(t[:, None], 2 * math.pi / 9) - np.array([0, 2 * math.pi / 9])), axis=1) > 0.05]
diff = np.abs(F_comb(t, Ns, bs, M) - F_fourier(t, Ns, bs, M, 200000)).max()
log("(1) comb formula vs Fourier series (200000 terms), away from the jumps: max diff = %.2e" % diff)

# (2) Example 9.2
log("(2) Example 9.2: N = (1,3), F = 0 on (4pi/3, 2pi)")
for M in (0.1, 0.25, 0.5, 1.0, 2.0):
    A0 = 2 * math.pi / (1 - math.exp(-2 * math.pi * M))
    A1 = (2 * math.pi / 3) * math.exp(4 * math.pi * M / 3) / (1 - math.exp(-2 * math.pi * M / 3))
    b0 = A1 / (2 * (A1 - A0)); b1 = 0.5 - b0
    tt = np.linspace(4 * math.pi / 3 + 1e-6, 2 * math.pi - 1e-6, 2001)
    Fv = F_comb(tt, [1, 3], [b0, b1], M)
    tt2 = np.linspace(1e-6, 4 * math.pi / 3 - 1e-6, 2001)
    Fo = F_comb(tt2, [1, 3], [b0, b1], M)
    s2 = sum((b0 if n % 3 else 0.5) * 8 * M * M / (M * M + n * n) for n in range(1, 200001))
    log("    M = %4.2f: A1/A0 = %.4f, weights m_n = %.6f (3 not | n), 0.5 (3 | n);  max|F| on arc = %.1e;  min|F| off arc = %.3f;  sigma^2 ~ %.4f"
        % (M, A1 / A0, b0, np.abs(Fv).max(), np.abs(Fo).min(), s2))

# (3) integer nested family
log("(3) integer nested family N_k = 3^k, k = 0..K, b_k pseudo-random in {-1,0,1}")
rng = np.random.default_rng(1)
K = 11
bs = list(rng.integers(-1, 2, size=K + 1)); bs[0] = 1
Ns = [3 ** k for k in range(K + 1)]
tt = np.linspace(0, 2 * math.pi, 3_000_001, endpoint=False)
Fv = F_comb(tt, Ns, bs, M=0.5)
log("    b = %s" % bs)
for eps in (0.1, 0.03, 0.01, 0.003, 0.001):
    log("    eps = %6.3f   measure{|F| < eps}/2pi = %.5f" % (eps, np.mean(np.abs(Fv) < eps)))
for c in (1.0, -1.0, 3.0):
    log("    level c = %4.1f: measure{|F - c| < 0.001}/2pi = %.5f" % (c, np.mean(np.abs(Fv - c) < 0.001)))
out.close()
