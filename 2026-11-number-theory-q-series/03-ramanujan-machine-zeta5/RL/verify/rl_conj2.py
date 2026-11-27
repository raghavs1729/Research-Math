"""Conjecture 2 of arXiv:2412.12361v2:  C[n^2/(k^2 (1-4n^2))] = (2/k)/(ln(k+1) - ln(k-1)),  k >= 1.
It is the classical continued fraction (DLMF 4.9.2)  ln((1+z)/(1-z)) = 2z/(1 - z^2/(3 - 4z^2/(5 - 9z^2/(7 - ...))))  at z = 1/k,
after the equivalence transformation r_n = 1/(2n+1).
(a) exact: the n-th approximants of C[f] and of (1/k) * (1 - (1/k^2)/(3 - (4/k^2)/(5 - ...))) coincide (Fractions, n <= 30, k = 1..6);
(b) numeric: k = 2..8, 40 digits;   (c) k = 1: the N-th approximant of C[n^2/(1-4n^2)] is exactly 1/H_{N+1} (harmonic number) -> 0."""
from fractions import Fraction as Fr
from mpmath import mp, mpf, log
def approx_C(k, N):
    v = Fr(1)
    for n in range(N, 0, -1): v = 1 + Fr(n * n, k * k * (1 - 4 * n * n)) / v
    return v
def approx_G(k, N):          # 1 - (1/k^2)/(3 - (4/k^2)/(5 - ... - (N^2/k^2)/(2N+1)))
    v = Fr(2 * N + 1)
    for n in range(N, 0, -1): v = (2 * n - 1) - Fr(n * n, k * k) / v
    return v
for k in range(1, 7):
    for N in range(1, 31):
        if k == 1 and approx_C(1, N) == 0: continue
        assert approx_C(k, N) == approx_G(k, N), (k, N)
print("(a) PASS: approximants of C[n^2/(k^2(1-4n^2))] = approximants of the denominator of Gauss's fraction at z = 1/k")
mp.dps = 40
for k in range(2, 9):
    v = mpf(1)
    for n in range(400, 0, -1): v = 1 + mpf(n * n) / (k * k * (1 - 4 * n * n)) / v
    assert abs(v - (mpf(2) / k) / (log(k + 1) - log(k - 1))) < mpf(10)**-35, k
print("(b) PASS: k = 2..8 to 35 digits")
H = Fr(1)
for N in range(1, 60):
    H += Fr(1, N + 1)                      # H_{N+1}
    assert approx_C(1, N) == 1 / H, N
print("(c) PASS: for k = 1 the N-th approximant is exactly 1/H_{N+1} (N < 60), so the value is 0 = 2/(ln 2 - ln 0)")
