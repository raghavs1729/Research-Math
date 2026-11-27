"""AUDIT a1. Independent numerical evaluation of the continued fraction (does NOT use ../../code/asym_match.py,
nor any statement of the proof).

Method: exact integer convergents p_n, q_n from the definition
    u_n = a_n u_{n-1} + b_n u_{n-2},  (p_{-1},p_0) = (1,a_0), (q_{-1},q_0) = (0,1),
for n <= N, then Richardson extrapolation in 1/n on the last k+1 exact convergents x_n = p_n/q_n:
    R_k(N) = sum_{j=0..k} (-1)^(k-j) (N-k+j)^k / (j! (k-j)!) * x_{N-k+j},
which is exact on every sequence c_0 + c_1/n + ... + c_k/n^k.  The quality of the extrapolation is measured
*internally* (differences between different (N,k)) before comparing with the claimed closed form.
"""
import sys
from fractions import Fraction as Fr
from math import factorial
from mpmath import mp, mpf, zeta, pi, log10, nstr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
mp.dps = 900
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
b = lambda n: -64 * n**10
KMAX = 220
pm, p, qm, q = 1, a(0), 0, 1
x = {}
for n in range(1, N + 1):
    pm, p = p, a(n) * p + b(n) * pm
    qm, q = q, a(n) * q + b(n) * qm
    if n >= N - 2 * KMAX - 400 or n in (10, 100, 1000):
        x[n] = mpf(p) / mpf(q)          # exact integers -> correctly rounded 900-digit quotient
def rich(Nt, k):
    s = mpf(0)
    for j in range(k + 1):
        n = Nt - k + j
        s += (-1) ** (k - j) * mpf(n) ** k / (factorial(j) * factorial(k - j)) * x[n]
    return s
target = 96 / (372 * zeta(5) - pi**5)
print("N =", N)
for n in (10, 100, 1000, N):
    print("  raw convergent n=%5d: x_n - target = %s,  n*(x_n-target) = %s" % (n, nstr(x[n] - target, 8), nstr(n * (x[n] - target), 10)))
print("internal consistency of Richardson (no reference to the target):")
best = None
for k in (20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220):
    r1, r2 = rich(N, k), rich(N - 300, k)
    d_int = abs(r1 - r2)
    print("  k=%3d  |R_k(N) - R_k(N-300)| = %s   |R_k(N) - target| = %s" % (k, nstr(d_int, 3), nstr(abs(r1 - target), 3)))
    if best is None or d_int < best[0]:
        best = (d_int, k, r1)
d_int, k, val = best
digs_int = int(-log10(d_int)) if d_int > 0 else mp.dps
print("best k =", k, "; internally justified digits ~", digs_int)
print("PCF  ~ ", nstr(val, min(digs_int, 120)))
print("claim  ", nstr(target, min(digs_int, 120)))
print("|PCF - 96/(372 zeta(5) - pi^5)| =", nstr(abs(val - target), 5))
# the constant c in x_n - PCF ~ c/n, extrapolated the same way from n (x_n - val)
