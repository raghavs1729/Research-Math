"""Referee 2: the x = 1 evaluations used in Theorem 5 (and Theorem 6 at x = 1), as q-series mod q^L, k = 1..kmax.
 (a) sum_n q^{(k+1)n^2} b^Omega_n(1) with b^Omega_n(1) = (-1)^n (q^{n^2-2n} + q^{n^2+2n}) (n>=1), 1 (n=0)
       == (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_inf                                    [Jacobi triple product step, D8]
 (b) (q^2;q^2)_inf C_k(1) == (q^{k+1}, q^{k+3}, q^{2k+4}; q^{2k+4})_inf C_1(1)       [J_k(1), from the multisums C_k, C_1]
 (c) (q;q)_inf B^{(k)}_{k+1}(1) == (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_inf           [Theorem 6 at x = 1, Bressoud]
Usage: python3 r2_x1.py [Nmax=400] [kmax=14]
"""
import sys
from flint import fmpz_poly
from r2_lib import mono, poch, InvCache

Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400
kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 14
L = Nmax + 1
cache = InvCache(L, 4 * int(L ** 0.5) + 10)


def prod3(e1, e2, m):
    return poch(e1, m, 'inf', L).mul_low(poch(e2, m, 'inf', L), L).mul_low(poch(m, m, 'inf', L), L)


def C1(k):  # C_k(1) by the s-chain DP
    smax = int(L ** 0.5) + 1
    f = {s: cache.Q(s).mul_low(mono(s * s, L), L) for s in range(smax + 1)}
    for _ in range(k - 1):
        f = {s2: sum((f[s].mul_low(cache.Q(s2 - s), L) for s in range(s2 + 1)), fmpz_poly()).mul_low(mono(s2 * s2, L), L)
             for s2 in range(smax + 1)}
    return sum(f.values(), fmpz_poly()).truncate(L)


def Bk1(k):  # B^{(k)}_{k+1}(1)
    nmax = int(L ** 0.5) + 1
    lin = lambda j, n: n if j in (k, k + 1) else 0
    f = {n: cache.Q(n).mul_low(mono(n * n + lin(k + 1, n), L), L) for n in range(nmax + 1)}
    for j in range(k, 0, -1):
        f = {n2: sum((f[n].mul_low(cache.q(n2 - n), L) for n in range(n2 + 1)), fmpz_poly()).mul_low(mono(n2 * n2 + lin(j, n2), L), L)
             for n2 in range(nmax + 1)}
    return sum(f.values(), fmpz_poly()).truncate(L)


ok_all = True
lines = ["Referee 2: x = 1 evaluations mod q^%d, k = 1..%d" % (L, kmax)]
qinf = poch(1, 1, 'inf', L)
Qinf = poch(2, 2, 'inf', L)
c1 = C1(1)
for k in range(1, kmax + 1):
    m = 2 * k + 4
    th = fmpz_poly([1])
    n = 1
    while (k + 2) * n * n - 2 * n < L:
        th += mono((k + 1) * n * n + n * n - 2 * n, L, (-1) ** n) + mono((k + 1) * n * n + n * n + 2 * n, L, (-1) ** n)
        n += 1
    a = (th.truncate(L) - prod3(k, k + 4, m)).truncate(L).is_zero()
    b = (Qinf.mul_low(C1(k), L) - prod3(k + 1, k + 3, m).mul_low(c1, L)).truncate(L).is_zero()
    c = (qinf.mul_low(Bk1(k), L) - prod3(k, k + 4, m)).truncate(L).is_zero()
    ok_all &= a and b and c
    lines.append("  k=%2d: (a) theta = (q^k,q^{k+4},q^{2k+4})_inf %s; (b) J_k(1) %s; (c) Bressoud B^{(k)}_{k+1}(1) %s"
                 % (k, "OK" if a else "FAIL", "OK" if b else "FAIL", "OK" if c else "FAIL"))
# negative control: (a) with B = -1 instead of -2
k = 4
th = fmpz_poly([1])
n = 1
while 6 * n * n < L:
    th += mono(6 * n * n - n, L, (-1) ** n) + mono(6 * n * n + n, L, (-1) ** n)
    n += 1
neg = not (th.truncate(L) - prod3(4, 8, 12)).truncate(L).is_zero()
lines.append("  negative control (theta with linear term -n instead of -2n, k=4) detected: %s" % neg)
ok_all &= neg
lines.append("RESULT: %s" % ("ALL OK" if ok_all else "FAILURE"))
print("\n".join(lines))
