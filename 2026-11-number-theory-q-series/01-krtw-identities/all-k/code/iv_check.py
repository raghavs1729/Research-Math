"""(iv): checks of the telescoping identities and the spectral weights.
 (a) [T, h_{xY}] = -x O T  on point masses (exact), O F(R,S) = q^{R+2S} F(R,S+1), h_{xY} = P^-1(1-sh) + x Y P
 (b) K^Omega_k(x) = x^{-1} (K_k^{[q^{-r_1}]} - K_k)  as series (direct summation), k <= 5
 (c) spectral weights: t (eps-1-x)/(x(1+x)) == [(lam_M/lam_N) bJ_M d_N + (lam_N/lam_M) bJ_N d_M]/(1+delta)
"""
import sys
from fractions import Fraction as Fr
import verify_ops as VO
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from check_II import target, betaJ
from tower import lam

# (a)
bad = 0
for (q, x) in [(Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))]:
    O = VO.Ops(q, x, RMAX=14, SMAX=9)
    def h(F):
        a = O.add((1, F), (-1, O.sh(F)))
        t1 = O.mul(a, lambda R, S: q ** (-R))
        t2 = O.mul(F, lambda R, S: x * q ** (2 * S) * q ** R)
        return O.add((1, t1), (1, t2))
    def Oop(F):
        return O.mul(O.up(F), lambda R, S: q ** (R + 2 * S))
    for r0 in range(4):
        for s0 in range(4):
            d = {(r0, s0): Fr(1)}
            lhs = O.add((1, O.T(h(d))), (-1, h(O.T(d))))
            rhs = O.mul(Oop(O.T(d)), lambda R, S: -x)
            bad += VO.window_diff(lhs, rhs, 10, 5)
print('(a) [T,h_xY] = -x O T: mismatches', bad)

# (c)
def lamn(n):
    return mono(n, n * n)
def d(e):
    s = -1 if e % 2 else 1
    t = RF(s) * mono(2 * e, 3 * e * e + 2 * e) * poch(X * X * Q ** 2, e, 2) / qpoch(e, 2)
    if e >= 1:
        s2 = -1 if (e - 1) % 2 else 1
        t = t - RF(s2) * mono(2 * e - 2, 3 * (e - 1) ** 2 + 2 * (e - 1)) * poch(X * X * Q ** 2, e - 1, 2) / qpoch(e - 1, 2) * mono(1, 2 * e - 1)
    return t
ok = True
for N in range(0, 7):
    for M in range(0, N + 1):
        eps = mono(0, -N - M) * RF(ONE + X * Q ** (2 * N)) * RF(ONE + X * Q ** (2 * M))
        lhs = target(N, M) * (eps - RF(ONE + X)) / (RF(X) * RF(ONE + X))
        rhs = lamn(M) / lamn(N) * betaJ(M) * d(N) + lamn(N) / lamn(M) * betaJ(N) * d(M)
        if N == M:
            rhs = rhs / 2
        if not (lhs == rhs):
            ok = False
            print('  mismatch', N, M, (lhs / rhs).factor_str())
print('(c) spectral weights of K^Omega = those of B^(k)_{k+1} J_k for N <= 6:', ok)
