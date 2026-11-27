"""Referee check of Thm 6.4 (expansion of u): u(R,S) = 1/(q)_R  vs  sum_{N>=M} u~_{N,M} phi_{N,M}(R,S),
with phi from the position-space solve (ref_eigen.eig_T) and u~ = V# (x^2q^{2L+2};Q)_{2M}(xq^{2M+1};q)_{2L}/(x^2Q;Q)_inf.
Only L <= R and M >= S contribute; the tail M > Mmax is O(x^{2(Mmax+1-S)})."""
import sys
from flint import fmpq
import mpmath as mp
from ref_common import Point, poch
from ref_eigen import eig_T, Vs_closed

mp.mp.dps = 80


def tomp(a):
    return mp.mpf(int(a.p)) / mp.mpf(int(a.q))


def run(qv, xv, Rbox=3, Sbox=2, Mmax=10):
    P = Point(qv, xv)
    q, x, Q = P.q, P.x, P.Q
    x2Qinf = mp.qp(tomp(x)**2 * tomp(Q), tomp(Q))
    acc = {(R, S): fmpq(0) for R in range(Rbox + 1) for S in range(Sbox + 1)}
    for L in range(Rbox + 1):
        for M in range(Mmax + 1):
            N = L + M
            phi = eig_T(P, N, M, Rbox)
            ut = Vs_closed(P, N, M) * poch(x**2 * q**(2 * L + 2), Q, 2 * M) * poch(x * q**(2 * M + 1), q, 2 * L)
            for (R, S) in acc:
                if S <= M:
                    acc[(R, S)] += ut * phi[(R, S)]
    worst = mp.mpf(0)
    for (R, S), v in acc.items():
        lhs = tomp(v) / x2Qinf
        rhs = tomp(P.iq(R))
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
        print(f"  (R,S)=({R},{S}): sum = {mp.nstr(lhs, 25)}   1/(q)_R = {mp.nstr(rhs, 25)}")
    return worst


if __name__ == "__main__":
    for (qv, xv) in ((fmpq(2, 7), fmpq(1, 50)), (fmpq(-3, 5), fmpq(1, 40))):
        w = run(qv, xv)
        print(f"(q,x)=({qv},{xv}), L<=3, M<=10: worst rel. err {mp.nstr(w, 3)}")
