"""Reversed-chain operator T#:  (T# G)(r,s) = sum_{R>=r, S<=s} w(R,S) G(R,S) / ((q)_{R-r} (q^2;q^2)_{s-S}),
w(R,S) = x^{R+2S} q^{(R+S)^2+S^2}.  Left eigenvectors of T are psi = w * (T#-eigenfunctions).
Basis Xi#_{e,rho}(R,S) = delta_{R,rho} / ((Q;Q)_{S-e} (x^2 q^{2 rho} Q;Q)_{S+e}),  Q = q^2.
Diagonal entry lam_{rho+e, e}.  Eigenvector of mode (N,M): top Xi#_{M, N-M}, supported on rho <= N-M.
(lam_{N,M} - lam_{rho+e,e}) c_{e,rho} = sum_{rho'>rho} sum_{e''} lam_{rho'+e'',e''} GQ^{(rho'-rho)}_{e'',e}(x^2 q^{2rho}) c_{e'',rho'} / (q)_{rho'-rho}
GQ^{(c)}_{d,e}(a) = coefficient of 1/((Q;Q)_{S-e}(aQ;Q)_{S+e}) in 1/((Q;Q)_{S-d}(aQ^{c+1};Q)_{S+d}).
"""
import sys
from functools import lru_cache
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import lam


def qbinQ(n, k):
    if k < 0 or k > n:
        return RF(0)
    return qpoch(n, 2) / (qpoch(k, 2) * qpoch(n - k, 2))


@lru_cache(maxsize=None)
def GQ(c, d, e, rho):
    """change of a in base Q, a = x^2 q^{2 rho}, a -> a Q^c"""
    n = e - d
    if n < 0 or n > c:
        return RF(0)
    a = X * X * Q ** (2 * rho)
    if d + e == 0:
        return RF(1)
    num = RF(ONE - a * Q ** (4 * e)) * poch(a * Q ** 2, d + e - 1, 2)   # (1-aQ^{2e}) (a;Q)_{d+e}/(1-a)
    den = poch(a * Q ** (2 * c + 2), d + e, 2)
    sgn = -1 if n % 2 else 1
    # a^n Q^{binom(n,2) + (e+d) n}
    am = RF(a ** n) * mono(0, 2 * (n * (n - 1) // 2 + (e + d) * n))
    return num / den * am * qbinQ(c, n) * sgn


def eig_sharp(N, M):
    """dict (e, rho) -> c, top (M, N-M) = 1"""
    L = N - M
    c = {(M, L): RF(1)}
    lm = lam(N, M)
    for rho in range(L - 1, -1, -1):
        acc = {}
        for (e2, rp), v in list(c.items()):
            if rp <= rho:
                continue
            dr = rp - rho
            pre = lam(rp + e2, e2) * v / qpoch(dr)
            for e in range(e2, e2 + dr + 1):
                g = GQ(dr, e2, e, rho)
                if g.iszero():
                    continue
                acc[e] = acc.get(e, RF(0)) + pre * g
        for e, v in acc.items():
            if v.iszero():
                continue
            c[(e, rho)] = v / (lm - lam(rho + e, e))
    return c


if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2])
    c = eig_sharp(N, M)
    for (e, rho), v in sorted(c.items(), key=lambda t: (-t[0][1], t[0][0])):
        print('rho=%d e=%d: %s' % (rho, e, v.factor_str()))
    s = sum((v for (e, rho), v in c.items() if rho == 0), RF(0))
    print('sum at rho=0:', s.factor_str())
