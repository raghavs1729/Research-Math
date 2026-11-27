"""Check (II'):  t_{L,0}(x q^{2M}) * sum_e gamma^{L,M}_e(x) == t_{L+M,M}(x)   (exact, rational functions)."""
import sys
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import tower, lam


def betaB(n):
    s = -1 if n % 2 else 1
    if n == 0:
        return RF(1)
    return RF(s) * mono(0, n * n) * poch(X * X, n, 2) * RF(ONE - X * Q ** (2 * n)) / (qpoch(n, 2) * RF(ONE - X))


def betaJ(m):
    s = -1 if m % 2 else 1
    if m == 0:
        return RF(1)
    return RF(s) * mono(0, m * m - m) * poch(X * X, m, 2) * RF(ONE - X * X * Q ** (4 * m)) / (qpoch(m, 2) * RF(ONE - X * X))


def target(N, M):
    t = betaB(N) * betaJ(M) + betaB(M) * betaJ(N)
    if N == M:
        t = t / 2
    return lam(N, M) * t


if __name__ == '__main__':
    LM = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    MM = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    for M in range(1, MM + 1):
        for L in range(0, LM + 1):
            g = tower(L, M)
            s = RF(0)
            for e, v in g.items():
                s = s + v
            lhs = target(L, 0).subs_x(X * Q ** (2 * M)) * s
            rhs = target(L + M, M)
            ok = (lhs == rhs)
            print('L=%d M=%d: %s' % (L, M, 'OK' if ok else 'FAIL'))
            if not ok:
                print('   ratio:', (lhs / rhs).factor_str())
            sys.stdout.flush()
