"""Eigenfunction phi_{N,M}(R,S) of T in A = R+S coordinates:
   phi_{N,M}(A-S, S) = sum_n kappa_n(S) phi_n(A; x),  n in [N-(M-S), N+(M-S)],
   kappa_n(S) = (xq;q)_{2S} gamma^{N-M, M-S}_{n-S}(x q^{2S}).
Prints kappa_n(S)/kappa_N(M) factored."""
import sys
from ratfun import RF, X, Q, ONE, mono, poch, qpoch
from tower import tower


def kappa(N, M):
    out = {}
    for S in range(M + 1):
        g = tower(N - M, M - S)
        xs = X * Q ** (2 * S)
        pre = poch(X * Q, 2 * S)
        for e, v in g.items():
            out[(e + S, S)] = pre * v.subs_x(xs)
    return out


if __name__ == '__main__':
    N = int(sys.argv[1]); M = int(sys.argv[2])
    k = kappa(N, M)
    top = k[(N, M)]
    for S in range(M, -1, -1):
        for n in range(N - (M - S), N + (M - S) + 1):
            v = k.get((n, S))
            if v is None:
                print('S=%d n=%d: 0' % (S, n)); continue
            print('S=%d n-N=%+d: %s' % (S, n - N, (v / top).factor_str()))
        print()
