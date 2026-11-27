"""eigenfunction phi_{N,M}(R,S) in position space (R,S coordinates): F(R,S) = sum_n kappa_n(S) phi_n(R+S; x)"""
import sys
from hrec import eig
from ratfun import RF, X, Q, ONE, poch, qpoch, mono

def phi(n, A):
    if A < n: return RF(0)
    return RF(1) / (qpoch(A - n) * poch(X * Q, A + n))

N = int(sys.argv[1]); M = int(sys.argv[2]); RMAX = int(sys.argv[3])
e = eig(N, M)
for S in range(M, -1, -1):
    for R in range(N - M, N - M + RMAX + 1):
        A = R + S
        v = sum((c * phi(n, A) for (n, B), c in e.items() if B == S), RF(0))
        # normalise by the S=M, R=L value
        print('S=%d R=%d: %s' % (S, R, v.factor_str()))
    print()
