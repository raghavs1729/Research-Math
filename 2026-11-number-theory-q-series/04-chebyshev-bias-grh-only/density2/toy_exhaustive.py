"""Exhaustive small-N check of the toy small-ball quantity ([num]).
For all m in {0,1}^N (N<=22): q_eta(m) = eta-quantile of |F_N| on a K-point grid.
Prints min over m for several eta, and the minimizer's pattern."""
import numpy as np, sys, itertools
from toy_smallball import basis

M = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
etas = [0.05, 0.10, 0.15, 0.20, 0.30]
for N in [10, 14, 18, 22]:
    K = 16*N
    t, B = basis(N, K, M)
    B = B.astype(np.float64)
    best = {e: (1e9, None) for e in etas}
    chunk = 1 << min(N, 15)
    bits = ((np.arange(chunk)[:, None] >> np.arange(min(N, 15))) & 1).astype(np.float64)
    nlow = min(N, 15)
    for hi in range(1 << (N - nlow)):
        mhi = np.array([(hi >> k) & 1 for k in range(N - nlow)], float)
        base = 1 + mhi @ B[nlow:] if N > nlow else 1 + 0*B[0]
        F = base[None, :] + bits @ B[:nlow]
        A = np.abs(F)
        for e in etas:
            kq = int(e*K)
            q = np.partition(A, kq, axis=1)[:, kq]
            i = np.argmin(q)
            if q[i] < best[e][0]:
                m = np.r_[bits[i], mhi]
                best[e] = (q[i], m.astype(int))
    print(f"N={N} K={K}: " + "  ".join(f"eta={e}: {best[e][0]:.3e}" for e in etas), flush=True)
    for e in [0.05, 0.15, 0.30]:
        print(f"    eta={e} minimizer m = {''.join(map(str, best[e][1]))}")
