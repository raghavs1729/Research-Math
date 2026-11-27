"""Beam search for the toy small-ball quantity ([num]).
Grow m one coordinate at a time; keep the `width` prefixes with the smallest
eta-quantile of |F_n| (exact on a K-point grid). Since appending a 0 leaves F
unchanged, tau*(n) is nonincreasing in n along the beam."""
import numpy as np, sys
from toy_smallball import basis

M = float(sys.argv[1]); eta = float(sys.argv[2]); Nmax = int(sys.argv[3]); width = int(sys.argv[4])
K = 8*Nmax
t, B = basis(Nmax, K, M)
kq = int(eta*K)
Fs = np.ones((1, K)); ms = np.zeros((1, 0), dtype=np.int8)
report = {16, 22, 32, 64, 128, 256, 512, 1024}
for n in range(Nmax):
    F1 = Fs + B[n][None, :]
    allF = np.vstack([Fs, F1])
    allm = np.vstack([np.hstack([ms, np.zeros((len(ms), 1), np.int8)]),
                      np.hstack([ms, np.ones((len(ms), 1), np.int8)])])
    q = np.partition(np.abs(allF), kq, axis=1)[:, kq]
    # dedupe by F (rounded) to keep diversity
    order = np.argsort(q)
    keep = []; seen = set()
    for i in order:
        key = tuple(np.round(allF[i, ::max(1, K//64)], 6))
        if key in seen: continue
        seen.add(key); keep.append(i)
        if len(keep) >= width: break
    Fs = allF[keep]; ms = allm[keep]
    if (n+1) in report:
        print(f"n={n+1:5d}  tau*={q[order[0]]:.4e}  sqrt(n)*tau*={np.sqrt(n+1)*q[order[0]]:.3f}  "
              f"best m prefix={''.join(map(str, allm[order[0]][:40]))}", flush=True)
