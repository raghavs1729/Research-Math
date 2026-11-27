"""Beam search, longer runs, and structure of the best sequence ([num])."""
import numpy as np, sys
from toy_smallball import basis
M = float(sys.argv[1]); eta = float(sys.argv[2]); Nmax = int(sys.argv[3]); width = int(sys.argv[4])
K = 8*Nmax
t, B = basis(Nmax, K, M)
kq = int(eta*K)
Fs = np.ones((1, K)); ms = np.zeros((1, 0), dtype=np.int8)
for n in range(Nmax):
    allF = np.vstack([Fs, Fs + B[n][None, :]])
    allm = np.vstack([np.hstack([ms, np.zeros((len(ms), 1), np.int8)]),
                      np.hstack([ms, np.ones((len(ms), 1), np.int8)])])
    q = np.partition(np.abs(allF), kq, axis=1)[:, kq]
    order = np.argsort(q); keep = []; seen = set()
    for i in order:
        key = tuple(np.round(allF[i, ::max(1, K//64)], 6))
        if key in seen: continue
        seen.add(key); keep.append(i)
        if len(keep) >= width: break
    Fs = allF[keep]; ms = allm[keep]
    if (n+1) & n == 0 and n >= 15:
        print(f"n={n+1:5d} tau*={q[order[0]]:.4e} sqrt(n)tau*={np.sqrt(n+1)*q[order[0]]:.3f}", flush=True)
F = Fs[0]; m = ms[0]
ones = np.nonzero(m)[0] + 1
print("positions of 1s:", ones.tolist()[:80])
tau = np.partition(np.abs(F), kq)[kq]
S = np.abs(F) <= tau
edges = np.nonzero(np.diff(np.r_[S[-1], S]).astype(int))[0]
print(f"tau={tau:.3e}, |S|/K={S.mean():.3f}, number of runs={len(edges)//2}")
# where is S? histogram over 32 bins of the circle
h, _ = np.histogram(t[S], bins=32, range=(0, 2*np.pi))
print("S histogram (32 bins over [0,2pi)):", h.tolist())
np.save(f"beam_best_eta{eta}_N{Nmax}.npy", m)
