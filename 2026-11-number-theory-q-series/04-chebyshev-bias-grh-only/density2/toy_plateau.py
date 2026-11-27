"""Plateau of the beam-search small-ball quantity vs eta ([num]).
Entropy heuristic predicts log2(1/tau_floor) ~ c/eta."""
import numpy as np, sys
from toy_smallball import basis
M = float(sys.argv[1]); Nmax = int(sys.argv[2]); width = int(sys.argv[3])
etas = [float(x) for x in sys.argv[4].split(',')]
K = 8*Nmax
t, B = basis(Nmax, K, M)
for eta in etas:
    kq = int(eta*K)
    Fs = np.ones((1, K)); hist = {}
    for n in range(Nmax):
        allF = np.vstack([Fs, Fs + B[n][None, :]])
        q = np.partition(np.abs(allF), kq, axis=1)[:, kq]
        order = np.argsort(q); keep = []; seen = set()
        for i in order:
            key = tuple(np.round(allF[i, ::max(1, K//64)], 6))
            if key in seen: continue
            seen.add(key); keep.append(i)
            if len(keep) >= width: break
        Fs = allF[keep]
        if (n+1) in (64, 128, 256, 512, 1024, 2048): hist[n+1] = q[order[0]]
    vals = " ".join(f"{k}:{v:.3e}" for k, v in hist.items())
    v = hist[max(hist)]
    print(f"eta={eta:.2f}  {vals}   log2(1/tau)={np.log2(1/v):.2f}  eta*log2(1/tau)={eta*np.log2(1/v):.2f}", flush=True)
