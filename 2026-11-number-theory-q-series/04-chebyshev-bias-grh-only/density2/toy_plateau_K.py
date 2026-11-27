"""Plateau vs eta for weights in {0,...,Kw} ([num]). Entropy heuristic:
eta*log2(1/tau_floor) ~ log2(Kw+1)."""
import numpy as np, sys
from toy_smallball import basis
M = float(sys.argv[1]); Nmax = int(sys.argv[2]); width = int(sys.argv[3]); Kw = int(sys.argv[4])
etas = [float(x) for x in sys.argv[5].split(',')]
K = 8*Nmax
t, B = basis(Nmax, K, M)
for eta in etas:
    kq = int(eta*K)
    Fs = np.ones((1, K)); hist = {}
    for n in range(Nmax):
        allF = np.vstack([Fs + c*B[n][None, :] for c in range(Kw+1)])
        q = np.partition(np.abs(allF), kq, axis=1)[:, kq]
        order = np.argsort(q); keep = []; seen = set()
        for i in order:
            key = tuple(np.round(allF[i, ::max(1, K//64)], 6))
            if key in seen: continue
            seen.add(key); keep.append(i)
            if len(keep) >= width: break
        Fs = allF[keep]
        if (n+1) in (64, 128, 256, 512, 1024): hist[n+1] = q[order[0]]
    v = hist[max(hist)]
    vals = " ".join(f"{k}:{x:.3e}" for k, x in hist.items())
    print(f"Kw={Kw} eta={eta:.2f}  {vals}  eta*log2(1/tau)={eta*np.log2(1/v):.2f}  (pred {np.log2(Kw+1):.2f})", flush=True)
