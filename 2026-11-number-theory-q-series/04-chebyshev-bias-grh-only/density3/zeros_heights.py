#!/usr/bin/env python3
"""zeros_heights.py -- [num] MAHLER.md, Heuristic 7.2: heights of the zeros of P_N(z) = z^N F_N(z) for cyclic
models at M = 1/2 (roots by numpy.roots).  delta_N = log(1/|C_N|) - sum_{|z|>1} log|z| (Jensen).
Output: zeros_heights_results.txt"""
import numpy as np, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
fout = open(os.path.join(HERE, "zeros_heights_results.txt"), "w")
M = 0.5
rng = np.random.default_rng(3)
for N in [100, 300, 1000]:
    for name in ["random01", "all1", "beam0.15"]:
        if name == "random01":
            w = rng.integers(0, 2, N+1).astype(float); w[0] = 0; w[N] = 1
        elif name == "all1":
            w = np.ones(N+1); w[0] = 0
        else:
            if N != 300: continue
            m01 = np.load(os.path.join(HERE, "beam_best_eta0.15_N512.npy")).astype(float)[:N]; w = np.r_[0.0, m01]
            top = np.nonzero(w)[0].max(); w = w[:top+1]
        Nn = len(w)-1
        n = np.arange(Nn+1); C = np.zeros(Nn+1, complex); C[1:] = w[1:]*2*M/(M+1j*n[1:]); C[0] = 1
        p = np.zeros(2*Nn+1, complex); p[Nn] = 1
        for k in range(1, Nn+1): p[Nn+k] = C[k]; p[Nn-k] = np.conj(C[k])
        r = np.roots(p[::-1])
        out = np.abs(r[np.abs(r) > 1+1e-9])
        on = np.sum(np.abs(np.abs(r)-1) <= 1e-9)
        h = (out-1)*Nn/math.log(Nn)
        S = np.sum(np.log(out)); J = math.log(1/abs(C[Nn]))
        line = (f"N={Nn:5d} {name:9s} #outer={len(out):5d} #on-circle={on:4d} sum log|z| (outer)={S:.3f} vs log(1/|C_N|)={J:.3f} -> delta={J-S:.3f}; "
              f"N(|z|-1)/logN: median={np.median(h):.2f} q10={np.quantile(h,0.1):.2f} q90={np.quantile(h,0.9):.2f}")
        print(line); fout.write(line + "\n")
fout.close()
