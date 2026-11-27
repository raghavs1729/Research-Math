"""deficit_beam.py -- [num] MAHLER.md, Numerical check 7.3.
Beam search maximising the Mahler deficit delta_N = -mean log|F_N| over 0/1 weights on the frequencies
NMIN..Nmax (cyclic toy, M = 1/2); the objective is evaluated on an offset grid of K = 16 Nmax points,
log clipped at 1e-12.  Usage: python3 deficit_beam.py Nmax width [NMIN]
Runs recorded in deficit_beam_results.txt:  "4096 40"  and  "4096 40 13" (support n >= 13, so that
sigma^2 <= sum_{n>=13} 8M^2/(M^2+n^2) = 0.1598)."""
import numpy as np, sys, math
M = 0.5; Nmax = int(sys.argv[1]); W = int(sys.argv[2]); K = 16*Nmax
t = 2*np.pi*(np.arange(K)+0.37)/K
n = np.arange(1, Nmax+1); a = 2*M/(M+1j*n)
def basis_row(k):  # k = frequency index (1-based)
    return 2*np.real(a[k-1]*np.exp(1j*k*t))
Fs = np.ones((1, K)); ms = np.zeros((1, 0), dtype=np.int8)
def obj(F): return -np.mean(np.log(np.maximum(np.abs(F), 1e-12)), axis=-1)
NMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 1
for k in range(1, Nmax+1):
    B = basis_row(k)
    allF = np.vstack([Fs, Fs + B[None, :]]) if k >= NMIN else np.vstack([Fs, Fs + 1e9])
    allm = np.vstack([np.hstack([ms, np.zeros((len(ms), 1), np.int8)]), np.hstack([ms, np.ones((len(ms), 1), np.int8)])])
    o = obj(allF)
    order = np.argsort(-o)[:W]
    Fs, ms = allF[order], allm[order]
    if k & (k-1) == 0 and k >= 8:
        F = Fs[0]; m = ms[0]
        top = np.nonzero(m)[0].max()+1 if m.any() else 0
        sig2 = float(np.sum(m*8*M*M/(M*M+n[:k]**2)))
        J = math.log(1/abs(a[top-1])) if top else float('nan')
        q = [float(np.mean(np.abs(F) < e)) for e in (1e-1, 1e-2, 1e-3)]
        print(f"N={k:5d} delta={o[order[0]]:.3f} J_top={J:.3f} delta/logN={o[order[0]]/math.log(k):.3f} sigma2={sig2:.3f} "
              f"lambda(|F|<1e-1,1e-2,1e-3)={q[0]:.3f},{q[1]:.4f},{q[2]:.5f}", flush=True)
import os
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"deficit_best_N{Nmax}_W{W}_n{NMIN}.npy"), ms[0])
