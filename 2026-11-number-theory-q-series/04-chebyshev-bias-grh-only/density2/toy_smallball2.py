"""Toy question, small-ball form, stronger search ([num]).
Simulated annealing on the eta-quantile of |F_N| over m in {0,1}^N, several starts
(relaxation-rounded, random, periodic patterns). Also records the structure of the
sublevel set (number of components at resolution 2pi/K)."""
import numpy as np, sys
from toy_smallball import basis, relax, quant

def anneal(B, m, eta, steps, T0, rng):
    N, K = B.shape
    F = 1 + m @ B
    q = quant(F, eta); best = (q, m.copy(), F.copy())
    for it in range(steps):
        T = T0*(1 - it/steps) + 1e-6
        k = rng.integers(N)
        Fn = F + (1-2*m[k])*B[k]
        qn = quant(Fn, eta)
        if qn < q or rng.random() < np.exp(-(qn-q)/T):
            F = Fn; m[k] = 1-m[k]; q = qn
            if q < best[0]:
                best = (q, m.copy(), F.copy())
    return best

def ncomp(F, tau):
    s = (np.abs(F) <= tau).astype(int)
    return int(np.sum(np.abs(np.diff(np.r_[s, s[0]])))//2)

if __name__ == '__main__':
    M = float(sys.argv[1]); eta = float(sys.argv[2])
    rng = np.random.default_rng(7)
    print(f"M={M} eta={eta}")
    for N in [32, 64, 128, 256]:
        K = 8*N
        t, B = basis(N, K, M)
        c = int(K*2.15/(2*np.pi)); L = int(eta*K)
        S0 = (np.arange(c-L//2, c+L//2)) % K
        w, Fr, S = relax(B, eta, 1.0, S0=S0)
        starts = [np.round(w)] + [(w > rng.random(N)).astype(float) for _ in range(3)] \
                 + [rng.integers(0, 2, N).astype(float) for _ in range(3)]
        res = []
        for m0 in starts:
            q, m, F = anneal(B, m0.copy(), eta, steps=60*N, T0=0.02, rng=rng)
            res.append((q, m, F))
        q, m, F = min(res, key=lambda r: r[0])
        print(f"N={N:4d} K={K} tau01={q:.4e} sqrtN*tau={np.sqrt(N)*q:.3f} "
              f"components={ncomp(F,q)} density_m={m.mean():.2f} relax={quant(Fr,eta):.2e}", flush=True)
