"""Control experiment ([num]): small-ball search at level c for 0/1 weights.
Parity heuristic: flat arcs are possible only at F = +-1 (odd levels); level 0 forbidden."""
import numpy as np, sys
from toy_smallball import basis, quant
from toy_smallball2 import anneal

def anneal_level(B, m, eta, c, steps, T0, rng):
    N, K = B.shape
    F = 1 + m @ B - c
    q = quant(F, eta); best = (q, m.copy())
    for it in range(steps):
        T = T0*(1 - it/steps) + 1e-6
        k = rng.integers(N)
        Fn = F + (1-2*m[k])*B[k]
        qn = quant(Fn, eta)
        if qn < q or rng.random() < np.exp(-(qn-q)/T):
            F = Fn; m[k] = 1-m[k]; q = qn
            if q < best[0]: best = (q, m.copy())
    return best

if __name__ == '__main__':
    M = 0.5; eta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.10
    rng = np.random.default_rng(3)
    for c in [0.0, 1.0, -1.0, 0.5, 2.0]:
        out = []
        for N in [32, 64, 128]:
            K = 8*N
            t, B = basis(N, K, M)
            best = 1e9
            # starts: random, and periodic patterns of small period
            starts = [rng.integers(0, 2, N).astype(float) for _ in range(4)]
            for P in [2, 3, 4, 5, 6]:
                for r in range(1, 2**P):
                    pat = np.array([(r >> (k % P)) & 1 for k in range(1, N+1)], float)
                    starts.append(pat)
            # pick the 6 best starts by quantile, anneal them
            sc = sorted(starts, key=lambda m: quant(1 + m @ B - c, eta))[:6]
            for m0 in sc:
                q, m = anneal_level(B, m0.copy(), eta, c, 40*N, 0.02, rng)
                best = min(best, q)
            out.append(best)
        print(f"level c={c:+.1f}: tau*(N=32,64,128) = " + "  ".join(f"{x:.2e}" for x in out), flush=True)
