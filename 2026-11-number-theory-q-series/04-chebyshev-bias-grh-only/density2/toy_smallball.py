"""Toy question, small-ball form (numerical, [num]).

F_N(t) = 1 + 2 Re sum_{n<=N} m_n a_n e^{int},  a_n = 2M/(M+in).
By compactness (PROGRESS.md 3.5): 0/1 weights admit F = 0 on a set of measure eta*2pi
iff  tau*_eta(N) := min_m  quantile_eta |F_N|  is O(N^{-1/2}) along N -> infinity.
Heuristic search: relaxation (real weights in [0,W]) on the adaptive sublevel set,
rounding, then greedy bit flips on sum_S F^2 with S = current eta-sublevel set.
Prints tau*(N) and sqrt(N)*tau*(N) for 0/1 and for the relaxation.
"""
import numpy as np, sys
from scipy.optimize import lsq_linear

def basis(N, K, M):
    t = 2*np.pi*np.arange(K)/K
    n = np.arange(1, N+1)
    a = 2*M/(M+1j*n)
    B = 2*np.real(a[:,None]*np.exp(1j*np.outer(n, t)))   # N x K
    return t, B

def quant(F, eta):
    return np.quantile(np.abs(F), eta)

def relax(B, eta, W, iters=6, S0=None):
    N, K = B.shape
    S = S0 if S0 is not None else np.arange(K)[:int(eta*K)]
    for _ in range(iters):
        res = lsq_linear(B[:,S].T, -np.ones(len(S)), bounds=(0, W), lsmr_tol='auto', max_iter=2000)
        w = res.x
        F = 1 + w @ B
        S = np.argsort(np.abs(F))[:int(eta*K)]
    return w, F, S

def greedy01(B, m, eta, sweeps=40, rng=None):
    N, K = B.shape
    F = 1 + m @ B
    best = (quant(F, eta), m.copy())
    for sw in range(sweeps):
        S = np.argsort(np.abs(F))[:int(eta*K)]
        BS = B[:, S]
        G = BS @ BS.T
        g = 2*BS @ F[S]                 # gradient-like: d/dm of sum F^2
        s = 1 - 2*m                      # +1 add, -1 remove
        improved = False
        for step in range(4*N):
            d = s*g + np.diag(G)
            k = np.argmin(d)
            if d[k] >= -1e-15:
                break
            # flip k
            g = g + 2*s[k]*G[:, k]
            F = F + s[k]*B[k]
            m[k] = 1 - m[k]; s[k] = -s[k]
            improved = True
        q = quant(F, eta)
        if q < best[0]:
            best = (q, m.copy())
        if not improved:
            # random kick
            idx = rng.choice(N, size=max(1, N//50), replace=False)
            for k in idx:
                F = F + (1-2*m[k])*B[k]; m[k] = 1-m[k]
    return best

if __name__ == '__main__':
    M = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    eta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
    rng = np.random.default_rng(1)
    print(f"M={M} eta={eta}")
    print(" N   tau_relax(W=1)  sqrtN*  tau_01   sqrtN*   N*tau_01")
    for N in [32, 64, 128, 256, 512]:
        K = 12*N
        t, B = basis(N, K, M)
        # start S: an arc around t=2.15 (where Prop 7.3 puts its flat arc)
        c = int(K*2.15/(2*np.pi)); L = int(eta*K)
        S0 = (np.arange(c-L//2, c+L//2)) % K
        w, Fr, S = relax(B, eta, 1.0, S0=S0)
        qr = quant(Fr, eta)
        bestq = 1e9
        for trial in range(4):
            m0 = (w > rng.random(N)).astype(float) if trial else np.round(w)
            q, m = greedy01(B, m0.copy(), eta, rng=rng)
            bestq = min(bestq, q)
        print(f"{N:4d}  {qr:.3e}  {np.sqrt(N)*qr:.3f}   {bestq:.3e}  {np.sqrt(N)*bestq:.3f}  {N*bestq:.3f}", flush=True)
