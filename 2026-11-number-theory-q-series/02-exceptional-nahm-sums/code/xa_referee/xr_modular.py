"""Asymptotic test of modularity / prefactor near q = 1.

For a weight-0 modular function f = q^{-c/24} F(q) with positive coefficients (c_eff = c):
   g(t) := log F(e^{-t}) - pi^2 c /(6t) + c t/24  =  const + O(exp(-4 pi^2 Delta / t)).
A wrong c in the prefactor gives a drift (c - c')t/24, a non-modular series gives power corrections.
We evaluate F(e^{-t}) by direct summation (float64, log-domain), fit g(t) on a grid by a polynomial in t and report.
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fractions import Fraction as Fr
from xr_lib import sw_cartan, coxeter, nonneg_points


def nahm_matrix(G, k):
    l = len(G)
    n = l * (k - 1)
    M = np.zeros((n, n))
    for i in range(l):
        for j in range(l):
            for s in range(1, k):
                for t in range(1, k):
                    M[i*(k-1)+s-1, j*(k-1)+t-1] = G[i][j] * (min(s, t) - s * t / k)
    return M


def logF(M, dvec, t, B, pts_cache={}):
    key = (M.tobytes(), tuple(dvec), B)
    if key not in pts_cache:
        P = np.array(nonneg_points(M, B), dtype=np.int64)
        Qv = 0.5 * np.einsum('ij,jk,ik->i', P, M, P)
        pts_cache[key] = (P, Qv)
    P, Qv = pts_cache[key]
    nmax = P.max()
    lp = {}
    for d in set(dvec):
        x = -np.log1p(-np.exp(-t * d * np.arange(1, nmax + 1)))
        lp[d] = np.concatenate([[0.0], np.cumsum(x)])
    E = -t * Qv
    for col, d in enumerate(dvec):
        E = E + lp[d][P[:, col]]
    m = E.max()
    return m + math.log(np.exp(E - m).sum()), E.max(), E[Qv > 0.9 * B].max() if (Qv > 0.9 * B).any() else -np.inf, len(P)


def test(label, G, dvec, k, c, ts, B, out=sys.stdout):
    """free fit log F(e^-t) = a/t + b0 + b1 t + b2 t^2 + b3 t^3, and constrained fits
    log F + c' t/24 = a/t + b0 for c' = c, c +- 1/10, c +- 1 (residual RMS)."""
    M = nahm_matrix(G, k)
    dd = [dvec[i] for i in range(len(G)) for _ in range(k - 1)]
    rows = [logF(M, dd, t, B) for t in ts]
    tt = np.array(ts, dtype=float); lf = np.array([r[0] for r in rows])
    V2 = np.column_stack([1/tt, np.ones_like(tt), tt, tt**2, tt**3])
    co2, *_ = np.linalg.lstsq(V2, lf, rcond=None)
    tail = max(r[2] - r[1] for r in rows)
    out.write(f"{label}: c_SW={c}; t in [{tt[0]:.3f},{tt[-1]:.3f}] ({len(tt)} pts); B={B}, #terms={rows[0][3]}, tail/max<=e^{tail:.0f}\n")
    out.write(f"   free fit: a={co2[0]:.10f} => c_eff=6a/pi^2={6*co2[0]/math.pi**2:.10f}; b1={co2[2]:.10f} (-c_SW/24={-float(c)/24:.10f}); b2={co2[3]:.1e}; b3={co2[4]:.1e}\n")
    V = np.column_stack([1/tt, np.ones_like(tt)])
    res = []
    for cp in [c, c + Fr(1, 10), c - Fr(1, 10), c + 1, c - 1]:
        y = lf + float(cp) * tt / 24
        co, *_ = np.linalg.lstsq(V, y, rcond=None)
        rms = float(np.sqrt(np.mean((V @ co - y) ** 2)))
        res.append(f"c'={str(cp)}: rms={rms:.1e}")
    out.write("   constrained fits (a/t+b0 after adding c't/24): " + "; ".join(res) + "\n")
    out.flush()
    return co2


def cSW(X, l, k):
    C, D = sw_cartan(X, l)
    h = coxeter(X, l)
    return Fr(sum(D) * (k - 1) * h, h + k)


if __name__ == '__main__':
    out = sys.stdout
    jobs = sys.argv[1:] or ['C3', 'B3', 'D5', 'T1A5']
    for j in jobs:
        if j == 'C3':
            X, l, k, ts, B = 'C', 3, 2, np.linspace(0.10, 0.30, 9), 600
        elif j == 'B3':
            X, l, k, ts, B = 'B', 3, 2, np.linspace(0.10, 0.30, 9), 700
        elif j == 'D5':
            X, l, k, ts, B = 'D', 5, 2, np.linspace(0.22, 0.45, 9), 220
        elif j == 'T1A5':
            X, l, k, ts, B = 'T', 1, 6, np.linspace(0.22, 0.45, 9), 220
        C, D = sw_cartan(X, l)
        G = [[C[i][j] * D[j] for j in range(l)] for i in range(l)]
        c = cSW(X, l, k)
        test(f"({X}_{l},A_{k-1})", G, D, k, c, ts, B, out=out)
