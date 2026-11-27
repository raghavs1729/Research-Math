"""Referee test 5 (independent re-implementation): the [num] claims of TOY.md section 5.

tau*_eta(N) := min over m in {0,1}^N of the eta-quantile of |F_{m,N}| on T, M = 1/2, where
F_{m,N}(t) = 1 + sum_{n<=N} m_n (4M/(M^2+n^2)) (M cos nt + n sin nt).
The eta-quantile on a K-point grid is the ceil(eta K)-th smallest value of |F|.

Part 1 (exhaustive, N <= 20): tau* on a coarse grid K = 16N, and the TRUE quantile of the
coarse-grid optimiser re-evaluated on a fine grid (K = 4096 N); also a full exhaustive search on
a fine grid (K = 256 N) for N <= 14, to see whether the optimiser changes.
Part 2 (beam search, own implementation, width 200, coarse grid K = 8N during the search):
N = 256 and N = 1024 for eta in {0.05, 0.25, 0.40}; the final sequence is re-evaluated on a
grid 64 times finer.  Usage: python3 t5_numerics_check.py [part1|part2|all]
"""
import math
import sys
import time

import numpy as np

M = 0.5
ETAS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
out = []


def log(s=""):
    print(s, flush=True)
    out.append(s)


def basis(N, K):
    t = 2 * np.pi * np.arange(K) / K
    n = np.arange(1, N + 1)[:, None]
    return (4 * M / (M * M + n * n)) * (M * np.cos(n * t) + n * np.sin(n * t))  # N x K


def quantiles(absvals, etas):
    """absvals: (rows, K). Returns (rows, len(etas)) of eta-quantiles."""
    K = absvals.shape[1]
    ks = [max(0, math.ceil(e * K) - 1) for e in etas]
    part = np.partition(absvals, sorted(set(ks)), axis=1)
    return np.stack([part[:, k] for k in ks], axis=1)


def eval_quantiles(m, K, etas):
    B = basis(len(m), K)
    F = 1 + np.asarray(m, dtype=float) @ B
    return quantiles(np.abs(F)[None, :], etas)[0]


def exhaustive(N, K, etas, low_bits=12):
    B = basis(N, K)
    lb = min(low_bits, N)
    # all sums over the low bits (bit i <-> n = i+1)
    L = np.zeros((1, K))
    for i in range(lb):
        L = np.concatenate([L, L + B[i]], axis=0)
    best = np.full(len(etas), np.inf)
    arg = [None] * len(etas)
    hb = N - lb
    for h in range(2 ** hb):
        base = 1.0 + sum(B[lb + j] for j in range(hb) if (h >> j) & 1) if hb else np.ones(K)
        q = quantiles(np.abs(L + base), etas)
        imin = np.argmin(q, axis=0)
        for e in range(len(etas)):
            v = q[imin[e], e]
            if v < best[e]:
                best[e] = v
                lowpat = imin[e]
                arg[e] = [((lowpat >> i) & 1) for i in range(lb)] + [((h >> j) & 1) for j in range(hb)]
    return best, arg


def part1():
    log("Part 1: exhaustive search, M = 1/2, weights {0,1}")
    log("  columns: eta | tau* on grid 16N | true quantile (grid 4096N) of that optimiser")
    for N in [10, 14, 18, 20]:
        t0 = time.time()
        best, arg = exhaustive(N, 16 * N, ETAS)
        fine = [eval_quantiles(arg[e], 4096 * N, [ETAS[e]])[0] for e in range(len(ETAS))]
        log(f"  N = {N} ({time.time() - t0:.0f}s)")
        for e, eta in enumerate(ETAS):
            log(f"    eta={eta:.2f}: {best[e]:.4e} | {fine[e]:.4e}   m = {''.join(map(str, arg[e]))}")
    log("  full exhaustive on a fine grid K = 256N (N = 14):")
    bestc, argc = exhaustive(14, 16 * 14, ETAS)
    bestf, argf = exhaustive(14, 256 * 14, ETAS)
    for e, eta in enumerate(ETAS):
        log(f"    eta={eta:.2f}: coarse-grid optimum {bestc[e]:.4e}, fine-grid optimum {bestf[e]:.4e}, same optimiser: {argc[e] == argf[e]}")


def beam(N, K, eta, width=200):
    B = basis(N, K)
    states_m = np.zeros((1, 0), dtype=np.int8)
    states_F = np.ones((1, K))
    kq = max(0, math.ceil(eta * K) - 1)
    for i in range(N):
        F0 = states_F
        F1 = states_F + B[i]
        allF = np.concatenate([F0, F1], axis=0)
        sc = np.partition(np.abs(allF), kq, axis=1)[:, kq]
        keep = np.argsort(sc, kind="stable")[:width]
        nb = states_F.shape[0]
        bits = (keep >= nb).astype(np.int8)
        parent = keep % nb
        states_m = np.concatenate([states_m[parent], bits[:, None]], axis=1)
        states_F = allF[keep]
        best_sc = sc[keep[0]]
    return states_m[0], best_sc


def part2():
    log("Part 2: beam search (own implementation, width 200), M = 1/2, weights {0,1}")
    log("  TOY.md table (beam, N = 1024): eta 0.40 -> 0.282, 0.25 -> 0.0615, 0.05 -> 0.0017")
    for N in [256, 1024]:
        for eta in [0.40, 0.25, 0.05]:
            t0 = time.time()
            m, sc = beam(N, 8 * N, eta)
            q_fine = eval_quantiles(m, 512 * N, [eta])[0]
            q16 = eval_quantiles(m, 16 * N, [eta])[0]
            even_only = bool(np.all(m[0::2] == 0))  # m[0] <-> n = 1, so m[0::2] are the odd n
            log(f"  N={N}, eta={eta:.2f}: beam score (grid 8N) {sc:.4e}; same m on grid 16N {q16:.4e}; "
                f"on grid 512N {q_fine:.4e}; ones = {int(m.sum())}, only even n: {even_only} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("part1", "all"):
        part1()
    if which in ("part2", "all"):
        part2()
    suffix = "" if which == "all" else "_" + which
    with open(__file__.replace(".py", f"{suffix}_output.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
