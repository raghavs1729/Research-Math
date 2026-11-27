"""Referee test 5b: exhaustive tau*_eta(N) for N = 22 (M = 1/2, weights {0,1}), grids 8N and 16N,
to compare with the exhaustive values recorded by the authors (PROGRESS.md section 6.3:
tau* ~ 6.0e-3, 1.37e-2, 3.2e-2, 6.0e-2, 0.16 for eta = .05, .10, .15, .20, .30, 'N <= 22').
Re-uses only the functions of t5_numerics_check.py (referee code)."""
import time

from t5_numerics_check import ETAS, eval_quantiles, exhaustive

out = []
for K_per_N in [8, 16]:
    N = 22
    t0 = time.time()
    best, arg = exhaustive(N, K_per_N * N, ETAS)
    out.append(f"N = {N}, grid {K_per_N}N ({time.time() - t0:.0f}s)")
    for e, eta in enumerate(ETAS):
        fine = eval_quantiles(arg[e], 4096 * N, [eta])[0]
        out.append(f"  eta={eta:.2f}: tau* = {best[e]:.4e}; true quantile of that optimiser (grid 4096N) = {fine:.4e};"
                   f" m = {''.join(map(str, arg[e]))}")
    print("\n".join(out[-(len(ETAS) + 1):]), flush=True)
with open(__file__.replace(".py", "_output.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
