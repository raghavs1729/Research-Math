#!/usr/bin/env python3
"""
explore_model_arc.py -- [num] exploration for DENSITY.md Remark 7.2 (no rigor).

Question: in the positive-weight cyclic model
    F(t) = 1 + sum_{n>=1} c_n (cos nt + (n/M) sin nt),  c_n >= 0,  sigma^2 = 2 sum c_n,
can F vanish on an arc J while sigma^2 < 2 ?  (Psi = 1 + sum c_n cos nt, F = Psi - Psi'/M.)
We solve the nonnegative least-squares problem  min int_J F^2  over c_1..c_N >= 0 for
growing N and report the residual and sum c_n.  A residual tending to 0 with sum c_n
bounded below 1 suggests that such models exist (so sigma^2 < 2 and positivity alone do
not exclude atoms); a residual bounded away from 0 would suggest an obstruction.
"""
import numpy as np
from scipy.optimize import nnls
M = 0.5
for (a, b) in ((1.0, 1.2), (1.0, 1.5), (2.0, 2.3), (0.3, 0.5)):
    t = np.linspace(a, b, 400)
    wq = np.full_like(t, (b - a) / len(t))
    print("arc J = [%.2f, %.2f], M = %.2f" % (a, b, M))
    for N in (25, 50, 100, 200, 400, 800):
        n = np.arange(1, N + 1)
        A = np.cos(np.outer(t, n)) + (n / M) * np.sin(np.outer(t, n))
        # weight columns by 1/n to condition; c_n = d_n / n
        As = A / n
        d, res = nnls(As * np.sqrt(wq)[:, None], -np.sqrt(wq), maxiter=20000)
        c = d / n
        F = 1 + A @ c
        print("   N=%4d  rms_J(F)=%.3e  max_J|F|=%.3e  sum c_n=%.4f (sigma^2=%.4f)  sum c_n^2 n^2=%.3f  #nonzero=%d"
              % (N, np.sqrt(np.sum(wq * F * F) / (b - a)), np.abs(F).max(), c.sum(), 2 * c.sum(),
                 np.sum(c * c * n * n), np.count_nonzero(c)))

# ---- second part: solve on a coarse grid, evaluate on a 10x finer grid; stability of the
# ---- minimiser; and the same with the extra constraint that weights are bounded (c_n <= C/n^2-type)
print()
print("fine-grid check (fit on 300 points, evaluate on 6000 points)")
from scipy.optimize import lsq_linear
for (a, b) in ((1.0, 1.2), (2.0, 2.3)):
    t = np.linspace(a, b, 300); tf = np.linspace(a, b, 6000)
    prev = None
    for N in (100, 200, 400, 800, 1600):
        n = np.arange(1, N + 1)
        A = (np.cos(np.outer(t, n)) + (n / M) * np.sin(np.outer(t, n))) / n
        d, res = nnls(A, -np.ones_like(t), maxiter=50000)
        c = d / n
        Ff = 1 + (np.cos(np.outer(tf, n)) + (n / M) * np.sin(np.outer(tf, n))) @ c
        # weights in the normalisation F = 1 + Re sum w_n 4M/(M+in) e^{int}:  c_n = w_n 4M^2/(M^2+n^2)
        w = c * (M * M + n * n) / (4 * M * M)
        chg = None if prev is None else float(np.abs(c[:len(prev)] - prev).sum())
        print("   J=[%.1f,%.1f] N=%4d  max|F| on fine grid = %.2e   sum c_n = %.4f   sum c_n^2 n^2 = %.4f   max weight w_n = %.3f   l1 change vs previous N: %s"
              % (a, b, N, np.abs(Ff).max(), c.sum(), np.sum(c * c * n * n), w.max(), "-" if chg is None else "%.3f" % chg))
        prev = c
