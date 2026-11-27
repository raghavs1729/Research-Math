#!/usr/bin/env python3
"""
jensen_newton_check.py -- [num] sanity checks for MAHLER.md, Lemma 2.1 (Jensen-Newton bound),
Lemma 2.2 (small-ball from the Mahler measure) and Example 3.7 (the model Gamma_0 = {b + n d}).

(1) Random Laurent polynomials P on T^2 and T^3 with random integer supports and complex
    coefficients: m(P) := int log|P| (grid quadrature) versus max over vertices v of the
    Newton polytope of log|c_v| (vertices by linear programming).  Lemma 2.1 says
    m(P) >= that maximum.
(2) The small-ball inequality  lambda{|P| < delta} log(1/delta) <= -m(P) + (1/2) log+ ||P||_2^2 + 1/(2e).
(3) The finite models of Example 3.7 on T^2: frequencies b + n d (0 <= n <= N), integer
    multiplicities, weight M = 1/2; F_N(z1, z2) = 1 + 2 Re( z1 * sum_n m_n r(b+nd) z2^n ),
    r(g) = 2M/(M + i g).  Lemma 2.1 with the vertex b gives m(F_N - c) >= log|m_0 r(b)|,
    uniformly in N.  The script prints m(F_N - c) for several N and c.
Output: jensen_newton_check_results.txt
"""
import math, os
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, "jensen_newton_check_results.txt"), "w")
def log(s=""):
    print(s); out.write(s + "\n"); out.flush()

def vertices(S):
    """indices i such that S[i] is not in the convex hull of the other points (LP feasibility)"""
    S = np.asarray(S, float); n, d = S.shape; V = []
    for i in range(n):
        others = np.delete(S, i, axis=0)
        A_eq = np.vstack([others.T, np.ones(n - 1)]); b_eq = np.r_[S[i], 1.0]
        res = linprog(np.zeros(n - 1), A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * (n - 1), method="highs")
        if res.status != 0:          # infeasible -> vertex
            V.append(i)
    return V

def grid_values(S, c, K, shift):
    """values of P = sum c_j z^{S_j} on the shifted grid (2 pi k/K + shift) of the torus T^d
    (S: integer points, may be negative; the shift avoids sampling the zero set exactly)"""
    d = S.shape[1]
    A = np.zeros((K,) * d, dtype=complex)
    for s, cj in zip(S, c):
        A[tuple(np.mod(s, K))] += cj * np.exp(1j * float(np.dot(s, shift)))
    return np.fft.ifftn(A) * K ** d

rng = np.random.default_rng(20261127)
log("(1)+(2) random Laurent polynomials: Lemma 2.1 and Lemma 2.2")
log("  d  #terms   m(P)      max_v log|c_v|   slack    | small-ball: max over delta of [lhs - rhs] (must be <= 0)")
worst_slack = 1e9; worst_sb = -1e9
for trial in range(24):
    d = 2 if trial < 16 else 3
    K = 256 if d == 2 else 48
    nt = int(rng.integers(3, 12))
    box = 6 if d == 2 else 3
    S = np.unique(rng.integers(-box, box + 1, size=(nt, d)), axis=0)
    if len(S) < 2:
        continue
    c = (rng.normal(size=len(S)) + 1j * rng.normal(size=len(S))) * np.exp(rng.normal(size=len(S)))
    V = vertices(S)
    vb = max(math.log(abs(c[i])) for i in V)
    P = grid_values(S, c, K, rng.uniform(0, 2 * np.pi / K, size=d)); a = np.abs(P)
    mP = float(np.mean(np.log(np.maximum(a, 1e-300))))
    L2 = float(np.mean(a ** 2))
    sb = -1e9
    for delta in [0.5, 0.2, 0.1, 0.03, 0.01, 1e-3]:
        eta = float(np.mean(a < delta))
        lhs = eta * math.log(1 / delta); rhs = -mP + 0.5 * max(0.0, math.log(L2)) + 1 / (2 * math.e)
        sb = max(sb, lhs - rhs)
    worst_slack = min(worst_slack, mP - vb); worst_sb = max(worst_sb, sb)
    log(f"  {d}  {len(S):3d}   {mP:+9.4f}   {vb:+9.4f}      {mP - vb:+8.4f} | {sb:+.4f}")
log(f"  minimum of m(P) - vertex bound over all trials = {worst_slack:+.4f} (Lemma 2.1: >= 0 up to quadrature error)")
log(f"  maximum of small-ball defect = {worst_sb:+.4f} (Lemma 2.2: <= 0)")
log()

log("(3) Example 3.7: frequencies b + n d, b = 6.0209, d = 2*sqrt(2) (Q-independent), M = 1/2, on T^2")
M = 0.5; b = 6.020948904697597; dd = 2 * math.sqrt(2)
r = lambda g: 2 * M / (M + 1j * g)
log(f"  vertex bound log|m_0 r(b)| (m_0 = 1) = {math.log(abs(r(b))):+.4f}")
for mult_kind in ["all ones", "random in {1,2,3}"]:
    for N in [4, 16, 64, 200]:
        n = np.arange(N + 1)
        mm = np.ones(N + 1) if mult_kind == "all ones" else rng.integers(1, 4, N + 1).astype(float)
        mm[0] = 1.0
        coef = mm * r(b + n * dd)
        K1, K2 = 512, max(256, 4 * N)
        A = np.zeros((K1, K2), dtype=complex)
        s1, s2 = math.pi / K1 * 0.7371, math.pi / K2 * 0.4129     # grid offsets (avoid singular lines)
        for k, ck in zip(n, coef):
            A[1, k % K2] += ck * np.exp(1j * (s1 + k * s2))       # z1^1 z2^n
        G = np.fft.ifft2(A) * K1 * K2
        Fv = 1 + 2 * G.real
        row = []
        for cl in [0.0, 1.0, 2.0]:
            row.append(float(np.mean(np.log(np.maximum(np.abs(Fv - cl), 1e-300)))))
        sig2 = float(np.sum(mm ** 2 * 8 * M * M / (M * M + (b + n * dd) ** 2)))
        log(f"  {mult_kind:18s} N={N:4d}  V={sig2:.3f}   m(F_N - c) for c = 0, 1, 2: " +
            ", ".join(f"{x:+.4f}" for x in row))
        if N == 200:
            for cl in [0.0, 1.0]:
                a = np.abs(Fv - cl)
                cst = math.log(1 / abs(r(b))) + 0.5 * max(0.0, math.log((1 - cl) ** 2 + sig2)) + 1 / (2 * math.e)
                log("      small ball c=%g: " % cl + ", ".join(
                    f"eps={e:.0e}: lambda={float(np.mean(a < e)):.2e} <= bound {cst / math.log(1 / e):.3f}"
                    for e in [1e-1, 1e-2, 1e-3]))
log("  (all values must be >= the vertex bound, up to quadrature error ~3e-4.  For c = 1,")
log("   m(F_N - 1) = m(h_N) with h_N(z) = sum_n m_n r(b + n d) z^n (since E log|2 cos| = 0), and this")
log("   equals log|h_N(0)| = the vertex bound exactly when h_N has no zeros in the open unit disk.)")
out.close()

# (4) explicit constant of Theorem 3.3(b) for chi_{-4}, M = 1/2, gamma_0 = gamma_1, under V(gamma_1),
#     assuming m_{gamma_1} = 1 and all zeros simple (then V_{1/2} = sigma^2_{1/2} = 0.155568, DENSITY.md section 9).
out = open(os.path.join(HERE, "jensen_newton_check_results.txt"), "a")
g1 = 6.020948904697597; V12 = 0.155568
C = math.log(abs(0.5 + 1j * g1) / (2 * 0.5)) + 0.5 * math.log(1 + V12) + 1 / (2 * math.e)
log()
log(f"(4) Theorem 3.3(b) constant for chi_-4, M = 1/2, gamma_0 = {g1:.6f}: log|M+i g|/(2M) = {math.log(abs(0.5 + 1j * g1)):.4f}, "
    f"(1/2)log(1+V) = {0.5 * math.log(1 + V12):.4f}, 1/(2e) = {1 / (2 * math.e):.4f}; total C = {C:.4f}")
log(f"    nu_(1/2)((-eps, eps)) <= C/log(1/eps): eps = 1e-10 -> {C / math.log(1e10):.4f}, eps = 1e-20 -> {C / math.log(1e20):.4f}")
out.close()
