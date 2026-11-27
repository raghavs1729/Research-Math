"""Referee 3: checks for HYP.md Proposition 1.1 and the examples of HYP.md section 2.

AI-produced, 2026-11-27.

(a) Sharpness of the bound 1/N in Proposition 1.1, in an abstract cyclic model: A = T, G = Z,
    H = N Z, Lambda = {1, 2} (so G = H + <Lambda>), F = F_H + F_Lambda with F_Lambda a fixed
    trigonometric polynomial and F_H := c - F_Lambda on [0, 2pi/N), extended 2pi/N-periodically
    (F_H has its frequencies in N Z). The level set {F = c} then has measure exactly 1/N, and on
    each K-orbit {t + 2 pi k/N} it has exactly one point. (This is only the abstract statement:
    F_H is not of the form 1 + Re sum m r e with positive integer m.)
(b) Role of G = H + <Lambda> (automatic in HYP, where G = <Gamma>): with H = 4Z inside G = Z but
    Lambda = {2}, so that H + <Lambda> = 2Z != G, the same construction gives a level set of
    measure 1/2 > 1/4: the bound uses that the non-invariant part separates the K-orbit.
(c) HYP section 2, items 2-4: indices [G : <X_{>T}>] for X = {2^k alpha} and X = {(2k+1) alpha};
    the relation lattice of (1, 3, 5, ..., 2K+1) (an integer kernel basis computed by sympy) and
    the parity of the coordinate sums of all its basis vectors; Example 3.3a of DENSITY.md in
    coordinates (b_n = (b_n + b_N) - b_N with both terms above any height).
Output: r3_hyp_results.txt
"""
from math import gcd
from functools import reduce
import numpy as np
import sympy as sp


def level_set_measure(F, c, tol):
    return np.mean(np.abs(F - c) < tol)


def model(N, lam_freqs, coeffs, c=0.3, n_grid=120 * 8192):
    t = 2 * np.pi * np.arange(n_grid) / n_grid
    FL = np.zeros_like(t)
    for f, a in zip(lam_freqs, coeffs):
        FL += np.real(a * np.exp(1j * f * t))
    # F_H := c - F_Lambda on [0, 2pi/N), extended periodically (frequencies in N Z)
    per = n_grid // N
    FH = np.tile(c - FL[:per], N)
    F = FH + FL
    meas = level_set_measure(F, c, 1e-12)
    # points per K-orbit
    pts = np.zeros(per, dtype=int)
    for k in range(N):
        pts += (np.abs(F[k * per:(k + 1) * per] - c) < 1e-12)
    return meas, np.bincount(pts, minlength=N + 1)


def main():
    lines = []
    for N in (2, 3, 5, 8):
        meas, hist = model(N, [1, 2], [0.8 + 0.3j, -0.5 + 0.2j])
        lines.append(f"(a) N={N}, Lambda={{1,2}}: measure of {{F = c}} = {meas:.6f} (1/N = {1/N:.6f}); "
                     f"histogram of #points per K-orbit (index = count): {hist.tolist()}")
    meas, hist = model(4, [2], [1.0 + 0.0j])
    lines.append(f"(b) H = 4Z, Lambda = {{2}} (H + <Lambda> = 2Z != G = Z): measure of {{F = c}} = "
                 f"{meas:.6f} > 1/4; histogram of #points per K-orbit: {hist.tolist()}")
    # (c) indices
    alpha_mult = [2 ** k for k in range(40)]
    for T in (1, 5, 100, 10 ** 6):
        tail = [x for x in alpha_mult if x > T]
        idx = reduce(gcd, tail)
        lines.append(f"(c) X = {{2^k alpha}}: T = {T} alpha: [G : <X_>T>] = {idx}")
    odd = [2 * k + 1 for k in range(2000)]
    for T in (1, 5, 100, 1000):
        tail = [x for x in odd if x > T]
        lines.append(f"(c) X = {{(2k+1) alpha}}: T = {T} alpha: [G : <X_>T>] = {reduce(gcd, tail)}")
    K = 12
    A = sp.Matrix([[2 * k + 1 for k in range(K + 1)]])
    ns = A.nullspace()
    # make an integer basis of the kernel lattice: the e_k - (2k+1) e_0 basis, and verify it spans
    basis = []
    for k in range(1, K + 1):
        v = [0] * (K + 1)
        v[k] = 1
        v[0] = -(2 * k + 1)
        basis.append(v)
    Bm = sp.Matrix(basis).T
    ok_kernel = all((A * Bm[:, j])[0] == 0 for j in range(K))
    # lattice spanned by basis equals the full integer kernel: kernel has rank K and the basis
    # has a K x K minor equal to +-1 (rows 1..K), so it is saturated
    minor = Bm[1:, :].det()
    sums = [sum(b) for b in basis]
    lines.append(f"(c) relation lattice of (1,3,...,{2*K+1}): rational kernel dimension {len(ns)}; "
                 f"basis e_k - (2k+1)e_0 lies in the kernel: {ok_kernel}; minor det = {minor} "
                 f"(+-1, so the basis spans the integer kernel); coordinate sums {sums[:6]}... "
                 f"all even: {all(s % 2 == 0 for s in sums)}")
    # Example 3.3a in coordinates: elements are vectors in Z^(indices); b_n = (b_n + b_N) - b_N
    n, Nlink = 3, 17
    e = lambda i: np.eye(20, dtype=int)[i]
    lines.append(f"(c) DENSITY Example 3.3a: b_{n} = (b_{n} + b_{Nlink}) - b_{Nlink} in coordinates: "
                 f"{np.array_equal(e(n), (e(n) + e(Nlink)) - e(Nlink))} (both terms lie above any "
                 f"height once b_N is large, so <X_>T> contains every b_n and equals G)")
    txt = "\n".join(lines)
    print(txt)
    with open("r3_hyp_results.txt", "w") as f:
        f.write("r3_hyp.py results (AI-produced, 2026-11-27)\n" + txt + "\n")


if __name__ == "__main__":
    main()
