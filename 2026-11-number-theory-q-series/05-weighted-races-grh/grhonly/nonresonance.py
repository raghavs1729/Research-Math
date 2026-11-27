#!/usr/bin/env python3
"""
nonresonance.py -- beyond Cantelli without LI: the "finite non-resonance"
bound of GRH-ONLY.md, Theorem 6.1 (degree-one multiplier), q = 4 and q = 3.

Idea.  LI is an infinite statement, but LI *up to a fixed order among the
certified zeros* is a finite one and can be decided by comparing certified
enclosures.  With F = 1 + H + W_K + W_t  (H: the first N1 zeros, W_K: the
other certified zeros up to T' = T - 2 gamma_{N1}, W_t: everything else) and
the test function  Psi = (1 - kappa F)^2 >= 1_{F <= 0},
kappa = kappa_0 (1 - sum_n b_n cos(theta_n)) >= 0  a multiplier depending on
the first N1 phases only, one gets

   nu_M((-inf,0])  <=  J(kappa_0, b)
     := (1-kappa_0)^2 + 2 (1-kappa_0) kappa_0 rho
        + kappa_0^2 [ sum (a_n-b_n)^2/2 + (1/4)(sum b_n^2)(sum a_n^2)
                      + (1/2)(sum a_n b_n)^2 - (3/8) sum a_n^2 b_n^2
                      + (1 + sum b_n^2/2) V_K + (1 + sum |b_n|)^2 V_t ],
   rho = sum a_n b_n / 2,   a_n = 4M/sqrt(M^2+gamma_n^2),

PROVIDED the non-resonance conditions (NR1)-(NR3) hold:

 (NR1) no integer vector n != 0 supported on the head with |n|_1 <= 4 has
       n . gamma = 0;
 (NR2) for |n|_1 <= 3,  |n . gamma|  is not an ordinate outside the head;
 (NR3) for 1 <= |n|_1 <= 2,  n . gamma  is not of the form  +-g_i +- g_j
       with g_i, g_j certified ordinates outside the head.

This script (i) verifies (NR1)-(NR3) in double precision with an explicit
margin (the certified enclosures have radius 1e-20 and the stored ordinates
25 significant digits; a margin > 1e-9 in IEEE doubles therefore proves
non-coincidence: the accumulated rounding error of the at most 4-term sums
of numbers < 1300 is < 1e-12), and (ii) minimises J numerically.  Any
(kappa_0, b) with sum |b_n| <= 1 gives a valid bound, so the minimisation
need not be certified; the final J is re-evaluated in mpmath with the
certified variance bounds taken from certified_bounds_results.txt (rounded UP).

Status of the output: SEARCH ONLY.  This script finds good (kappa_0, b) and
prints double-precision margins for (NR1)-(NR3); it carries no rigor.  The
certificate for Theorem 6.2 (exact-integer verification of (NR1)-(NR3) and an
Arb enclosure of J at an exactly rational feasible point) is produced by
nonresonance_certified.py, which imports only the float search from here.
"""
import itertools
import os
import sys

import numpy as np
import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")


def load(fn):
    with open(os.path.join(ZDIR, fn)) as f:
        return [mp.mpf(l.strip()) for l in f if l.strip()]


def int_vectors(N, order):
    """all nonzero integer vectors of length N with |n|_1 <= order, up to sign"""
    seen = set()
    for k in range(1, order + 1):
        for idx in itertools.combinations_with_replacement(range(N), k):
            for signs in itertools.product((1, -1), repeat=k):
                v = [0] * N
                for i, s in zip(idx, signs):
                    v[i] += s
                v = tuple(v)
                if any(v):
                    # normalise sign
                    first = next(x for x in v if x)
                    if first < 0:
                        v = tuple(-x for x in v)
                    if sum(abs(x) for x in v) <= order:
                        seen.add(v)
    return sorted(seen)


def check_nonresonance(g, N1, T, log):
    """returns the minimal margins of (NR1), (NR2), (NR3)."""
    gf = np.array([float(x) for x in g])
    head = gf[:N1]
    rest = gf[N1:]
    # (NR1)
    m1 = np.inf
    arg1 = None
    for v in int_vectors(N1, 4):
        s = abs(float(np.dot(v, head)))
        if s < m1:
            m1, arg1 = s, v
    # (NR2): |n.gamma| for |n|_1 <= 3 against ordinates outside the head.
    #        All such values are <= 3*gamma_{N1} < T, so only listed zeros matter.
    assert 3 * head[-1] < T
    m2 = np.inf
    arg2 = None
    for v in int_vectors(N1, 3):
        s = abs(float(np.dot(v, head)))
        j = np.searchsorted(rest, s)
        for jj in (j - 1, j):
            if 0 <= jj < len(rest):
                d = abs(rest[jj] - s)
                if d < m2:
                    m2, arg2 = d, (v, float(rest[jj]))
    # (NR3): h = n.gamma, 1 <= |n|_1 <= 2, against  +-g_i +- g_j  (i, j outside head)
    R = rest
    sums = (R[:, None] + R[None, :])[np.triu_indices(len(R))]          # g_i + g_j, i <= j
    diffs = (R[None, :] - R[:, None])[np.triu_indices(len(R), k=1)]    # g_j - g_i > 0, i < j
    D = np.sort(np.concatenate([sums, diffs]))
    m3 = np.inf
    arg3 = None
    for v in int_vectors(N1, 2):
        s = abs(float(np.dot(v, head)))
        j = np.searchsorted(D, s)
        for jj in (j - 1, j):
            if 0 <= jj < len(D):
                d = abs(D[jj] - s)
                if d < m3:
                    m3, arg3 = d, (v, float(D[jj]))
    log("    (NR1) min |n.gamma|, 1<=|n|_1<=4            : %.3e   at n = %s" % (m1, arg1))
    log("    (NR2) min dist(|n.gamma|, other ordinates)  : %.3e   at %s" % (m2, arg2))
    log("    (NR3) min dist(n.gamma, +-g_i+-g_j)         : %.3e   at %s  (|D| = %d)"
        % (m3, arg3, len(D)))
    return m1, m2, m3


def J_value(k0, b, a, VK, Vt):
    """the bound J(kappa_0, b) (works for floats, numpy arrays, or mpf lists)"""
    Sab = sum(x * y for x, y in zip(a, b))
    Saa = sum(x * x for x in a)
    Sbb = sum(y * y for y in b)
    Sa2b2 = sum(x * x * y * y for x, y in zip(a, b))
    Sd = sum((x - y) ** 2 for x, y in zip(a, b))
    Sabs = sum(abs(y) for y in b)
    rho = Sab / 2
    br = Sd / 2 + Saa * Sbb / 4 + Sab ** 2 / 2 - 3 * Sa2b2 / 8 \
        + (1 + Sbb / 2) * VK + (1 + Sabs) ** 2 * Vt
    return (1 - k0) ** 2 + 2 * (1 - k0) * k0 * rho + k0 ** 2 * br


def minimise(a, VK, Vt, iters=400):
    """crude coordinate search; any feasible point is a valid bound"""
    a = [float(x) for x in a]
    N = len(a)
    V = sum(x * x for x in a) / 2 + VK + Vt
    k0 = 1 / (1 + V)
    b = [0.0] * N
    best = J_value(k0, b, a, VK, Vt)
    step = 0.1
    for it in range(iters):
        improved = False
        for i in range(-1, N):
            for sgn in (1, -1):
                if i == -1:
                    k1, b1 = k0 + sgn * step, b
                    if k1 <= 0:
                        continue
                else:
                    b1 = list(b)
                    b1[i] += sgn * step
                    k1 = k0
                    if b1[i] < 0 or sum(abs(y) for y in b1) > 1 - 1e-9:
                        continue
                val = J_value(k1, b1, a, VK, Vt)
                if val < best - 1e-15:
                    best, k0, b, improved = val, k1, b1, True
        if not improved:
            step /= 2
            if step < 1e-7:
                break
    return best, k0, b


def run(name, q, zfile, M, Vplus_A, Vplus_B, out):
    def log(s):
        print(s)
        out.write(s + "\n")

    g = load(zfile)
    T = float(g[-1]) + 0.25
    M = mp.mpf(M)
    amp = [4 * M / mp.sqrt(M * M + x * x) for x in g]
    var_each = [x * x / 2 for x in amp]
    log("=" * 78)
    log("%s, M = %s (m = %s): %d certified simple zeros, T = %.6f" % (name, M, M - mp.mpf(1) / 2, len(g), T))
    log("  Cantelli with V^+ (tier A) = %.9f : nu((-inf,0]) <= %.6f, density >= %.6f"
        % (Vplus_A, Vplus_A / (1 + Vplus_A), 1 / (1 + Vplus_A)))
    results = {}
    for N1 in (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20):
        gN = float(g[N1 - 1])
        Tprime = T - 2 * gN
        if 3 * gN >= T or Tprime <= gN:
            continue
        idxK = [i for i in range(N1, len(g)) if float(g[i]) <= Tprime]
        Vh = sum(var_each[:N1])
        VK = sum(var_each[i] for i in idxK)
        log("  head N1 = %d (gamma_N1 = %.4f), T' = %.3f, #K = %d, V_head = %.6f, V_K = %.6f"
            % (N1, gN, Tprime, len(idxK), float(Vh), float(VK)))
        m1, m2, m3 = check_nonresonance(g, N1, T, log)
        ok = min(m1, m2, m3) > 1e-9
        log("    non-resonance verified with margin > 1e-9: %s" % ok)
        if not ok:
            continue
        for tier, Vp in (("A", Vplus_A), ("B", Vplus_B)):
            # V_t = everything not in head or K; round UP by 1e-12
            Vt = mp.mpf(Vp) - Vh - VK + mp.mpf("1e-12")
            best, k0, b = minimise(amp[:N1], float(VK), float(Vt))
            # re-evaluate in mpmath at the found point
            bmp = [mp.mpf(repr(y)) for y in b]
            assert sum(abs(y) for y in bmp) <= 1
            Jmp = J_value(mp.mpf(repr(k0)), bmp, amp[:N1], VK, Vt)
            log("    tier %s: V_t = %.6f  ->  J = %s   density >= %s   (kappa_0 = %.5f, b = %s)"
                % (tier, float(Vt), mp.nstr(Jmp, 9), mp.nstr(1 - Jmp, 9), k0,
                   "[" + ", ".join("%.4f" % y for y in b) + "]"))
            results[(N1, tier)] = (Jmp, k0, b)
    for tier in ("A", "B"):
        cand = [(v[0], k) for k, v in results.items() if k[1] == tier]
        if cand:
            Jb, kk = min(cand)
            log("  BEST tier %s: N1 = %d, nu((-inf,0]) <= %s, lower log-density{D_m>0} >= %s"
                % (tier, kk[0], mp.nstr(Jb, 8), mp.nstr(1 - Jb, 8)))
    return results


def main():
    outfn = os.path.join(HERE, "nonresonance_results.txt")
    with open(outfn, "w") as out:
        # certified variance bounds, copied from certified_bounds_results.txt and rounded UP
        run("chi_4", 4, "chi4_zeros.txt", "0.5", 0.169630835, 0.155573631, out)
        run("chi_3", 3, "chi3_zeros.txt", "0.5", 0.125451412, 0.113234208, out)
        run("chi_4", 4, "chi4_zeros.txt", "0.1", 0.00680295734, 0.00624066911, out)
        run("chi_4", 4, "chi4_zeros.txt", "0.25", 0.0424941433, 0.0389798419, out)
        run("chi_4", 4, "chi4_zeros.txt", "0.75", 0.380388092, 0.348759386, out)
        run("chi_4", 4, "chi4_zeros.txt", "1.0", 0.673112007, 0.616883206, out)
        run("chi_3", 3, "chi3_zeros.txt", "0.25", 0.0313951130, 0.0283408118, out)
        run("chi_3", 3, "chi3_zeros.txt", "1.0", 0.499764841, 0.450896037, out)

if __name__ == "__main__":
    main()
