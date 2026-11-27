#!/usr/bin/env python3
"""
sieve_checks.py -- numerical sanity checks for GRH-ONLY.md against an exact
sieve to X = 10^7 (q = 4 and q = 3).  Unconditional, double precision;
nothing here is part of a proof.

 (1) the exact resolvent identity of Theorem 2.1(a) for the theta-form
     E^theta_M(x) = 2M x^{-M} R_m(x), R_m(x) = sum_{p<=x} chi*(p) p^m log p;
 (2) the exact unlogged identity (Remark 2.3) with the factor y/(y-u);
 (3) the pointwise O(1/y) form of Theorem 2.1(b) for E_M(x) = 2M x^{-M} log x D_m(x);
 (4) E_M - E^theta_M = O(1/log x)  (eq. (1.1));
 (5) the absolutely convergent formula of Theorem 2.1(c) against the zero list;
 (6) empirical time averages over y = log x in [log 10^3, log 10^7]:
     variance vs sigma_M^2, and monotonicity in M of convex functionals
     (Theorem 3.1, Corollary 3.5) -- finite-range averages, so only indicative.

Runtime ~ 1 minute on one core.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, "..", "..", "prototype", "explore")
X = 10 ** 7


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def chistar(p, q):
    """chi*(p) = -chi(p): +1 on nonresidues, -1 on residues, 0 at p | q"""
    if q == 4:
        r = p % 4
        return np.where(r == 3, 1.0, np.where(r == 1, -1.0, 0.0))
    r = p % 3
    return np.where(r == 2, 1.0, np.where(r == 1, -1.0, 0.0))


def main():
    t0 = time.time()
    out = open(os.path.join(HERE, "sieve_checks_results.txt"), "w")

    def P(s=""):
        print(s)
        out.write(s + "\n")

    p = primes_upto(X).astype(np.float64)
    lp = np.log(p)
    P("primes up to %d: %d" % (X, len(p)))

    for q, zfile in ((4, "chi4_zeros.txt"), (3, "chi3_zeros.txt")):
        cs = chistar(p.astype(np.int64), q)
        zeros = np.array([float(l) for l in open(os.path.join(ZDIR, zfile)) if l.strip()])
        P("")
        P("=" * 76)
        P("q = %d" % q)

        def R(m):                      # R_m at the primes (value just after the jump)
            return np.cumsum(cs * p ** m * lp)

        def D(m):
            return np.cumsum(cs * p ** m)

        # ---------- (1),(2): exact identities ----------
        P("(1)-(2) exact identities (relative residuals; should be ~1e-12):")
        for (m, m2) in ((0.0, 1.0), (-0.25, 0.5), (0.5, 2.0)):
            M, M2, kap = m + 0.5, m2 + 0.5, m2 - m
            for name, S, S2 in (("theta", R(m), R(m2)), ("unlogged", D(m), D(m2))):
                # I(x) = kappa * int_2^x t^{-kappa-1} S2(t) dt, S2 a step function
                pk = p ** (-kap)
                seg = S2[:-1] * (pk[:-1] - pk[1:])           # kappa * int over [p_i, p_{i+1})
                I = np.concatenate([[0.0], np.cumsum(seg)])   # value at x = p_i (left limit = right)
                lhs = S - pk * S2
                idx = [10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, len(p) - 1]
                res = max(abs(lhs[i] - I[i]) / max(abs(I[i]), 1e-300) for i in idx)
                P("    (m,m') = (%5.2f,%4.1f) %-8s  max rel. residual of S_m - x^-k S_m' = k*int : %.2e"
                  % (m, m2, name, res))
        # In normalised form the theta identity reads
        #   E^th_M(x) = (M/M') E^th_M'(x) + (1-M/M') * M x^{-M} int_2^x t^{M-1} E^th_M'(t) dt,
        # which is the same identity multiplied by 2M x^{-M}:  (1 - M/M') * M * 2M'/kappa = 2M.

        # ---------- (3),(4): pointwise O(1/y) statements ----------
        P("(3) Delta(y) = E_M(e^y) - int E_M'(e^{y-u}) mu(du)  (unlogged, no y/(y-u) factor);")
        P("(4) E_M - E^theta_M;   both should be O(1/y):")
        for (m, m2) in ((0.0, 1.0), (-0.25, 0.5)):
            M, M2, kap = m + 0.5, m2 + 0.5, m2 - m
            Dm, Dm2, Rm = D(m), D(m2), R(m)
            # int_2^x t^{-kap-1} log t * D_{m'}(t) dt, antiderivative -t^{-kap}(log t/kap + 1/kap^2)
            A = -p ** (-kap) * (lp / kap + 1 / kap ** 2)
            seg = Dm2[:-1] * (A[1:] - A[:-1])
            I = np.concatenate([[0.0], np.cumsum(seg)])
            EM = 2 * M * p ** (-M) * lp * Dm
            EM2 = 2 * M2 * p ** (-M2) * lp * Dm2
            conv = (M / M2) * EM2 + (1 - M / M2) * M * p ** (-M) * 2 * M2 * I
            Eth = 2 * M * p ** (-M) * Rm
            rows = []
            for xx in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
                i = np.searchsorted(p, xx, side="right") - 1
                y = lp[i]
                rows.append("x=%.0e: y*Delta=%+.3f  y*(E-Eth)=%+.3f" % (xx, y * (EM[i] - conv[i]), y * (EM[i] - Eth[i])))
            j0 = np.searchsorted(p, 1000.0)
            P("    (m,m')=(%5.2f,%4.1f): sup_{10^3<=x<=10^7} y|Delta| = %.3f,  sup y|E-Eth| = %.3f"
              % (m, m2, np.max(lp[j0:] * np.abs(EM[j0:] - conv[j0:])), np.max(lp[j0:] * np.abs(EM[j0:] - Eth[j0:]))))
            for r in rows:
                P("        " + r)

        # ---------- (5): the absolutely convergent formula ----------
        P("(5) G(y) = 2MM' x^{-M} int_2^x t^{-kappa-1} R_m'(t) dt  vs  1 + Re sum 4MM' x^{i gamma}/((M+ig)(M'+ig)):")
        for (m, m2) in ((0.0, 1.0),):
            M, M2, kap = m + 0.5, m2 + 0.5, m2 - m
            Rm2 = R(m2)
            pk = p ** (-kap)
            seg = Rm2[:-1] * (pk[:-1] - pk[1:]) / kap
            I = np.concatenate([[0.0], np.cumsum(seg)])
            G = 2 * M * M2 * p ** (-M) * I
            tailb = np.sum(4 * M * M2 / zeros ** 2) * 0 + 4 * M * M2 * (  # crude tail: sum_{g>T} 1/g^2
                (0.1560296498896449 if q == 4 else 0.1134021568527586) / 2 - np.sum(1 / zeros ** 2))
            for xx in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
                i = np.searchsorted(p, xx, side="right") - 1
                y = lp[i]
                zsum = 1 + np.sum((4 * M * M2 * np.exp(1j * zeros * y) / ((M + 1j * zeros) * (M2 + 1j * zeros))).real)
                P("    x=%.0e: G = %.5f, zero formula (%d zeros) = %.5f, diff = %+.5f  (omitted-zero bound %.4f; rest is the o(1))"
                  % (xx, G[i], len(zeros), zsum, G[i] - zsum, tailb))

        # ---------- (6): empirical averages ----------
        P("(6) time averages over y in [log 1e3, log 1e7] (uniform in y; theta-form g_M and unlogged f_M):")
        Y0, Y1 = np.log(1e3), np.log(1e7)
        ys = np.linspace(Y0, Y1, 400001)
        xs = np.exp(ys)
        idx = np.searchsorted(p, xs, side="right") - 1
        Ms = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.5]
        P("    %5s | %9s %9s | %9s %9s %9s %9s | %8s" %
          ("M", "mean", "var", "E|g-1|", "E g^-", "E(g-.5)^-", "E(g-1)^4", "frac g>0"))
        prev = None
        mono = True
        for M in Ms:
            m = M - 0.5
            g = 2 * M * xs ** (-M) * R(m)[idx]
            f = 2 * M * xs ** (-M) * ys * D(m)[idx]
            st = (g.mean(), g.var(), np.abs(g - 1).mean(), np.maximum(-g, 0).mean(),
                  np.maximum(0.5 - g, 0).mean(), ((g - 1) ** 4).mean(), (g > 0).mean())
            stf = (f.mean(), f.var(), np.abs(f - 1).mean(), np.maximum(-f, 0).mean(),
                   np.maximum(0.5 - f, 0).mean(), ((f - 1) ** 4).mean(), (f > 0).mean())
            P("    %5.2f | %9.5f %9.5f | %9.5f %9.6f %9.6f %9.5f | %8.5f   (theta)" % ((M,) + st))
            P("    %5s | %9.5f %9.5f | %9.5f %9.6f %9.6f %9.5f | %8.5f   (unlogged)" % (("",) + stf))
            cur = (st[1], st[2], st[3], st[4], st[5])
            if prev is not None and any(c < pv - 1e-12 for c, pv in zip(cur, prev)):
                mono = False
            prev = cur
        P("    convex functionals (var, E|g-1|, E g^-, E(g-.5)^-, E(g-1)^4) nondecreasing in M on this window: %s" % mono)
        P("    [closed-form variances sigma_M^2 (if all zeros simple), q=4: M=.25: 0.03898, .5: 0.15557, 1: 0.61686, 1.5: 1.36856;")
        P("     q=3: M=.25: 0.02834, .5: 0.11323, 1: 0.45088, 1.5: 1.00701.  Finite window => agreement only to ~10-20%.]")

    P("\nelapsed %.1f s" % (time.time() - t0))
    out.close()


if __name__ == "__main__":
    main()
