#!/usr/bin/env python3
"""
model_certified.py -- certified construction (Arb) for DENSITY.md Proposition 7.3:
at the classical weight M = 1/2 there are c_n >= 0 with sum c_n < 1 (sigma^2 < 2) such that
    F(t) = 1 + sum_{n>=1} c_n (cos nt + (n/M) sin nt)
vanishes identically on the arc J = [2, 2.3] of the circle.

Method (all inequalities are checked in ball arithmetic, 200 bits):
 1. a float nonnegative least-squares solution c^(N), N = 100, of F = 0 on J
    (scipy; stored in model_M0.5_coeffs.txt; the doubles are used as EXACT numbers);
 2. certified bounds E0 >= sup_J |F_N|, E1 >= sup_J |F_N'| by Taylor expansion of order 3
    around K midpoints, remainder bounded with L4 = sum c_n n^4 (1 + n/M) >= sup |F_N''''|;
 3. the residual e = -F_N is removed exactly on J by Psi_c(t) = -M e^{Mt} int_a^t e^{-Ms} e(s) ds
    (so Psi_c - Psi_c'/M = e on J), extended to a C^{1,1} even function supported in
    [a-delta, b+delta] and its mirror image, with
       ||Psi_c_ext''||_{L1(0,pi)} <= W2 := (b-a) Q2 + 2 (9 Q1 + 6 Q0/delta),
       Q0 = M e^{M(b-a)} (b-a) E0,  Q1 = M (Q0 + E0),  Q2 = M (Q1 + E1);
    its cosine coefficients are bounded by (2/pi) W2 / n^2;
 4. they are dominated by eps * (cosine coefficients of P(t) = (1-|t|/r)_+^2), r < a,
    which are (4/(pi r n^2)) (1 - sin(nr)/(nr)) >= 4 e_r/(pi r n^2), e_r = 1 - sin(r)/r
    (valid since sin(x)/x <= sin(r)/r for all x >= r when sin(r)/r >= 1/pi), with
    eps := r W2 / (2 e_r);  P vanishes on J;
 5. Psi := (Psi_N + Psi_c_ext + eps P)/Z has nonnegative cosine coefficients, mean 1,
    satisfies Psi' = M Psi on J, and  sum c_n = Psi(0) - 1 <= (1 + s_N + eps)/Zmin - 1,
    Zmin = 1 - (1/pi)(b - a + 2 delta)(Q0 + Q1 delta).
Output: model_certified_results.txt.  Runtime: a few seconds.
"""
import os, sys
from fractions import Fraction
import numpy as np
from flint import arb, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
M = Fraction(1, 2); a = Fraction(2); b = Fraction(23, 10); r = Fraction(19, 10); delta = Fraction(1, 20)
N = 100; CF = os.path.join(HERE, "model_M0.5_coeffs.txt")

def fit():
    from scipy.optimize import nnls
    t = np.linspace(float(a), float(b), 300); n = np.arange(1, N + 1)
    A = (np.cos(np.outer(t, n)) + (n / float(M)) * np.sin(np.outer(t, n))) / n
    d, _ = nnls(A, -np.ones_like(t), maxiter=50000)
    c = d / n
    with open(CF, "w") as f:
        for k, v in zip(n, c):
            if v > 0: f.write("%d %s\n" % (k, float(v).hex()))

def q(fr): return arb(fr.numerator) / arb(fr.denominator)

def main():
    if not os.path.exists(CF): fit()
    cs = [(int(l.split()[0]), float.fromhex(l.split()[1])) for l in open(CF) if l.strip()]
    assert all(v > 0 for _, v in cs)
    out = open(os.path.join(HERE, "model_certified_results.txt"), "w")
    def log(s=""):
        print(s); out.write(s + "\n")
    Ma = q(M); ca = [(k, arb(v)) for k, v in cs]          # arb(double) is exact
    sN = sum((v for _, v in ca), arb(0))
    L4 = sum((v * k ** 4 * (1 + k / Ma) for k, v in ca), arb(0))
    K = 20000; h = (b - a) / K; rho = q(h) / 2
    E0 = arb(0); E1 = arb(0)
    for i in range(K):
        t = q(a + (2 * i + 1) * h / 2)
        F0 = arb(1); F1 = arb(0); F2 = arb(0); F3 = arb(0)
        for k, v in ca:
            S, C = (k * t).sin_cos()
            kM = k / Ma
            F0 += v * (C + kM * S)
            F1 += v * k * (-S + kM * C)
            F2 += v * k * k * (-C - kM * S)
            F3 += v * k ** 3 * (S - kM * C)
        b0 = abs(F0) + abs(F1) * rho + abs(F2) * rho ** 2 / 2 + abs(F3) * rho ** 3 / 6 + L4 * rho ** 4 / 24
        b1 = abs(F1) + abs(F2) * rho + abs(F3) * rho ** 2 / 2 + L4 * rho ** 3 / 6
        E0 = E0.max(b0)
        E1 = E1.max(b1)
    ba = q(b - a); da = q(delta); ra = q(r)
    Q0 = Ma * (Ma * ba).exp() * ba * E0
    Q1 = Ma * (Q0 + E0)
    Q2 = Ma * (Q1 + E1)
    W2 = ba * Q2 + 2 * (9 * Q1 + 6 * Q0 / da)
    sr = ra.sin() / ra
    assert sr > 1 / arb.pi(), "need sin(r)/r >= 1/pi"
    assert ra < q(a) - da and q(b) + da < arb.pi() and ra > 0
    er = 1 - sr
    eps = ra * W2 / (2 * er)
    Zmin = 1 - (ba + 2 * da) * (Q0 + Q1 * da) / arb.pi()
    assert Zmin > 0
    sfin = (1 + sN + eps) / Zmin - 1
    log("M = %s, J = [%s, %s], r = %s, delta = %s, N = %d, %d nonzero coefficients" % (M, a, b, r, delta, N, len(cs)))
    log("sum c_n^(N)            = %s" % sN.str(15))
    log("L4                     = %s" % L4.str(8))
    log("E0 >= sup_J |F_N|      : %s" % E0.str(6))
    log("E1 >= sup_J |F_N'|     : %s" % E1.str(6))
    log("W2                     = %s" % W2.str(6))
    log("eps                    = %s" % eps.str(6))
    log("Zmin                   = %s" % Zmin.str(15))
    log("certified: sum c_n <= %s ;  sigma^2 = 2 sum c_n <= %s" % (sfin.str(12), (2 * sfin).str(12)))
    ok = sfin < q(Fraction(2717, 10000))
    assert ok
    log("certified comparison  sum c_n < 0.2717  (hence sigma^2 < 0.5434 < 2): %s" % bool(ok))
    wmax = max(float(v) * (float(M) ** 2 + k * k) / (4 * float(M) ** 2) for k, v in cs)
    log("[num] largest weight w_n = c_n (M^2+n^2)/(4M^2) among the N fitted ones: %.4f" % wmax)
    log("atom of the law of F at 0: mass >= |J|/(2 pi) = %.5f" % (0.3 / (2 * np.pi)))
    out.close()

if __name__ == "__main__":
    main()
