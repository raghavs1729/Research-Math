#!/usr/bin/env python3
"""certify_analytic_leg.py -- ball-arithmetic certificate for the analytic leg

    Proposition.  Assume GRH(chi_q), q in {3,4}.  Then delta_q'(M) < 0 for all
    M >= 300.

(paper/weighted-races.tex, Lemma lem:analytic-leg and the chain it rests on).
The mathematical argument is written out in monotonicity/ANALYTIC-LEG-AUDIT.md;
this script certifies, in Arb (python-flint) ball arithmetic, every NUMERICAL
fact that argument uses, and runs a number of independent consistency checks.

Sections (each line of output ends in PASS or FAIL; exit status 1 on any FAIL):

 [A] constants of the chain: Bessel constants (j_{0,1} > 2.4, 0.446, 0.733,
     0.447, c_3, (J2) spot sweep), 0.9973, -zeta'/zeta(2) < 0.57,
     sum Lambda(n) log n / n^2 < 0.89, the closed-form integrals behind
     64 D_+ + 512 l_0 + 1165, I_1', I_2'.
 [A'] the constants of the proof of thm:monotone-large (0.189, 0.18, 1.133,
     0.9510, 0.90, 0.517); not needed for the Proposition.
 [B] the inequality TAIL_an(300) < MAIN_an(300) for q = 3, 4.
 [C] the three logarithmic-derivative constants kappa_1, kappa_2, kappa_3 at
     M = 300 (see the audit, Step 6), i.e. the finitely many numerical facts
     the monotonicity argument needs:  kappa_1 > kappa_3 and kappa_2 > kappa_3.
     Also the paper's intermediate constants (-0.15, -0.09, 0.001, 0.007,
     l_{Y_1} >= 4.4, l_{Y_2} >= 5.1).
 [D] INDEPENDENT of [C]: a cell sweep of [300, 1e7] (ratio 1.01).  On every
     cell (i) TAIL_an/MAIN_an < 1 by monotone enclosures, and (ii) the EXACT
     logarithmic derivatives of T_1/MAIN_an and T_2/MAIN_an, evaluated on the
     whole cell as a ball, are <= -0.14 resp. <= -0.08.
 [E] closed forms: the explicit bounds of lem:Dprop(iii), lem:Nlow used for
     D, sigma^2, S_4 are compared with Arb evaluations of the exact closed
     forms at sample weights (sanity check of the Binet bookkeeping).
 [F] zero lists (sanity only; NOT used by the proof): n_1(M) <= N(M),
     n_2(M) <= N(2M), N(T) <= N^+(T) against the certified ordinate lists.
 [G] cross-check against the direct certificate of delta' (cert/results_q*_
     stretch.txt): the analytic bound  delta'(M) <= -(MAIN_an - TAIL_an)(M)
     must not contradict the certified enclosures lo <= delta'(M) <= hi.

Run from the repository root:   python3 monotonicity/analytic-leg/certify_analytic_leg.py
"""
from __future__ import annotations

import math
import os
import sys

from flint import arb, arb_series, ctx

ctx.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PI = arb.pi()

FAILS: list[str] = []


def A(x) -> arb:
    """Exact decimal ball."""
    return arb(str(x))


def check(name: str, cond: bool, detail: str = "") -> None:
    ok = bool(cond)
    print(f"    {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def s(x: arb, n: int = 8) -> str:
    return x.str(n, radius=False)


# --------------------------------------------------------------------------
# the explicit functions of lem:analytic-leg
# --------------------------------------------------------------------------
R1 = A("0.447") / A("0.733")
R2 = A("0.733")
C3 = 4 * arb(5).sqrt().sqrt() / PI
I1P = arb(16) ** A("-2.5") / A("2.5")
I2P = arb(16) ** A("-2.5") * (A("2.5") * arb(16).log() + 1) / A("6.25")
SQ3 = arb(3).sqrt()


def ell(q: int, M: arb) -> arb:
    return (arb(q) * M / (2 * PI)).log()


def n_i(q: int, M: arb, i: int) -> arb:
    Y = i * M / SQ3
    return Y * (ell(q, Y) - 2) / 16 - A("0.5")


def dn_i(q: int, M: arb, i: int) -> arb:
    """d n_i / dM = i (l_{Y_i} - 1) / (16 sqrt 3)."""
    Y = i * M / SQ3
    return i * (ell(q, Y) - 1) / (16 * SQ3)


def Dplus(q, M):
    return 2 * ell(q, M) + 2 + 8 / M


def Dminus(q, M):
    return 2 * ell(q, M) + 2 - 8 / M


def ell0(q, M):
    return (arb(q) * (2 * M + 2) / (2 * PI)).log() + A("1.14")


def beta_plus(q, M):
    l = ell(q, M)
    return M * l + A("1.55") + (32 * M * (l - 1) + 100) / 1088


def G(beta: arb) -> arb:
    """int_0^{1/8} t^2 exp(-beta t^2) dt."""
    rb = beta.sqrt()
    return PI.sqrt() / (4 * beta * rb) * (rb / 8).erf() - (-beta / 64).exp() / (16 * beta)


def H(beta: arb) -> arb:
    """int_0^{1/8} t^4 exp(-beta t^2) dt = (3/(2 beta)) G - 8^{-3} e^{-beta/64}/(2 beta)."""
    return 3 * G(beta) / (2 * beta) - (-beta / 64).exp() / (1024 * beta)


K3 = A("0.9973") / PI / 2


def MAIN_an(q, M, Mb=None):
    """MAIN_an(M); with Mb given: lower bound on the cell [M, Mb]
    (D_- increasing, beta_+ increasing, G decreasing)."""
    return K3 * Dminus(q, M) * G(beta_plus(q, M if Mb is None else Mb))


def P1(q, M):
    return 64 * Dplus(q, M) + 512 * ell0(q, M) + 1165


def P2(q, M):
    return Dplus(q, M) * I1P / 2 + 4 * (ell0(q, M) * I1P + I2P)


def T1(q, M, Mb=None):
    """first summand of TAIL_an; with Mb: upper bound on the cell [M, Mb]
    (P_1 increasing, n_i increasing)."""
    return P1(q, M if Mb is None else Mb) / (A("0.447") * PI) \
        * R1 ** n_i(q, M, 1) * R2 ** n_i(q, M, 2)


def T2(q, M, Mb=None):
    return C3 ** 9 / PI * P2(q, M if Mb is None else Mb) * R2 ** (n_i(q, M, 2) - 10)


def dlog_T1(q, M):
    dP = 64 * (2 / M - 8 / (M * M)) + 512 / (M + 1)
    return dP / P1(q, M) + dn_i(q, M, 1) * R1.log() + dn_i(q, M, 2) * R2.log()


def dlog_T2(q, M):
    dP = (1 / M - 4 / (M * M)) * I1P + 4 * I1P / (M + 1)
    return dP / P2(q, M) + dn_i(q, M, 2) * R2.log()


def dlog_MAIN(q, M):
    l = ell(q, M)
    b = beta_plus(q, M)
    db = l * 35 / 34 + 1
    dD = 2 / M + 8 / (M * M)
    return dD / Dminus(q, M) - db * H(b) / G(b)


# --------------------------------------------------------------------------
def j0_abs_upper(a: float, b: float) -> arb:
    """Rigorous upper bound for sup_{[a,b]} |J0|: |J0(mid)| + (b-a)/2, by the
    mean value theorem with |J0'| = |J1| <= 1 (J1(z) = (1/pi) int_0^pi
    cos(theta - z sin theta) d theta).  (Arb's bessel_j on a wide input ball
    is valid but very pessimistic, so wide balls are avoided.)"""
    mid = (arb(a) + arb(b)) / 2
    return abs(mid.bessel_j(0)) + (arb(b) - arb(a)) / 2


def sup_abs_j0(lo: float, hi: float, target: arb) -> bool:
    """Prove sup_{[lo,hi]} |J0| < target by adaptive bisection."""
    stack = [(lo, hi)]
    while stack:
        a, b = stack.pop()
        if j0_abs_upper(a, b) < target:
            continue
        if b - a < 1e-9:
            return False
        m = (a + b) / 2
        stack.append((a, m))
        stack.append((m, b))
    return True


def section_A() -> None:
    print("[A] constants of the chain")
    # lem:phi(iii)
    j = A("2.4").bessel_j(0)
    check("J0(2.4) > 0.0024  (j_{0,1} > 2.4 with J0 decreasing on [0, 2 sqrt 2])",
          j > A("0.0024"), f"J0(2.4) = {s(j)}")
    # J1 > 0 on (0, 2 sqrt 2]: brute-force confirmation on [1e-6, 2.8285]
    ok = True
    a = 1e-6
    while a < 2.8285:
        b = min(a * 1.02 + 1e-7, 2.8285)
        # J1 >= J1(mid) - (b-a)/2 on [a,b], since |J1'| = |J0 - J1/z| <= 1
        # (J1' = (J0 - J2)/2 and |J0|, |J2| <= 1)
        if not ((arb(a) + arb(b)) / 2).bessel_j(1) - (arb(b) - arb(a)) / 2 > 0:
            ok = False
        a = b
    check("J1 > 0 on [1e-6, 2.8285]  (J0 strictly decreasing on (0, 2 sqrt 2])", ok)
    # lem:J0env(ii)
    z8 = arb(8).sqrt()
    P8 = z8.bessel_j(0) ** 2 + z8.bessel_j(1) ** 2
    check("sqrt(P(2 sqrt 2)) < 0.446", P8.sqrt() < A("0.446"), f"= {s(P8.sqrt())}")
    check("J0(2 sqrt 2) in [-0.19655,-0.19654], J1(2 sqrt 2) in [0.40019,0.40020]",
          A("-0.19655") < z8.bessel_j(0) < A("-0.19654")
          and A("0.40019") < z8.bessel_j(1) < A("0.40020"),
          f"J0 = {s(z8.bessel_j(0))}, J1 = {s(z8.bessel_j(1))}")
    # lem:J0env(iii)
    thr1 = 12 / (5 * arb(5).sqrt())
    check("12/(5 sqrt 5) > 1.0733", thr1 > A("1.0733"), f"= {s(thr1)}")
    j = A("1.0733").bessel_j(0)
    check("0.446 < J0(1.0733) < 0.733", A("0.446") < j < A("0.733"), f"= {s(j)}")
    thr2 = 6 * arb(2).sqrt() / 5
    check("2 sqrt 2 * 3/5 > 1.697", thr2 > A("1.697"), f"= {s(thr2)}")
    j = A("1.697").bessel_j(0)
    check("J0(1.697) < 0.410  (so sup_{z>=1.697}|J0| <= 0.446 < 0.447)",
          abs(j) < A("0.410"), f"= {s(j)}")
    # independent brute force of the two suprema: [z0, 64] by bisection, and
    # (J2): |J0(z)| <= (8/pi) z^{-1/2} <= 0.3184 for z >= 64.
    check("(8/pi)/sqrt(64) < 0.447", 8 / PI / 8 < A("0.447"), f"= {s(8 / PI / 8)}")
    check("brute force: sup_{[1.0733, 64]} |J0| < 0.733",
          sup_abs_j0(1.0733, 64.0, A("0.733")))
    check("brute force: sup_{[1.697, 64]} |J0| < 0.447",
          sup_abs_j0(1.697, 64.0, A("0.447")))
    # (J2) sweep on [1, 400] (the bound is proved in lem:J0facts; confirmation)
    ok = True
    a = 1.0
    while a < 400.0:
        b = a + 1 / 32
        if not j0_abs_upper(a, b) * arb(b).sqrt() < 8 / PI:
            ok = False
        a = b
    check("(J2) confirmation: |J0(z)| sqrt z < 8/pi on [1, 400]", ok)
    check("c_3 = 4*5^(1/4)/pi < 1.9041 and c_3/sqrt(16) < 0.48",
          C3 < A("1.9041") and C3 / 4 < A("0.48"), f"c_3 = {s(C3)}")
    check("1 - 1/384 > 0.9973", 1 - arb(1) / 384 > A("0.9973"))
    # (J3) and lem:phi(ii) via positivity of the Rayleigh coefficients (lem:phi(i)):
    # (-log J0(z) - z^2/4)/z^4 and phi(z)/z are power series in z^2 with positive
    # coefficients, hence increasing; their suprema on (0, z1] are the values at z1.
    for z1 in ("0.5", "1"):
        z = A(z1)
        v = (-z.bessel_j(0).log() - z * z / 4) / z ** 4
        check(f"(J3): sup_(0,{z1}] (-log J0(z) - z^2/4)/z^4 = value at {z1} <= 1/17",
              v <= arb(1) / 17, f"= {s(v)}  (1/17 = 0.0588)")
    z = A("0.5")
    v = z.bessel_j(1) / z.bessel_j(0) / z
    check("lem:phi(ii): sup_(0,1/2] phi(z)/z = phi(1/2)/(1/2) <= 0.517 (and >= 1/2)",
          A("0.5") <= v <= A("0.517"), f"= {s(v)}")
    # zeta constants
    zs = arb_series([2, 1], prec=4).zeta()
    z0, z1, z2 = zs[0], zs[1], 2 * zs[2]
    lam1 = -z1 / z0                      # sum Lambda(n) n^-2
    lam2 = z2 / z0 - (z1 / z0) ** 2      # sum Lambda(n) log n n^-2
    check("sum Lambda(n)/n^2 = -zeta'/zeta(2) < 0.57", lam1 < A("0.57"), f"= {s(lam1, 10)}")
    check("sum Lambda(n) log n/n^2 = (zeta'/zeta)'(2) < 0.89", lam2 < A("0.89"), f"= {s(lam2, 10)}")
    # the integrals behind 64 D_+ + 512 l_0 + 1165
    i_t = (arb(256) - A("0.36")) / 2                      # int_{3/5}^{16} t dt
    i_tl = 128 * arb(16).log() - arb(255) / 4             # int_1^16 t log t dt
    check("int_{3/5}^{16} t dt <= 127.9  (so /2 <= 64 and *4 <= 512)",
          i_t <= A("127.9") and i_t / 2 < 64 and 4 * i_t < 512, f"= {s(i_t)}")
    check("4 * int_1^16 t log t dt <= 1165", 4 * i_tl < 1165, f"int = {s(i_tl)}")
    check("I_1', I_2' closed forms positive",
          I1P > 0 and I2P > 0, f"I_1' = {s(I1P)}, I_2' = {s(I2P)}")
    # nonvanishing at the centre (unconditional): L(s,chi) = q^-s sum_a chi(a) zeta(s, a/q)
    from flint import acb
    h = acb("0.5")
    L4 = (h.zeta(acb(1) / 4) - h.zeta(acb(3) / 4)).real / 2
    L3 = (h.zeta(acb(1) / 3) - h.zeta(acb(2) / 3)).real / arb(3).sqrt()
    check("L(1/2, chi_4) > 0 and L(1/2, chi_3) > 0  (no hypothesis needed at s = 1/2)",
          L4 > 0 and L3 > 0, f"L(1/2,chi_4) = {s(L4, 10)}, L(1/2,chi_3) = {s(L3, 10)}")
    # trivial algebra used in the audit
    check("64*17 = 1088; 35/34 < 1.03", 64 * 17 == 1088 and arb(35) / 34 < A("1.03"))


def section_A2() -> None:
    print("[A'] constants of thm:monotone-large / lem:phi(ii) (not needed for the Proposition)")
    s2p = (2 * PI).sqrt()
    c = A("0.9973") * (arb(33) / 32) ** A("-1.5") / (2 * s2p)
    check("0.9973 (33/32)^(-3/2) / (2 sqrt(2 pi)) > 0.189 and (33/32)^(-3/2) > 0.9548",
          c > A("0.189") and (arb(33) / 32) ** A("-1.5") > A("0.9548"), f"= {s(c)}")
    check("0.189*0.98 - 0.005 >= 0.18", A("0.189") * A("0.98") - A("0.005") >= A("0.18"))
    u = A("0.517") / (4 * PI.sqrt()) * (arb(17) / 8) ** A("1.5")
    check("0.517/(4 sqrt pi) (17/8)^(3/2) <= 0.226 and 0.226 * 2 sqrt(2 pi) <= 1.133",
          u <= A("0.226") and A("0.226") * 2 * s2p <= A("1.133"),
          f"= {s(u)}, {s(A('0.226') * 2 * s2p)}")
    # The paper's sandwich step says "(4/sqrt5)(1/8) = 0.447, factor <= J0(0.447) <= 0.9510".
    # That arithmetic is wrong: (4/sqrt5)/8 = 1/(2 sqrt5) = 0.2236 (0.447 is 1/sqrt5).
    # Corrected constant: J0(1/(2 sqrt 5)) <= 0.9876 (the rho of lem:bigt).  See audit, R3.
    w = 4 / (8 * arb(5).sqrt())
    print(f"      NOTE (audit R3): (4/sqrt 5)/8 = {s(w)}, not 0.447 as printed in the sandwich step; "
          f"J0 at the true threshold = {s(w.bessel_j(0))} (paper: 0.9510 = J0(0.447))")
    check("corrected sandwich constant: (4/sqrt 5)/8 < 0.447 and J0(1/(2 sqrt 5)) <= 0.9876 < 1",
          w < A("0.447") and w.bessel_j(0) <= A("0.9876") and w.bessel_j(0) > A("0.9510"),
          f"J0(0.447) = {s(A('0.447').bessel_j(0))}")
    check("erf(2.65) >= 0.9998, sqrt(450)/8 >= 2.65, 0.9973*0.9548*0.996 >= 0.90",
          A("2.65").erf() >= A("0.9998") and arb(450).sqrt() / 8 >= A("2.65")
          and A("0.9973") * A("0.9548") * A("0.996") >= A("0.90"))
    b = arb(450)
    check("e^(-beta/64) sqrt(beta)/(4 sqrt pi) <= 0.003 at beta = 450 (decreasing for beta > 32)",
          (-b / 64).exp() * b.sqrt() / (4 * PI.sqrt()) <= A("0.003"))
    check("1/2 + 1/64 + 1/1440 < 0.517", arb(1) / 2 + arb(1) / 64 + arb(1) / 1440 < A("0.517"))


def section_B() -> dict:
    print("[B] TAIL_an(300) < MAIN_an(300)")
    out = {}
    M = arb(300)
    for q, claim in ((4, "5e-4"), (3, "3e-3")):
        t1, t2, mn = T1(q, M), T2(q, M), MAIN_an(q, M)
        ratio = (t1 + t2) / mn
        print(f"      q={q}: n_1 = {s(n_i(q, M, 1))}, n_2 = {s(n_i(q, M, 2))}, "
              f"D_- = {s(Dminus(q, M))}, D_+ = {s(Dplus(q, M))}, beta_+ = {s(beta_plus(q, M))}")
        print(f"            MAIN_an = {s(mn, 12)}  T_1 = {s(t1, 12)}  T_2 = {s(t2, 12)}")
        check(f"q={q}: TAIL_an(300)/MAIN_an(300) < {claim} < 1",
              ratio < A(claim), f"ratio = {s(ratio, 12)}")
        check(f"q={q}: n_2(299) - 10 >= 8  (>= 18 ordinates below 2(M-1): lem:diffint, t >= 16 step)",
              n_i(q, arb(299), 2) - 10 >= 8, f"n_2(299) = {s(n_i(q, arb(299), 2))}")
        # the script's weakened variant n_i - 9/2
        w = arb(9) / 2
        t1w = t1 * R1 ** (-w) * R2 ** (-w)
        t2w = t2 * R2 ** (-w)
        claimw = {4: "1.9e-3", 3: "1.2e-2"}[q]
        check(f"q={q}: weakened variant (n_i - 9/2) ratio < {claimw}",
              (t1w + t2w) / mn < A(claimw), f"= {s((t1w + t2w) / mn, 6)}")
        out[q] = ratio
    return out


def section_C() -> dict:
    print("[C] logarithmic-derivative constants at M = 300 (worst case of M >= 300)")
    M = arb(300)
    kap = {}
    for q in (3, 4):
        l = ell(q, M)
        k1 = -(640 / (M * P1(q, M)) + dn_i(q, M, 1) * R1.log() + dn_i(q, M, 2) * R2.log())
        k2 = -(5 / (M * (Dplus(q, M) / 2 + 4 * ell0(q, M))) + dn_i(q, M, 2) * R2.log())
        k3 = 3 / (2 * M) * (arb(35) / 34 + 1 / l)
        kap[q] = (k1, k2, k3)
        print(f"      q={q}: kappa_1 = {s(k1, 10)}  kappa_2 = {s(k2, 10)}  kappa_3 = {s(k3, 10)}")
        check(f"q={q}: kappa_1 - kappa_3 >= 0.14 > 0", k1 - k3 >= A("0.14"), f"= {s(k1 - k3, 10)}")
        check(f"q={q}: kappa_2 - kappa_3 >= 0.08 > 0", k2 - k3 >= A("0.08"), f"= {s(k2 - k3, 10)}")
        # the paper's intermediate constants
        lY1, lY2 = ell(q, M / SQ3), ell(q, 2 * M / SQ3)
        check(f"q={q}: l_(Y_1) >= 4.4, l_(Y_2) >= 5.1, l > 1",
              lY1 >= A("4.4") and lY2 >= A("5.1") and l > 1, f"{s(lY1)}, {s(lY2)}, l = {s(l)}")
        paper1 = -(R1.log() * -1) * (lY1 - 1) / A("27.8") - (-R2.log()) * (lY2 - 1) / A("13.9") \
            + 640 / (M * P1(q, M))
        paper2 = -(-R2.log()) * (lY2 - 1) / A("13.9") + A("0.001")
        check(f"q={q}: paper's display for the first summand <= -0.15", paper1 <= A("-0.15"),
              f"= {s(paper1, 8)}")
        check(f"q={q}: paper's display for the second summand <= -0.09", paper2 <= A("-0.09"),
              f"= {s(paper2, 8)}")
        check(f"q={q}: 5/(M(D_+/2 + 4 l_0)) <= 0.001",
              5 / (M * (Dplus(q, M) / 2 + 4 * ell0(q, M))) <= A("0.001"))
        check(f"q={q}: (3/2)(1.03 l + 1)/(M l) <= 0.007",
              A("1.5") * (A("1.03") * l + 1) / (M * l) <= A("0.007"),
              f"= {s(A('1.5') * (A('1.03') * l + 1) / (M * l), 8)}")
        check(f"q={q}: 16 sqrt 3 <= 27.8 and 8 sqrt 3 <= 13.9",
              16 * SQ3 <= A("27.8") and 8 * SQ3 <= A("13.9"))
        # exact log-derivatives at M = 300, for the record
        e1 = dlog_T1(q, M) - dlog_MAIN(q, M)
        e2 = dlog_T2(q, M) - dlog_MAIN(q, M)
        print(f"            exact at M=300: (log T_1/MAIN)' = {s(e1, 8)}, "
              f"(log T_2/MAIN)' = {s(e2, 8)}, (log MAIN)' = {s(dlog_MAIN(q, M), 8)}, "
              f"T_2/T_1 = {s(T2(q, M) / T1(q, M), 4)}")
        check(f"q={q}: exact derivatives at 300 are <= -(kappa_i - kappa_3)",
              e1 <= -(k1 - k3) and e2 <= -(k2 - k3))
    return kap


def section_D(kap: dict, top: float = 1.0e7, ratio: float = 1.01) -> None:
    print(f"[D] independent cell sweep of [300, {top:g}], ratio {ratio}")
    for q in (3, 4):
        grid = [300.0]
        while grid[-1] < top:
            grid.append(min(grid[-1] * ratio, top))
        k1, k2, k3 = kap[q]
        worst_ratio = arb(0)
        worst_d1 = arb(-1)
        worst_d2 = arb(-1)
        ok_ratio = ok_d = ok_int = True
        prev = None
        for a, b in zip(grid[:-1], grid[1:]):
            Ma, Mb = arb(a), arb(b)
            # (i) monotone enclosure of the ratio on the cell
            r = (T1(q, Ma, Mb) + T2(q, Ma, Mb)) / MAIN_an(q, Ma, Mb)
            if not r < 1:
                ok_ratio = False
            if r.upper() > worst_ratio.upper():
                worst_ratio = r
            # (ii) exact log-derivatives on the whole cell as a ball
            Mball = Ma.union(Mb)
            dm = dlog_MAIN(q, Mball)
            d1 = dlog_T1(q, Mball) - dm
            d2 = dlog_T2(q, Mball) - dm
            if not (d1 <= A("-0.14") and d2 <= A("-0.08")):
                ok_d = False
            if d1.upper() > worst_d1.upper():
                worst_d1 = d1
            if d2.upper() > worst_d2.upper():
                worst_d2 = d2
            # (iii) integrated form of the kappa bounds between grid points
            cur = ((T1(q, Mb) / MAIN_an(q, Mb)).log(), (T2(q, Mb) / MAIN_an(q, Mb)).log())
            if prev is None:
                prev = ((T1(q, Ma) / MAIN_an(q, Ma)).log(), (T2(q, Ma) / MAIN_an(q, Ma)).log())
            if not (cur[0] - prev[0] <= -(k1 - k3) * (Mb - Ma)
                    and cur[1] - prev[1] <= -(k2 - k3) * (Mb - Ma)):
                ok_int = False
            prev = cur
        check(f"q={q}: TAIL_an/MAIN_an < 1 on all {len(grid)-1} cells of [300, {top:g}]",
              ok_ratio, f"largest cell bound {s(worst_ratio, 6)} (first cell)")
        check(f"q={q}: exact (log T_1/MAIN_an)' <= -0.14 and (log T_2/MAIN_an)' <= -0.08 on every cell",
              ok_d, f"largest cell values {worst_d1.upper().str(6, radius=False)}, "
                    f"{worst_d2.upper().str(6, radius=False)}")
        check(f"q={q}: log(T_i/MAIN_an) drops by >= (kappa_i - kappa_3)(M_b - M_a) across every cell",
              ok_int)


# --------------------------------------------------------------------------
SERIES_N = 96


def _prime_power_base(n: int):
    for p in range(2, int(math.isqrt(n)) + 1):
        if n % p == 0:
            k = n
            while k % p == 0:
                k //= p
            return p if k == 1 else None
    return n


def chi(q: int, n: int) -> int:
    r = n % q
    if q == 4:
        return 1 if r == 1 else (-1 if r == 3 else 0)
    return 1 if r == 1 else (-1 if r == 2 else 0)


def sds(M: arb, q: int):
    """Exact closed forms (sigma^2, D, S4) at weight M >= 10 as balls.
    xi'/xi(s) = (1/2)log(q/pi) + (1/2)psi((s+1)/2) + L'/L(s), s = M + 1/2."""
    sv = M + A("0.5")
    x = (sv + 1) / 2
    ll = arb(0)
    ll1 = arb(0)
    for n in range(2, SERIES_N + 1):
        c = chi(q, n)
        p = _prime_power_base(n)
        if c == 0 or p is None:
            continue
        ns = (-sv * arb(n).log()).exp()
        ll -= c * arb(p).log() * ns
        ll1 += c * arb(p).log() * arb(n).log() * ns
    N = arb(SERIES_N)
    u = sv - 1
    ln = N.log()
    pref = ((1 - sv) * ln).exp()
    tail1 = pref * (ln / u + 1 / (u * u))
    tail2 = pref * (ln * ln / u + 2 * ln / (u * u) + 2 / (u * u * u))
    ll += arb(0, tail1.upper())
    ll1 += arb(0, tail2.upper())
    trigamma = x.zeta_hurwitz2() if hasattr(x, "zeta_hurwitz2") else None
    if trigamma is None:
        from flint import acb
        trigamma = acb(2).zeta(acb(x)).real
    d1 = (arb(q) / PI).log() / 2 + x.digamma() / 2 + ll
    d2 = trigamma / 4 + ll1
    return 4 * M * d1, 4 * d1 + 4 * M * d2, 16 * 4 * M * d1 - 64 * M * M * d2


def section_E() -> None:
    print("[E] explicit closed-form bounds vs exact closed forms (sample weights)")
    pts = ["10", "12.5", "20", "34.7", "50", "100", "173.21", "299", "300", "320", "346.42",
           "350", "600", "1000", "5000", "20000", "100000"]
    for q in (3, 4):
        ok10 = ok20 = oklow = True
        worstD = arb(0)
        for p in pts:
            M = A(p)
            s2, D, S4 = sds(M, q)
            l = ell(q, M)
            dev = (D - 2 * l - 2) * M
            if abs(dev).upper() > abs(worstD).upper():
                worstD = dev
            if not (Dminus(q, M) <= D <= Dplus(q, M) and 2 * M * l <= s2 <= 2 * M * l + A("3.1")):
                ok10 = False
            if M >= 20:
                if not (S4 <= 32 * M * (l - 1) + 100 and S4 >= 32 * M * l - 32 * M - 34
                        and s2 <= 2 * M * l + A("1.4")):
                    ok20 = False
                if not ((S4 / 8 - s2) / 32 >= M * (l - 2) / 16 - A("0.5")):
                    oklow = False
        check(f"q={q}: D_- <= D <= D_+ and 2Ml <= sigma^2 <= 2Ml + 3.1 at {len(pts)} weights in [10, 1e5]",
              ok10, f"max |M (D - 2l - 2)| = {abs(worstD).upper().str(5, radius=False)} (claim <= 8)")
        check(f"q={q}: S_4 <= 32M(l-1)+100, S_4 >= 32Ml-32M-34, sigma^2 <= 2Ml+1.4 for M >= 20", ok20)
        check(f"q={q}: (S_4/8 - sigma^2)/32 >= Y(l_Y - 2)/16 - 1/2 for Y >= 20", oklow)
        # beta <= beta_+ and D >= D_- at 300
        M = arb(300)
        s2, D, S4 = sds(M, q)
        print(f"      q={q}, M=300: sigma^2 = {s(s2, 12)}, D = {s(D, 12)}, S_4 = {s(S4, 12)}; "
              f"beta = {s(s2 / 2 + S4 / 1088, 10)} <= beta_+ = {s(beta_plus(q, M), 10)}")


def load_zeros(q: int):
    path = os.path.join(ROOT, "prototype", "explore", {4: "chi4_zeros.txt", 3: "chi3_zeros.txt"}[q])
    return [float(x) for x in open(path, encoding="utf-8").read().split()]


def section_F() -> None:
    print("[F] zero lists (sanity only; not used by the proof)")
    import bisect
    for q in (3, 4):
        g = load_zeros(q)
        gmax = g[-1]
        ok1 = ok2 = ok3 = True
        slack1 = slack2 = float("inf")
        M = 35.0
        while 2 * M < gmax:
            n1 = float(n_i(q, arb(M), 1).upper())
            n2 = float(n_i(q, arb(M), 2).upper())
            N1 = bisect.bisect_right(g, M)
            N2 = bisect.bisect_right(g, 2 * M)
            ok1 &= n1 <= N1
            ok2 &= n2 <= N2
            slack1 = min(slack1, N1 - n1)
            slack2 = min(slack2, N2 - n2)
            M += 0.25
        T = 1.5
        while T < gmax:
            Np = 0.5 * T * math.log(q * (T + 2) / (2 * math.pi)) + 0.57 * T
            ok3 &= bisect.bisect_right(g, T) <= Np
            T += 0.25
        check(f"q={q}: n_1(M) <= N(M) for M in [35, {gmax/2:.1f}]", ok1, f"min slack {slack1:.1f}")
        check(f"q={q}: n_2(M) <= N(2M) for M in [35, {gmax/2:.1f}]", ok2, f"min slack {slack2:.1f}")
        check(f"q={q}: N(T) <= N^+(T) for T in [1.5, {gmax:.1f}]", ok3)
        M = 300.0
        print(f"      q={q}: N(300) = {bisect.bisect_right(g, 300.0)}, N(600) = {bisect.bisect_right(g, 600.0)} "
              f"vs n_1 = {float(n_i(q, arb(M), 1).mid()):.2f}, n_2 = {float(n_i(q, arb(M), 2).mid()):.2f}")


def load_cert(q: int):
    path = os.path.join(ROOT, "monotonicity", "cert", f"results_q{q}_stretch.txt")
    cells = []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or "PASS" not in line:
            continue
        f = line.split()
        d = dict(x.split("=", 1) for x in f[2:] if "=" in x)
        cells.append((f[0], f[1], d["hi"], d["lo"]))
    cells.sort(key=lambda c: float(c[0]))
    return cells


def section_G() -> None:
    print("[G] cross-check against the direct certificate of delta'")
    for q, pts in ((4, (300, 310, 320)), (3, (300, 310, 320, 335, 350))):
        cells = load_cert(q)
        top = float(cells[-1][1])
        # point checks
        for Mv in pts:
            M = arb(Mv)
            U = -(MAIN_an(q, M) - T1(q, M) - T2(q, M))        # delta'(M) <= U
            s2, D, _ = sds(M, q)
            asym = -D / (2 * (2 * PI).sqrt() * s2 * s2.sqrt())
            for a, b, hi, lo in cells:
                if float(a) <= Mv <= float(b):
                    check(f"q={q} M={Mv}: cell [{a},{b}] certified [{lo}, {hi}]; analytic bound "
                          f"U = {s(U, 7)} >= lo (consistent)", U >= A(lo),
                          f"U/lo = {s(U / A(lo), 5)}, U/hi = {s(U / A(hi), 5)}, "
                          f"-D/(2 sqrt(2pi) sigma^3) = {s(asym, 7)}")
        # sharp tiny cells (sharp_cells.txt, produced by sharp_cells.sh)
        sp = os.path.join(HERE, "sharp_cells.txt")
        if os.path.exists(sp):
            for line in open(sp, encoding="utf-8"):
                f = line.split()
                if len(f) < 5 or f[0] != f"q={q}" or "PASS" not in line:
                    continue
                d = dict(x.split("=", 1) for x in f[3:] if "=" in x)
                Ma, Mb = A(f[1]), A(f[2])
                # on the cell:  U(M) = -(MAIN_an - TAIL_an)(M) >= -K3 D_-(Mb) G(beta_+(Ma))
                Umin = -K3 * Dminus(q, Mb) * G(beta_plus(q, Ma))
                check(f"q={q}: sharp cell [{f[1]},{f[2]}] certified delta' in [{d['lo']}, {d['hi']}]; "
                      f"STRONG: hi <= min_cell U (the analytic bound is implied by the direct certificate)",
                      A(d["hi"]) <= Umin,
                      f"min_cell U = {s(Umin, 7)}; U/hi = {s(Umin / A(d['hi']), 5)}")
        else:
            print("      (sharp_cells.txt not found: run sharp_cells.sh for the strong check)")
        # every cell of the overlap [300, top]
        okc = True
        nstrong = ncell = 0
        worst = arb(0)
        for a, b, hi, lo in cells:
            if float(b) <= 300:
                continue
            Ma, Mb = A(a), A(b)
            if Ma < 300:
                Ma = arb(300)
            # sup over the cell of (MAIN_an - TAIL_an) <= K3 D_-(Mb) G(beta_+(Ma))
            sup_claim = K3 * Dminus(q, Mb) * G(beta_plus(q, Ma))
            ncell += 1
            if not sup_claim <= -A(lo):
                okc = False
            if sup_claim <= -A(hi):
                nstrong += 1
            rr = sup_claim / (-A(lo))
            if rr.upper() > worst.upper():
                worst = rr
        check(f"q={q}: on all {ncell} certified cells meeting [300, {top:g}]: "
              f"sup_cell (MAIN_an - TAIL_an) <= -lo  (analytic claim never exceeds certified |delta'|_max)",
              okc, f"max ratio {s(worst, 5)}; on {nstrong} cells even <= -hi")


def main() -> None:
    section_A()
    section_A2()
    section_B()
    kap = section_C()
    section_D(kap)
    section_E()
    section_F()
    section_G()
    if FAILS:
        print(f"FAILED: {len(FAILS)} check(s): " + "; ".join(FAILS))
        sys.exit(1)
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
