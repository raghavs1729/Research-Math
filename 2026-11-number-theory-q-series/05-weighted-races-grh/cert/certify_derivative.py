#!/usr/bin/env python3
"""Fail-closed Arb certificate that delta'(M) < 0 on whole M-cells.

WHAT IS PROVED
==============
Let chi be the real primitive character mod q, q in {3, 4}, M = m + 1/2 > 0,
c_g(M) = 4M/sqrt(M^2 + g^2) (= 2 a_gamma in paper/weighted-races.tex) and

    Phi_M(t) = prod_{g>0} J_0(c_g(M) t),
    delta(M) = 1/2 + (1/pi) int_0^inf (sin t / t) Phi_M(t) dt,

the product running over the positive ordinates g of L(s, chi).  For each
cell [M_a, M_b] this program prints an interval [lo, hi] and the word PASS
only if it has verified, in outward-rounded ball arithmetic, that

    lo <= delta'(M) <= hi < 0        for EVERY real M in [M_a, M_b].

HYPOTHESES (not proved here)
============================
(H1) GRH(chi) and LI(chi): they make delta(M) the logarithmic density of the
     weighted race (Theorem thm:interp of the paper).  Everything this file
     does concerns the function delta(M) defined by the displayed formula;
     for that function only GRH (real ordinates) and the two items below are
     used.
(H2) The variance identity prop:varid / eq:varidq and its derivative
     lem:Dprop(i) (both proved in the paper under GRH):
        sigma^2(M) := sum_g c_g^2/2      = 4M (log Lambda)'(M + 1/2),
        D(M)       := sum_g c_g c_g'     = 4 (log Lambda)'(M+1/2)
                                           + 4M (log Lambda)''(M+1/2),
     and lem:diffint(a) (pointwise formula eq:dPhi for d/dM Phi).  The
     differentiation under the integral sign, lem:diffint(c), is NOT taken
     from the paper: its hypothesis (8 ordinates below 2 M_a) fails for
     small M.  It is re-derived in METHOD.md section 4 from the integrable
     dominant of the large-t remainder below, which this program verifies
     for every cell (it requires K >= 8 decaying Bessel factors).
(H3) The zero files chi4_zeros.txt / chi3_zeros.txt: every listed decimal is
     within 1e-20 of an ordinate, and the list is complete up to its last
     entry (certified separately by prototype/explore/certify_arb.py --K 0).
(H4) Correctness of FLINT/Arb ball arithmetic (python-flint).

METHOD (details and proofs of every inequality: METHOD.md)
==========================================================
Split the zeros into a head H (the first N listed ordinates) and a tail
(all other zeros, listed or not; they satisfy g > g_N).  Phi = P * Psi.

 * 0 <= t <= T1, where c_{g_N}(M_b) T1 <= z0 = 1:
       d/dM Phi = Psi * (dP - P * B_T),
       G e^{-C4 S4 t^4} <= Psi <= G := exp(-V_T t^2/2),
       (t^2/2) D_T <= B_T <= (t^2/2) D_T (1 + 2 kappa c_{g_N}^2 t^2),
   with V_T = sigma^2 - sum_H c^2/2 and D_T = D - sum_H c c' EXACT closed
   forms, both nondecreasing in M (so they are enclosed over the cell by
   their endpoint values), S4 <= 2 c_{g_N}^2 V_T.  The resulting integrand is
   evaluated on boxes [t_a, t_b] x [M_a, M_b] in ball arithmetic, and
       int_{t_a}^{t_b} F dt  in  (t_b - t_a) * F([t_a,t_b] x [M_a,M_b]).
 * t >= T1:  |d/dM Phi| <= (t^2/2) D(M) * max_g prod_{g' in S\\{g}} |J_0|
   with |J_1(z)| <= z/2 and |J_0(x)| <= x^{-1/2} (Watson 13.74), giving the
   closed-form remainder E2.

Nothing is sampled: no value of M or t is skipped, and no floating-point
number enters a certified quantity (floats are used only to steer the
adaptive subdivision, which cannot affect validity).  Any non-finite ball,
failed positivity or failed hypothesis check raises.

USAGE
=====
  certify_derivative.py cell  --q 4 --Ma 3.6 --Mb 3.636       one cell
  certify_derivative.py range --q 4 --lo 3.6 --hi 17.4 --ratio 1.005 \
                               --tol 0.02 --out results_q4.txt --workers 3
  certify_derivative.py aux   --out aux_results.txt
  certify_derivative.py selftest
"""

from __future__ import annotations

import argparse
import heapq
import math
import os
import sys
import time
from decimal import ROUND_CEILING, ROUND_FLOOR, Context, Decimal
from fractions import Fraction

from flint import acb, arb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
ZERO_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "prototype", "explore"))
ZERO_FILES = {4: "chi4_zeros.txt", 3: "chi3_zeros.txt"}

PREC_HI = 128          # closed forms, zero balls, head sums
PREC_LO = 64           # the integrand (balls there are wide anyway)
NSER = 20000           # prime powers n <= NSER in the Dirichlet series
SMIN_SERIES = 3        # the series is only used for s = M + 1/2 >= 3
Z0 = arb(1)            # tail-model range: c_g t <= Z0 for every tail zero
E2_TARGET = 1e-10      # steering target for the large-t remainder (float)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def fr(x) -> Fraction:
    """Exact rational from a decimal string / int / Fraction."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, float):
        raise TypeError("floats are not accepted as certified inputs")
    return Fraction(str(x))


def ex(x: Fraction) -> arb:
    """Ball containing the rational x (exact when x is dyadic and short)."""
    return arb(x.numerator) / arb(x.denominator)


def hull(a: arb, b: arb) -> arb:
    """Ball containing both balls and everything between them."""
    return a.union(b)


def rad_ball(r: arb) -> arb:
    """[0 +/- r_upper]."""
    u = r.upper()
    if not u.is_finite():
        raise ArithmeticError("non-finite radius")
    return arb(0, u)


def need_finite(x: arb, label: str) -> None:
    if not x.is_finite():
        raise ArithmeticError(f"non-finite ball: {label}")


def fl(x: arb) -> float:
    return float(x.mid())


def _dfmt(x: float, digits: int, rounding, sci: bool) -> str:
    """Decimal string of the double x, rounded in the given DIRECTION (the
    conversion double -> Decimal is exact, so the printed number lies on the
    stated side of x)."""
    d = Decimal(x)
    if not d.is_finite():
        return str(x)
    if sci:
        c = Context(prec=digits + 1, rounding=rounding)
        return f"{c.create_decimal(d):.{digits}e}"
    return str(d.quantize(Decimal(1).scaleb(-digits), rounding=rounding))


def e_up(x: float, digits: int = 6) -> str:
    return _dfmt(x, digits, ROUND_CEILING, True)


def e_dn(x: float, digits: int = 6) -> str:
    return _dfmt(x, digits, ROUND_FLOOR, True)


def f_up(x: float, digits: int = 12) -> str:
    return _dfmt(x, digits, ROUND_CEILING, False)


def f_dn(x: float, digits: int = 12) -> str:
    return _dfmt(x, digits, ROUND_FLOOR, False)


def up(x: arb) -> float:
    """Float >= upper endpoint (display/steering only)."""
    return math.nextafter(float(x.upper()), math.inf)


def dn(x: arb) -> float:
    return math.nextafter(float(x.lower()), -math.inf)


# --------------------------------------------------------------------------
# zeros
# --------------------------------------------------------------------------
_ZCACHE: dict[int, list[arb]] = {}


def load_zeros(q: int) -> list[arb]:
    if q in _ZCACHE:
        return _ZCACHE[q]
    old = ctx.prec
    ctx.prec = PREC_HI
    try:
        path = os.path.join(ZERO_DIR, ZERO_FILES[q])
        with open(path, encoding="ascii") as fh:
            vals = [arb(line.strip(), "1e-20") for line in fh if line.strip()]
        if len(vals) < 100:
            raise ArithmeticError("zero file too short")
        if not vals[0].lower() > 1:
            raise ArithmeticError("first ordinate not > 1")
        for a, b in zip(vals, vals[1:]):
            if not a.upper() < b.lower():
                raise ArithmeticError("zero balls not strictly ordered")
    finally:
        ctx.prec = old
    _ZCACHE[q] = vals
    return vals


# --------------------------------------------------------------------------
# closed forms  (log Lambda)'(s), (log Lambda)''(s)   for chi odd mod q
# --------------------------------------------------------------------------
def chi_val(q: int, n: int) -> int:
    if q == 4:
        return {1: 1, 3: -1}.get(n % 4, 0)
    if q == 3:
        return {1: 1, 2: -1}.get(n % 3, 0)
    raise ValueError("q must be 3 or 4")


_PP: dict[int, list[tuple[int, int, int]]] = {}


def prime_powers(q: int) -> list[tuple[int, int, int]]:
    """(n, p, chi(n)) for prime powers n = p^k <= NSER with chi(n) != 0."""
    if q in _PP:
        return _PP[q]
    sieve = bytearray([1]) * (NSER + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(NSER ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    out = []
    for p in range(2, NSER + 1):
        if sieve[p]:
            n = p
            while n <= NSER:
                c = chi_val(q, n)
                if c:
                    out.append((n, p, c))
                n *= p
    out.sort()
    _PP[q] = out
    return out


def loglam_derivs_series(q: int, s: arb) -> tuple[arb, arb]:
    """Enclosures of (log Lambda)'(s) and (log Lambda)''(s), real s >= 3.

    Lambda(s) = (q/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi)   (chi odd), so
      (log Lambda)'  = log(q/pi)/2 + psi((s+1)/2)/2 + (L'/L)(s),
      (log Lambda)'' = psi'((s+1)/2)/4 + (L'/L)'(s),
      (L'/L)(s)  = - sum_n chi(n) Lambda(n) n^{-s},
      (L'/L)'(s) =   sum_n chi(n) Lambda(n) log(n) n^{-s},
    absolutely convergent for s > 1.  Tails beyond N = NSER:
      sum_{n>N} log(n)   n^{-s} <= N^{1-s} (L/a + 1/a^2),
      sum_{n>N} log(n)^2 n^{-s} <= N^{1-s} (L^2/a + 2L/a^2 + 2/a^3),
    L = log N, a = s - 1 (integral test; both integrands decrease on
    [N, inf) because N >= e).
    """
    if not s.lower() >= SMIN_SERIES:
        raise ArithmeticError("Dirichlet-series closed form needs s >= 3")
    s1 = arb(0)
    s2 = arb(0)
    logs: dict[int, arb] = {}
    for n, p, c in prime_powers(q):
        lp = logs.get(p)
        if lp is None:
            lp = logs[p] = arb(p).log()
        ln = arb(n).log()
        term = lp * (-s * ln).exp()
        if c > 0:
            s1 += term
            s2 += term * ln
        else:
            s1 -= term
            s2 -= term * ln
    N = arb(NSER)
    L = N.log()
    a = s - 1
    pw = (-a * L).exp()
    tail1 = pw * (L / a + 1 / (a * a))
    tail2 = pw * (L * L / a + 2 * L / (a * a) + 2 / (a * a * a))
    lpl = -s1 + rad_ball(tail1)
    dlpl = s2 + rad_ball(tail2)
    x = (s + 1) / 2
    d1 = (arb(q) / arb.pi()).log() / 2 + x.digamma() / 2 + lpl
    d2 = arb(2).zeta(x) / 4 + dlpl          # psi'(x) = zeta(2, x) (Hurwitz)
    need_finite(d1, "(log Lambda)'")
    need_finite(d2, "(log Lambda)''")
    return d1, d2


def loglam_d1_complex_step(q: int, s: arb) -> arb:
    """(log Lambda)'(s) for real s >= 1 by the complex step of
    race_density_arb.sigma2_arb: L'(s) = Im L(s+ih)/h + remainder, the
    remainder bounded by Cauchy on |z-s| <= R = 1/4 with
    |L(z)| <= |z|/Re z <= (s+R)/(s-R) (Abel summation, partial sums of chi
    bounded by 1 for q = 3, 4)."""
    if not s.lower() >= 1:
        raise ArithmeticError("complex step used only for s >= 1")
    old = ctx.prec
    ctx.prec = 384          # the step h = 2^-60 costs ~120 bits (more at s = 1)
    try:
        chi = dirichlet_char(q, q - 1)
        h = arb(2) ** -60
        R = arb(1) / 4
        Ls = chi.l(acb(s)).real
        if not Ls.lower() > 0:
            raise ArithmeticError("L(s) not certified positive")
        Lh = chi.l(acb(s, h))
        dL = Lh.imag / h
        ML = (s + R) / (s - R)
        rem = ML * h ** 2 / (R ** 3 * (1 - (h / R) ** 2))
        dL = dL + rad_ball(rem)
        x = (s + 1) / 2
        out = (arb(q) / arb.pi()).log() / 2 + x.digamma() / 2 + dL / Ls
    finally:
        ctx.prec = old
    need_finite(out, "(log Lambda)' complex step")
    return out


def sigma2_any(q: int, M: arb) -> arb:
    """sigma^2(M) = 4M (log Lambda)'(M+1/2), any exact M >= 1/2."""
    s = M + arb(1) / 2
    if s.lower() >= SMIN_SERIES + 1:
        d1, _ = loglam_derivs_series(q, s)
    else:
        d1 = loglam_d1_complex_step(q, s)
    out = 4 * M * d1
    if not out.lower() > 0:
        raise ArithmeticError("sigma^2 not certified positive")
    return out


def head_sums(zeros: list[arb], M: arb) -> tuple[arb, arb]:
    """sum_H c^2/2 = sum 8M^2/(M^2+g^2)  and  sum_H c c' = sum 16 M g^2/(M^2+g^2)^2."""
    v = arb(0)
    d = arb(0)
    M2 = M * M
    for g in zeros:
        g2 = g * g
        den = M2 + g2
        v += 8 * M2 / den
        d += 16 * M * g2 / (den * den)
    return v, d


def tail_VD(q: int, zeros: list[arb], M: arb) -> tuple[arb, arb, arb, arb]:
    """(V_T, D_T, sigma^2, D) at a thin ball M (precision PREC_HI)."""
    s = M + arb(1) / 2
    d1, d2 = loglam_derivs_series(q, s)
    sig2 = 4 * M * d1
    D = 4 * d1 + 4 * M * d2
    hv, hd = head_sums(zeros, M)
    VT = sig2 - hv
    DT = D - hd
    if not VT.lower() > 0:
        raise ArithmeticError(f"tail variance not certified positive: {VT}")
    if not DT.lower() > 0:
        raise ArithmeticError(f"tail variance derivative not certified positive: {DT}")
    return VT, DT, sig2, D


# --------------------------------------------------------------------------
# Bessel enclosures on (possibly wide) balls
# --------------------------------------------------------------------------
ONE_BALL = None  # set in setup()


def bessel01(x: arb) -> tuple[arb, arb]:
    """Enclosures of J_0 and J_1 over the ball x (x > 0 throughout).

    With x0 = mid(x), r = rad(x):  Taylor with Lagrange remainder,
        |J_0(x) - J_0(x0) + J_1(x0)(x-x0)| <= r^2/2 sup|J_0''|,
        |J_1(x) - J_1(x0) - J_1'(x0)(x-x0)| <= r^2/2 sup|J_1''|,
    J_0'' = (J_2 - J_0)/2, J_1'' = (J_3 - 3 J_1)/4, |J_n| <= 1 on the real
    line, J_1'(x0) = J_0(x0) - J_1(x0)/x0.  The result is intersected with
    the trivial bound |J_n| <= 1 and, for J_0, with Watson's envelope
    |J_0(x)| <= x^{-1/2}.
    """
    x0 = arb(x.mid())
    r = arb(x.rad())
    if not x0 > 0:
        raise ArithmeticError("bessel01: midpoint not positive")
    j0 = x0.bessel_j(0)
    j1 = x0.bessel_j(1)
    d1 = j0 - j1 / x0
    hr2 = r * r / 2
    e0 = abs(j1) * r + hr2
    e1 = abs(d1) * r + hr2
    J0 = j0 + rad_ball(e0)
    J1 = j1 + rad_ball(e1)
    if not J0.abs_upper() <= 1:
        J0 = J0.intersection(ONE_BALL)
    if not J1.abs_upper() <= 1:
        J1 = J1.intersection(ONE_BALL)
    xl = x.lower()
    if xl > 1:
        w = 1 / xl.sqrt()
        if not J0.abs_upper() <= w:
            J0 = J0.intersection(rad_ball(w))
    return J0, J1


# --------------------------------------------------------------------------
# one M-cell
# --------------------------------------------------------------------------
class Cell:
    """All M-dependent data for a cell [Ma, Mb] (exact rationals)."""

    def __init__(self, q: int, Ma: Fraction, Mb: Fraction, nhead: int | None = None):
        if not (0 < Ma < Mb):
            raise ValueError("need 0 < Ma < Mb")
        self.q, self.Ma, self.Mb = q, Ma, Mb
        allz = load_zeros(q)
        N = len(allz) if nhead is None else nhead
        if not (20 <= N <= len(allz)):
            raise ValueError("bad head size")
        self.N = N
        zeros = allz[:N]
        gN = zeros[-1]

        ctx.prec = PREC_HI
        a, b = ex(Ma), ex(Mb)
        M = hull(a, b)
        self.M = M
        # monotonicity of D_T on the cell needs g^2 > 3 M^2 for tail zeros
        if not gN.lower() ** 2 > 3 * b.upper() ** 2:
            raise ArithmeticError("tail zeros not above sqrt(3)*Mb")
        VTa, DTa, s2a, Da = tail_VD(q, zeros, a)
        VTb, DTb, s2b, Db = tail_VD(q, zeros, b)
        if not (VTa.upper() <= VTb.upper() and DTa.upper() <= DTb.upper()):
            raise ArithmeticError("endpoint tail quantities not ordered")
        self.V = hull(VTa, VTb)          # V_T(M) for all M in the cell
        self.DT = hull(DTa, DTb)         # D_T(M) for all M in the cell
        self.sig2 = hull(s2a, s2b)       # sigma^2 increasing in M
        self.info = dict(VT=(dn(VTa), up(VTb)), DT=(dn(DTa), up(DTb)),
                         sig2=(dn(s2a), up(s2b)))

        # coefficients over the cell (M occurs once in each expression)
        self.cs = []
        self.cps = []
        for g in zeros:
            c = 4 / (1 + (g / M) ** 2).sqrt()
            cp = (4 / g) / ((1 + (M / g) ** 2) ** 3).sqrt()
            self.cs.append(c)
            self.cps.append(cp)
        # c at endpoints, for the Watson remainder
        self.cs_a = [4 / (1 + (g / a) ** 2).sqrt() for g in zeros]
        self.c1_b = 4 / (1 + (zeros[0] / b) ** 2).sqrt()
        cG = 4 / (1 + (gN / b) ** 2).sqrt()       # c_{g_N}(Mb) >= c_g(M), tail g
        self.cG = arb(cG.upper())
        # D(M) <= sum_H c c' + D_T(Mb)
        hd = arb(0)
        for c, cp in zip(self.cs, self.cps):
            hd += c * cp
        self.Dmax = arb(hd.upper()) + arb(DTb.upper())
        self.D_info = (dn(Da), up(Db), up(self.Dmax))

        # tail-model constants at z0 (monotone ratios, see METHOD.md)
        j0z, j1z = Z0.bessel_j(0), Z0.bessel_j(1)
        if not j0z.lower() > 0:
            raise ArithmeticError("J0(z0) not positive")
        C4 = (-(j0z.log()) - Z0 ** 2 / 4) / Z0 ** 4
        kap = (j1z / j0z - Z0 / 2) / Z0 ** 3
        if not (C4.lower() > 0 and kap.lower() > 0):
            raise ArithmeticError("tail-model constants not positive")
        S4max = self.cG ** 2 * 2 * arb(self.V.upper())
        self.C4S4 = arb((C4 * S4max).upper())
        self.eta2 = arb((2 * kap * self.cG ** 2).upper())
        self.T1max = Z0 / self.cG                     # ball; use lower()
        self.zeros = zeros
        ctx.prec = PREC_LO

    # ---- large-t remainder ------------------------------------------------
    def remainder(self, T1: Fraction) -> tuple[arb, int]:
        """E2 >= | int_{T1}^inf (sin t/t) dPhi/dM dt | for all M in the cell."""
        T = ex(T1)
        if not (self.cG * T).upper() <= Z0.lower():
            raise ArithmeticError("T1 outside the tail-model range")
        prod = arb(1)
        K = 0
        for ca in self.cs_a:
            x = ca * T
            if not x.lower() >= 1:
                break                       # c_g decreasing in g: stop at first failure
            prod *= 1 / x.sqrt()            # (c_g(Ma) T1)^{-1/2} <= 1
            K += 1
        if K < 8:
            return arb("inf"), K
        if not (self.c1_b * T).lower() >= 1:
            raise ArithmeticError("c_max T1 < 1")
        # int_{T1}^inf t^{1-(K-1)/2} dt = T1^{2-(K-1)/2} / ((K-5)/2)
        # bound = (Dmax/2) * prod_S (c T1)^{-1/2} * (c_max T1)^{1/2} * T1^2 * 2/(K-5)
        E2 = self.Dmax / 2 * prod * (self.c1_b * T).sqrt() * T * T * 2 / (K - 5)
        need_finite(E2, "E2")
        return arb(E2.upper()), K

    def choose_T1(self) -> tuple[Fraction, arb, int]:
        tmax = float(self.T1max.lower())
        tmaxq = Fraction(math.floor(tmax * 1024 * (1 - 1e-9)), 1024)
        ts = 1.0 / math.sqrt(self.info["sig2"][0])
        cand = []
        x = max(3.0 * ts, 0.25)
        while x < tmax:
            cand.append(Fraction(math.ceil(x * 1024), 1024))
            x *= 1.25
        cand = [c for c in cand if c < tmaxq] + [tmaxq]
        best = None
        for T1 in cand:
            E2, K = self.remainder(T1)
            if E2.is_finite() and up(E2) <= E2_TARGET:
                return T1, E2, K
            best = (T1, E2, K)
        if best is None or not best[1].is_finite():
            raise ArithmeticError("no admissible T1 with a finite remainder")
        return best

    # ---- integrand on a t-box ----------------------------------------------
    def integrand(self, ta: arb, tb: arb) -> tuple[arb, arb]:
        """Ball containing F(t, M) = sinc(t) * dPhi_M/dM (t) for all t in
        [ta, tb] and all M in the cell; second value: upper bound for the
        part of |F| due to the tail model (diagnostic only)."""
        t = hull(ta, tb)
        n = self.N
        J0 = [None] * n
        J1 = [None] * n
        cs = self.cs
        for i in range(n):
            J0[i], J1[i] = bessel01(cs[i] * t)
        pre = [None] * (n + 1)
        p = arb(1)
        pre[0] = p
        for i in range(n):
            p = p * J0[i]
            pre[i + 1] = p
        P = p
        suf = arb(1)
        acc = arb(0)
        cps = self.cps
        for i in range(n - 1, -1, -1):
            acc += cps[i] * J1[i] * pre[i] * suf
            suf = suf * J0[i]
        dP = -t * acc
        t2 = t * t
        G = (-self.V * t2 / 2).exp()
        em = (-self.C4S4 * t2 * t2).exp()
        one = arb(1)
        Psi = G * hull(em, one)
        BT0 = t2 / 2 * self.DT
        BT = BT0 * hull(one, 1 + self.eta2 * t2)
        sc = t.sinc()
        F = sc * Psi * (dP - P * BT)
        need_finite(F, "integrand")
        model = abs(sc) * G * ((1 - em) * abs(dP - P * BT) + abs(P) * abs(BT0) * self.eta2 * t2)
        return F, arb(model.upper())

    # ---- adaptive interval Riemann sum on [0, T1] ---------------------------
    def integrate(self, T1: Fraction, tol: float, max_cells: int):
        ts = 1.0 / math.sqrt(self.info["sig2"][0])
        hres = Fraction(max(1, int(0.25 * ts * 4096)), 4096)
        # initial breakpoints: uniform on [0, 6 ts], geometric beyond
        pts = [Fraction(0)]
        h0 = Fraction(max(1, int(ts / 4 * 4096)), 4096)
        x = h0
        while x < T1 and x <= 6 * ts:
            pts.append(x)
            x += h0
        while x < T1:
            pts.append(x)
            x = Fraction(math.ceil(float(x) * 1.2 * 4096), 4096)
        pts.append(T1)

        heap = []
        done = []
        serial = 0
        tot_mid = 0.0
        tot_rad = 0.0
        ncell = 0

        def evalcell(a: Fraction, b: Fraction):
            F, model = self.integrand(ex(a), ex(b))
            w = ex(b - a)
            return (a, b, F * w, model * w, F)

        def push(item, saturated: bool):
            nonlocal serial, tot_mid, tot_rad, ncell
            c = item[2]
            r = float(c.rad())
            tot_mid += fl(c)
            tot_rad += r
            ncell += 1
            if saturated or r == 0.0:
                done.append(item)
            else:
                heapq.heappush(heap, (-r, serial, item))
                serial += 1

        for a, b in zip(pts, pts[1:]):
            push(evalcell(a, b), False)

        while heap:
            if tot_rad <= tol * abs(tot_mid) or ncell >= max_cells:
                break
            negr, _, item = heapq.heappop(heap)
            a, b, c, _, F = item
            tot_mid -= fl(c)
            tot_rad -= -negr
            ncell -= 1
            mid = (a + b) / 2
            left = evalcell(a, mid)
            right = evalcell(mid, b)
            rsum = float(left[2].rad()) + float(right[2].rad())
            resolved = (b - a) <= hres and float(F.rad()) < 0.5 * abs(fl(F))
            sat = resolved and rsum > 0.8 * (-negr)
            push(left, sat)
            push(right, sat)

        leaves = done + [it for _, _, it in heap]
        leaves.sort(key=lambda it: it[0])
        # verify the leaves tile [0, T1] exactly
        pos = Fraction(0)
        for it in leaves:
            if it[0] != pos:
                raise ArithmeticError("t-cells do not tile [0, T1]")
            pos = it[1]
        if pos != T1:
            raise ArithmeticError("t-cells do not reach T1")
        total = arb(0)
        model = arb(0)
        for it in leaves:
            total += it[2]
            model += it[3]
        need_finite(total, "integral")
        return total, arb(model.upper()), len(leaves)


def certify_cell(q: int, Ma, Mb, nhead=None, tol=0.05, max_cells=6000) -> dict:
    """Enclosure of delta'(M) valid for every M in [Ma, Mb]."""
    t0 = time.time()
    setup()
    Ma, Mb = fr(Ma), fr(Mb)
    cell = Cell(q, Ma, Mb, nhead)
    T1, E2, K = cell.choose_T1()
    I1, model, ncells = cell.integrate(T1, tol, max_cells)
    pi = arb.pi()
    enc = (I1 + rad_ball(E2)) / pi
    need_finite(enc, "delta' enclosure")
    lo, hi = dn(enc), up(enc)
    ok = bool(enc.upper() < 0) and hi < 0.0
    return dict(q=q, Ma=Ma, Mb=Mb, lo=lo, hi=hi, mid=fl(enc),
                E2=up(E2 / pi), Emodel=up(model / pi),
                rad=up(arb(enc.rad())), T1=float(T1), K=K, ncells=ncells,
                N=cell.N, VT=cell.info["VT"], DT=cell.info["DT"],
                sig2=cell.info["sig2"], D=cell.D_info,
                secs=time.time() - t0, ok=ok)


def setup() -> None:
    global ONE_BALL
    ctx.prec = PREC_LO
    ONE_BALL = arb(0, 1)


# --------------------------------------------------------------------------
# range driver
# --------------------------------------------------------------------------
def dec(x: Fraction) -> str:
    """Exact decimal string of a rational with power-of-ten denominator."""
    s = x.limit_denominator(10 ** 12)
    if s != x:
        raise ValueError("endpoint is not a short decimal")
    k = 0
    while (x * 10 ** k).denominator != 1:
        k += 1
        if k > 12:
            raise ValueError("endpoint is not a short decimal")
    n = int(x * 10 ** k)
    if k == 0:
        return str(n)
    txt = f"{n:0{k + 1}d}"
    return txt[:-k] + "." + txt[-k:]


def make_grid(lo: Fraction, hi: Fraction, ratio: float, digits: int = 4) -> list[Fraction]:
    pts = [lo]
    scale = 10 ** digits
    while True:
        nxt = Fraction(math.floor(float(pts[-1]) * ratio * scale), scale)
        if nxt <= pts[-1]:
            raise ValueError("ratio too small for the endpoint resolution")
        if nxt >= hi or float(hi - nxt) < 0.3 * float(nxt - pts[-1]):
            pts.append(hi)
            return pts
        pts.append(nxt)


def fmt_line(r: dict) -> str:
    return (f"{dec(r['Ma'])} {dec(r['Mb'])} hi={e_up(r['hi'])} lo={e_dn(r['lo'])} "
            f"E2={e_up(r['E2'], 2)} Emodel={e_up(r['Emodel'], 2)} rad={e_up(r['rad'], 2)} "
            f"T1={r['T1']:.4f} K={r['K']} N={r['N']} tcells={r['ncells']} "
            f"VT=[{e_dn(r['VT'][0], 5)},{e_up(r['VT'][1], 5)}] "
            f"DT=[{e_dn(r['DT'][0], 5)},{e_up(r['DT'][1], 5)}] "
            f"secs={r['secs']:.1f} {'PASS' if r['ok'] else 'FAIL'}")


def _work(job):
    q, Ma, Mb, nhead, tol, max_cells, depth = job
    try:
        r = certify_cell(q, Ma, Mb, nhead, tol, max_cells)
        return job, r, None
    except Exception as e:  # fail closed: report, never pass
        return job, None, f"{type(e).__name__}: {e}"


def parse_results(path: str) -> list[tuple[Fraction, Fraction, float, float, bool]]:
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding="ascii") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            try:
                a, b = Fraction(parts[0]), Fraction(parts[1])
                hi = float(parts[2].split("=")[1])
                lo = float(parts[3].split("=")[1])
            except (ValueError, IndexError):
                continue
            out.append((a, b, hi, lo, parts[-1] == "PASS"))
    return out


def run_range(q, lo, hi, ratio, out, workers, nhead, tol, max_cells, max_depth=6):
    import multiprocessing as mp
    lo, hi = fr(lo), fr(hi)
    grid = make_grid(lo, hi, ratio)
    done = {(a, b) for a, b, _, _, ok in parse_results(out) if ok}
    jobs = [(q, a, b, nhead, tol, max_cells, 0) for a, b in zip(grid, grid[1:])]

    def covered(a, b):
        # already certified by cells in the file that tile [a, b] exactly
        pos = a
        cells = sorted(c for c in done if a <= c[0] and c[1] <= b)
        for c in cells:
            if c[0] == pos:
                pos = c[1]
        return pos == b

    jobs = [j for j in jobs if not covered(j[1], j[2])]
    new = not os.path.exists(out)
    fh = open(out, "a", encoding="ascii")
    if new:
        fh.write(f"# certify_derivative.py range q={q} [{dec(lo)}, {dec(hi)}] ratio={ratio} "
                 f"tol={tol} nhead={nhead}\n"
                 "# columns: M_a M_b hi lo (lo <= delta'(M) <= hi for ALL M in [M_a,M_b]); "
                 "E2 = large-t remainder/pi; Emodel = tail-model part/pi; rad = enclosure radius\n")
        fh.flush()
    print(f"q={q}: {len(jobs)} cells to do", flush=True)
    if not jobs:
        fh.close()
        return
    t0 = time.time()
    with mp.Pool(workers) as pool:
        pending = list(jobs)
        while pending:
            retry = []
            for job, r, err in pool.imap_unordered(_work, pending):
                if r is not None and r["ok"]:
                    fh.write(fmt_line(r) + "\n")
                    fh.flush()
                    continue
                depth = job[6]
                msg = err if r is None else fmt_line(r)
                fh.write(f"# not certified at depth {depth}: {dec(job[1])} {dec(job[2])} :: {msg}\n")
                fh.flush()
                if depth < max_depth:
                    a, b = job[1], job[2]
                    mid = (a + b) / 2
                    retry.append(job[:1] + (a, mid) + job[3:6] + (depth + 1,))
                    retry.append(job[:1] + (mid, b) + job[3:6] + (depth + 1,))
                else:
                    fh.write(f"{dec(job[1])} {dec(job[2])} hi=nan lo=nan FAIL\n")
                    fh.flush()
            pending = retry
    fh.close()
    print(f"q={q}: finished in {time.time() - t0:.0f} s", flush=True)
    summarize(out, lo, hi)


def summarize(out: str, lo=None, hi=None) -> None:
    rows = parse_results(out)
    ok = sorted((a, b, h, l) for a, b, h, l, good in rows if good)
    bad = [(a, b) for a, b, _, _, good in rows if not good]
    # union of certified cells
    comps = []
    for a, b, h, l in ok:
        if comps and a <= comps[-1][1]:
            comps[-1][1] = max(comps[-1][1], b)
        else:
            comps.append([a, b])
    worst = max(ok, key=lambda r: r[2]) if ok else None
    lines = ["# SUMMARY " + time.strftime("%Y-%m-%d %H:%M:%S")]
    lines.append("# certified cells: %d; failed cells: %d" % (len(ok), len(bad)))
    lines.append("# union of certified cells: " +
                 " U ".join(f"[{dec(a)}, {dec(b)}]" for a, b in comps))
    if worst:
        lines.append(f"# worst (closest to 0) certified upper bound: delta' <= {e_up(worst[2])} "
                     f"on [{dec(worst[0])}, {dec(worst[1])}]")
    if lo is not None:
        full = len(comps) == 1 and comps[0][0] <= lo and comps[0][1] >= hi and not bad
        lines.append(f"# requested range [{dec(lo)}, {dec(hi)}] fully certified: {full}")
    with open(out, "a", encoding="ascii") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)


# --------------------------------------------------------------------------
# auxiliary certificates (a), (b)
# --------------------------------------------------------------------------
def Q_enclosure(q: int, M: Fraction) -> arb:
    """Enclosure of Q(M) = sum_g 4M^3/(M^2+g^2)^{3/2} (M >= 2.5).

    Head: all listed zeros.  Tail: the zeros beyond the list, all with
    g > g_N.  Put a_g := 4M^2/(M^2+g^2); the tail term is a_g^{3/2}/2 and
        A1 := sum_tail a_g   = V_T/2            (exact, variance identity),
        A2 := sum_tail a_g^2 = 2 V_T - M D_T    (exact: D_T = (2/M) V_T - A2/M),
        a_g <= amax := 4M^2/(M^2+g_N^2).
    Then
        A2 / (2 sqrt(amax))  <=  tail  <=  min( sqrt(amax) A1, sqrt(A1 A2) ) / 2,
    (a^2 <= sqrt(amax) a^{3/2};  a^{3/2} <= sqrt(amax) a;  Cauchy-Schwarz)."""
    ctx.prec = PREC_HI
    zeros = load_zeros(q)
    Mb = ex(M)
    VT, DT, _, _ = tail_VD(q, zeros, Mb)
    head = arb(0)
    for g in zeros:
        head += 4 * Mb ** 3 / ((Mb * Mb + g * g) ** 3).sqrt()
    A1 = VT / 2
    A2 = 2 * VT - Mb * DT
    if not A2.lower() > 0:
        raise ArithmeticError("tail fourth moment not certified positive")
    gN = zeros[-1]
    samax = (4 * Mb * Mb / (Mb * Mb + gN * gN)).sqrt()
    up1 = samax * A1 / 2
    up2 = (A1 * A2).sqrt() / 2
    tail_hi = up2 if up2.upper() < up1.upper() else up1
    tail_lo = A2 / (2 * samax)
    if not tail_lo.lower() <= tail_hi.upper():
        raise ArithmeticError("inconsistent tail bounds for Q")
    return head + hull(arb(tail_lo.lower()), arb(tail_hi.upper()))


def Lambda_enclosure(q: int, M: Fraction, Mp: Fraction) -> tuple[arb, arb, arb]:
    """Enclosure of Lambda(M,M') = 1/(M M') - sum_g 4/(|M+ig||M'+ig|), M < M'.

    Tail (zeros beyond the list): with W_M := sum_{g>T} 4/(M^2+g^2)
    = V_T(M)/(2M^2),  W_{M'} <= tail <= sqrt(W_M W_{M'})  (termwise
    monotonicity in M for the lower bound, Cauchy-Schwarz for the upper)."""
    ctx.prec = PREC_HI
    if not M < Mp:
        raise ValueError("need M < M'")
    zeros = load_zeros(q)
    a, b = ex(M), ex(Mp)
    head = arb(0)
    for g in zeros:
        g2 = g * g
        head += 4 / ((a * a + g2) * (b * b + g2)).sqrt()
    W = []
    sig = []
    for x in (a, b):
        s2 = sigma2_any(q, x)
        hv, _ = head_sums(zeros, x)
        VT = s2 - hv
        if not VT.lower() > 0:
            raise ArithmeticError("tail variance not positive")
        W.append(VT / (2 * x * x))
        sig.append(s2)
    tail = hull(arb(W[1].lower()), arb((W[0] * W[1]).sqrt().upper()))
    lam = 1 / (a * b) - head - tail
    need_finite(lam, "Lambda")
    return lam, sig[0], sig[1]


def bisect_threshold(f, lo: Fraction, hi: Fraction, digits: int) -> tuple[Fraction, Fraction]:
    """f(M) returns +1 if certified 'below', -1 if certified 'above', 0 if
    undecided; f is monotone.  Returns (largest grid M certified below,
    smallest grid M certified above) on the grid 10^-digits."""
    scale = 10 ** digits
    a, b = int(lo * scale), int(hi * scale)
    if f(Fraction(a, scale)) != 1 or f(Fraction(b, scale)) != -1:
        raise ArithmeticError("bisection bracket not certified")
    a2, b2 = a, b
    while b2 - a2 > 1:
        mid = (a2 + b2) // 2
        s = f(Fraction(mid, scale))
        if s == 1:
            a2 = mid
        else:
            b2 = mid
    below = a2
    a2, b2 = below, b
    while b2 - a2 > 1:
        mid = (a2 + b2) // 2
        s = f(Fraction(mid, scale))
        if s == -1:
            b2 = mid
        else:
            a2 = mid
    return Fraction(below, scale), Fraction(b2, scale)


def run_aux(out: str) -> None:
    ctx.prec = PREC_HI
    L = []
    P = L.append
    P("# aux_results.txt -- produced by certify_derivative.py aux")
    P("# All enclosures are Arb balls (outward rounded).  Hypotheses: GRH(chi) (for the")
    P("# variance identity) and the certified zero lists; see METHOD.md section 8.")
    P("")
    P("== (a)  Q_chi(M) = sum_{gamma>0} 4 M^3 / (M^2+gamma^2)^{3/2} ==")
    for q, Mreq, br in ((4, "3.68", ("3.5", "3.9")), (3, "4.47", ("4.3", "4.7"))):
        Qv = Q_enclosure(q, fr(Mreq))
        ok = bool(Qv.upper() < 1)
        P(f"q={q}: Q({Mreq}) in [{f_dn(dn(Qv))}, {f_up(up(Qv))}]   Q <= 1 certified: {ok}")
        if not ok:
            raise ArithmeticError("Q <= 1 not certified at the requested point")

        def side(M, q=q):
            v = Q_enclosure(q, M)
            if v.upper() < 1:
                return 1
            if v.lower() > 1:
                return -1
            return 0

        below, above = bisect_threshold(side, fr(br[0]), fr(br[1]), 6)
        vb, va = Q_enclosure(q, below), Q_enclosure(q, above)
        P(f"q={q}: largest 1e-6-grid M with Q(M) < 1 certified : M = {dec(below)}   "
          f"Q in [{f_dn(dn(vb))}, {f_up(up(vb))}]")
        P(f"q={q}: smallest 1e-6-grid M with Q(M) > 1 certified: M = {dec(above)}   "
          f"Q in [{f_dn(dn(va))}, {f_up(up(va))}]")
        P(f"q={q}: hence (Q strictly increasing)  {dec(below)} < M_0(chi_{q}) < {dec(above)}")
    P("")
    P("== (b1)  M* := sup{M : sigma_M^2 <= 2}  (sigma^2 strictly increasing in M) ==")
    for q, br in ((4, ("1.5", "2.0")), (3, ("1.8", "2.4"))):
        def side(M, q=q):
            v = sigma2_any(q, ex(M))
            if v.upper() <= 2:
                return 1
            if v.lower() > 2:
                return -1
            return 0

        below, above = bisect_threshold(side, fr(br[0]), fr(br[1]), 8)
        vb, va = sigma2_any(q, ex(below)), sigma2_any(q, ex(above))
        P(f"q={q}: sigma^2({dec(below)}) in [{f_dn(dn(vb), 15)}, {f_up(up(vb), 15)}]  <= 2 certified")
        P(f"q={q}: sigma^2({dec(above)}) in [{f_dn(dn(va), 15)}, {f_up(up(va), 15)}]  >  2 certified")
        P(f"q={q}: hence  {dec(below)} <= M* < {dec(above)}   (m* = M* - 1/2)")
    P("")
    P("== (b2)  Lambda(M,M') = 1/(M M') - sum_{gamma>0} 4/(|M+i gamma| |M'+i gamma|),  M=m+1/2 ==")
    P("# columns: q m m'  Lambda_lower  Lambda_upper  sign   sigma_M*sigma_M' (upper)   [Cauchy-Schwarz test: < 2]")
    for q in (4, 3):
        for m in range(4):
            for mp_ in range(m + 1, 4):
                M, Mp = Fraction(2 * m + 1, 2), Fraction(2 * mp_ + 1, 2)
                lam, s1, s2 = Lambda_enclosure(q, M, Mp)
                sgn = "+" if lam.lower() > 0 else ("-" if lam.upper() < 0 else "?")
                ss = (s1 * s2).sqrt()
                P(f"q={q} m={m} m'={mp_}  [{f_dn(dn(lam))}, {f_up(up(lam))}]  sign {sgn}   "
                  f"sigma*sigma' <= {f_up(up(ss), 9)}  ({'< 2' if ss.upper() < 2 else 'not < 2'})")
    P("")
    P("== sigma_m^2 at m = 0..3 (for reference) ==")
    for q in (4, 3):
        for m in range(4):
            v = sigma2_any(q, ex(Fraction(2 * m + 1, 2)))
            P(f"q={q} m={m}: sigma^2 in [{f_dn(dn(v), 15)}, {f_up(up(v), 15)}]")
    txt = "\n".join(L) + "\n"
    with open(out, "w", encoding="ascii") as fh:
        fh.write(txt)
    print(txt)


# --------------------------------------------------------------------------
# self-test of the closed forms (cross-checks; not part of the certificate)
# --------------------------------------------------------------------------
def selftest() -> None:
    ctx.prec = PREC_HI
    from mpmath import mp, mpf
    mp.dps = 40
    for q in (4, 3):
        chi = [0, 1, 0, -1] if q == 4 else [0, 1, -1]
        zeros = load_zeros(q)
        for Ms in ("3.6", "8.5", "17.4", "100.5"):
            M = ex(fr(Ms))
            s = M + arb(1) / 2
            d1, d2 = loglam_derivs_series(q, s)
            d1c = loglam_d1_complex_step(q, s)
            if not d1.overlaps(d1c):
                raise ArithmeticError("series vs complex step disagree")
            sm = mpf(Ms) + mpf(1) / 2

            def logLam(x):
                return (mp.log(mpf(q) / mp.pi) * (x + 1) / 2 + mp.loggamma((x + 1) / 2)
                        + mp.log(mp.dirichlet(x, chi)))

            r1 = mp.diff(logLam, sm, 1)
            r2 = mp.diff(logLam, sm, 2)
            e1 = abs(mpf(d1.mid().str(40, radius=False)) - r1)
            e2 = abs(mpf(d2.mid().str(40, radius=False)) - r2)
            # the balls must contain the independent mpmath values (up to
            # mpmath's own numerical-differentiation error, allowed 1e-25)
            if e1 > float(d1.rad()) + 1e-25 or e2 > float(d2.rad()) + 1e-25:
                raise ArithmeticError("closed form disagrees with mpmath")
            VT, DT, s2, D = tail_VD(q, zeros, M)
            print(f"q={q} M={Ms}: sigma^2={s2.str(18)} D={D.str(14)} VT={VT.str(10)} "
                  f"DT={DT.str(8)} |d1-mp|={float(e1):.1e} |d2-mp|={float(e2):.1e}")
    print("selftest OK")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cell")
    c.add_argument("--q", type=int, choices=[3, 4], required=True)
    c.add_argument("--Ma", required=True)
    c.add_argument("--Mb", required=True)
    c.add_argument("--nhead", type=int)
    c.add_argument("--tol", type=float, default=0.05)
    c.add_argument("--max-cells", type=int, default=6000)
    r = sub.add_parser("range")
    r.add_argument("--q", type=int, choices=[3, 4], required=True)
    r.add_argument("--lo", required=True)
    r.add_argument("--hi", required=True)
    r.add_argument("--ratio", type=float, default=1.01)
    r.add_argument("--out", required=True)
    r.add_argument("--workers", type=int, default=3)
    r.add_argument("--nhead", type=int)
    r.add_argument("--tol", type=float, default=0.05)
    r.add_argument("--max-cells", type=int, default=6000)
    s = sub.add_parser("summarize")
    s.add_argument("--out", required=True)
    s.add_argument("--lo")
    s.add_argument("--hi")
    a = sub.add_parser("aux")
    a.add_argument("--out", required=True)
    sub.add_parser("selftest")
    args = ap.parse_args()
    setup()
    if args.cmd == "cell":
        res = certify_cell(args.q, args.Ma, args.Mb, args.nhead, args.tol, args.max_cells)
        print(fmt_line(res))
        if not res["ok"]:
            sys.exit(1)
    elif args.cmd == "range":
        if args.workers > 3:
            ap.error("at most 3 workers")
        run_range(args.q, args.lo, args.hi, args.ratio, args.out, args.workers,
                  args.nhead, args.tol, args.max_cells)
    elif args.cmd == "summarize":
        summarize(args.out, fr(args.lo) if args.lo else None, fr(args.hi) if args.hi else None)
    elif args.cmd == "aux":
        run_aux(args.out)
    elif args.cmd == "selftest":
        selftest()


if __name__ == "__main__":
    main()
