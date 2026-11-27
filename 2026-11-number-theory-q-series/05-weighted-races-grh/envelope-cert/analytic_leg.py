#!/usr/bin/env python3
"""Arb evaluation of the explicit constants of Theorem E (zero-data-free analytic leg):

    for EVERY real primitive character chi (GRH(chi), L(1/2,chi) != 0) and every M >= M_E,
        - M delta'(M)  >=  G(sigma^2(M)) / (sqrt(2 pi) sigma(M))  >  0 .

The proof (ANALYTIC-LEG.md) reduces the statement to four inequalities between elementary
functions at the single point M_E; this script evaluates them in ball arithmetic, for the two
parity families  (a = 1, q >= 3)  and  (a = 0, q >= 5).  No zero of any L-function is used.

 (E2)  1/(2(s+1)) > Z2(s)  at s = M_E + 1/2          [=> h'(M) > 0, S_4 <= 16 sigma^2, all M >= M_E]
 (E3)  sig0^2 := min_families 4 M_E [ log(q_min/pi)/2 + psi((s+a)/2)/2 - Z1(s) ]   (lower bound of
       sigma_chi^2(M) for all chi, all M >= M_E)
 (E4)  N_chi(kappa M_E) >= (Y/2) [ (1 - 2Y/Gam) log(q_min/pi)/2 + psi((Y+1/2+a)/2)/2
                                   - (2Y/Gam) psi((Gam+1/2+a)/2)/2 - Z1(Y+1/2) - (2Y/Gam) Z1(Gam+1/2) ] >= 9,
       Gam = kappa M_E
 (E5)  G(sig0^2) > 0,
       G(x) = 1 - (1/2 + 48 beta)/x - (1/2)(1 - 8/x)^{-3/2}
              - sqrt(2 pi x) [ rho^{x/8}/(pi x tau) + (9/(7 pi)) T_c rho^{(x-72)/8} ]
              - (x^{3/2}/sqrt(2 pi)) (4/3) T_c^3 rho^{(x-80)/8},
       rho = 0.446, tau = z_*/4, T_c = (8/pi)^2 sqrt(1+kappa^2)/4, beta = beta(z_*).
 Z1(s) = L 2^{-s} + 2^{1-s}(L/a + 1/a^2),  Z2(s) = L^2 2^{-s} + 2^{1-s}(L^2/a + 2L/a^2 + 2/a^3),
 L = log 2, a = s - 1   (bounds for sum_{n>=2} log n n^{-s}, sum_{n>=2} log^2 n n^{-s}; valid for
 s > 2/log 2 = 2.885..., the script asks for s >= 4, more than is needed).

PRINTING (REFEREE-5 R1): every number printed after ">=" is a lower endpoint rounded DOWN, every
number printed after "<=" an upper endpoint rounded UP (p_dn / p_up).  Decisions use the balls.

Usage: python3 analytic_leg.py            (prints the table for the stated parameter triples)
"""
import math
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal
from fractions import Fraction

from flint import arb, ctx

ctx.prec = 128
RHO = arb(446) / 1000
ZS = 2 * (-(RHO.log())).sqrt()
TAU = ZS / 4
L2 = arb(2).log()


def Z1(s):
    a = s - 1
    return L2 * arb(2) ** (-s) + 2 * arb(2) ** (-s) * (L2 / a + 1 / (a * a))


def Z2(s):
    a = s - 1
    return L2 * L2 * arb(2) ** (-s) + 2 * arb(2) ** (-s) * (L2 * L2 / a + 2 * L2 / (a * a) + 2 / a ** 3)


def beta_star():
    J = ZS.bessel_j(0)
    assert J.lower() > 0
    return (-(J.log()) - ZS * ZS / 4) / ZS ** 4


def G(x, kappa):
    Tc = (8 / arb.pi()) ** 2 * (1 + arb(kappa) ** 2).sqrt() / 4
    beta = beta_star()
    two_pi = 2 * arb.pi()
    e1 = (two_pi * x).sqrt() * (RHO ** (x / 8) / (arb.pi() * x * TAU)
                                + arb(9) / (7 * arb.pi()) * Tc * RHO ** ((x - 72) / 8))
    e3 = x ** arb("1.5") / two_pi.sqrt() * arb(4) / 3 * Tc ** 3 * RHO ** ((x - 80) / 8)
    return 1 - (arb(1) / 2 + 48 * beta) / x - (1 - 8 / x) ** arb("-1.5") / 2 - e1 - e3, Tc, beta


def p_dn(x: arb, digits: int) -> str:
    """Decimal string <= the lower endpoint of the ball x (directed rounding: display only).
    The double below the endpoint is converted to Decimal exactly, then rounded down."""
    d = Decimal(math.nextafter(float(x.lower()), -math.inf))
    return str(d.quantize(Decimal(1).scaleb(-digits), rounding=ROUND_FLOOR))


def p_up(x: arb, digits: int) -> str:
    """Decimal string >= the upper endpoint of the ball x (directed rounding: display only)."""
    d = Decimal(math.nextafter(float(x.upper()), math.inf))
    return str(d.quantize(Decimal(1).scaleb(-digits), rounding=ROUND_CEILING))


FAMILIES = ((1, 3), (0, 5))        # (parity a, smallest conductor of that parity)


def check(ME: Fraction, kappa: int, Y: Fraction, verbose=True):
    M = arb(ME.numerator) / ME.denominator
    s = M + arb(1) / 2
    Yb = arb(Y.numerator) / Y.denominator
    Gam = kappa * M
    ok = True
    # hypotheses on the parameters
    ok &= bool(s >= 4) and bool(Yb + arb(1) / 2 >= 4) and bool(2 * Yb <= Gam) and bool(Yb <= Gam)
    # (E2)
    e2 = 1 / (2 * (s + 1)) - Z2(s)
    ok &= bool(e2.lower() > 0)
    # (E3), (E4)
    sig0 = None
    cnt = None
    for a, qmin in FAMILIES:
        lq = (arb(qmin) / arb.pi()).log() / 2
        hl = lq + ((s + a) / 2).digamma() / 2 - Z1(s)
        sg = 4 * M * hl
        sig0 = sg if sig0 is None else sig0.min(sg)
        w = 2 * Yb / Gam
        c = Yb / 2 * ((1 - w) * lq + ((Yb + arb(1) / 2 + a) / 2).digamma() / 2
                      - w * ((Gam + arb(1) / 2 + a) / 2).digamma() / 2
                      - Z1(Yb + arb(1) / 2) - w * Z1(Gam + arb(1) / 2))
        cnt = c if cnt is None else cnt.min(c)
    x0 = sig0.lower()
    ok &= bool(x0 > 80) and bool(cnt.lower() >= 9)
    g, Tc, beta = G(x0, kappa)
    ok &= bool(g.lower() > 0)
    if verbose:
        # lower endpoints rounded DOWN, upper endpoints rounded UP (REFEREE-5 R1)
        print(f"M_E={float(ME):g} kappa={kappa} Y={float(Y):g}:  h' margin (E2) >= {p_dn(e2, 7)};  "
              f"sig0^2 >= {p_dn(x0, 5)};  N(kappa M_E) >= {p_dn(cnt, 4)};  T_c <= {p_up(Tc, 6)};  "
              f"beta(z_*) <= {p_up(beta, 7)};  G(sig0^2) >= {p_dn(g, 7)}   "
              f"{'ALL CERTIFIED' if ok else 'NOT CERTIFIED'}")
    return ok


if __name__ == "__main__":
    print("# Theorem E constants (Arb, prec 128).  z_* =", ZS.str(12), " tau = z_*/4 =", TAU.str(12))
    for ME, kappa, Y in ((300, 4, 100), (100, 4, 40), (60, 4, 48), (50, 4, 40), (45, 5, 40), (40, 5, 40)):
        check(Fraction(ME), kappa, Fraction(Y))
    # smallest M_E on the grid 0.5 for which some (kappa, Y) passes
    best = None
    for twoM in range(100, 60, -1):
        ME = Fraction(twoM, 2)
        found = None
        for kappa in (3, 4, 5, 6, 7):
            for Yn in range(16, 80, 2):
                Y = Fraction(Yn)
                if 2 * Y > kappa * ME:
                    continue
                if check(ME, kappa, Y, verbose=False):
                    found = (kappa, Y)
                    break
            if found:
                break
        if not found:
            break
        best = (ME, found)
    if best:
        print("# smallest M_E on the half-integer grid with a passing (kappa, Y):")
        check(best[0], best[1][0], best[1][1])
