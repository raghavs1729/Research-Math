#!/usr/bin/env python3
"""
verify_A.py -- numerical verification of every identity listed in SURVEY-A.md
(Ramanujan-type hypergeometric series for 1/pi^2, 1/pi^3, 1/pi^4 and relatives).

Self-contained: needs only mpmath and sympy.  Run:  python3 verify_A.py > verify_A.out.txt
Working precision mp.dps = 80 (110 for the Taylor-coefficient checks).

Sections
  O  : identities that are OPEN / status unclear as far as the survey determined
  P  : identities that are PROVED in the literature (do not re-prove)
  SC : supercongruences (exact rational arithmetic, primes 7..43)
  T  : Taylor-expansion ("extended series") conjectures of Cohen-Guillera, via Cauchy integral
  W  : WZ pairs (F,G) transcribed in the survey; checks F(n+1,k)-F(n,k) = G(n,k+1)-G(n,k)

Notation: (a)_n Pochhammer.  hsum(z, num, den, P, n0) = sum_{n>=n0} z^n prod(num)_n/prod(den)_n P(n).
"""
import sys, time
from fractions import Fraction as Fr
from mpmath import (mp, mpf, mpc, mpmathify, pi, sqrt, log, asin, zeta, gamma, exp, sin, cos,
                    quad, linspace, harmonic, factorial, nstr)
from sympy import bernoulli as sbern, Rational, primerange
from sympy.functions.combinatorial.numbers import kronecker_symbol

mp.dps = 80
H = mpf(1) / 2


def Q(*a):
    return [mpf(p) / q for (p, q) in a]


def hsum(z, num, den, P, n0=0, tol=None):
    """sum_{n>=n0} z^n prod (a)_n / prod (b)_n * P(n), direct summation (|z|<1)."""
    if tol is None:
        tol = mpf(10) ** (-mp.dps - 15)
    z = mpmathify(z)
    t = mpf(1)
    s = mpf(0)
    n = 0
    while True:
        if n >= n0:
            term = t * P(mpf(n))
            s += term
            if n > n0 + 5 and abs(term) < tol and abs(t) < tol:
                break
        r = z
        for a in num:
            r *= (a + n)
        for b in den:
            r /= (b + n)
        t *= r
        n += 1
        if n > 10 ** 6:
            raise RuntimeError("series too slow")
    return s


def Sshift(a, ss, P, x):
    """Cohen-Guillera shifted sum S(H,a,P;x) = sum_{n>=0} P(n+x) H_{n+x} / (sign(a)^n |a|^(n+x)),
    H_{n+x} = prod_s Gamma(s+n+x)/Gamma(s) / Gamma(1+n+x)^d   (d = number of s)."""
    a = mpmathify(a)
    x = mpmathify(x)
    d = len(ss)
    pref = mpf(1)
    for s in ss:
        pref *= gamma(s + x) / gamma(s)
    pref /= gamma(1 + x) ** d
    pref /= abs(a) ** x
    return pref * hsum(1 / a, [s + x for s in ss], [1 + x] * d, lambda n: P(n + x))


def barnes(z, ss, P, T=62):
    """Analytic continuation to z<0, |z|>1 of sum_{n>=0} prod (s_i)_n/(1)_n^d P(n) z^n
    = (1/2 pi i) int_{c-i oo}^{c+i oo} prod (s_i)_s/(1)_s^d P(s) Gamma(1+s)Gamma(-s) (-z)^s ds,
    -min(s_i) < c < 0  (Guillera, arXiv:1012.2681, eq. (5))."""
    d = len(ss)
    c = -min(ss) / 2
    lg = log(-mpmathify(z))
    gs = [gamma(s) for s in ss]

    def f(t):
        s = mpc(c, t)
        v = mpf(1)
        for si, g in zip(ss, gs):
            v *= gamma(si + s) / g
        v /= gamma(1 + s) ** d
        return v * P(s) * (-pi / sin(pi * s)) * exp(s * lg)

    return quad(f, linspace(-T, T, 2 * T + 1)) / (2 * pi)


def kron(D, n):
    return int(kronecker_symbol(D, n))


def LD(D, s):
    """Dirichlet L(chi_D, s) for a positive fundamental discriminant D (D=1: zeta)."""
    if D == 1:
        return zeta(s)
    return sum(kron(D, a) * zeta(s, mpf(a) / D) for a in range(1, D)) / mpf(D) ** s


RESULTS = []


def report(ident, lhs, rhs, note=""):
    err = abs(lhs - rhs)
    scale = max(abs(rhs), mpf(10) ** -30)
    digits = int(-mp.log10(err / scale)) if err != 0 else mp.dps
    RESULTS.append((ident, digits))
    print("%-10s |lhs-rhs| = %-12s  digits ~ %3d   %s" % (ident, nstr(err, 5), digits, note))
    sys.stdout.flush()


# --------------------------------------------------------------------------------------
# parameter sets
s_2_3_4 = Q((1, 2), (1, 3), (2, 3), (1, 4), (3, 4))
s_2_4_6 = Q((1, 2), (1, 4), (3, 4), (1, 6), (5, 6))
s_2_3_6 = Q((1, 2), (1, 3), (2, 3), (1, 6), (5, 6))
s_2_8 = Q((1, 2), (1, 8), (3, 8), (5, 8), (7, 8))
s_2_5 = Q((1, 2), (1, 5), (2, 5), (3, 5), (4, 5))
s_222_3 = [H, H, H] + Q((1, 3), (2, 3))
s_222_4 = [H, H, H] + Q((1, 4), (3, 4))
one5 = [1] * 5
sq5 = sqrt(5)
phi5 = ((sq5 - 1) / 2) ** 5          # = (5 sqrt5 - 11)/2 = golden-ratio^(-5)
Phi5 = ((sq5 + 1) / 2) ** 5          # golden ratio^5 = 1/phi5

print("=" * 100)
print("SECTION O: open / status-unclear identities (mp.dps = %d)" % mp.dps)
print("=" * 100)

# O1-O5: Guillera / Almkvist-Guillera rational 1/pi^2 series still unproved
report("O1-G252", hsum(mpf(-1) / 48, s_2_3_4, one5, lambda n: 252 * n * n + 63 * n + 5), 48 / pi ** 2,
       "sum (-1/48)^n (1/2)(1/3)(2/3)(1/4)(3/4)/(1)^5 (252n^2+63n+5) = 48/pi^2")
report("O2-G1640", hsum(mpf(-1) / 1024, s_2_4_6, one5, lambda n: 1640 * n * n + 278 * n + 15),
       256 / (sqrt(3) * pi ** 2), "... (-1/2^10)^n (1/2)(1/4)(3/4)(1/6)(5/6) (1640n^2+278n+15) = 256/(sqrt3 pi^2)")
report("O3-G5418", hsum(mpf(-1) / 80 ** 3, s_2_3_6, one5, lambda n: 5418 * n * n + 693 * n + 29),
       128 * sq5 / pi ** 2, "... (-1/80^3)^n (1/2)(1/3)(2/3)(1/6)(5/6) (5418n^2+693n+29) = 128 sqrt5/pi^2")
report("O4-G1920", hsum(mpf(1) / 7 ** 4, s_2_8, one5, lambda n: 1920 * n * n + 304 * n + 15),
       56 * sqrt(7) / pi ** 2, "... (1/7^4)^n (1/2)(1/8)(3/8)(5/8)(7/8) (1920n^2+304n+15) = 56 sqrt7/pi^2")
report("O5-AG532", hsum((mpf(3) / 5) ** 6, s_2_3_6, one5, lambda n: 532 * n * n + 126 * n + 9),
       375 / (4 * pi ** 2), "... (3/5)^(6n) (1/2)(1/3)(2/3)(1/6)(5/6) (532n^2+126n+9) = 375/(4 pi^2)")

# O6: golden-ratio (irrational z) 1/pi^2 series
report("O6-GOLD", hsum((3 * phi5) ** 3, s_222_3, one5,
                       lambda n: (32 - 216 * phi5) * n * n + (18 - 162 * phi5) * n + (3 - 30 * phi5)),
       3 / pi ** 2, "z=(3 phi5)^3, phi5=((sqrt5-1)/2)^5")
# same series in Zudilin's binomial normalisation (arXiv:2602.06679 eq (8))
report("O6b-GOLD", hsum(27 * phi5 ** 3, s_222_3, one5,
                        lambda n: 20 * (61 - 27 * sq5) * n * n + 9 * (101 - 45 * sq5) * n + 3 * (56 - 25 * sq5)),
       3 / pi ** 2, "Zudilin 2026 form (same identity)")

# O8: upside-down companion of the golden 'divergent' series
L53 = LD(5, 3)
report("O8-GOLDud", hsum((-phi5 / 3) ** 3, one5, s_222_3,
                         lambda n: ((2408 + 216 * phi5) * n * n - (1800 + 162 * phi5) * n + (333 + 30 * phi5)) / n ** 5,
                         n0=1),
       mpf(1125) / 4 * sq5 * L53 - 448 * zeta(3), "= (1125/4) sqrt5 L_5(3) - 448 zeta(3)")

# O9: convergent 1/pi companion of the golden divergent series (Guillera 1610.04839 eq (18))
report("O9-GOLD1pi", hsum((-1 / (3 * Phi5)) ** 3, [H] * 5, [1, 1, 1, mpf(1) / 6, mpf(5) / 6],
                          lambda n: ((2408 + 216 * phi5) * n * n + (608 + 54 * phi5) * n + (35 + 3 * phi5)) / (6 * n + 1)),
       3 * Phi5 ** 1.5 / pi, "= 3 Phi^(3/2)/pi, Phi = golden ratio^5")

# O11, O12: Guillera's Ramanujan-Orr type 1/pi relatives (arXiv:1501.06413 (21)-(22) and last display)
report("O11-RO10", hsum(mpf(1) / 64, Q((1, 10), (3, 10), (7, 10), (9, 10)), [H, 1, 1, 1],
                        lambda n: (2100 * n * n + 1160 * n + 63) / (2 * n + 1)), 200 / pi,
       "sum 2^(-6n) (1/10)(3/10)(7/10)(9/10)/((1/2)(1)^3) (2100n^2+1160n+63)/(2n+1) = 200/pi")
report("O12-S35", hsum((mpf(3) / 5) ** 6, Q((1, 3), (2, 3), (1, 6), (5, 6)), [H, 1, 1, 1],
                       lambda n: (133 * n * n + 79 * n + 6) / (2 * n + 1)), 625 / (32 * pi),
       "sum (3/5)^(6n) (1/3)(2/3)(1/6)(5/6)/((1/2)(1)^3) (133n^2+79n+6)/(2n+1) = 625/(32 pi)")

# O13-O17: shifted (x = 1/2) evaluations, Cohen-Guillera table (sec. 6.3), for the unproved series
CG = {  # id: (s-list, a (=1/z), P, value of sum * pi^2, value of pi^3 S(1/2)/(sum*pi^2))
    5: (s_2_3_4, -48, lambda n: 252 * n * n + 63 * n + 5, mpf(48), 4 * log(mpf(2) ** 5 / 3 ** 3)),
    6: (s_2_3_6, -mpf(4) ** 6 / 3 ** 6, lambda n: 1930 * n * n + 549 * n + 45, mpf(384), 2 * log(2)),
    7: (s_2_3_6, -80 ** 3, lambda n: 5418 * n * n + 693 * n + 29, 128 * sq5, log(mpf(2) ** 74 * 5 ** 5 / mpf(3) ** 54)),
    8: (s_2_3_6, mpf(5) ** 6 / 3 ** 6, lambda n: 532 * n * n + 126 * n + 9, mpf(375) / 4, 2 * asin(mpf(164833) / 5 ** 8)),
    9: (s_2_4_6, -1024, lambda n: 1640 * n * n + 278 * n + 15, 256 / sqrt(3), log(mpf(3) ** 9 / 2 ** 14)),
    10: (s_2_8, 7 ** 4, lambda n: 1920 * n * n + 304 * n + 15, 56 * sqrt(7), 2 * asin(mpf(2241857) / 2 ** 25)),
    1: ([H] * 5, -4, lambda n: 20 * n * n + 8 * n + 1, mpf(8), 14 * zeta(3) / pi ** 2),
    2: ([H] * 5, -1024, lambda n: 820 * n * n + 180 * n + 13, mpf(128), 2 * zeta(3) / pi ** 2),
    3: (s_222_3, mpf(64) / 27, lambda n: 74 * n * n + 27 * n + 3, mpf(48), 2 * pi / 3),
    4: (s_222_4, 16, lambda n: 120 * n * n + 34 * n + 3, mpf(32), pi / 3),
}
for lab, i in (("O13", 5), ("O14", 7), ("O15", 8), ("O16", 9), ("O17", 10), ("O18", 6)):
    ss, a, P, tot, val = CG[i]
    report("%s-CG%dh" % (lab, i), Sshift(a, ss, P, H), tot * val / pi ** 3,
           "Cohen-Guillera #%d shifted by x=1/2: S(1/2) = (sum*pi^2) * V / pi^3" % i)
S12 = ([H] * 5 + Q((1, 3), (2, 3), (1, 4), (3, 4)), -mpf(256) / 27,
       lambda n: 4528 * n ** 4 + 3180 * n ** 3 + 972 * n * n + 147 * n + 9)
S13 = ([H] * 7 + Q((1, 4), (3, 4)), 4096, lambda n: 43680 * n ** 4 + 20632 * n ** 3 + 4340 * n * n + 466 * n + 21)
report("O19-CG12h", Sshift(S12[1], S12[0], S12[2], H), 9216 * zeta(3) / pi ** 7, "Zhao 1/pi^4 series shifted by 1/2 = 9216 zeta(3)/pi^7")
report("O20-CG13h", Sshift(S13[1], S13[0], S13[2], H), mpf(2048) / 15 / pi ** 4, "Cullen series shifted by 1/2 = (2048/15)/pi^4")

# O21: Guillera 2025 note 'On a formula for zeta(5) conjectured by Almkvist'
report("O21-Gz5", hsum(mpf(1) / 64, [1] * 7, [H] * 7,
                       lambda n: 4 * (84 * n ** 3 - 110 * n * n + 48 * n - 7) / n ** 8
                       + 56 * (7 * n * n - 5 * n + 1) * (3 * n - 1) * (harmonic(2 * n) - harmonic(n)) / n ** 7, n0=1),
       4 * pi ** 4 * log(2) + 186 * zeta(5), "= 4 pi^4 log2 + 186 zeta(5)")

# O7: golden 'divergent' series, by Barnes integral (analytic continuation)
t0 = time.time()
v = barnes(-27 * Phi5 ** 3, s_222_3,
           lambda n: (2408 + 216 * phi5) * n * n + (1800 + 162 * phi5) * n + (333 + 30 * phi5))
report("O7-GOLDdiv", v, 36 / pi ** 2, "analytic continuation at z=-27 Phi^3 (Barnes integral) '=' 36/pi^2  [%.0fs]" % (time.time() - t0))

print()
print("=" * 100)
print("SECTION P: proved identities (reference values; do not re-prove)")
print("=" * 100)
report("P1-G20", hsum(mpf(-1) / 4, [H] * 5, one5, lambda n: 20 * n * n + 8 * n + 1), 8 / pi ** 2, "(1/2)^5, z=-1/4: 8/pi^2")
report("P2-G820", hsum(mpf(-1) / 1024, [H] * 5, one5, lambda n: 820 * n * n + 180 * n + 13), 128 / pi ** 2, "(1/2)^5, z=-1/2^10: 128/pi^2")
report("P3-G120", hsum(mpf(1) / 16, s_222_4, one5, lambda n: 120 * n * n + 34 * n + 3), 32 / pi ** 2, "(1/2)^3(1/4)(3/4), z=1/16: 32/pi^2")
report("P4-G74", hsum(mpf(27) / 64, s_222_3, one5, lambda n: 74 * n * n + 27 * n + 3), 48 / pi ** 2, "(1/2)^3(1/3)(2/3), z=27/64: 48/pi^2")
report("P5-G1930", hsum(-mpf(3) ** 6 / 4 ** 6, s_2_3_6, one5, lambda n: 1930 * n * n + 549 * n + 45), 384 / pi ** 2,
       "(1/2)(1/3)(2/3)(1/6)(5/6), z=-(3/4)^6: 384/pi^2 (Au 2023)")
report("P6-GOUR", hsum(mpf(1) / 64, [H] * 7, [1] * 7, lambda n: 168 * n ** 3 + 76 * n * n + 14 * n + 1), 32 / pi ** 3,
       "Gourevitch: (1/2)^7, z=1/64: 32/pi^3 (Au 2023)")
report("P7-CULLEN", hsum(mpf(1) / 4096, [H] * 7 + Q((1, 4), (3, 4)), [1] * 9,
                         lambda n: 43680 * n ** 4 + 20632 * n ** 3 + 4340 * n * n + 466 * n + 21), 2048 / pi ** 4,
       "Cullen: (1/2)^7(1/4)(3/4), z=1/2^12: 2048/pi^4 (Au 2023)")
report("P8-ZHAO", hsum(-mpf(27) / 256, [H] * 5 + Q((1, 3), (2, 3), (1, 4), (3, 4)), [1] * 9,
                       lambda n: 4528 * n ** 4 + 3180 * n ** 3 + 972 * n * n + 147 * n + 9), 768 / pi ** 4,
       "Zhao: (1/2)^5(1/3)(2/3)(1/4)(3/4), z=-27/256: 768/pi^4 (Au 2023)")
report("P9-SUN12", hsum(mpf(4) / 27, [H] * 7, [mpf(1) / 6, mpf(5) / 6, 1, 1, 1, 1, 1],
                        lambda n: (92 * n ** 3 + 54 * n * n + 12 * n + 1) / (6 * n + 1)), 12 / pi ** 2,
       "Sun: (4/27)^n (1/2)^7/((1/6)(5/6)(1)^5) (92n^3+54n^2+12n+1)/(6n+1) = 12/pi^2 (Au 2023)")


def shalev():
    s = mpf(0)
    n = 0
    while True:
        t = ((-1) ** n * factorial(2 * n) ** 10 * factorial(5 * n) / (factorial(n) ** 15 * factorial(10 * n))
             * 3 * (5532 * n ** 4 + 5464 * n ** 3 + 2173 * n ** 2 + 416 * n + 32) / mpf((10 * n + 1) * (10 * n + 3)))
        s += t
        if abs(t) < mpf(10) ** -95:
            return s
        n += 1


report("P10-SHALEV", shalev(), 240 / pi ** 2,
       "sum (-1)^n (2n)!^10 (5n)!/(n!^15 (10n)!) 3(5532n^4+5464n^3+2173n^2+416n+32)/((10n+1)(10n+3)) = 240/pi^2")

# upside-down (dual) series, all proved
UD = [
    ("U14", mpf(-1) / 4, one5, [H] * 5, lambda n: (10 * n * n - 6 * n + 1) / n ** 5, -28 * zeta(3), "=-28 zeta(3)"),
    ("U15", mpf(-1) / 1024, one5, [H] * 5, lambda n: (205 * n * n - 160 * n + 32) / n ** 5, -2 * zeta(3), "=-2 zeta(3) (Amdeberhan-Zeilberger)"),
    ("U16", mpf(-1) / 27, one5, s_222_3, lambda n: (28 * n * n - 18 * n + 3) / n ** 5, -14 * zeta(3), "=-14 zeta(3)"),
    ("U17", mpf(-256) / 3125, one5, s_2_5, lambda n: (483 * n * n - 245 * n + 30) / n ** 5, -896 * zeta(3), "=-896 zeta(3) (Au Ex. IX)"),
    ("U18", mpf(-16) / 27, one5, s_2_3_4, lambda n: (172 * n * n - 75 * n + 9) / n ** 5, -1792 * zeta(3), "=-1792 zeta(3) (Chu-Zhang)"),
    ("U19", mpf(1) / 64, [1] * 7, [H] * 7, lambda n: (21 * n ** 3 - 22 * n * n + 8 * n - 1) / n ** 7, pi ** 4 / 8, "=pi^4/8 (Au Ex. XI)"),
    ("U20", mpf(4) / 27, [1] * 7, [H] * 5 + Q((1, 3), (2, 3)), lambda n: (92 * n ** 3 - 84 * n * n + 27 * n - 3) / n ** 7, 8 * pi ** 4, "=8 pi^4 (Au Ex. IV)"),
    ("U21", mpf(-1024) / 3125, [1] * 9, [H] * 5 + Q((1, 5), (2, 5), (3, 5), (4, 5)),
     lambda n: (5532 * n ** 4 - 5600 * n ** 3 + 2275 * n * n - 425 * n + 30) / n ** 9, -380928 * zeta(5), "=-380928 zeta(5) (Zhao; Au Ex. X)"),
    ("U22", mpf(16) / 27, [1] * 3, [H] + Q((1, 3), (2, 3)), lambda n: (11 * n - 3) / n ** 3, 8 * pi ** 2, "=8 pi^2 (Guillera 1012.2681 (14))"),
]
for (lab, z, nu, de, P, rhs, note) in UD:
    report(lab, hsum(z, nu, de, P, n0=1), rhs, "upside-down " + note)
report("U18b-CZ42", hsum(mpf(-16) / 27, one5, Q((3, 2), (4, 3), (5, 3), (5, 4), (7, 4)), lambda n: 172 * n * n + 269 * n + 106),
       63 * zeta(3), "Chu-Zhang Ex. 42 form: sum (-16/27)^n (1)^5/((3/2)(4/3)(5/3)(5/4)(7/4)) (172n^2+269n+106) = 63 zeta(3)")
report("U16b-Chu14", hsum(mpf(-1) / 27, one5, Q((3, 2), (3, 2), (3, 2), (4, 3), (5, 3)), lambda n: 28 * n * n + 38 * n + 13),
       mpf(21) / 2 * zeta(3), "Chu Ex. 14 form: = (21/2) zeta(3)")
report("P11-G18", hsum(mpf(27) / 64, [1, 1, 1, mpf(1) / 6, mpf(5) / 6], [H] * 5,
                       lambda n: (74 * n * n + 101 * n + 35) * (6 * n + 1) / (2 * n + 1) ** 5), 16 * pi ** 2 / 3,
       "Guillera 1003.1915 (18): sum (27/64)^n (1)^3(1/6)(5/6)/(1/2)^5 (74n^2+101n+35)(6n+1)/(2n+1)^5 = 16 pi^2/3")
report("P17-RO7", hsum(mpf(1) / 7 ** 4, Q((1, 8), (3, 8), (5, 8), (7, 8)), [H, 1, 1, 1],
                      lambda n: (1920 * n * n + 1072 * n + 55) / (2 * n + 1)), 196 * sqrt(7) / (3 * pi),
       "sum 7^(-4n) (1/8)(3/8)(5/8)(7/8)/((1/2)(1)^3) (1920n^2+1072n+55)/(2n+1) = 196 sqrt7/(3 pi) (Guillera NZJM 2017 per DPVZ)")
# P18: Guillera 1610.04839 eq (17) = Sun conjecture; proved by Hou-Sun, arXiv:2604.15172 (2026)
report("P18-HS", hsum(mpf(-1) / 27, [H] * 5, [1, 1, 1, mpf(1) / 6, mpf(5) / 6],
                       lambda n: (28 * n * n + 10 * n + 1) / (6 * n + 1)), 3 / pi,
       "sum (-1/27)^n (1/2)^5/((1)^3(1/6)(5/6)) (28n^2+10n+1)/(6n+1) = 3/pi")

# P19: Guillera's k-generalisation of Gourevitch's formula (arXiv:2503.00570 sec. 4.5), tested at k = 3/10
def gour_gen(k):
    k = mpmathify(k)
    pk = lambda a: gamma(a + k) / gamma(a)
    pref = pk(H) * pk(mpf(1) / 6) ** 2 * pk(mpf(5) / 6) ** 2 / pk(mpf(1)) ** 5 * (mpf(729) / 1024) ** k
    return pref * hsum(mpf(1) / 64, [H + 3 * k] * 2 + [H - k] * 2 + [H + k] * 3, [1, 1] + [1 + 2 * k] * 2 + [1 + k] * 3,
                       lambda n: (6 * n + 6 * k + 1) * (28 * n * n + 8 * n + 1 + 12 * k * k + 56 * n * k + 8 * k))


report("P19-GGOUR", gour_gen(mpf(3) / 10), 32 / pi ** 3,
       "sum_n (1/2+3k)_n^2(1/2-k)_n^2(1/2+k)_n^3/((1)_n^2(1+2k)_n^2(1+k)_n^3) (1/2)_k(1/6)_k^2(5/6)_k^2/(1)_k^5 64^-n (729/1024)^k (6n+6k+1)(28n^2+8n+1+12k^2+56nk+8k) = 32/pi^3 at k=0.3")
for lab, i in (("P12", 1), ("P13", 2), ("P14", 3), ("P15", 4)):
    ss, a, P, tot, val = CG[i]
    report("%s-CG%dh" % (lab, i), Sshift(a, ss, P, H), tot * val / pi ** 3, "Cohen-Guillera #%d shifted by x=1/2" % i)
report("P16-CG11h", Sshift(64, [H] * 7, lambda n: 168 * n ** 3 + 76 * n * n + 14 * n + 1, H), 8 / pi ** 3,
       "Gourevitch shifted by 1/2 = 8/pi^3 (equivalent to U19)")

# 'divergent' series (analytic continuation by Barnes integral), all proved
DIV = [
    ("D14", -4, [H] * 5, lambda n: 10 * n * n + 6 * n + 1, 4 / pi ** 2, "(1/2)^5, z=-4, 10n^2+6n+1 '=' 4/pi^2"),
    ("D15", -1024, [H] * 5, lambda n: 205 * n * n + 160 * n + 32, 16 / pi ** 2, "(1/2)^5, z=-2^10, 205n^2+160n+32 '=' 16/pi^2"),
    ("D16", -27, s_222_3, lambda n: 28 * n * n + 18 * n + 3, 6 / pi ** 2, "(1/2)^3(1/3)(2/3), z=-27 '=' 6/pi^2"),
    ("D17", -mpf(3125) / 256, s_2_5, lambda n: 483 * n * n + 245 * n + 30, 80 / pi ** 2, "(1/2)(k/5), z=-5^5/2^8 '=' 80/pi^2"),
    ("D18", -mpf(27) / 16, s_2_3_4, lambda n: 172 * n * n + 75 * n + 9, 48 / pi ** 2, "(1/2)(1/3)(2/3)(1/4)(3/4), z=-27/16 '=' 48/pi^2"),
    ("D21", -mpf(3125) / 1024, [H] * 5 + Q((1, 5), (2, 5), (3, 5), (4, 5)),
     lambda n: 5532 * n ** 4 + 5600 * n ** 3 + 2275 * n * n + 425 * n + 30, 1280 / pi ** 4, "(1/2)^5(k/5), z=-5^5/2^10 '=' 1280/pi^4"),
]
for (lab, z, ss, P, rhs, note) in DIV:
    t0 = time.time()
    report(lab, barnes(z, ss, P), rhs, "Barnes integral: " + note + " [%.0fs]" % (time.time() - t0))

print()
print("=" * 100)
print("SECTION SC: supercongruences for the OPEN 1/pi^2 series (exact arithmetic)")
print("  S_p = sum_{n<p} P(n) prod(s)_n/(1)_n^5 z^n == P(0) (D/p) p^2                      (mod p^5)   [Zudilin type]")
print("  S_p == P(0) (D/p) p^2 + C5 * L(chi_D, 4-p) * p^5                                   (mod p^6)   [Zhao type; C5 from Cohen-Guillera]")
print("=" * 100)


def Lneg(D, m):
    """L(chi_D, 1-m) = -B_{m,chi}/m (exact rational); D positive fundamental discriminant."""
    if D == 1:
        return -Rational(sbern(m)) / m
    tot = sum(kron(D, a) * sbern(m, Rational(a, D)) for a in range(1, D + 1))
    return -Rational(D) ** (m - 1) * tot / m


def trunc(ss, a, P, p):
    t = Fr(1)
    s = Fr(0)
    d = len(ss)
    for n in range(p):
        s += t * P(n)
        for si in ss:
            t *= (si + n)
        t /= Fr(n + 1) ** d
        t /= a
    return s


def modp(fr, M):
    return fr.numerator * pow(fr.denominator, -1, M) % M


SC = {
    "SC1-G252": ([Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 4), Fr(3, 4)], Fr(-48), lambda n: 252 * n * n + 63 * n + 5, 1, Fr(-91, 9)),
    "SC2-G1640": ([Fr(1, 2), Fr(1, 4), Fr(3, 4), Fr(1, 6), Fr(5, 6)], Fr(-1024), lambda n: 1640 * n * n + 278 * n + 15, 12, Fr(-15, 2)),
    "SC3-G5418": ([Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 6), Fr(5, 6)], Fr(-80 ** 3), lambda n: 5418 * n * n + 693 * n + 29, 5, Fr(-35, 216)),
    "SC4-G1920": ([Fr(1, 2), Fr(1, 8), Fr(3, 8), Fr(5, 8), Fr(7, 8)], Fr(7 ** 4), lambda n: 1920 * n * n + 304 * n + 15, 28, Fr(35, 8)),
    "SC5-AG532": ([Fr(1, 2), Fr(1, 3), Fr(2, 3), Fr(1, 6), Fr(5, 6)], Fr(5 ** 6, 3 ** 6), lambda n: 532 * n * n + 126 * n + 9, 1, Fr(1953, 50)),
}
for name, (ss, a, P, D, C5) in SC.items():
    ok5, ok6, used = True, True, []
    for p in primerange(7, 45):
        if a.numerator % p == 0 or a.denominator % p == 0 or D % p == 0:
            continue
        used.append(p)
        S = trunc(ss, a, P, p)
        main = P(0) * kron(D, p) * p * p
        ok5 &= (modp(S, p ** 5) - main) % p ** 5 == 0
        L = Lneg(D, p - 3)
        corr = Fr(int(L.p), int(L.q)) * C5 * p ** 5
        ok6 &= (modp(S - corr, p ** 6) - main) % p ** 6 == 0
    print("%-10s D=%-2d C5=%-8s primes %s : mod p^5 %s ; mod p^6 (with L-term) %s"
          % (name, D, C5, used, "ALL HOLD" if ok5 else "FAIL", "ALL HOLD" if ok6 else "FAIL"))
    RESULTS.append((name, 99 if (ok5 and ok6) else 0))
    sys.stdout.flush()

print()
print("=" * 100)
print("SECTION T: Cohen-Guillera Taylor expansions of the shifted OPEN series (mp.dps=110, Cauchy integral)")
print("  S(x) = a0*D*[ sqrt(D)/pi^2 + 0 x - a2 D sqrt(D) x^2 + 0 x^3 + a4 D^3 L(D,2) x^4 - a5 a0 D^4 L(D,3) x^5 + O(x^6) ]")
print("=" * 100)
mp.dps = 110


def taylor_coeffs(f, N, r, M=64):
    vals = [f(r * exp(2j * pi * k / M)) for k in range(M)]
    return [sum(vals[k] * exp(-2j * pi * k * m / M) for k in range(M)) / M / r ** m for m in range(N + 1)]


TAY = {  # id: (cg number, D, a0, a2, a4, a5)
    "T1-G252": (5, 1, mpf(48), mpf(3) / 2, mpf(157) / 4, mpf(91) / 9),
    "T2-G5418": (7, 5, mpf(128) / 5, mpf(3) / 2, mpf(887) / 32, mpf(21) / 8),
    "T3-AG532": (8, 1, mpf(375) / 4, mpf(4) / 3, mpf(40), mpf(6944) / 1125),
    "T4-G1640": (9, 12, mpf(32) / 9, mpf(7) / 24, mpf(757) / 1152, mpf(3) / 16),
    "T5-G1920": (10, 28, mpf(1), mpf(1) / 7, mpf(31) / 336, mpf(1) / 21),
}
PARS = {5: ([(1, 2), (1, 3), (2, 3), (1, 4), (3, 4)], -48, (252, 63, 5)),
        7: ([(1, 2), (1, 3), (2, 3), (1, 6), (5, 6)], -80 ** 3, (5418, 693, 29)),
        8: ([(1, 2), (1, 3), (2, 3), (1, 6), (5, 6)], Fr(5 ** 6, 3 ** 6), (532, 126, 9)),
        9: ([(1, 2), (1, 4), (3, 4), (1, 6), (5, 6)], -1024, (1640, 278, 15)),
        10: ([(1, 2), (1, 8), (3, 8), (5, 8), (7, 8)], 7 ** 4, (1920, 304, 15))}
for name, (i, D, a0, a2, a4, a5) in TAY.items():
    sl, a, (c2, c1, c0) = PARS[i]
    ss = [mpf(p) / q for (p, q) in sl]
    a = mpf(a.numerator) / a.denominator if isinstance(a, Fr) else mpf(a)
    c = taylor_coeffs(lambda x: Sshift(a, ss, lambda n: c2 * n * n + c1 * n + c0, x), 5, mpf(1) / 100)
    sD = sqrt(mpf(D))
    pred = [a0 * D * sD / pi ** 2, 0, -a0 * D * a2 * D * sD, 0, a0 * D * a4 * D ** 3 * LD(D, 2), -a0 * D * a5 * a0 * D ** 4 * LD(D, 3)]
    worst = min(int(-mp.log10(abs(c[m] - pred[m]) / max(abs(pred[m]), 1))) for m in range(6))
    print("%-10s coefficients x^0..x^5 agree; worst agreement ~ %d digits" % (name, worst))
    RESULTS.append((name, worst))
    sys.stdout.flush()
mp.dps = 80

print()
print("=" * 100)
print("SECTION W: transcribed WZ pairs; max relative defect of F(n+1,k)-F(n,k) - (G(n,k+1)-G(n,k)) at sample (n,k)")
print("=" * 100)
mp.dps = 65
f = lambda x: gamma(mpmathify(x) + 1)
po = lambda a, n: gamma(mpmathify(a) + n) / gamma(mpmathify(a))
bi = lambda a, b: f(a) / (f(b) * f(a - b))
sg = lambda n: cos(pi * n)
PTS = ((1, 0), (1, 1), (2, 1), (3, 2), (4, 1), (2, 0), (5, 3))


def wz(F, G, pts=PTS, sign=1):
    m = mpf(0)
    for (n, k) in pts:
        n = mpf(n)
        k = mpf(k)
        lhs = F(n + 1, k) - F(n, k)
        rhs = sign * (G(n, k + 1) - G(n, k))
        sc = max(abs(F(n + 1, k)), abs(F(n, k)), abs(G(n, k + 1)), abs(G(n, k)))
        m = max(m, abs(lhs - rhs) / sc)
    return m


def wrep(name, d, note):
    print("%-5s defect %-10s %s" % (name, nstr(d, 3), note))
    RESULTS.append((name, int(-mp.log10(d)) if d != 0 else 65))


# W1: Guillera, Ramanujan J. 15 (2008) = arXiv:1104.0396, Identity 8, pair for eq. (13)  [20n^2+8n+1]
B = lambda n, k: sg(n) / (mpf(2) ** (12 * n) * mpf(2) ** (8 * k)) * f(2 * k) ** 4 * f(2 * n) ** 5 / (f(n + k) ** 4 * f(k) ** 4 * f(n) ** 6)
wrep("W1", wz(lambda n, k: B(n, k) * 8 * n * (2 * n + 4 * k + 1),
              lambda n, k: B(n, k) * (20 * n * n + 8 * n + 1 + 24 * k * n + 8 * k * k + 4 * k)), "1104.0396 Id.8 (13): 20n^2+8n+1")
# W2: same paper, Identity 8, pair for eq. (14)
B = lambda n, k: sg(k) * sg(n) / (mpf(2) ** (12 * n) * mpf(2) ** (6 * k)) * f(2 * n) ** 4 * f(2 * n + 2 * k) * f(2 * k) ** 2 * f(n - k - H) / (
    f(n) ** 7 * f(n + k) ** 3 * f(k) ** 2 * f(n - H))
wrep("W2", wz(lambda n, k: B(n, k) * 32 * n ** 3 / (2 * n - 2 * k - 1),
              lambda n, k: B(n, k) * (20 * n * n + 12 * k * n + 8 * n + 2 * k + 1)), "1104.0396 Id.8 (14): 20n^2+8n+1")
# W3: Identity 9, pair for eq. (16)  [820n^2+180n+13]
B = lambda n, k: sg(n) / (mpf(2) ** (20 * n) * mpf(2) ** (8 * k)) * f(2 * n + 2 * k) ** 4 * f(2 * n) ** 5 / (
    f(2 * n + k) ** 4 * f(n + k) ** 4 * f(n) ** 6)
G3 = lambda n, k: B(n, k) * ((2 * n + 2 * k + 1) ** 4 / (2 * n + k + 1) ** 4 * (820 * n * n + 180 * n + 13 + 8 * k * k + 20 * k + 72 * n * k)
                             - (296 * n * k ** 3 + 1056 * n * n * k * k + 1280 * n ** 3 * k + 528 * n ** 4 + 800 * n ** 3 + 1344 * n * n * k
                                + 608 * n * k * k + 28 * k ** 3 + 408 * n * n + 384 * n * k + 40 * k * k + 72 * n + 16 * k + 1)
                             * 32 * n * k / (2 * n + k + 1) ** 4)
wrep("W3", wz(lambda n, k: B(n, k) * 128 * n * (6 * n + 4 * k + 1), G3), "1104.0396 Id.9 (16): 820n^2+180n+13")
# W4: Identity 9, pair for eq. (17)
B = lambda n, k: sg(k) * sg(n) / (mpf(2) ** (20 * n) * mpf(2) ** (6 * k)) * f(2 * n) ** 4 * f(2 * n + 2 * k) ** 3 * f(n - k - H) / (
    f(n) ** 7 * f(n + k) ** 3 * f(2 * n + k) ** 2 * f(n - H))
G4 = lambda n, k: B(n, k) * (820 * n * n + 180 * n + 13 + k / (2 * n + k + 1) ** 2 * (
    1312 * n ** 3 + 1340 * n * n * k + 336 * n * k * k + 1456 * n * n + 828 * n * k + 40 * k * k + 472 * n + 79 * k + 36))
wrep("W4", wz(lambda n, k: B(n, k) * 2048 * n ** 3 / (2 * n - 2 * k - 1), G4), "1104.0396 Id.9 (17): 820n^2+180n+13")
# W5: Identity 10, pair 1  [120n^2+34n+3]
B = lambda n, k: 1 / (mpf(2) ** (16 * n) * mpf(2) ** (8 * k)) * f(2 * k) ** 2 * f(2 * n + 2 * k) ** 2 * f(4 * n) * f(2 * n) ** 2 / (
    f(2 * n + k) ** 2 * f(n + k) ** 4 * f(k) ** 2 * f(n) ** 4)
G5 = lambda n, k: B(n, k) * (120 * n * n + 34 * n + 3 + 32 * k * k + 128 * k * n + 16 * k
                             + k * (32 * n ** 3 + 8 * n * n * k + 16 * n + 6 * k * n + 40 * n * n + k + 2) / (2 * n + k + 1) ** 2)
# NOTE: the pdftotext layout of the paper shows extra denominators (4n+4k+1) in G; that reading FAILS the WZ equation.
# The form above (no (4n+4k+1) denominators) was determined numerically from F and G(n,0) and satisfies it exactly.
wrep("W5", wz(lambda n, k: B(n, k) * 32 * n * (4 * n + 4 * k + 1), G5), "1104.0396 Id.10 pair 1 (numerically validated reading): 120n^2+34n+3")
# W6: Identity 10, pair 2
B = lambda n, k: sg(k) / (mpf(2) ** (16 * n) * mpf(2) ** (6 * k)) * f(2 * n) ** 2 * f(4 * n) * f(2 * k) ** 3 * f(2 * n - k - H) / (
    f(n) ** 6 * f(n + k) ** 2 * f(k) ** 3 * f(2 * n - H))
wrep("W6", wz(lambda n, k: B(n, k) * 512 * n ** 3 / (4 * n - 2 * k - 1),
              lambda n, k: B(n, k) * (120 * n * n + 84 * k * n + 34 * n + 10 * k + 3)), "1104.0396 Id.10 pair 2: 120n^2+34n+3")
# W7: Guillera, Ramanujan J. 26 (2011) = arXiv:1003.1915, eqs (10),(11)
U = lambda n, k: po(H, n) ** 3 * po(1 + k / 2, n) * po(H + k / 2, n) / (po(1, n) ** 3 * po(1 + k, n) ** 2) * po(H, k) ** 2 / po(1, k) ** 2
P5 = ((1, 1), (2, 1), (3, 2), (1, 2), (2, 3))
wrep("W7", wz(lambda n, k: U(n, k) * 96 * n ** 3 / (2 * n + k),
              lambda n, k: U(n, k) * 12 * k * (8 * n * n + 6 * k * n + 2 * n + k) / (2 * n + k), pts=P5), "1003.1915 (10)-(11)")
# W8: same paper eqs (15),(16)  [74n^2+27n+3]
U = lambda n, k: (mpf(27) / 64) ** n * po(H, n) ** 3 * po(H + k, n) ** 2 * po(1 + k / 3, n) * po((1 + k) / 3, n) * po((2 + k) / 3, n) / (
    po(1, n) ** 3 * po(1 + k, n) * po(1 + k / 2, n) ** 2 * po(H + k / 2, n) ** 2) * po(H, k) ** 2 / po(1, k) ** 2
Pp = lambda n, k: ((2 * n + 1) * (296 * n ** 3 + 164 * n * n + 26 * n + 1) + (480 * n ** 3 + 360 * n * n + 78 * n + 5) * k
                   + (176 * n * n + 80 * n + 8) * k * k + (24 * n + 4) * k ** 3)
wrep("W8", wz(lambda n, k: U(n, k) * 96 * n ** 3 / (3 * n + k),
              lambda n, k: U(n, k) * (n * (2 * n + 1) ** 2 * (74 * n * n + 27 * n + 3) + k * Pp(n, k)) / ((n + k / 3) * (2 * n + k + 1) ** 2),
              pts=P5 + ((1, 0), (2, 0))), "1003.1915 (15)-(16): 74n^2+27n+3")
# W9: Guillera, Adv. Appl. Math. 29 (2002) = arXiv:math/0503345, first pair (1/pi, 20n+3 <-> 42n+5)
C = lambda n, k: sg(n) * sg(k) / (mpf(2) ** (10 * n) * mpf(2) ** (2 * k)) * bi(2 * k, k) ** 2 * bi(2 * n, n) ** 2 * bi(4 * n - 2 * k, 2 * n - k) / (
    bi(2 * n, k) * bi(n + k, n))
P6 = ((1, 0), (1, 1), (2, 1), (3, 2), (4, 1), (2, 0))
wrep("W9", wz(lambda n, k: 64 * C(n, k) * n ** 2 / (4 * n - 2 * k - 1), lambda n, k: C(n, k) * (20 * n + 2 * k + 3), pts=P6),
     "math/0503345 pair 1")
# W10: same paper, second pair (120n^2+34n+3 <-> 820n^2+180n+13)
C = lambda n, k: sg(k) / (mpf(2) ** (16 * n) * mpf(2) ** (4 * k)) * bi(2 * k, k) ** 3 * bi(2 * n, n) ** 4 * bi(4 * n - 2 * k, 2 * n - k) / (
    bi(2 * n, k) * bi(n + k, n) ** 2)
wrep("W10", wz(lambda n, k: 512 * C(n, k) * n ** 3 / (4 * n - 2 * k - 1),
               lambda n, k: C(n, k) * (120 * n * n + 84 * n * k + 34 * n + 10 * k + 3), pts=P6), "math/0503345 pair 2")
# W11: Guillera, arXiv:1012.2681 Example 1, eq. (12).  NOTE sign: holds with F -> -F relative to the printed F.
U = lambda n, k: f(2 * n) ** 2 * f(2 * n + k) * f(2 * k) / (f(n) ** 4 * f(k) * f(n + k) ** 2) / (mpf(16) ** n * mpf(4) ** k)
wrep("W11", wz(lambda n, k: U(n, k) * 2 * n * n / (2 * n + k), lambda n, k: U(n, k) * (6 * n * n + 2 * n + k + 4 * n * k) / (2 * n + k),
               pts=P5, sign=-1), "1012.2681 Ex.1 (12): holds as F(n+1,k)-F(n,k) = -(G(n,k+1)-G(n,k)) for the printed F,G")
# W12: Guillera, arXiv:2503.00570 sec. 4.4.  NOTE: needs z=16/27 (paper prints 4/27) and the same opposite sign.
z = mpf(16) / 27
U = lambda n, k: po(1, n) * po(1 + k, n) * po(1 + 2 * k, n) * po(1 + k / 2, n) * po(H + k / 2, n) / (
    po(H, n) * po(H + k, n) * po(1 + 2 * k / 3, n) * po((1 + 2 * k) / 3, n) * po((2 + 2 * k) / 3, n)) * po(1, k) ** 2 / po(H, k) ** 2
S = lambda n, k: n * (2 * n + 2 * k + mpf(3) / 2) / ((1 + 2 * k) * (2 * n + 2 * k + 1) * (3 * n + 2 * k + 1) * (3 * n + 2 * k + 2))
Rr = lambda n, k: (22 * n ** 3 + (34 * k + 49) * n * n + (12 * k * k + 48 * k + 35) * n + 8 * (k + 1) ** 2) / (
    2 * (1 + 2 * n) * (3 + 2 * k + 3 * n) * (2 * n + 2 * k + 1) * (3 * n + 2 * k + 1) * (3 * n + 2 * k + 2))
wrep("W12", wz(lambda n, k: U(n, k) * S(n, k) * z ** n, lambda n, k: U(n, k) * Rr(n, k) * z ** n, pts=P5 + ((1, 0), (2, 0)), sign=-1),
     "2503.00570 sec.4.4 with z=16/27: F(n+1,k)-F(n,k) = -(G(n,k+1)-G(n,k))")
mp.dps = 80
for kk in (mpf(0), mpf(3) / 2):
    v = hsum(z, [1, 1 + kk, 1 + 2 * kk, 1 + kk / 2, H + kk / 2], [H, H + kk, 1 + 2 * kk / 3, (1 + 2 * kk) / 3, (2 + 2 * kk) / 3],
             lambda n: Rr(n, kk)) * (gamma(1 + kk) ** 2 / gamma(H + kk) ** 2 * gamma(H) ** 2)
    report("W12sum k=%s" % nstr(kk, 2), v, pi ** 2 / 8, "sum_n G(n,k) = pi^2/8 (2503.00570 sec 4.4, z=16/27)")

print()
print("=" * 100)
bad = [r for r in RESULTS if r[1] < 50]
print("SUMMARY: %d checks; minimum digits of agreement = %d; checks below 50 digits: %s"
      % (len(RESULTS), min(r[1] for r in RESULTS), bad if bad else "none"))
