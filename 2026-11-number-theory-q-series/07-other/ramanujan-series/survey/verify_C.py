#!/usr/bin/env python3
# verify_C.py -- numerical verification (mpmath, mp.dps = 80) of the identities listed in SURVEY-C.md
# Self-contained.  Run:  python3 verify_C.py > verify_C.out.txt
#
# Each entry is a dict:
#   id      : label used in SURVEY-C.md
#   src     : source (paper, equation number)
#   status  : 'OPEN' / 'OPEN?' / 'PROVED' (sanity/calibration) as assessed in the survey
#   start   : first summation index
#   rate    : |ratio| of consecutive terms (geometric convergence rate), used to choose #terms
#   term    : lambda n -> n-th summand (mpf)
#   rhs     : lambda -> conjectured closed form (mpf)
#
# Notation:  H(m,s) = sum_{k=1}^m k^-s ;  U(n,s) = sum_{i=0}^{n-1} (i+1/6)^-s ;
#            T(n,s) = 2^s H(2n,s) - H(n,s) = sum_{i=0}^{n-1} (i+1/2)^-s ;  rf = rising factorial.
import sys
from mpmath import (mp, mpf, mpc, pi, zeta, log, sqrt, catalan, polylog, rf, binomial,
                    factorial, gamma, nstr, floor, log10, fabs, im)

mp.dps = 80
WORK = 100                 # internal working precision (digits); results compared at 80
mp.dps = WORK

# ---------------------------------------------------------------- harmonic numbers (cached)
_Hc = {}
def H(m, s=1):
    m = int(m)
    L = _Hc.setdefault(s, [mpf(0)])
    while len(L) <= m:
        k = len(L)
        L.append(L[-1] + mpf(1) / mpf(k) ** s)
    return L[m]

_Uc = {}
def U(n, s=1):
    n = int(n)
    L = _Uc.setdefault(s, [mpf(0)])
    while len(L) <= n:
        i = len(L) - 1
        L.append(L[-1] + mpf(1) / (mpf(i) + mpf(1) / 6) ** s)
    return L[n]

def T(n, s=1):
    return mpf(2) ** s * H(2 * n, s) - H(n, s)

def memo(f):
    c = {}
    def g(n):
        if n not in c:
            c[n] = f(n)
        return c[n]
    return g

# ---------------------------------------------------------------- constants
def cvz(a, n=170):
    """Cohen-Villegas-Zagier acceleration of sum_{k>=0} (-1)^k a(k)."""
    d = (3 + sqrt(mpf(8))) ** n
    d = (d + 1 / d) / 2
    b = mpf(-1); c = -d; s = mpf(0)
    for k in range(n):
        c = b - c
        s += c * a(k)
        b = (k + n) * (k - n) * b / ((k + mpf(1) / 2) * (k + 1))
    return s / d

L3_2 = (zeta(2, mpf(1) / 3) - zeta(2, mpf(2) / 3)) / 9                       # L_{-3}(2)
L4_4 = (zeta(4, mpf(1) / 4) - zeta(4, mpf(3) / 4)) / 256                     # L_{-4}(4) = beta(4)
L8_2 = (zeta(2, mpf(1) / 8) + zeta(2, mpf(3) / 8)
        - zeta(2, mpf(5) / 8) - zeta(2, mpf(7) / 8)) / 64                    # L_{-8}(2)
G = catalan
ln2 = log(2); ln3 = log(3)
Li4h = polylog(4, mpf(1) / 2); Li5h = polylog(5, mpf(1) / 2)
ImLi3 = im(polylog(3, mpc(mpf(1) / 2, mpf(1) / 2)))
# zeta(bar5,1) := sum_{n>m>=1} (-1)^n / (n^5 m)   (Au's convention, bar = alternating sign on that index)
Zb51 = -cvz(lambda k: H(k) / mpf(k + 1) ** 5)
# Im Li_{4,1}(i,1)  = sum_{j>=0} (-1)^j H_{2j}/(2j+1)^4 ;  Li_{s1,s2}(a,b)=sum_{n1>n2>=1} a^n1 b^n2/(n1^s1 n2^s2)
ImLi41_i_1 = cvz(lambda j: H(2 * j) / mpf(2 * j + 1) ** 4)
# Im Li_{4,1}(i,-1) = sum_{j>=0} (-1)^j A_{2j}/(2j+1)^4,  A_N = sum_{m<=N} (-1)^m/m = H_{N/2} - H_N (N even)
ImLi41_i_m1 = cvz(lambda j: (H(j) - H(2 * j)) / mpf(2 * j + 1) ** 4)
CONST_CHECK = [
    ("zeta(bar5,1) CVZ n=170 vs n=220", Zb51, -cvz(lambda k: H(k) / mpf(k + 1) ** 5, 220)),
    ("Im Li41(i,1) CVZ n=170 vs n=220", ImLi41_i_1, cvz(lambda j: H(2 * j) / mpf(2 * j + 1) ** 4, 220)),
    ("Im Li41(i,-1) CVZ n=170 vs n=220", ImLi41_i_m1, cvz(lambda j: (H(j) - H(2 * j)) / mpf(2 * j + 1) ** 4, 220)),
]

h = lambda x: mpf(1) / 2 if x == 'h' else mpf(x)
q = lambda a, b: mpf(a) / mpf(b)

E = []          # the list of entries
def add(**kw):
    E.append(kw)

# =====================================================================================
# A.  K. C. Au, arXiv:2602.08721v4 (9 Nov 2026), Section 4 -- identities marked "?=" there
# =====================================================================================
# ---- Example VI:  a_n = (1/4)^n (1)_n^3/(1/2)_n^3,  n >= 1
a6 = memo(lambda n: mpf(4) ** (-n) * rf(1, n) ** 3 / rf(q(1, 2), n) ** 3)
S = "Au 2602.08721v4 "
add(id="Au-VI-(4.2)", src=S + "(4.2)", status="PROVED (Li-Chu 2023; calibration)", start=1, rate=0.25,
    term=lambda n: a6(n) * (3 * n - 1) * (n * H(n) - 1) / mpf(n) ** 4,
    rhs=lambda: -2 * pi ** 2 / 3 * ln2 + 14 * zeta(3) / 3)
add(id="Au-VI-(4.6)", src=S + "(4.6)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(17 * n - 6) / (6 * mpf(n) ** 5) - (3 * n - 1) * H(n) * (n * H(n) - 2) / mpf(n) ** 4),
    rhs=lambda: -128 * Li4h / 9 + 533 * pi ** 4 / 6480 - q(16, 27) * ln2 ** 4 - q(8, 27) * pi ** 2 * ln2 ** 2)
add(id="Au-VI-(4.7)", src=S + "(4.7)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(7 * n - 2) / (4 * mpf(n) ** 5) + (4 * n - 1) * H(n) / (2 * mpf(n) ** 4)
                            - (3 * n - 1) * (n * H(n) - 1) * H(2 * n) / mpf(n) ** 4),
    rhs=lambda: -112 * Li4h / 9 + 629 * pi ** 4 / 12960 - q(14, 27) * ln2 ** 4 + q(2, 27) * pi ** 2 * ln2 ** 2)
add(id="Au-VI-(4.8)", src=S + "(4.8)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(19 * n - 6) / (24 * mpf(n) ** 5) + (4 * n - 1) * H(2 * n) / mpf(n) ** 4
                            - (3 * n - 1) * H(2 * n) ** 2 / mpf(n) ** 3),
    rhs=lambda: -80 * Li4h / 9 - 17 * pi ** 4 / 5184 - q(10, 27) * ln2 ** 4 + q(4, 27) * pi ** 2 * ln2 ** 2)
add(id="Au-VI-(4.10)", src=S + "(4.10)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(13 * n - 5) / (3 * mpf(n) ** 6) + (5 * n - 2) * H(n) / (2 * mpf(n) ** 5)
                            - (3 * n - 1) * (-3 * H(n, 2) + 3 * n * H(n) * H(n, 2) - 2 * n * H(n, 3)) / (3 * mpf(n) ** 4)),
    rhs=lambda: -17 * pi ** 2 * zeta(3) / 12 + 341 * zeta(5) / 24 + pi ** 4 * ln2 / 36)
add(id="Au-VI-(4.11)", src=S + "(4.11)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(29 * n - 13) / (6 * mpf(n) ** 6) + (n - 2) * H(n) / (4 * mpf(n) ** 5)
                            + (3 * n - 1) * (5 * n * H(n, 3) + 6 * H(2 * n, 2) - 6 * n * H(n) * H(2 * n, 2)) / (3 * mpf(n) ** 4)),
    rhs=lambda: -19 * pi ** 2 * zeta(3) / 8 + 713 * zeta(5) / 48 + 13 * pi ** 4 * ln2 / 72)
add(id="Au-VI-(4.12)", src=S + "(4.12)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * ((19 * n - 7) / (3 * mpf(n) ** 6) - (17 * n - 6) * H(n) / (2 * mpf(n) ** 5)
                            - (3 * n - 1) * (-9 * H(n) ** 2 + 3 * n * H(n) ** 3 + 4 * n * H(n, 3)) / (3 * mpf(n) ** 4)),
    rhs=lambda: (-512 * Li5h / 9 + 67 * pi ** 2 * zeta(3) / 36 + 1705 * zeta(5) / 72 + 64 * ln2 ** 5 / 135
                 + q(32, 81) * pi ** 2 * ln2 ** 3 - 533 * pi ** 4 * ln2 / 1620))
add(id="Au-VI-(4.13)", src=S + "(4.13)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * ((200 * n - 71) / (12 * mpf(n) ** 6) - (19 * n - 6) * H(2 * n) / (4 * mpf(n) ** 5)
                            + 3 * (4 * n - 1) * H(2 * n) ** 2 / mpf(n) ** 4
                            - (3 * n - 1) * (6 * H(2 * n) ** 3 + 17 * H(n, 3)) / (3 * mpf(n) ** 3)),
    rhs=lambda: (-320 * Li5h / 9 + 545 * pi ** 2 * zeta(3) / 72 - 24149 * zeta(5) / 288 + 8 * ln2 ** 5 / 27
                 - q(16, 81) * pi ** 2 * ln2 ** 3 + 17 * pi ** 4 * ln2 / 1296))
add(id="Au-VI-(4.14)", src=S + "(4.14)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(7 * n - 3) / (2 * mpf(n) ** 7) + (5 * n - 2) * H(n, 2) / mpf(n) ** 5
                            - (3 * n - 1) * (2 * H(n, 2) ** 2 - H(n, 4)) / (2 * mpf(n) ** 3)),
    rhs=lambda: -pi ** 6 / 1440)
add(id="Au-VI-(4.15)", src=S + "(4.15)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * ((9 * n - 1) / (4 * mpf(n) ** 7) + (n - 2) * H(n, 2) / (4 * mpf(n) ** 5)
                            + (5 * n - 2) * H(2 * n, 2) / mpf(n) ** 5
                            - (3 * n - 1) * (3 * H(n, 4) + 8 * H(n, 2) * H(2 * n, 2)) / (4 * mpf(n) ** 3)),
    rhs=lambda: -17 * pi ** 6 / 2880)
add(id="Au-VI-(4.16)", src=S + "(4.16)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * ((73 * n - 21) / (8 * mpf(n) ** 7) + (n - 2) * H(2 * n, 2) / mpf(n) ** 5
                            - (3 * n - 1) * (21 * H(n, 4) + 32 * H(2 * n, 2) ** 2 + 32 * H(2 * n, 4)) / (8 * mpf(n) ** 3)),
    rhs=lambda: -337 * pi ** 6 / 5760)
add(id="Au-VI-(4.17)", src=S + "(4.17)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * (-(3 * n - 1) * (n * H(n) - 1)
                            * (-8 + 7 * n ** 3 * H(n, 3) + 8 * n ** 3 * H(2 * n, 3)) / (8 * mpf(n) ** 7)),
    rhs=lambda: -7 * zeta(3) ** 2 / 3 + q(2, 3) * pi ** 2 * zeta(3) * ln2 - pi ** 6 / 224)
add(id="Au-VI-(4.18)", src=S + "(4.18)", status="OPEN", start=1, rate=0.25,
    term=lambda n: a6(n) * ((299 * n - 129) / (36 * mpf(n) ** 7) - (29 * n - 13) * H(n) / (3 * mpf(n) ** 6)
                            + (n - 2) * H(n) ** 2 / (4 * mpf(n) ** 5) - (17 * n - 6) * H(2 * n, 2) / (3 * mpf(n) ** 5)
                            + (3 * n - 1) * (-40 * H(n, 3) + 40 * n * H(n) * H(n, 3) + 3 * n * H(n, 4)
                                             + 48 * H(n) * H(2 * n, 2) - 24 * n * H(n) ** 2 * H(2 * n, 2)) / (12 * mpf(n) ** 4)),
    rhs=lambda: (q(184, 9) * Zb51 - 179 * zeta(3) ** 2 / 6 + q(19, 3) * pi ** 2 * zeta(3) * ln2
                 - q(713, 18) * zeta(5) * ln2 + 827 * pi ** 6 / 25920 - q(13, 54) * pi ** 4 * ln2 ** 2))

# ---- Example VII:  a_n = 64^-n (1/2)_n^3/(1)_n^3,  n >= 0
a7 = memo(lambda n: mpf(64) ** (-n) * rf(q(1, 2), n) ** 3 / rf(1, n) ** 3)
p0_7 = lambda n: 42 * n + 5
add(id="Au-VII-(4.19)", src=S + "(4.19)", status="PROVED (Sun Conj.29 case; Zhou 2025; calibration)", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (1 + q(1, 7) * (H(2 * n) - H(n)) * p0_7(n)),
    rhs=lambda: 32 * ln2 / (7 * pi))
add(id="Au-VII-(4.20)", src=S + "(4.20)", status="OPEN? (Au marks ?=, but says provable by WZ and modular methods [Zhou])", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (mpf(31 + 162 * n) / (3 + 18 * n) + (-H(n) / 2 - H(3 * n) + 2 * H(6 * n)) * p0_7(n)),
    rhs=lambda: 48 * ln2 / pi)
add(id="Au-VII-(4.21)", src=S + "(4.21)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (mpf(35 + 162 * n) / (1 + 6 * n) + (-H(n) / 2 + 2 * U(n, 1)) * p0_7(n)),
    rhs=lambda: 32 / sqrt(mpf(3)) + 80 * ln2 / pi)
add(id="Au-VII-(4.24)", src=S + "(4.24)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (mpf(43 + 396 * n + 972 * n * n) / (9 * (1 + 2 * n) * (1 + 6 * n) ** 2)
                            + (-H(3 * n, 2) / 4 + H(6 * n, 2)) * p0_7(n)),
    rhs=lambda: 14 * pi / 9)
def _t426(n):
    X = H(n) + 2 * H(3 * n) - 4 * H(6 * n)
    return a7(n) * (mpf((6 * n + 5) * (78 * n + 17)) / (9 * (2 * n + 1) * (6 * n + 1) ** 2)
                    + mpf(p0_7(n)) / 4 * X ** 2 - (162 * n + 31) * X / (18 * n + 3))
add(id="Au-VII-(4.26)", src=S + "(4.26)", status="OPEN", start=0, rate=1 / 64, term=_t426,
    rhs=lambda: 144 * ln2 ** 2 / pi - 34 * pi / 9)
def _t427(n):
    p0 = p0_7(n); p2 = 35 + 162 * n; p5 = 77 + 366 * n
    u1 = U(n, 1); u2 = U(n, 2); h1 = H(n); h2n = H(2 * n)
    return a7(n) * (14 * u1 - mpf(8) / (6 * n + 1) ** 2 - mpf(23) / (6 * n + 3)
                    + q(1, 6) * (3 * h1 ** 2 + 12 * h2n * u1 - 3 * h1 * (h2n + 4 * u1) - 2 * u2) * p0
                    + h2n * p2 / (1 + 6 * n) - h1 * p5 / (2 + 12 * n))
add(id="Au-VII-(4.27)", src=S + "(4.27)", status="OPEN", start=0, rate=1 / 64, term=_t427,
    rhs=lambda: -120 * L3_2 / pi - 104 * pi / 9 + 160 * ln2 ** 2 / pi + 64 * ln2 / sqrt(mpf(3)))
add(id="Au-VII-(4.28)", src=S + "(4.28)", status="OPEN (would follow from two conjectures of Sun)", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (-mpf(9) / (1 + 2 * n) ** 2 + q(1, 3) * (-2 * H(n, 3) + 17 * H(2 * n, 3)) * p0_7(n)),
    rhs=lambda: 80 * zeta(3) / pi - 128 * G / 3)
add(id="Au-VII-(4.29)", src=S + "(4.29)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * ((-479 - 1275 * (1 + 2 * n) * H(n) + 1275 * (1 + 2 * n) * H(2 * n)
                             + (357 + 1428 * n + 1428 * n * n) * H(2 * n, 2)) / mpf(1 + 2 * n) ** 2
                            + (-4 * H(n, 3) + 51 * (H(2 * n) - H(n)) * H(2 * n, 2)) * p0_7(n)),
    rhs=lambda: -2432 * G - 64 * zeta(3) / pi + 816 * pi * ln2)
add(id="Au-VII-(4.30)", src=S + "(4.30)", status="OPEN (also Sun arXiv:2307.03086 (5.14))", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (mpf(25) / (9 * mpf(1 + 2 * n) ** 3) + (-H(n, 4) / 16 + H(2 * n, 4)) * p0_7(n)),
    rhs=lambda: 5 * pi ** 3 / 54)
add(id="Au-VII-(4.31)", src=S + "(4.31)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (4 * (-5620 + 621 * (1 + 2 * n) ** 2 * H(n, 2)) / (27 * mpf(1 + 2 * n) ** 3)
                            + (H(n, 2) ** 2 + H(n, 4)) * p0_7(n) / 2),
    rhs=lambda: -10874 * pi ** 3 / 405)
add(id="Au-VII-(4.32)", src=S + "(4.32)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (25 * H(n, 2) / (1 + 2 * n)
                            + 4 * (-3103 + 621 * (1 + 2 * n) ** 2 * H(2 * n, 2)) / (27 * mpf(1 + 2 * n) ** 3)
                            + (H(n, 4) / 4 + H(n, 2) * H(2 * n, 2)) * p0_7(n)),
    rhs=lambda: -6002 * pi ** 3 / 405)
add(id="Au-VII-(4.33)", src=S + "(4.33)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (5 * (-673 + 270 * (1 + 2 * n) ** 2 * H(2 * n, 2)) / (27 * mpf(1 + 2 * n) ** 3)
                            + (H(n, 4) / 16 + H(2 * n, 2) ** 2) * p0_7(n)),
    rhs=lambda: -3251 * pi ** 3 / 810)
add(id="Au-VII-(4.34)", src=S + "(4.34)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (9 * p0_7(n) * (16 * (H(n) - H(2 * n)) * (2 * H(n, 3) - 17 * H(2 * n, 3)) - H(n, 4))
                            + 3888 * (H(n) - H(2 * n)) / mpf(2 * n + 1) ** 2 - 2016 * H(n, 3) + 17136 * H(2 * n, 3)
                            + mpf(9392) / mpf(2 * n + 1) ** 3),
    rhs=lambda: 73728 * ImLi3 + 69120 * zeta(3) * ln2 / pi - 22904 * pi ** 3 / 15 - 2304 * pi * ln2 ** 2)
add(id="Au-VII-(4.35)", src=S + "(4.35)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (H(n, 2) * ((17 * H(2 * n, 3) - 2 * H(n, 3)) * p0_7(n) - mpf(27) / (2 * n + 1) ** 2)
                            + 4 * (5 * p0_7(n) * H(n, 5) + (391 * H(2 * n, 3) - 46 * H(n, 3)) / (2 * n + 1)
                                   - 159 * p0_7(n) * H(2 * n, 5) - mpf(4) / mpf(2 * n + 1) ** 4)),
    rhs=lambda: 1536 * L4_4 + 440 * pi * zeta(3) - 9856 * zeta(5) / pi)
add(id="Au-VII-(4.36)", src=S + "(4.36)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * ((H(2 * n) - H(n)) * ((144 * H(2 * n, 4) - 9 * H(n, 4)) * p0_7(n) + mpf(400) / mpf(2 * n + 1) ** 3)
                            - 63 * H(n, 4) + 6 * (42 * n + 5) * H(n, 5) - 8064 * n * H(2 * n, 5)
                            + 16 * (63 * H(2 * n, 4) - 60 * H(2 * n, 5) - mpf(25) / mpf(2 * n + 1) ** 4)),
    rhs=lambda: q(80, 3) * pi ** 3 * ln2 - 2976 * zeta(5) / pi)
add(id="Au-VII-(4.37)", src=S + "(4.37)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (-50 * H(n, 3) / (2 * n + 1) + 425 * H(2 * n, 3) / (2 * n + 1)
                            + p0_7(n) * (6 * H(n, 5) - 2 * H(n, 3) * H(2 * n, 2) + 17 * H(2 * n, 2) * H(2 * n, 3) - 191 * H(2 * n, 5))
                            - 27 * H(2 * n, 2) / mpf(2 * n + 1) ** 2 - mpf(2) / mpf(2 * n + 1) ** 4),
    rhs=lambda: 512 * L4_4 + 120 * pi * zeta(3) - 2960 * zeta(5) / pi)
add(id="Au-VII-(4.38)", src=S + "(4.38)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (p0_7(n) * (-27 * H(n, 2) * H(n, 4) - 172 * H(n, 6) + 432 * H(n, 2) * H(2 * n, 4) + 11008 * H(2 * n, 6))
                            + 1200 * H(n, 2) / mpf(2 * n + 1) ** 3 - 2484 * (H(n, 4) - 16 * H(2 * n, 4)) / (2 * n + 1)
                            - mpf(192) / mpf(2 * n + 1) ** 5),
    rhs=lambda: 44 * pi ** 5 / 15)
add(id="Au-VII-(4.39)", src=S + "(4.39)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (p0_7(n) * (27 * H(n, 2) ** 3 + 81 * H(n, 4) * H(n, 2) + 13022 * H(n, 6) - 829952 * H(2 * n, 6))
                            + 7452 * H(n, 2) ** 2 / (2 * n + 1) - 134880 * H(n, 2) / mpf(2 * n + 1) ** 3
                            + 7452 * H(n, 4) / (2 * n + 1) + mpf(428928) / mpf(2 * n + 1) ** 5),
    rhs=lambda: 120706 * pi ** 5 / 105)
add(id="Au-VII-(4.40)", src=S + "(4.40)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (p0_7(n) * (-155 * H(n, 6) - 108 * H(n, 4) * H(2 * n, 2) + 1728 * H(2 * n, 2) * H(2 * n, 4) + 9920 * H(2 * n, 6))
                            - 2700 * H(n, 4) / (2 * n + 1) + 4800 * H(2 * n, 2) / mpf(2 * n + 1) ** 3
                            + 43200 * H(2 * n, 4) / (2 * n + 1) - mpf(1152) / mpf(2 * n + 1) ** 5),
    rhs=lambda: mpf(0))
add(id="Au-VII-(4.41)", src=S + "(4.41)", status="OPEN", start=0, rate=1 / 64,
    term=lambda n: a7(n) * (p0_7(n) * (4 * H(n, 3) ** 2 - 68 * H(2 * n, 3) * H(n, 3) + 289 * H(2 * n, 3) ** 2
                                       + 2 * H(n, 6) - 127 * H(2 * n, 6))
                            + 108 * H(n, 3) / mpf(2 * n + 1) ** 2 - 918 * H(2 * n, 3) / mpf(2 * n + 1) ** 2
                            - mpf(54) / mpf(2 * n + 1) ** 5),
    rhs=lambda: (3072 * ImLi41_i_1 + 1024 * ImLi41_i_m1 + 2048 * L4_4 * ln2 + 3600 * zeta(3) ** 2 / pi
                 - 9434 * pi ** 5 / 945))

# ---- Example VIII:  a_n = 9^-n (1/2)_n(1/4)_n(3/4)_n/(1)_n^3, p0 = 1+8n, n >= 0
a8 = memo(lambda n: mpf(9) ** (-n) * rf(q(1, 2), n) * rf(q(1, 4), n) * rf(q(3, 4), n) / rf(1, n) ** 3)
r3 = sqrt(mpf(3)); r2 = sqrt(mpf(2))
add(id="Au-VIII-1", src=S + "Ex. VIII, 1st, corrected placement of p0 (= Sun, Conj. 14 of arXiv:2210.07238)", status="OPEN? (Sun's conjecture; belongs to Sun's lists)", start=0, rate=1 / 9,
    # NOTE: as printed in v4 the summand reads a_n*p0*(1+(1/6)(-4H_n+3H_2n)); that fails numerically (1.1108 vs 1.0191).
    # The form below, a_n*(1+(1/6)(-4H_n+3H_2n)*p0), is the one that holds (misplaced p0 in the paper).
    term=lambda n: a8(n) * (1 + q(1, 6) * (-4 * H(n) + 3 * H(2 * n)) * (1 + 8 * n)),
    rhs=lambda: 8 * ln2 / (r3 * pi))
add(id="Au-VIII-2", src=S + "Ex. VIII, 2nd", status="PROVED (Sun Conj.29 case; Zhou; calibration)", start=0, rate=1 / 9,
    term=lambda n: a8(n) * (2 + (1 + 8 * n) * (H(4 * n) - H(n))),
    rhs=lambda: 4 * r3 * ln2 / pi + r3 * ln3 / pi)
add(id="Au-VIII-4", src=S + "Ex. VIII, 4th", status="OPEN", start=0, rate=1 / 9,
    term=lambda n: a8(n) * (-4 * H(n) + 3 * H(2 * n) + q(1, 24) * (32 * H(n) ** 2 - 48 * H(n) * H(2 * n)
                            + 18 * H(2 * n) ** 2 + 5 * H(n, 2)) * (1 + 8 * n)),
    rhs=lambda: 32 * ln2 ** 2 / (r3 * pi) - 37 * pi / (24 * r3))
add(id="Au-VIII-5", src=S + "Ex. VIII, 5th", status="OPEN", start=0, rate=1 / 9,
    term=lambda n: a8(n) * (-H(n) + H(4 * n) + q(1, 16) * (4 * H(n) ** 2 - 8 * H(n) * H(4 * n) + 4 * H(4 * n) ** 2
                            + H(n, 2) - 4 * H(4 * n, 2)) * (1 + 8 * n)),
    rhs=lambda: -3 * r3 * pi / 16 + 2 * r3 * ln2 ** 2 / pi + r3 * ln3 ** 2 / (8 * pi) + r3 * ln3 * ln2 / pi)
add(id="Au-VIII-6", src=S + "Ex. VIII, 6th", status="OPEN", start=0, rate=1 / 9,
    term=lambda n: a8(n) * (-5 * H(n, 2) + 18 * H(2 * n, 2) + q(1, 6) * (20 * H(n) * H(n, 2) - 15 * H(2 * n) * H(n, 2)
                            + 14 * H(n, 3) - 72 * H(n) * H(2 * n, 2) + 54 * H(2 * n) * H(2 * n, 2) - 108 * H(2 * n, 3)) * (1 + 8 * n)),
    rhs=lambda: 45 * L3_2 / 2 - 94 * zeta(3) / (r3 * pi) + 4 * pi * ln2 / (3 * r3))
add(id="Au-VIII-7", src=S + "Ex. VIII, 7th", status="OPEN", start=0, rate=1 / 9,
    term=lambda n: a8(n) * (1 + 8 * n) * (50 * H(n, 2) ** 2 - 360 * H(2 * n, 2) * H(n, 2) + 648 * H(2 * n, 2) ** 2
                                          + 41 * H(n, 4) - 648 * H(2 * n, 4)),
    rhs=lambda: -2 * pi ** 3 / (15 * r3))

# ---- Example IX:  a_n = (-27/512)^n (1/2)_n(1/6)_n(5/6)_n/(1)_n^3,  n >= 0
a9 = memo(lambda n: mpf(-27) ** n / mpf(512) ** n * rf(q(1, 2), n) * rf(q(1, 6), n) * rf(q(5, 6), n) / rf(1, n) ** 3)
P0 = lambda n: 15 + 154 * n
P1 = lambda n: 61 + 110 * n
P2 = lambda n: 1 + 6 * n
add(id="Au-IX-1", src=S + "Ex. IX, 1st", status="OPEN", start=0, rate=27 / 512,
    term=lambda n: a9(n) * ((H(2 * n) - H(n)) * P0(n) + mpf(P1(n)) / (3 * (1 + 2 * n))),
    rhs=lambda: 64 * r2 * ln2 / pi)
add(id="Au-IX-2", src=S + "Ex. IX, 2nd", status="PROVED (Sun Conj.29 case; Zhou; calibration)", start=0, rate=27 / 512,
    term=lambda n: a9(n) * (q(1, 3) + q(1, 154) * (-H(n) + T(3 * n, 1)) * P0(n)),
    rhs=lambda: 80 * r2 * ln2 / (77 * pi))
add(id="Au-IX-3", src=S + "Ex. IX, 3rd", status="OPEN", start=0, rate=27 / 512,
    term=lambda n: a9(n) * (q(1, 15) * (-4 * H(n, 2) + 15 * H(2 * n, 2)) * P0(n) - mpf(P2(n)) / (1 + 2 * n) ** 2),
    rhs=lambda: -16 * r2 * pi / 45)
add(id="Au-IX-4", src=S + "Ex. IX, 4th", status="OPEN", start=0, rate=27 / 512,
    term=lambda n: a9(n) * (q(1, 18) * (9 * H(n) ** 2 - 18 * H(n) * H(2 * n) + 9 * H(2 * n) ** 2 + H(n, 2)) * P0(n)
                            - (H(n) - H(2 * n)) * P1(n) / (3 * (1 + 2 * n)) - mpf(P2(n)) / (2 * (1 + 2 * n) ** 2)),
    rhs=lambda: 64 * r2 * ln2 ** 2 / pi - 88 * r2 * pi / 27)
add(id="Au-IX-5", src=S + "Ex. IX, 5th", status="OPEN", start=0, rate=27 / 512,
    term=lambda n: a9(n) * (-2 * (H(n) - T(3 * n, 1)) + q(1, 154) * (3 * H(n) ** 2 + H(n, 2) - 6 * H(n) * T(3 * n, 1)
                            + 3 * T(3 * n, 1) ** 2 - 3 * T(3 * n, 2)) * P0(n)),
    rhs=lambda: 1200 * r2 * ln2 ** 2 / (77 * pi) - 16 * r2 * pi / 21)
add(id="Au-IX-6", src=S + "Ex. IX, 6th", status="OPEN", start=0, rate=27 / 512,
    term=lambda n: a9(n) * (q(1, 3) * (77 * H(2 * n, 2) + mpf(6 * n - 1) / mpf(2 * n + 1) ** 3)
                            + (H(n) - T(3 * n, 1)) * P2(n) / (2 * mpf(1 + 2 * n) ** 2) - 308 * H(n, 2) / 45
                            + q(1, 90) * (12 * H(n) * H(n, 2) + 8 * H(n, 3) - 45 * H(n) * H(2 * n, 2) - 60 * H(2 * n, 3)
                                          - 12 * H(n, 2) * T(3 * n, 1) + 45 * H(2 * n, 2) * T(3 * n, 1)) * P0(n)),
    rhs=lambda: 512 * L8_2 / 45 - 832 * r2 * zeta(3) / (45 * pi) - q(8, 9) * r2 * pi * ln2)

# ---- Example X:  a_n = (-2^10/5^5)^n (1)_n^9/((1/2)_n^5 (1/5)_n(2/5)_n(3/5)_n(4/5)_n),  n >= 1
aX = memo(lambda n: mpf(-1024) ** n / mpf(3125) ** n * rf(1, n) ** 9
          / (rf(q(1, 2), n) ** 5 * rf(q(1, 5), n) * rf(q(2, 5), n) * rf(q(3, 5), n) * rf(q(4, 5), n)))
X0 = lambda n: 30 - 425 * n + 2275 * n ** 2 - 5600 * n ** 3 + 5532 * n ** 4
X1 = lambda n: 54 - 680 * n + 3185 * n ** 2 - 6720 * n ** 3 + 5532 * n ** 4
X2 = lambda n: 510 - 7225 * n + 37835 * n ** 2 - 87620 * n ** 3 + 76632 * n ** 4
add(id="Au-X-0", src="Zhao's zeta(5) series (Au 2312.14051 Ex. X)", status="PROVED (Au; calibration)", start=1, rate=1024 / 3125,
    term=lambda n: aX(n) * X0(n) / mpf(n) ** 9, rhs=lambda: -380928 * zeta(5))
add(id="Au-X-1", src=S + "Ex. X, 1st", status="OPEN (the one of six not covered by Shalev 2610.34482)", start=1, rate=1024 / 3125,
    term=lambda n: aX(n) * (X1(n) / mpf(n) ** 10 + X0(n) * (-3 * H(n) + 2 * H(2 * n) + H(5 * n)) / mpf(n) ** 9),
    rhs=lambda: 786432 * Zb51 - 294912 * zeta(3) ** 2 - 1523712 * zeta(5) * ln2 + 32768 * pi ** 6 / 105)
add(id="Au-X-2", src=S + "Ex. X, 2nd", status="PROVED (Shalev 2610.34482, (1.4); calibration)", start=1, rate=1024 / 3125,
    term=lambda n: aX(n) * (X0(n) * (11 * H(2 * n, 2) - 7 * H(n, 2)) / mpf(n) ** 9 + X2(n) / (4 * mpf(n) ** 11)),
    rhs=lambda: -3901440 * zeta(7))

# (further entries B..., appended below)
# ---- Sun's conjecture (4.1) of "New series involving binomial coefficients (I)", quoted with "?=" in
#      Au 2602.08721v4 after Example IV ("evidently does not follow from the formulas above ... We encourage
#      interested readers to seek other WZ-pairs that might prove it").  NOTE: as displayed by Au with
#      a_n = (-1/27)^n (1)_n^2/((1/3)_n(2/3)_n) the identity fails numerically; it holds with the extra factor 1/n^2:
add(id="Sun-(4.1)/Au-IV", src=S + "remark after (3.32) in Ex. IV (Sun, NUJM Biquarterly 2024, (4.1)); factor 1/n^2 restored",
    status="OPEN (Au: not reachable by his three WZ pairs)", start=1, rate=1 / 27,
    term=lambda n: (mpf(-1) ** n / (mpf(n) ** 2 * binomial(2 * n, n) * binomial(3 * n, n))
                    * (mpf(7 * n - 2) / (2 * n - 1) * H(n - 1) - mpf(1) / (6 * n))),
    rhs=lambda: pi ** 2 * ln2 / 6 - 11 * zeta(3) / 12)

# =====================================================================================
# B.  Rational Ramanujan-type series for 1/pi^2 from the Cohen-Guillera table (arXiv:2101.12592, Sect. 5)
#     that were NOT seen proved in this sweep (the 1/pi^k families are surveyed by another agent; listed
#     here because Au (arXiv:2312.14051, Sect. 1 and 3.4) names CG#10 as out of reach of the WZ seeds).
# =====================================================================================
def _hyp5(ss, z):
    return memo(lambda n: mpf(z) ** n * rf(ss[0], n) * rf(ss[1], n) * rf(ss[2], n) * rf(ss[3], n) * rf(ss[4], n) / rf(1, n) ** 5)
_cg10 = _hyp5([q(1, 2), q(1, 8), q(3, 8), q(5, 8), q(7, 8)], q(1, 7 ** 4))
add(id="CG-10", src="Cohen-Guillera 2101.12592 table #10; Au 2312.14051 p.2 ('?=')", status="OPEN (Au: not reachable by listed WZ seeds)",
    start=0, rate=1 / 7 ** 4,
    term=lambda n: _cg10(n) * (1920 * n * n + 304 * n + 15), rhs=lambda: 56 * sqrt(mpf(7)) / pi ** 2)
_cg7 = _hyp5([q(1, 2), q(1, 3), q(2, 3), q(1, 6), q(5, 6)], q(-1, 80 ** 3))
add(id="CG-7", src="Cohen-Guillera 2101.12592 table #7", status="OPEN? (no proof seen in sweep; see survey A)",
    start=0, rate=1 / 80 ** 3,
    term=lambda n: _cg7(n) * (5418 * n * n + 693 * n + 29), rhs=lambda: 128 * sqrt(mpf(5)) / pi ** 2)
_cg8 = _hyp5([q(1, 2), q(1, 3), q(2, 3), q(1, 6), q(5, 6)], q(3 ** 6, 5 ** 6))
add(id="CG-8", src="Cohen-Guillera 2101.12592 table #8", status="OPEN? (no proof seen in sweep; see survey A)",
    start=0, rate=(3 / 5) ** 6,
    term=lambda n: _cg8(n) * (532 * n * n + 126 * n + 9), rhs=lambda: mpf(375) / (4 * pi ** 2))
_cg9 = _hyp5([q(1, 2), q(1, 4), q(3, 4), q(1, 6), q(5, 6)], q(-1, 2 ** 10))
add(id="CG-9", src="Cohen-Guillera 2101.12592 table #9", status="OPEN? (no proof seen in sweep; see survey A)",
    start=0, rate=1 / 2 ** 10,
    term=lambda n: _cg9(n) * (1640 * n * n + 278 * n + 15), rhs=lambda: 256 / (sqrt(mpf(3)) * pi ** 2))
_cg5 = _hyp5([q(1, 2), q(1, 3), q(2, 3), q(1, 4), q(3, 4)], q(-1, 48))
add(id="CG-5", src="Cohen-Guillera 2101.12592 table #5", status="OPEN? (no proof seen in sweep; see survey A)",
    start=0, rate=1 / 48,
    term=lambda n: _cg5(n) * (252 * n * n + 63 * n + 5), rhs=lambda: 48 / pi ** 2)

# =====================================================================================
# C.  J. Zuniga's unproven fast series for log 2 (MO 466853; arXiv:2506.08245 eq. (23))
# =====================================================================================
add(id="Zuniga-log2", src="Zuniga, MO 466853 (2024) and arXiv:2506.08245 eq.(23)", status="OPEN (stated unproven, June 2025)",
    start=0, rate=1 / 221184,
    term=lambda n: (mpf(1) / 576 / mpf(8) ** n
                    * (686430 * n ** 3 + 1317033 * n ** 2 + 798173 * n + 153712)
                    / mpf((12 * n + 1) * (12 * n + 5) * (12 * n + 7) * (12 * n + 11))
                    * factorial(3 * n) * factorial(4 * n) * factorial(6 * n) / (factorial(n) * factorial(12 * n))),
    rhs=lambda: ln2)
# same series in the form of arXiv:2506.08245 (23): sum_{n>=1} P(n)/(3n(2n-1)(3n-1)(3n-2)) (2^-13 3^-3)^n
#     (1)_n(1/2)_n(1/3)_n(2/3)_n / ((1/12)_n(5/12)_n(7/12)_n(11/12)_n),  P = 686430n^3-742257n^2+223397n-13858
add(id="Zuniga-log2-b", src="arXiv:2506.08245 eq.(23) (as printed)", status="OPEN (same identity, other normalisation)",
    start=1, rate=1 / 221184,
    term=lambda n: (mpf(686430 * n ** 3 - 742257 * n ** 2 + 223397 * n - 13858) / (3 * n * (2 * n - 1) * (3 * n - 1) * (3 * n - 2))
                    / mpf(221184) ** n * rf(1, n) * rf(q(1, 2), n) * rf(q(1, 3), n) * rf(q(2, 3), n)
                    / (rf(q(1, 12), n) * rf(q(5, 12), n) * rf(q(7, 12), n) * rf(q(11, 12), n))),
    rhs=lambda: ln2)

# =====================================================================================
# D.  Deyi Chen, MathOverflow 486961 (2025) and 507693 (2026)
# =====================================================================================
add(id="Chen-486961-(3)", src="D. Chen, MO 486961 eq.(3)", status="OPEN? (no proof posted; half-shift of Zhao 1/pi^4 pair expected)",
    start=0, rate=27 / 256,
    term=lambda k: ((-1) ** k * binomial(6 * k, k) * binomial(5 * k, k) * binomial(4 * k, k) ** 2 * binomial(3 * k, k)
                    * (4 * k + 1) * (6 * k + 1) * (2264 * k ** 4 + 6118 * k ** 3 + 6267 * k ** 2 + 2884 * k + 503)
                    / (384 * mpf(2 * k + 1) ** 9 * binomial(2 * k, k) ** 10)),
    rhs=lambda: zeta(3))
add(id="Chen-486961-(5)", src="D. Chen, MO 486961 eq.(5)", status="OPEN? (no proof posted; related to Sun Conj.4.6 / Hou-Sun 2604.15172)",
    start=1, rate=1 / 27,
    term=lambda k: ((-1) ** k * (28 * k * k + 10 * k + 1) * binomial(2 * k, k) ** 2
                    / ((6 * k + 1) * mpf(2 * k + 1) ** 2 * mpf(k) ** 3 * binomial(3 * k, k) * binomial(6 * k, 3 * k))),
    rhs=lambda: 24 - 20 * zeta(3))
add(id="Chen-507693-(2)", src="D. Chen, MO 507693 eq.(2)", status="OPEN (unanswered)",
    start=0, rate=1 / 4,
    term=lambda k: (H(2 * k) - H(k)) * binomial(2 * k, k) ** 3 / mpf(256) ** k,
    rhs=lambda: (log(256) - sqrt(mpf(3)) * pi) / 6 * sqrt(mpf(3)) * gamma(q(1, 3)) ** 6 / (4 * mpf(2) ** q(2, 3) * pi ** 4))
# Chen's Conjecture 1 (sum f'(k) = -2 Im(tau) pi sum f(k)) in "vanishing" form, two instances given there:
add(id="Chen-507693-C1a", src="D. Chen, MO 507693 Conj.1, f=C(4k,2k)C(2k,k)^2/396^(4k), lambda=-sqrt(58)pi", status="OPEN (unanswered)",
    start=0, rate=256 / 396 ** 4,
    term=lambda k: (4 * (H(4 * k) - H(k)) - 4 * log(396) + sqrt(mpf(58)) * pi)
                   * binomial(4 * k, 2 * k) * binomial(2 * k, k) ** 2 / mpf(396) ** (4 * k),
    rhs=lambda: mpf(0))
add(id="Chen-507693-C1b", src="D. Chen, MO 507693 Conj.1, f=C(6k,3k)C(3k,k)C(2k,k)/255^(3k), lambda=-2sqrt(7)pi", status="OPEN (unanswered)",
    start=0, rate=1728 / 255 ** 3,
    term=lambda k: (3 * (2 * H(6 * k) - H(3 * k) - H(k)) - 3 * log(255) + 2 * sqrt(mpf(7)) * pi)
                   * binomial(6 * k, 3 * k) * binomial(3 * k, k) * binomial(2 * k, k) / mpf(255) ** (3 * k),
    rhs=lambda: mpf(0))

# =====================================================================================
# E.  Guillera's conjecture: Taylor coefficient x^(2m+1) of the "x-extended" Ramanujan-like series for 1/pi^m
#     equals r*L(chi,m+1), r rational (arXiv:2312.16827 Sect. 6; MO 473746, unanswered; cf. Cohen-Guillera Sect. 6).
#     S(x) = sum_{n>=0} prod_i (s_i)_{n+x}/(1)_{n+x} * P(n+x) * (-1)^n |z0|^(n+x),   (s)_y := Gamma(s+y)/Gamma(s).
#     Coefficients extracted by a discrete Cauchy integral on |x| = r (N = 128 nodes).
# =====================================================================================
from mpmath import expjpi, re as _re
def _taylor_coeff(Sf, r, N, k):
    return sum(Sf(r * expjpi(mpf(2) * j / N)) * expjpi(-mpf(2) * j * k / N) for j in range(N)) / N / r ** k
def _make_S(ss, poly, absz, sign, nt):
    def Sf(x):
        t = absz ** x
        for s in ss:
            t *= gamma(s + x) / (gamma(s) * gamma(1 + x))
        tot = 0
        for n in range(nt):
            y = n + x
            tot += t * poly(y)
            rr = sign * absz
            for s in ss:
                rr *= (s + y) / (1 + y)
            t *= rr
        return tot
    return Sf
def _L(chi, m, s):
    return sum(chi[a] * zeta(s, mpf(a) / m) for a in range(1, m)) / mpf(m) ** s
def _guillera1():
    Sf = _make_S([q(1, 2), q(1, 4), q(3, 4)], lambda y: 65 * y + 8, q(16, 63) ** 2, -1, 115)
    c3 = _re(_taylor_coeff(Sf, q(1, 16), 128, 3))
    return c3, q(2205, 2) * _L([0, 1, 1, -1, 1, -1, -1], 7, 2)
def _guillera2():
    Sf = _make_S([q(1, 2), q(1, 3), q(2, 3), q(1, 6), q(5, 6)], lambda y: 5418 * y * y + 693 * y + 29, q(1, 80 ** 3), -1, 32)
    c5 = _re(_taylor_coeff(Sf, q(1, 24), 128, 5))
    return c5, -5376000 * _L([0, 1, -1, -1, 1], 5, 3)
add(id="Guillera-L(-7,2)", src="Guillera MO 473746 ex.1 / arXiv:2312.16827 Sect.6; rational 2205/2 found here by PSLQ",
    status="OPEN (general conjecture; r identified numerically)", custom=_guillera1)
add(id="Guillera-L(5,3)", src="Guillera MO 473746 ex.2 / arXiv:2312.16827 Sect.6; rational -5376000 found here",
    status="OPEN (general conjecture; r identified numerically)", custom=_guillera2)


# ---------------------------------------------------------------- driver
def nterms(rate, digits=WORK + 8):
    import math
    return int(digits * math.log(10) / -math.log(rate)) + 25

def run():
    out = []
    print("verify_C.py   mp.dps(work) = %d, agreement reported in decimal digits (capped at 80)" % WORK)
    print("constant self-checks:")
    for name, a, b in CONST_CHECK:
        d = fabs(a - b)
        print("  %-40s |diff| = %s" % (name, nstr(d, 3)))
    print("-" * 110)
    print("%-16s %-7s %-6s %-28s %s" % ("id", "digits", "terms", "status", "lhs (25 digits)"))
    nbad = 0
    for e in E:
        if "custom" in e:
            N = 0
            s, r = e["custom"]()
        else:
            N = e.get("N") or nterms(e["rate"])
            s = mpf(0)
            for n in range(e["start"], e["start"] + N):
                s += e["term"](n)
            r = e["rhs"]()
        diff = fabs(s - r)
        scale = max(fabs(r), fabs(s), mpf(1))
        if diff == 0:
            dig = 80
        else:
            dig = int(min(80, floor(-log10(diff / scale))))
        ok = dig >= 50
        if not ok:
            nbad += 1
        print("%-16s %-7d %-6d %-28s %s   %s" % (e["id"], dig, N, e["status"][:28], nstr(s, 25), "" if ok else "<<< MISMATCH"))
    print("-" * 110)
    print("entries: %d, mismatches (<50 digits): %d" % (len(E), nbad))

if __name__ == "__main__":
    run()
