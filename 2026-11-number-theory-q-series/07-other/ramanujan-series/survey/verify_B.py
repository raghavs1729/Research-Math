#!/usr/bin/env python3
"""verify_B.py -- numerical verification (mpmath, >= 50 digits) of the identities
listed as OPEN in SURVEY-B.md (Zhi-Wei Sun's conjectural hypergeometric series,
with or without harmonic numbers).

Each entry of ENTRIES is a dict with keys
    id      : short id used in SURVEY-B.md
    summand : python function of the integer k (returns an mpf)
    start   : first summation index
    rhs     : function of no arguments returning the conjectured value (mpf)
    src     : short source string
Importable:  from verify_B import ENTRIES, lhs, check
Run:         python3 verify_B.py            (prints id, #terms, |lhs-rhs|)
"""
from mpmath import mp, mpf, binomial, pi, zeta, log, sqrt, catalan, nsum, inf

mp.dps = 80
GUARD_DPS = 95          # working precision while summing
TOL_EXP = 88            # stop when 6 consecutive terms are below 10^-TOL_EXP


# ---------------------------------------------------------------- helpers
def C(n, k):
    return binomial(n, k)


class _Harm:
    """H(n, m) = sum_{0<j<=n} 1/j^m by exact incremental accumulation (cached)."""

    def __init__(self):
        self.tab = {}

    def reset(self):
        self.tab = {}

    def __call__(self, n, m=1):
        t = self.tab.setdefault(m, [mpf(0)])
        while len(t) <= n:
            j = len(t)
            t.append(t[-1] + mpf(1) / mpf(j) ** m)
        return t[n]


H = _Harm()


def L3(s):
    """L(s, (-3/.)) = sum_{n>=0} (1/(3n+1)^s - 1/(3n+2)^s)"""
    return (zeta(s, mpf(1) / 3) - zeta(s, mpf(2) / 3)) / mpf(3) ** s


def beta4():
    """Dirichlet beta(4) = sum (-1)^k/(2k+1)^4"""
    return (zeta(4, mpf(1) / 4) - zeta(4, mpf(3) / 4)) / mpf(4) ** 4


def K():
    return L3(2)


G = lambda: +catalan

# kernels used repeatedly
def ram23(k):   # binom(2k,k)^2 binom(3k,k)
    return C(2 * k, k) ** 2 * C(3 * k, k)


def ram24(k):   # binom(2k,k)^2 binom(4k,2k)
    return C(2 * k, k) ** 2 * C(4 * k, 2 * k)


def ram236(k):  # binom(2k,k) binom(3k,k) binom(6k,3k)
    return C(2 * k, k) * C(3 * k, k) * C(6 * k, 3 * k)


def d6(k):      # (2k+1)(6k+1)(6k+5) binom(6k,3k)
    return (2 * k + 1) * (6 * k + 1) * (6 * k + 5) * C(6 * k, 3 * k)


def c344(k):    # binom(3k,k) binom(4k,k)^2
    return C(3 * k, k) * C(4 * k, k) ** 2


def P88(k): return 88 * k**3 + 108 * k**2 + 36 * k + 3
def P368(k): return 368 * k**3 + 400 * k**2 + 118 * k + 9
def P896(k): return 896 * k**3 + 992 * k**2 + 296 * k + 21
def P74(k): return 74 * k**3 - 84 * k**2 + 29 * k - 3
def P77(k): return 77 * k**3 - 86 * k**2 + 29 * k - 3
def P360(k): return 360 * k**3 + 612 * k**2 + 230 * k + 15
def P1344(k): return 1344 * k**3 + 944 * k**2 + 156 * k + 9
def P28(k): return 28 * k**2 + 10 * k + 1
def P21(k): return 21 * k**3 - 22 * k**2 + 8 * k - 1
def Q92(k): return 92 * k**3 - 84 * k**2 + 27 * k - 3
def R43680(k): return 43680 * k**4 + 20632 * k**3 + 4340 * k**2 + 466 * k + 21
def P344(k): return 344 * k**3 - 386 * k**2 + 115 * k - 9
def Q1952(k): return 1952 * k**3 - 1732 * k**2 + 315 * k + 12


def q(a, b=1):
    return mpf(a) / mpf(b)


E = []


def add(id_, start, summand, rhs, src):
    E.append(dict(id=id_, start=start, summand=summand, rhs=rhs, src=src))


# ======================================================================
# Section A: hypergeometric, no harmonic numbers
# ======================================================================
add("A1", 1,
    lambda k: C(2*k, k) * mpf(16)**k * (1656*k**4 - 2604*k**3 + 1158*k**2 - 229*k + 15)
    / (mpf(k)**2 * (2*k - 1)**2 * (6*k - 1) * (6*k - 5) * C(6*k, 3*k) * C(3*k, k)),
    lambda: mpf(0), "MO 506695 (2026-01-06)")
add("A2", 1,
    lambda k: (-1)**(k - 1) * (448*k**4 - 104*k**3 + 24*k**2 - 8*k + 1) * C(2*k, k)**2
    / (mpf(2*k - 1)**4 * mpf(k)**3 * C(6*k, 3*k) * C(3*k, k)),
    lambda: 20 * zeta(3), "MO 487361 eq.(3) (2025-02-07)")
add("A3", 1,
    lambda k: mpf(10*k - 3) / ((2*k - 1) * mpf(k)**2 * mpf(3)**k * C(2*k, k) * C(3*k, k)),
    lambda: K() / 2, "II v10 Conj 5.1(i) eq.1")
add("A4", 1,
    lambda k: (4*k - 1) * mpf(3)**k / ((2*k - 1) * mpf(k)**2 * C(2*k, k) * C(3*k, k)),
    lambda: 2 * K(), "II v10 Conj 5.1(i) eq.2")
add("A5", 1,
    lambda k: (8*k - 3) * C(4*k, 2*k) / (mpf(k) * (4*k - 1) * mpf(9)**k * C(2*k, k)**2),
    lambda: sqrt(3) * pi / 18, "II v10 Conj 5.9(i) eq.1; MO 456779")
add("A6", 1,
    lambda k: (13*k**2 - 15*k + 2) * mpf(9)**(k - 1) / (mpf(k) * (3*k - 1) * (3*k - 2) * C(4*k, k)),
    lambda: pi / sqrt(3), "II v10 Conj 4.17(i)")
add("A7", 0,
    lambda k: P360(k) * C(3*k, k) * C(6*k, 3*k)**2 / ((3*k + 1) * (3*k + 2) * mpf(2)**(15*k)),
    lambda: 32 * sqrt(2) / pi, "II v10 Conj 5.25(i) eq.1 (=(5.20))")
add("A8", 0,
    lambda k: P1344(k) * C(2*k, k)**2 * C(4*k, 2*k)**3 / ((2*k + 1) * mpf(2)**(22*k)),
    lambda: 64 * sqrt(2) / pi**2, "II v10 Conj 5.33(i) eq.1 (=(5.31))")
add("A9", 1,
    lambda k: P77(k) * mpf(-32)**k / (mpf(k)**3 * (3*k - 1) * (3*k - 2) * c344(k)),
    lambda: -6 * G(), "II v10 Conj 4.33(i) eq.1")
add("A10", 0,
    lambda k: (3*k + 2)**2 * mpf(16)**k / d6(k),
    lambda: pi / (2 * sqrt(3)), "III v6 Conj 3.2")

# ======================================================================
# Section B: hypergeometric with harmonic numbers
# ======================================================================
# --- K1: binom(2k,k)^2 binom(3k,k) / m^k  (Ramanujan 1/pi kernels), 2026
add("B1", 0, lambda k: ram23(k) / mpf(-192)**k * ((5*k + 1) * (H(2*k) - H(k)) + q(3*k + 1, 2*k + 1)),
    lambda: 4 * log(3) / (sqrt(3) * pi), "2603.29973v3 Conj 3.1; MO 507592")
add("B2", 0, lambda k: ram23(k) / mpf(-192)**k * ((5*k + 1) * (H(3*k) - H(k)) + q(4*k + 3, 6*k + 3)),
    lambda: 4 / (3 * sqrt(3) * pi) * log(q(64, 3)), "2603.29973v3 Conj 3.1")
add("B3", 0, lambda k: ram23(k) / mpf(-12)**(3*k) * ((51*k + 7) * (H(2*k) + H(k)) - q(3 * (13*k + 7), 2*k + 1)),
    lambda: -36 * sqrt(3) * log(3) / pi, "2603.29973v3 Conj 3.2")
# NOTE: as printed (arXiv:2603.29973v3 and MO 507592 eq.(4)) the rational term is "-(180k+93)/(2k+1)";
# that version is numerically FALSE.  PSLQ gives the corrected sign "+", used here.
add("B4", 0, lambda k: ram23(k) / mpf(-12)**(3*k) * ((51*k + 7) * (3*H(3*k) - 7*H(k)) + q(180*k + 93, 2*k + 1)),
    lambda: 36 * sqrt(3) * log(108) / pi, "2603.29973v3 Conj 3.2 (sign of rational term corrected)")
# NOTE: as printed (arXiv:2603.29973v3 and MO 507592 eq.(5)) the rational term is "-(27k+13)/(6k+3)";
# that version is numerically FALSE.  PSLQ gives the corrected sign "+", used here.
add("B5", 0, lambda k: ram23(k) / mpf(-8640)**k * ((9*k + 1) * (H(2*k) - 2*H(k)) + q(27*k + 13, 6*k + 3)),
    lambda: 48 * log(3) / (sqrt(15) * pi), "2603.29973v3 Conj 3.3 (sign of rational term corrected)")
add("B6", 0, lambda k: ram23(k) / mpf(-8640)**k * ((9*k + 1) * (3*H(3*k) - H(k)) + q(1, 6*k + 3)),
    lambda: 12 / (sqrt(15) * pi) * log(q(320, 243)), "2603.29973v3 Conj 3.3")
# --- K2: binom(2k,k)^2 binom(4k,2k) / (-12288)^k
add("B7", 0, lambda k: ram24(k) / mpf(-12288)**k * ((28*k + 3) * (H(2*k) - H(k)) + q(20*k + 9, 4*k + 2)),
    lambda: 4 * log(432) / (sqrt(3) * pi), "2603.29973v3 Conj 3.4")
add("B8", 0, lambda k: ram24(k) / mpf(-12288)**k * ((28*k + 3) * (H(4*k) - H(2*k)) + q(8*k + 5, 4*k + 2)),
    lambda: 8 / (sqrt(3) * pi) * log(q(16, 3)), "2603.29973v3 Conj 3.4")
add("B9", 0, lambda k: ram24(k) / mpf(-12288)**k
    * ((28*k + 3) * (36*H(4*k, 2) - 9*H(2*k, 2) - H(k, 2)) + q(12, 4*k + 1)),
    lambda: 16 * pi / (3 * sqrt(3)), "2603.29973v3 Conj 3.5; MO 507600 (3)")
# --- K3: binom(2k,k) binom(3k,k) binom(6k,3k) / m^k
add("B10", 0, lambda k: ram236(k) / mpf(-15)**(3*k) * ((63*k + 8) * (3*H(2*k) - 5*H(k)) + q(198*k + 83, 4*k + 2)),
    lambda: 30 * sqrt(15) * log(3) / pi, "2603.29973v3 Conj 3.6")
add("B11", 0, lambda k: ram236(k) / mpf(-96)**(3*k) * ((342*k + 25) * (3*H(2*k) - 7*H(k)) + q(1122*k + 571, 2*k + 1)),
    lambda: 64 * sqrt(6) * log(mpf(2)**7 * 3**6) / pi, "2603.29973v3 Conj 3.7")
# --- K4: binom(2k,k)^3 / 4096^k
add("B12", 0, lambda k: C(2*k, k)**3 / mpf(4096)**k
    * (9 * (42*k + 5) * (H(2*k, 4) - H(k, 4) / 16) + mpf(25) / mpf(2*k + 1)**3),
    lambda: 5 * pi**3 / 6, "II v10 Conj 5.13(i) (=I (5.14)); Au 2602.08721v4 '?='")
# --- K5: 1/(k^3 binom(2k,k)^3)
add("B13", 1, lambda k: ((21*k - 8) * (2*H(3*k - 1) - H(k - 1)) - q(4 * (27*k - 11), 3 * (3*k - 1)))
    / (mpf(k)**3 * C(2*k, k)**3),
    lambda: 3 * zeta(3), "2603.29973v3 Conj 3.8; MO 507600 (1)")
# --- K6: z^k/(k^3 binom(2k,k)^2 binom(3k,k)),  z = 8, 64, -27
def z23(k): return mpf(k)**3 * ram23(k)
add("B14", 1, lambda k: mpf(8)**k * ((10*k - 3) * (H(2*k - 1) - H(k - 1)) - 1) / z23(k),
    lambda: 7 * zeta(3) / 2, "2210.07238v9 Conj 23(i) eq.1 (=(63))")
add("B15", 1, lambda k: mpf(8)**k * ((10*k - 3) * (H(3*k - 1) - H(k - 1)) - q(8, 3)) / z23(k),
    lambda: (2 * pi**2 * log(2) + 7 * zeta(3)) / 4, "2210.07238v9 Conj 23(i) eq.2")
add("B16", 1, lambda k: mpf(8)**k * ((50*k - 15) * H(k - 1, 2) + q(4, k)) / z23(k),
    lambda: pi**4 / 24, "II v10 Conj 5.15(i)")
add("B17", 1, lambda k: mpf(64)**(k - 1) * ((11*k - 3) * (2*H(2*k - 1) + H(k - 1)) - 4) / z23(k),
    lambda: 7 * zeta(3) / 2, "2210.07238v9 Conj 26(i) eq.1 (=(71))")
add("B18", 1, lambda k: mpf(64)**(k - 1) * ((11*k - 3) * (3*H(3*k - 1) - 6*H(k - 1)) - 7) / z23(k),
    lambda: (6 * pi**2 * log(2) - 21 * zeta(3)) / 8, "2210.07238v9 Conj 26(i) eq.2")
add("B19", 1, lambda k: mpf(64)**k * ((55*k - 15) * H(k - 1, 2) + q(8, k)) / z23(k),
    lambda: 8 * pi**4 / 3, "II v10 Conj 5.17(i)")
add("B20", 1, lambda k: mpf(-27)**k * ((15*k - 4) * (H(2*k - 1) - H(k - 1)) + q(1, 2*k - 1)) / z23(k),
    lambda: -4 * pi**3 / (3 * sqrt(3)), "2210.07238v9 Conj 25(i)")
add("B21", 1, lambda k: mpf(-27)**(k - 1) * ((15*k - 4) * (9*H(3*k - 1, 2) - 5*H(k - 1, 2)) - q(18, 3*k - 1)) / z23(k),
    lambda: 9 * L3(4), "2603.29973v3 Conj 3.9; MO 507600 (2)")
# --- K7: 81^k/(k^3 binom(2k,k)^2 binom(4k,2k))
def z24(k): return mpf(k)**3 * ram24(k)
add("B22", 1, lambda k: mpf(81)**k * ((35*k - 8) * (H(4*k - 1) - H(k - 1)) - q(35, 4)) / z24(k),
    lambda: 12 * pi**2 * log(3) + 39 * zeta(3), "2210.07238v9 Conj 27(i) (=(75))")
add("B23", 1, lambda k: mpf(81)**k * ((350*k - 80) * H(k - 1, 2) + q(27, k)) / z24(k),
    lambda: 4 * pi**4, "II v10 Conj 5.19(i)")
# --- K8: 1/(k^3 (2k-1) binom(2k,k) binom(3k,k)^2), (-1)^k/(k^3 (2k-1) binom(2k,k) binom(4k,2k)^2)
add("B24", 1, lambda k: (17 * (145*k**2 - 104*k + 18) * H(k - 1, 3) + q(28 * (2*k - 1), k**2))
    / (mpf(k)**3 * (2*k - 1) * C(2*k, k) * C(3*k, k)**2),
    lambda: 528 * zeta(5) - 46 * pi**2 * zeta(3), "II v10 Conj 5.21(i)")
add("B25", 1, lambda k: (-1)**k * (59 * (410*k**2 - 197*k + 24) * H(k - 1, 3) - q(62 * (2*k - 1), k**2))
    / (mpf(k)**3 * (2*k - 1) * C(2*k, k) * C(4*k, 2*k)**2),
    lambda: 6 * (57 * pi**2 * zeta(3) - 652 * zeta(5)), "II v10 Conj 5.23(i)")
# --- K9: 3^{-k} resp. 3^k / ((2k-1) k^2 binom(2k,k) binom(3k,k))
def z9(k): return (2*k - 1) * mpf(k)**2 * C(2*k, k) * C(3*k, k)
add("B26", 1, lambda k: ((10*k - 3) * (3*H(3*k - 1) + 4*H(2*k - 1) - 6*H(k - 1)) - 12) / (z9(k) * mpf(3)**k),
    lambda: 4 * pi**3 / (27 * sqrt(3)), "II v10 Conj 5.1(i) eq.3")
add("B27", 1, lambda k: mpf(3)**k * ((4*k - 1) * (3*H(3*k - 1) - 3*H(2*k - 1) - H(k - 1)) + q(6*k, 2*k - 1)) / z9(k),
    lambda: 8 * pi**3 / (27 * sqrt(3)), "II v10 Conj 5.1(i) eq.4")
# --- K10: (-16)^k / ((2k-1) k^2 binom(2k,k) binom(4k,2k))
def z10(k): return mpf(k)**2 * C(2*k, k) * C(4*k, 2*k)
add("B28", 1, lambda k: mpf(-16)**k / z10(k) * (q(5*k - 1, 2*k - 1) * H(k - 1) - q(1, 8*k)),
    lambda: pi**2 * log(2) - 21 * zeta(3) / 4, "II v10 Conj 5.3(i) first")
add("B29", 1, lambda k: mpf(-16)**k / z10(k) * (q(5*k - 1, 2*k - 1) * H(2*k - 1) + q(1, 12*k)),
    lambda: pi**2 * log(2) / 3 - 35 * zeta(3) / 6, "II v10 Conj 5.3(i) second")
add("B30", 1, lambda k: mpf(-16)**k / ((2*k - 1) * z10(k)) * ((5*k - 1) * H(4*k - 1) - q(46*k - 5, 24*k)),
    lambda: pi**2 * log(2) / 6 - 77 * zeta(3) / 12, "II v10 Conj 5.3(i) third")
add("B31", 1, lambda k: mpf(-16)**k / ((2*k - 1) * z10(k)) * ((5*k - 1) * H(k - 1, 2) - q(4 * (4*k - 1), (2*k - 1)**2)),
    lambda: pi**4 / 6, "II v10 Conj 5.5(i) first")
add("B32", 1, lambda k: mpf(-16)**k / ((2*k - 1) * z10(k)) * ((5*k - 1) * H(k - 1, 3) + q(12 * (4*k - 1), (2*k - 1)**3)),
    lambda: -4 * pi**2 * zeta(3), "II v10 Conj 5.5(i) second")
# --- K11: (-256)^k / (k^2 (2k-1)(4k-1)(4k-3) binom(3k,k) binom(6k,3k))
def z11(k): return mpf(k)**2 * (2*k - 1) * (4*k - 1) * (4*k - 3) * C(3*k, k) * C(6*k, 3*k)
add("B33", 1, lambda k: mpf(-256)**k * (P344(k) * (H(2*k - 1) - 2*H(k - 1))
                                        - q((4*k - 1) * (58*k**2 + 181*k - 66), 26*k)) / z11(k),
    lambda: 1680 * zeta(3) / 13 - 32 * pi**2 * log(2), "II v10 Conj 5.7 first")
add("B34", 1, lambda k: mpf(-256)**k * (P344(k) * (2*H(6*k - 1) - H(3*k - 1) - 3*H(k - 1))
                                        - q(Q1952(k), 13*k)) / z11(k),
    lambda: 2464 * zeta(3) / 13 - 64 * pi**2 * log(2), "II v10 Conj 5.7 second")
# --- K12: binom(4k,2k) / (k (4k-1) 9^k binom(2k,k)^2)
def z12(k): return mpf(k) * (4*k - 1) * mpf(9)**k * C(2*k, k)**2
add("B35", 1, lambda k: C(4*k, 2*k) * ((8*k - 3) * (5*H(2*k - 1) - 4*H(k - 1)) - 6) / z12(k),
    lambda: 3 * K() / 2, "II v10 Conj 5.9(i) eq.2")
add("B36", 1, lambda k: C(4*k, 2*k) * (8*k - 3) * (2*H(2*k - 1, 2) - 5*H(k - 1, 2)) / z12(k),
    lambda: pi**3 / (36 * sqrt(3)), "II v10 Conj 5.9(i) eq.3")
# --- K13: binom(3k,k) binom(4k,k)^2 / ((3k+1)(3k+2) m^k)   and inverse (Zeilberger-type)
def z13(k, m): return (3*k + 1) * (3*k + 2) * mpf(m)**k
add("B37", 0, lambda k: c344(k) * (P88(k) * (2*H(2*k) - 3*H(k)) + 40*k**2 + 56*k + 14) / z13(k, 1024),
    lambda: 32 * log(2) / pi, "II v10 Conj 4.25 first")
add("B38", 0, lambda k: c344(k) * (P88(k) * (2*H(4*k) - 3*H(2*k)) + 40*k**2 + 32*k + 8) / z13(k, 1024),
    lambda: 16 * log(2) / pi, "II v10 Conj 4.25 second")
add("B39", 0, lambda k: c344(k) * (P368(k) * (H(2*k) + H(k)) - 128*k**2 - 136*k - 31) / z13(k, 4096),
    lambda: -64 * log(2) / pi, "II v10 Conj 4.27 first")
add("B40", 0, lambda k: c344(k) * (P368(k) * (2*H(4*k) - 5*H(2*k)) + 192*k**2 + 212*k + 46) / z13(k, 4096),
    lambda: 96 * log(2) / pi, "II v10 Conj 4.27 second")
add("B41", 0, lambda k: c344(k) * (P896(k) * (4*H(4*k) - 4*H(2*k) + 3*H(k)) - 576*k**2 - 544*k - 110) / z13(k, -2**14),
    lambda: -256 * log(2) / pi, "II v10 Conj 4.29")
def z16(k): return mpf(k)**3 * (3*k - 1) * (3*k - 2) * c344(k)
add("B42", 1, lambda k: mpf(256)**k * (P74(k) * (3*H(2*k - 1) - 4*H(k - 1)) - 57*k**2 + 33*k - 5) / z16(k),
    lambda: 24 * pi**2 * log(2) - 84 * zeta(3), "II v10 Conj 4.31(i) first")
add("B43", 1, lambda k: mpf(256)**k * (P74(k) * (4*H(4*k - 1) - 7*H(2*k - 1)) + 51*k**2 - 15*k - 1) / z16(k),
    lambda: 24 * pi**2 * log(2) - 84 * zeta(3), "II v10 Conj 4.31(i) second")
add("B44", 1, lambda k: mpf(256)**k * (P74(k) * (4*H(4*k - 1, 2) - H(2*k - 1, 2) - 4*H(k - 1, 2)) + 4*k - 4) / z16(k),
    lambda: 2 * pi**4, "II v10 Conj 4.31(i) third")
add("B45", 1, lambda k: mpf(-32)**(k - 1) * (P77(k) * (4*H(4*k - 1, 2) - H(2*k - 1, 2) - 4*H(k - 1, 2)) - 8*k + 4) / z16(k),
    lambda: 3 * beta4() / 4, "II v10 Conj 4.33(i) eq.2")
# --- K15: binom(3k,k) binom(6k,3k)^2 / ((3k+1)(3k+2) 2^{15k})
add("B46", 0, lambda k: C(3*k, k) * C(6*k, 3*k)**2 * (P360(k) * (H(2*k) - H(k)) - 180*k**2 + 36*k + 23)
    / ((3*k + 1) * (3*k + 2) * mpf(2)**(15*k)),
    lambda: 48 * sqrt(2) * log(2) / pi, "II v10 Conj 5.25(i) eq.2")
# --- K16: binom(2k,k)^2 binom(4k,2k)^3 / ((2k+1) 2^{22k})
def k16(k): return C(2*k, k)**2 * C(4*k, 2*k)**3 / ((2*k + 1) * mpf(2)**(22*k))
add("B47", 0, lambda k: k16(k) * ((3*H(4*k) - 2*H(2*k) - H(k)) * P1344(k) + 672*k**2 + 404*k + 34 + q(1, 4*k + 2)),
    lambda: 352 * sqrt(2) * log(2) / pi**2, "II v10 Conj 5.33(i) eq.2")
add("B48", 0, lambda k: k16(k) * ((32*H(4*k, 2) - 8*H(2*k, 2) - H(k, 2)) * P1344(k) + 32 * (12*k + 5)),
    lambda: 352 * sqrt(2) / 3, "II v10 Conj 5.33(i) eq.3")
# --- K17: binom(2k,k)^5 / ((6k+1)(-64)^k binom(3k,k) binom(6k,3k))
def k17(k): return C(2*k, k)**5 / ((6*k + 1) * mpf(-64)**k * C(3*k, k) * C(6*k, 3*k))
add("B49", 0, lambda k: k17(k) * (P28(k) * (10*H(2*k, 2) - 3*H(k, 2)) + 2),
    lambda: pi / 2, "II v10 Conj 5.37(i)")
add("B50", 0, lambda k: k17(k) * (P28(k) * (28*H(2*k, 3) - 3*H(k, 3)) + q(8, 2*k + 1)),
    lambda: 75 * zeta(3) / pi - 24 * G(), "II v10 Conj 5.39")
# --- K18: 256^k / (k^7 binom(2k,k)^7)
def z18(k): return mpf(k)**7 * C(2*k, k)**7
add("B51", 1, lambda k: mpf(256)**k * P21(k) * (64*H(2*k - 1, 6) - 65*H(k - 1, 6)) / z18(k),
    lambda: 31 * pi**10 / 3780, "II v10 Conj 5.41(i); MO 486518")
add("B52", 1, lambda k: mpf(256)**k * (14 * P21(k) * (H(2*k - 1) - H(k - 1)) - 63*k**2 + 44*k - 8) / z18(k),
    lambda: 93 * zeta(5) / 2 + pi**4 * log(2), "2603.29973v3 Conj 3.11")
# --- K19: 4096^k / (k^7 binom(2k,k)^6 binom(3k,k))
def z19(k): return mpf(k)**7 * C(2*k, k)**6 * C(3*k, k)
add("B53", 1, lambda k: mpf(4096)**k * (Q92(k) * (5*H(2*k - 1) - 2*H(k - 1)) - 94*k**2 + 57*k - 9) / z19(k),
    lambda: 2976 * zeta(5), "II v10 Conj 5.45 first; MO 487145")
add("B54", 1, lambda k: mpf(4096)**k * (Q92(k) * (H(3*k - 1) - 3*H(k - 1)) - q(88 * k**2, 3) + 18*k - 3) / z19(k),
    lambda: 16 * (2 * pi**4 * log(2) - 93 * zeta(5)), "II v10 Conj 5.45 second")
add("B55", 1, lambda k: mpf(4096)**k * (Q92(k) * (8*H(2*k - 1, 2) - 7*H(k - 1, 2)) - 44*k + 12) / z19(k),
    lambda: 16 * pi**6 / 3, "II v10 Conj 5.45 third")
add("B56", 1, lambda k: mpf(4096)**k * (Q92(k) * (64*H(2*k - 1, 4) + 13*H(k - 1, 4)) + q(32, k)) / z19(k),
    lambda: 248 * pi**8 / 45, "II v10 Conj 5.45 fourth")
# --- K20: binom(2k,k)^8 binom(4k,2k) / 2^{32k}
def k20(k): return C(2*k, k)**8 * C(4*k, 2*k) / mpf(2)**(32*k)
add("B57", 0, lambda k: k20(k) * (R43680(k) * (17*H(2*k, 4) - H(k, 4)) + q(1796*k + 193, 2*k + 1)),
    lambda: q(8704, 45), "II v10 Conj 5.47 first")
add("B58", 0, lambda k: k20(k) * (R43680(k) * (127*H(2*k, 6) - 2*H(k, 6)) + q(1323 * (4*k + 1), (2*k + 1)**3)),
    lambda: 126976 * pi**2 / 945, "II v10 Conj 5.47 second")
add("B59", 0, lambda k: k20(k) * (R43680(k) * (4354*H(2*k, 8) - 17*H(k, 8)) - q(4182 * (4*k + 1), (2*k + 1)**5)),
    lambda: -197632 * pi**4 / 4725, "II v10 Conj 5.47 third")
# --- K21: 1/pi^2 kernels, first-order harmonic numbers (2210.07238)
add("B60", 0, lambda k: C(2*k, k)**3 * C(3*k, k) * C(4*k, 2*k) / mpf(-24**4)**k
    * ((252*k**2 + 63*k + 5) * (4*H(4*k) + 3*H(3*k) - 7*H(k)) + 504*k + 63),
    lambda: 192 * log(24) / pi**2, "2210.07238v9 Conj 59(i)")
def k21(k): return C(2*k, k)**2 * C(3*k, k)**2 * C(6*k, 3*k)
add("B61", 0, lambda k: k21(k) / mpf(10)**(6*k) * (3 * (532*k**2 + 126*k + 9) * (H(6*k) - H(k)) + 532*k + 63),
    lambda: 1125 * log(10) / (4 * pi**2), "2210.07238v9 Conj 60(i)")
add("B62", 0, lambda k: k21(k) / mpf(-2**18)**k * (6 * (1930*k**2 + 549*k + 45) * (H(6*k) - H(k)) + 3860*k + 549),
    lambda: 6912 * log(2) / pi**2, "2210.07238v9 Conj 61(i)")
add("B63", 0, lambda k: k21(k) / mpf(-2**18 * 3**6 * 5**3)**k
    * (2 * (5418*k**2 + 693*k + 29) * (H(6*k) - H(k)) + 3612*k + 231),
    lambda: 128 * sqrt(5) / pi**2 * log(2**6 * 3**2 * 5), "2210.07238v9 Conj 62(i)")
# --- K22: z^k / ((2k+1)(6k+1)(6k+5) binom(6k,3k))   (III, Section 3)
add("B64", 0, lambda k: mpf(8)**k * ((63*k**2 + 78*k + 22) * H(3*k) + 12*k + 8) / d6(k),
    lambda: 6 * G() - 3 * pi * log(2) / 2, "III v6 Conj 3.1")
add("B65", 0, lambda k: (3*k + 2)**2 * mpf(16)**k * H(k) / d6(k),
    lambda: 5 * K() / 2 - pi * log(432) / (6 * sqrt(3)), "III v6 Conj 3.2 second")
add("B66", 0, lambda k: mpf(16)**k * ((3*k + 2)**2 * (3*H(3*k) - 2*H(2*k)) + q(2*k * (3*k + 2), 2*k + 1)) / d6(k),
    lambda: 5 * K() / 4 + pi / (6 * sqrt(3)) * log(q(27, 256)), "III v6 Conj 3.2 third")
add("B67", 0, lambda k: ((21*k**2 + 27*k + 8) * (H(2*k) + 2*H(k)) + q((3*k + 1) * (5*k + 4), 6*k + 3)) / d6(k),
    lambda: 14 * K() - 16 * log(3) * pi / (3 * sqrt(3)), "III v6 Conj 3.3 first")
add("B68", 0, lambda k: ((21*k**2 + 27*k + 8) * H(3*k) + 2*k + q(13, 9)) / d6(k),
    lambda: 8 * K() / 3 - 8 * log(3) * pi / (9 * sqrt(3)), "III v6 Conj 3.3 second")
add("B69", 0, lambda k: mpf(27)**k * ((111*k**2 + 133*k + 36) * (H(2*k) - 3*H(k)) - q((3*k + 1) * (3*k + 4), 2*k + 1)) / d6(k),
    lambda: 16 * log(3) * pi / sqrt(3) - 48 * K(), "III v6 Conj 3.4 first")
add("B70", 0, lambda k: mpf(27)**k * ((111*k**2 + 133*k + 36) * H(3*k) + 30*k + 19) / d6(k),
    lambda: 12 * K(), "III v6 Conj 3.4 second")
# --- K23: (-1)^k / (k^2 binom(2k,k) binom(3k,k))    [Sun (I) eq. (4.1); Au 2602.08721v4 '?=']
add("B71", 1, lambda k: (-1)**k / (mpf(k)**2 * C(2*k, k) * C(3*k, k)) * (q(7*k - 2, 2*k - 1) * H(k - 1) - q(1, 6*k)),
    lambda: pi**2 * log(2) / 6 - 11 * zeta(3) / 12, "Sun (I) (4.1) as quoted in Au 2602.08721v4")
# --- K24: nonlinear in harmonic numbers, 1/(k^2 binom(2k,k))
add("B72", 1, lambda k: H(k - 1, 2) / (mpf(k)**2 * C(2*k, k)) * (H(2*k) + 26*H(k) - q(8, k)),
    lambda: (34 * pi**2 * zeta(3) - 287 * zeta(5)) / 54, "2603.29973v3 Conj 3.10 (TeX has '-' misprinted for '='; value fixed by PSLQ)")

# --- K25: z^k / (k (3k-1)(3k-2) binom(4k,k))   (II, Section 4; single binomial coefficient)
def z25(k): return mpf(k) * (3*k - 1) * (3*k - 2) * C(4*k, k)
def P95(k): return 95*k**2 - 84*k + 16
def P13(k): return 13*k**2 - 15*k + 2
def P77q(k): return 77*k**2 - 53*k + 10
def P415(k): return 415*k**2 - 343*k + 62
add("B73", 1, lambda k: q(9, 8)**(k - 1) * (P95(k) * (H(4*k - 1) - H(k - 1)) - 25*k + 12) / z25(k),
    lambda: log(3) * pi / sqrt(3) + 15 * K() / 4, "II v10 Conj 4.15(ii)")
add("B74", 1, lambda k: mpf(9)**(k - 1) * (P13(k) * (3*H(2*k - 1) - 4*H(k - 1)) + 12 * (4*k - 1)) / z25(k),
    lambda: 27 * K() / 4, "II v10 Conj 4.19 eq.1")
add("B75", 1, lambda k: mpf(9)**(k - 1) * (P13(k) * (8*H(4*k - 1) - 9*H(2*k - 1)) - 160*k + 84) / z25(k),
    lambda: 4 * log(3) * pi / sqrt(3) - 93 * K() / 4, "II v10 Conj 4.19 eq.2")
add("B76", 1, lambda k: (P77q(k) * (H(2*k - 1) - H(k - 1)) - q(209*k**2 - 113*k + 18, 2*k)) / (z25(k) * mpf(-2)**k),
    lambda: 3 * log(2)**2, "II v10 Conj 4.21 eq.1")
add("B77", 1, lambda k: (P77q(k) * (2*H(4*k - 1) - H(2*k - 1) - H(k - 1)) + 11*k - 27 + q(6, k)) / (z25(k) * mpf(-2)**k),
    lambda: -pi**2 / 2, "II v10 Conj 4.21 eq.2")
add("B78", 1, lambda k: (P415(k) * (H(4*k - 1) - H(k - 1)) - q(581*k**2 - 229*k - 6, 4*k)) / (z25(k) * mpf(-8)**k),
    lambda: -pi**2 / 4, "II v10 Conj 4.23")

# --- K26: binom(4k,k) / ((3k+1)(3k+2) 16^k)  with second-order harmonic numbers
add("B79", 0, lambda k: C(4*k, k) * ((11*k**2 + 8*k + 1) * (H(2*k, 2) + q(3, 4) * H(k, 2)) + q(4*k + 1, 2*k + 1))
    / ((3*k + 1) * (3*k + 2) * mpf(16)**k),
    lambda: pi**2 / 6, "II v10 Conj 4.3(i) eq.4 (the only part of Conj 4.3 not in Jiang-Sun 2601.18645)")

ENTRIES = E


# ---------------------------------------------------------------- summation
def lhs(entry, maxterms=20000):
    """Direct summation; returns (value, number_of_terms_used)."""
    old = mp.dps
    mp.dps = GUARD_DPS
    H.reset()
    try:
        s = mpf(0)
        small = 0
        k = entry["start"]
        n = 0
        tol = mpf(10) ** (-TOL_EXP)
        while n < maxterms:
            t = entry["summand"](k)
            s += t
            n += 1
            k += 1
            if abs(t) < tol:
                small += 1
                if small >= 6 and n > 12:
                    break
            else:
                small = 0
        return +s, n
    finally:
        mp.dps = old


def check(entry):
    v, n = lhs(entry)
    old = mp.dps
    mp.dps = GUARD_DPS
    r = entry["rhs"]()
    d = abs(v - r)
    mp.dps = old
    return v, r, d, n


if __name__ == "__main__":
    import sys
    only = set(sys.argv[1:])
    bad = 0
    print("# id   terms   |lhs-rhs|      digits   source")
    for e in ENTRIES:
        if only and e["id"] not in only:
            continue
        v, r, d, n = check(e)
        if d == 0:
            digs = GUARD_DPS
        else:
            scale = max(abs(r), mpf(1))
            digs = int(-mp.log10(d / scale))
        flag = "" if digs >= 50 else "   <-- FAIL (<50 digits)"
        if digs < 50:
            bad += 1
        print("%-4s %6d   %s   %3d   %s%s" % (e["id"], n, mp.nstr(d, 5), digs, e["src"], flag))
    print("# entries: %d   failures: %d   (mp.dps=%d, working dps=%d)" % (len(ENTRIES), bad, mp.dps, GUARD_DPS))
