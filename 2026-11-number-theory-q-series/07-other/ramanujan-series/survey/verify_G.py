#!/usr/bin/env python3
"""verify_G.py -- Guillera's conjectured values of half-shifted Ramanujan series
(J. Guillera, "Ramanujan series with a shift", arXiv:1802.02023v4, J. Aust. Math. Soc. 107 (2019) 367-380,
Tables 1-6: levels 2, 3, 1; the paper proves only level 4).  One of them (level 2, z = 1/9) is Sun's
Conjecture 5.9(i) = our Theorem 2.

For a Ramanujan series   sum_n  [(1/2)_n (s1)_n (s2)_n / n!^3] lam (a + b n) z^n = 1/pi   (a, b integers,
lam the algebraic normalisation, computed here numerically from this very identity) the shifted series is

    sum_{n>=0} [(1)_n (s1+1/2)_n (s2+1/2)_n / (3/2)_n^3] (2a + b + 2b n) z^n  =  2 V / (lam sqrt|z|),

V = Guillera's tabulated value of G_l(1/2, q) (resp. -i G_l for z < 0).
Entry format as in verify_B.py (id, start, summand, rhs, src).  Run: python3 verify_G.py > verify_G.out.txt
"""
from mpmath import mp, mpf, pi, log, sqrt, atan, rf, nsum, inf, factorial

mp.dps = 80
LEV = {2: (mpf(1) / 4, mpf(3) / 4), 3: (mpf(1) / 3, mpf(2) / 3), 1: (mpf(1) / 6, mpf(5) / 6)}
ln = log
ROWS = [
    # id, level, z, a, b, V
    ('G1', 2, mpf(-1) / 4, 3, 20, lambda: ln(2)),
    ('G2', 2, -mpf(256) / 3969, 8, 65, lambda: ln(88 + 13 * sqrt(7)) - 4 * ln(3)),
    ('G3', 2, -mpf(1) / 48, 3, 28, lambda: mpf(3) / 2 * ln(3) - 2 * ln(2)),
    ('G4', 2, -mpf(1) / 324, 23, 260, lambda: 2 * ln(3) - 3 * ln(2)),
    ('G5', 2, -mpf(1) / 25920, 41, 644, lambda: 9 * ln(2) - 2 * ln(3) - mpf(5) / 2 * ln(5)),
    ('G6', 2, -mpf(1) / 777924, 1123, 21460, lambda: ln(2) + 10 * ln(3) - 6 * ln(7)),
    ('G7', 2, mpf(32) / 81, 2, 14, lambda: pi / 2 - 2 * atan(1 / (2 * sqrt(2)))),
    ('G8', 2, mpf(1) / 9, 1, 8, lambda: pi / 6),
    ('G9', 2, mpf(1) / 81, 4, 40, lambda: pi / 2 - 4 * atan(1 / (2 * sqrt(2)))),   # the table prints '+': numerically false; '-' is correct
    ('G10', 2, mpf(1) / 2401, 27, 360, lambda: -pi / 6 + 4 * atan(1 / (4 * sqrt(3)))),
    ('G11', 2, mpf(1) / 9801, 19, 280, lambda: -pi / 2 + 4 * atan(7 / (5 * sqrt(11)))),
    ('G12', 2, mpf(1) / 96059601, 4412, 105560, lambda: 13 * pi / 2 - 16 * atan(1 / sqrt(2)) - 24 * atan(sqrt(2) / 3)),
    ('G13', 3, -mpf(9) / 16, 1, 5, lambda: sqrt(3) / 4 * (3 * ln(3) - 2 * ln(2))),
    ('G14', 3, -mpf(1) / 16, 7, 51, lambda: 3 * sqrt(3) / 4 * (2 * ln(2) - ln(3))),
    ('G15', 3, -mpf(1) / 80, 1, 9, lambda: sqrt(3) / 4 * (9 * ln(3) - 2 * ln(2) - 5 * ln(5))),
    ('G16', 3, -mpf(1) / 1024, 106, 1230, lambda: 3 * sqrt(3) / 4 * (8 * ln(2) - 5 * ln(3))),
    ('G17', 3, -mpf(1) / 3024, 26, 330, lambda: sqrt(3) / 4 * (7 * ln(7) - 10 * ln(2) - 6 * ln(3))),
    ('G18', 3, -mpf(1) / 250000, 827, 14151, lambda: 3 * sqrt(3) / 4 * (6 * ln(5) - 6 * ln(2) - 5 * ln(3))),
    ('G19', 3, mpf(1) / 2, 1, 6, lambda: sqrt(3) / 4 * (3 * pi - 12 * atan(sqrt(2) / 2))),
    ('G20', 3, mpf(2) / 27, 2, 15, lambda: sqrt(3) / 4 * (5 * pi - 24 * atan(sqrt(2) / 2))),
    ('G21', 3, mpf(4) / 125, 4, 33, lambda: sqrt(3) / 4 * (-3 * pi + 12 * atan(sqrt(5) / 2))),
    ('G22', 1, -mpf(64) / 125, 8, 63, lambda: 3 * sqrt(3) / 8 * ln(mpf(27) / 5)),
    ('G23', 1, -mpf(27) / 512, 15, 154, lambda: 3 * sqrt(3) / 8 * ln(2)),
    ('G24', 1, -mpf(1) / 512, 25, 342, lambda: 3 * sqrt(3) / 8 * ln(mpf(32) / 27)),
    ('G25', 1, -mpf(9) / 64000, 279, 4554, lambda: 3 * sqrt(3) / 8 * ln(mpf(135) / 128)),
    ('G26', 1, -mpf(1) / 512000, 263, 5418, lambda: 3 * sqrt(3) / 8 * ln(mpf(4 * 3**9) / 5**7)),
    ('G27', 1, -mpf(1) / 440**3, 10177, 261702, lambda: 3 * sqrt(3) / 8 * ln(mpf(2**13 * 11**5) / (27 * 5**11))),
    ('G28', 1, -mpf(1) / 53360**3, 13591409, 545140134,
     lambda: 3 * sqrt(3) / 8 * ln(mpf(3**21 * 5**13 * 29**5) / (mpf(2)**38 * 23**11))),
    ('G29', 1, mpf(27) / 125, 3, 28, lambda: 3 * sqrt(3) / 8 * (pi - 4 * atan(mpf(1) / 2))),
    ('G30', 1, mpf(4) / 125, 1, 11, lambda: 3 * sqrt(3) / 8 * (-pi + 8 * atan(mpf(1) / 2))),
    ('G31', 1, mpf(8) / 1331, 5, 63, lambda: 3 * sqrt(3) / 8 * (3 * pi - 4 * atan(sqrt(2) / 3) - 12 * atan(sqrt(2) / 2))),
    ('G32', 1, mpf(64) / 614125, 8, 133, lambda: 3 * sqrt(3) / 8 * (3 * pi - 16 * atan(mpf(1) / 2) - 8 * atan(mpf(1) / 4))),
]


def lam(level, z, a, b):
    s1, s2 = LEV[level]
    S = nsum(lambda n: rf(mpf(1) / 2, n) * rf(s1, n) * rf(s2, n) / factorial(n) ** 3 * (a + b * n) * z ** n, [0, inf])
    return 1 / (pi * S)


E = []


def add(id_, start, summand, rhs, src):
    E.append(dict(id=id_, start=start, summand=summand, rhs=rhs, src=src))


def _mk(level, z, a, b, V):
    s1, s2 = LEV[level]
    half = mpf(1) / 2
    def summand(k):
        return rf(mpf(1), k) * rf(s1 + half, k) * rf(s2 + half, k) / rf(mpf(3) / 2, k) ** 3 * (2 * a + b + 2 * b * k) * z ** k
    def rhs():
        return 2 * V() / (lam(level, z, a, b) * sqrt(abs(z)))
    return summand, rhs


for (i, lev, z, a, b, V) in ROWS:
    sm, rh = _mk(lev, z, a, b, V)
    add(i, 0, sm, rh, "Guillera 1802.02023 level %d, z = %s, (%d + %d n)" % (lev, str(z)[:16], a, b))
ENTRIES = E

if __name__ == "__main__":
    bad = 0
    for e in ENTRIES:
        S = nsum(e['summand'], [0, inf])
        r = e['rhs']()
        d = abs(S - r)
        digs = 80 if d == 0 else int(-mp.log10(d / max(abs(r), 1)))
        if digs < 50:
            bad += 1
        print("%-4s |lhs-rhs| = %s  digits %3d  rhs = %s   %s%s" % (e['id'], mp.nstr(d, 5), digs, mp.nstr(r, 20), e['src'],
                                                                 "" if digs >= 50 else "   <-- FAIL"))
    print("# entries: %d   failures: %d" % (len(ENTRIES), bad))
