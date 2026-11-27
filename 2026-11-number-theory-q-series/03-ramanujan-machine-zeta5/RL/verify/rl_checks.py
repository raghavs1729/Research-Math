"""Checks for ../THEOREM.md  (Table 3 of arXiv:2412.12361v2, "The Ramanujan Library").
C[f] := 1 + f_1/(1 + f_2/(1 + ...)),  R2 := sqrt(pi e/2) erf(1/sqrt 2) = 1F1(1; 3/2; 1/2) = sum_j 1/(2j+1)!!.
 (1) the eight printed rows numerically (40 digits): six agree, two do not (misprints) -- corrected values by PSLQ;
 (2) coefficient identity behind  J_n + 2(n+1) J_{n-1} = 2(n+k) J_{n-2}   (sympy);
 (3) m_n = (-1/2)^n J_n/(n+1)!,  J_n = sum_j Gamma(n+k+2+j)/(2^j j! Gamma(n+j+7/2)),  solves
     u_n = u_{n-1} + f_n u_{n-2},  f_n = (n+k)/(2n(n+1))   (numerically, several k, n = 1..12);
 (4) Theorem RL:  C[(n+k)/(2n(n+1))] = 1 + ((k+1)/5) M(k+2; 7/2; 1/2)/M(k+1; 5/2; 1/2)   (k = 0..6 and non-integer k);
 (5) contiguity (exact, Fractions): for integer k >= 0 the right side is (alpha R2 + beta)/(gamma R2 + delta);
     table for k = 0..8, compared with (1) and (4)."""
from mpmath import mp, mpf, pi, e, erf, sqrt, hyp1f1, gamma, nsum, inf, pslq, nstr, factorial
from fractions import Fraction as Fr
import sympy as sp
mp.dps = 40
def C(f, N=6000):
    v = mpf(1)
    for n in range(N, 0, -1): v = 1 + f(n) / v
    return v
R = sqrt(pi * e / 2) * erf(1 / sqrt(2))
assert abs(R - hyp1f1(1, mpf(3) / 2, mpf(1) / 2)) < mpf(10)**-35
rows = [("C[(1-2n)/(4n(n+1))]", lambda n: mpf(1 - 2 * n) / (4 * n * (n + 1)), (R + 1) / (2 * R), "(R2+1)/(2R2)"),
        ("C[(1-2n)/(4(n+1)(n+2))]", lambda n: mpf(1 - 2 * n) / (4 * (n + 1) * (n + 2)), (2 * R + 1) / 4, "(2R2+1)/4"),
        ("C[(1-2n)/(4(n+2)(n+3))]", lambda n: mpf(1 - 2 * n) / (4 * (n + 2) * (n + 3)), (R + 1) / (6 * R - 1), "(R2+1)/(6R2-1)"),
        ("C[n/(2n(n+1))]", lambda n: mpf(n) / (2 * n * (n + 1)), 1 / (2 * R - 1), "1/(2R2-1)"),
        ("C[(n+1)/(2n(n+1))]", lambda n: mpf(n + 1) / (2 * n * (n + 1)), R, "R2"),
        ("C[(n+2)/(2n(n+1))]", lambda n: mpf(n + 2) / (2 * n * (n + 1)), (2 * R + 1) / (R + 1), "(2R2+1)/(R2+1)")]
print("(1) R2 = %s" % nstr(R, 30))
bad = []
for name, f, val, txt in rows:
    c = C(f); ok = abs(c - val) < mpf(10)**-30
    rel = pslq([c * R, c, R, mpf(1)], maxcoeff=1000, maxsteps=10000)
    print("    %-26s = %s   printed %-16s %s   relation [cR, c, R, 1]: %s" % (name, nstr(c, 25), txt, "OK" if ok else "FALSE (printed value %s)" % nstr(val, 12), rel))
    if not ok: bad.append(name)
assert bad == ["C[(1-2n)/(4(n+2)(n+3))]", "C[n/(2n(n+1))]"]
c3 = C(rows[2][1]); c4 = C(rows[3][1])
assert abs(c3 - (R + 1) / (6 * R - 6)) < mpf(10)**-30 and abs(c4 - 1 / (2 * R - 2)) < mpf(10)**-30
print("    corrected: C[(1-2n)/(4(n+2)(n+3))] = (R2+1)/(6R2-6),  C[n/(2n(n+1))] = 1/(2R2-2)   (30 digits)")
# (2)
J, n, k = sp.symbols('J n k')
lhs = 4 * J * (J - 1) + 4 * (n + 1) * J - 2 * (n + k)
rhs = -2 * (n + k + J) + 4 * J * (n + J + sp.Rational(1, 2))
assert sp.expand(lhs - rhs) == 0
print("(2) PASS: coefficient identity")
# (3), (4)
def Jn(nn, kk):       # = sum_j Gamma(n+k+2+j)/(2^j j! Gamma(n+j+7/2)), summed directly (positive terms, ratio ~ 1/(2j))
    tot = mpf(0); j = 0; term = gamma(nn + kk + 2) / gamma(nn + mpf(7) / 2)
    while True:
        tot += term
        term = term * (nn + kk + 2 + j) / (2 * (j + 1) * (nn + j + mpf(7) / 2)); j += 1
        if term < mpf(10)**(-mp.dps - 5) * tot: break
    return tot
mp.dps = 30
for kk in (0, 1, 2, 5, mpf(3) / 10, mpf(3) / 2, mpf(5) / 2, -mpf(1) / 2):
    m = {nn: (-mpf(1) / 2)**nn * Jn(nn, kk) / factorial(nn + 1) for nn in range(-1, 13)}
    worst = max(abs(m[nn] - m[nn - 1] - (nn + mpf(kk)) / (2 * nn * (nn + 1)) * m[nn - 2]) / abs(m[nn]) for nn in range(1, 13))
    cf = C(lambda t: (t + mpf(kk)) / (2 * t * (t + 1)))
    th = 1 + (mpf(kk) + 1) / 5 * hyp1f1(kk + 2, mpf(7) / 2, mpf(1) / 2) / hyp1f1(kk + 1, mpf(5) / 2, mpf(1) / 2)
    alt = 1 - m[0] / m[-1]
    print("(3,4) k = %-4s recurrence residual %s ;  C = %s ;  C - formula = %s ;  C - (1 - m_0/m_{-1}) = %s" % (nstr(kk, 3), nstr(worst, 2), nstr(cf, 20), nstr(cf - th, 2), nstr(cf - alt, 2)))
    assert worst < mpf(10)**-20 and abs(cf - th) < mpf(10)**-22 and abs(cf - alt) < mpf(10)**-22
# (5) contiguity:  M(a) := M(a;3/2;1/2) = x_a + y_a R2,  M(0) = 1, M(1) = R2,  a M(a+1) = (3/2 - a) M(a-1) + (2a - 1) M(a);
#     M(a;5/2;1/2) = 3 (M(a) - M(a-1));   M(a;7/2;1/2) = 5 (M(a;5/2) - M(a-1;5/2)).
KM = 10
M = {0: (Fr(1), Fr(0)), 1: (Fr(0), Fr(1))}
for a in range(1, KM + 3):
    M[a + 1] = tuple(((Fr(3, 2) - a) * M[a - 1][i] + (2 * a - 1) * M[a][i]) / a for i in range(2))
sub = lambda u, v: (u[0] - v[0], u[1] - v[1])
N5 = {a: tuple(3 * t for t in sub(M[a], M[a - 1])) for a in range(1, KM + 3)}
O7 = {a: tuple(5 * t for t in sub(N5[a], N5[a - 1])) for a in range(2, KM + 3)}
mp.dps = 40
for a in range(1, 6):
    assert abs(hyp1f1(a, mpf(5) / 2, mpf(1) / 2) - (mpf(N5[a][0].numerator) / N5[a][0].denominator + mpf(N5[a][1].numerator) / N5[a][1].denominator * R)) < mpf(10)**-30
assert abs(C(lambda t: (t + mpf(3) / 2) / (2 * t * (t + 1))) - mpf(3) / 2) < mpf(10)**-35 and abs(C(lambda t: (t + mpf(5) / 2) / (2 * t * (t + 1))) - mpf(5) / 3) < mpf(10)**-35
print("    remark: C = 3/2 for k = 3/2 and C = 5/3 for k = 5/2 (35 digits)")
print("(5) C[(n+k)/(2n(n+1))] = (alpha R2 + beta)/(gamma R2 + delta):")
from math import gcd
from functools import reduce
for kk in range(0, 9):
    num = (N5[kk + 1][0] + Fr(kk + 1, 5) * O7[kk + 2][0], N5[kk + 1][1] + Fr(kk + 1, 5) * O7[kk + 2][1])     # beta, alpha
    den = N5[kk + 1]
    L = reduce(lambda u, v: u * v // gcd(u, v), [t.denominator for t in num + den])
    ints = [int(t * L) for t in (num[1], num[0], den[1], den[0])]
    g = reduce(gcd, ints); ints = [t // g for t in ints]
    if ints[2] < 0 or (ints[2] == 0 and ints[3] < 0): ints = [-t for t in ints]
    al, be, ga, de = ints
    val = (al * R + be) / (ga * R + de); cf = C(lambda t: mpf(t + kk) / (2 * t * (t + 1)))
    print("    k = %d: (%d R2 + %d)/(%d R2 + %d) = %s ;  C - value = %s" % (kk, al, be, ga, de, nstr(val, 25), nstr(cf - val, 2)))
    assert abs(cf - val) < mpf(10)**-30
print("PASS")
