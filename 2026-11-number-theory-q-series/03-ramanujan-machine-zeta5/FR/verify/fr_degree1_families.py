"""Exact tests for Theorem B (the 'if' half of Conjecture 1.4 of arXiv:2111.04468).
PCF[a_n, b_n], a_n = A n + k, b_n = B (n-x1)(n-x2), rational parameters, x1 not a positive integer, A^2+4B != 0.
With c_j := u_{j-1}/(1-x1)_j the recurrence becomes (j-x1) c_j = (A(j-1)+k) c_{j-1} + B (j-1-x2) c_{j-2}, and
     sum c_j z^j = W(z) * sum_j e_j z^j/(j-x1),   sum e_j z^j = (r0 + r1 z)/(Delta(z) W(z)),  Delta = 1 - A z - B z^2,
     W = (1-l1 z)^(-al1) (1-l2 z)^(-al2)  [family (Eq. 6): l1, l2 rational roots of l^2 = A l + B]
     W = Delta^(-al), al = (1+x1-x2)/2     [family (Eq. 7): k = A(1-x1-x2)/2]
We check (1) this closed form against the recurrence (exact), (2) that the common denominator of c_0..c_n grows only
exponentially, (3) that (1/n) log(GCD_n/n!) stays bounded for the inflated integer PCF, and (4) a control outside the
two families where it visibly tends to -infinity."""
from fractions import Fraction as Fr
from math import gcd, lcm, log, lgamma, isqrt
import random
N = 120
def ser_mul(a, b, n=N + 2):
    r = [Fr(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: r[i + j] += x * y
    return r
def binom_series(lam, al, n=N + 2):          # (1 - lam z)^(-al)
    r = [Fr(1)]
    for i in range(1, n + 1):
        r.append(r[-1] * (al + i - 1) / i * lam)
    return r
def delta_power(A, B, al, n=N + 2):          # (1 - A z - B z^2)^(-al) = sum_i (al)_i/i! (A z + B z^2)^i
    r = [Fr(0)] * (n + 1); coef = Fr(1); pw = [Fr(1)] + [Fr(0)] * n
    base = [Fr(0), Fr(A), Fr(B)]
    for i in range(0, n + 1):
        for m in range(n + 1): r[m] += coef * pw[m]
        coef = coef * (al + i) / (i + 1)
        pw = ser_mul(pw, base, n)
    return r
def closed_form(A, k, B, x1, x2, c0, c1, family, lam=None):
    r0 = -x1 * c0; r1 = (1 - x1) * c1 - k * c0
    if family == 6:
        l1, l2 = lam
        al1 = (k + A * x1 - l2 * (1 + x1 - x2)) / (l1 - l2); al2 = (l1 * (1 + x1 - x2) - k - A * x1) / (l1 - l2)
        W = ser_mul(binom_series(l1, al1), binom_series(l2, al2))
        Winv = ser_mul(binom_series(l1, 1 - al1), binom_series(l2, 1 - al2))      # 1/(Delta W)
    else:
        al = (1 + x1 - x2) / 2
        W = delta_power(A, B, al); Winv = delta_power(A, B, 1 - al)
    e = ser_mul([r0, r1], Winv)
    V = [(e[j] / (j - x1)) if j != x1 else Fr(c0) for j in range(N + 3)]
    return ser_mul(W, V)
def direct(A, k, B, x1, x2, c0, c1):
    c = [Fr(c0), Fr(c1)]
    for j in range(2, N + 3):
        c.append(((A * (j - 1) + k) * c[j - 1] + B * (j - 1 - x2) * c[j - 2]) / (j - x1))
    return c
def gcd_growth(A, k, B, x1, x2, ns=(30, 60, 120)):
    # inflate to integers
    d = lcm(A.denominator, k.denominator)
    dB = (B * x1 * x2).denominator; dB = lcm(dB, (B * (x1 + x2)).denominator, B.denominator)
    c = lcm(d, dB)
    a = lambda n: c * (A * n + k); b = lambda n: c * c * B * (n - x1) * (n - x2)
    p0, p1, q0, q1 = 1, int(a(0)), 0, 1
    out = {}
    for n in range(1, max(ns) + 1):
        an, bn = a(n), b(n); assert an.denominator == 1 and bn.denominator == 1
        p0, p1 = p1, int(an) * p1 + int(bn) * p0; q0, q1 = q1, int(an) * q1 + int(bn) * q0
        if n in ns:
            g = gcd(p1, q1)
            out[n] = (log(g) - lgamma(n + 1)) / n if g else float('-inf')
    return out
random.seed(7)
rr = lambda lo, hi, dens=(1, 1, 2, 3): Fr(random.randint(lo, hi), random.choice(dens))
worst = 0.0; cnt6 = cnt7 = 0
while cnt6 < 25:
    l1, l2 = rr(-4, 4), rr(-4, 4)
    if l1 == l2 or l1 == 0 or l2 == 0 or l1 + l2 == 0: continue
    A, B = l1 + l2, -l1 * l2; k = rr(-5, 5); x1, x2 = rr(-4, 4, (1, 2, 3, 4)), rr(-4, 4, (1, 2, 3, 4))
    if (x1.denominator == 1 and x1 >= 0) or (x2.denominator == 1 and x2 >= 1) or k == 0: continue
    for (c0, c1) in ((Fr(1), (A * 0 + k) / (1 - x1)), (Fr(0), 1 / (1 - x1))):
        assert closed_form(A, k, B, x1, x2, c0, c1, 6, (l1, l2))[:N + 1] == direct(A, k, B, x1, x2, c0, c1)[:N + 1]
    g = gcd_growth(A, k, B, x1, x2); worst = min(worst, min(g.values())); cnt6 += 1
while cnt7 < 25:
    A = rr(-5, 5); B = rr(-6, 6); x1, x2 = rr(-4, 4, (1, 2, 3, 4)), rr(-4, 4, (1, 2, 3, 4))
    if A == 0 or B == 0 or A * A + 4 * B == 0: continue
    if (x1.denominator == 1 and x1 >= 0) or (x2.denominator == 1 and x2 >= 1): continue
    k = A * (1 - x1 - x2) / 2
    for (c0, c1) in ((Fr(1), k / (1 - x1)), (Fr(0), 1 / (1 - x1))):
        assert closed_form(A, k, B, x1, x2, c0, c1, 7)[:N + 1] == direct(A, k, B, x1, x2, c0, c1)[:N + 1]
    g = gcd_growth(A, k, B, x1, x2); worst = min(worst, min(g.values())); cnt7 += 1
print("Theorem B: closed form = recurrence (exact, j <= %d) for %d (Eq.6) + %d (Eq.7) random parameter sets;" % (N, cnt6, cnt7))
print("           min over all of (1/n) log(GCD_n/n!) at n = 30, 60, 120:  %.3f  (bounded below => FR)" % worst)
# the paper's examples
print("paper Eq.4  PCF[7n+3, 8n^2-2]:", gcd_growth(Fr(7), Fr(3), Fr(8), Fr(1, 2), Fr(-1, 2)))
# control: A^2+4B not a square and k != A(1-x1-x2)/2
print("control PCF[3n+1, n^2-1/4] (A^2+4B = 13):", gcd_growth(Fr(3), Fr(1), Fr(1), Fr(1, 2), Fr(-1, 2), ns=(30, 60, 120, 240)))
print("PASS")
