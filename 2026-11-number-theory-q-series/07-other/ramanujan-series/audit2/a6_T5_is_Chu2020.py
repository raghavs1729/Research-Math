#!/usr/bin/env python3
"""
audit2 / a6: Theorem 5 (Sun, arXiv:2307.03086v10, Conj. 4.33(i), dated 2025-02-28) is, term by term, the identity

    2G = sum_{n>=0} (-1/32)^n (1)_n (1)_n (1)_{3n} / ((3/2)_n (5/2)_{2n} (5/2)_{2n}) * (17 + 88n + 145n^2 + 77n^3)/9

printed as equation (7.10) (Example 7.3) in
    W. Chu, Infinite series formulae related to Gauss and Bailey 2F1(1/2)-sums,
    Acta Math. Sci. 40B (2020), no. 2, 293-315, doi:10.1007/s10473-020-0201-y.
(Chu's W_k(1,1,1) = (17+88k+145k^2+77k^3)/((1+2k)(3+4k)^2); Theorem 7.2, eq. (7.6), a = c = 1.)

Here: exact check that  Chu's n-th term  ==  -(1/3) * (Sun's k-th term at k = n+1),  for n < 400,
plus the polynomial identity P(n+1) = 77n^3+145n^2+88n+17 for P(k) = 77k^3-86k^2+29k-3,
plus a symbolic check of the term ratio (so the equality holds for all n), plus a 60-digit numerical check
of Chu's series against 2G.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import flint
from mpmath import mp, mpf, catalan

def poch(a, m):
    r = Fr(1)
    for i in range(m): r *= a + i
    return r
def chu(n):
    return Fr(-1, 32)**n * factorial(n)**2 * factorial(3*n) / (poch(Fr(3, 2), n) * poch(Fr(5, 2), 2*n)**2) * Fr(17 + 88*n + 145*n**2 + 77*n**3, 9)
def sun(k):
    return Fr((77*k**3 - 86*k**2 + 29*k - 3) * (-32)**k, k**3*(3*k-1)*(3*k-2)*comb(3*k, k)*comb(4*k, k)**2)

ok1 = all(chu(n) == -sun(n+1)/3 for n in range(400))
print('Chu (7.10) n-th term == -(1/3) * Sun Conj. 4.33(i) term at k = n+1, for n = 0..399 (exact):', ok1)

ctx = flint.fmpq_mpoly_ctx.get(['n'], 'lex'); n, = ctx.gens()
P = lambda k: 77*k**3 - 86*k**2 + 29*k - 3
ok2 = (P(n+1) == 77*n**3 + 145*n**2 + 88*n + 17)
print('P(n+1) = 77n^3+145n^2+88n+17:', ok2)
# term ratios as rational functions of n (cross-multiplied): chu(n+1)/chu(n) == sun(n+2)/sun(n+1)
half = flint.fmpq(1, 2)
C = lambda m: 77*m**3 + 145*m**2 + 88*m + 17
chu_num = -(n+1)**2*(3*n+1)*(3*n+2)*(3*n+3)*C(n+1)
chu_den = 32*(n+3*half)*((2*n+5*half)*(2*n+7*half))**2*C(n)
# sun: t(k+1)/t(k) with k = n+1;  C(3k,k)->(3k+1)(3k+2)(3k+3)/((k+1)(2k+1)(2k+2)),  C(4k,k)->(4k+1)(4k+2)(4k+3)(4k+4)/((k+1)(3k+1)(3k+2)(3k+3))
k = n + 1
sun_num = -32*P(k+1)*k**3*(3*k-1)*(3*k-2) * (k+1)*(2*k+1)*(2*k+2) * ((k+1)*(3*k+1)*(3*k+2)*(3*k+3))**2
sun_den = P(k)*(k+1)**3*(3*k+2)*(3*k+1) * (3*k+1)*(3*k+2)*(3*k+3) * ((4*k+1)*(4*k+2)*(4*k+3)*(4*k+4))**2
ok3 = (chu_num*sun_den == sun_num*chu_den)
print('term ratios agree identically in n (flint, cross-multiplied):', ok3)

mp.dps = 70
s = mpf(0)
for m in range(60):
    c = chu(m); s += mpf(c.numerator)/c.denominator
ok4 = abs(s - 2*catalan) < mpf(10)**(-60)
print('sum of Chu (7.10) (60 terms) - 2G =', mp.nstr(s - 2*catalan, 5), '->', ok4)
print('CONCLUSION: T5 == Chu 2020, eq. (7.10), reindexed k = n+1 and multiplied by -3:', ok1 and ok2 and ok3 and ok4)
