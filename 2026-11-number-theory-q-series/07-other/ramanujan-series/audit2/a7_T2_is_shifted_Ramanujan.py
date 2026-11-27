#!/usr/bin/env python3
"""
audit2 / a7: Theorem 2 (Sun, Conj. 5.9(i), 2023) is the Ramanujan series
        sum_{n>=0} (1+8n) (1/2)_n (1/4)_n (3/4)_n / (n!^3 9^n) = 2 sqrt3/pi
"with the shift 1/2" in the sense of J. Guillera, "Ramanujan series with a shift", arXiv:1802.02023
(J. Aust. Math. Soc. 107 (2019) 367-380): Table 1 lists the level-2 series (q = e^{-pi sqrt6}, a = 1/(2 sqrt3),
b = 8/(2 sqrt3), z = 1/9) and Table 2 ("Some conjectured values of G_2(1/2,q)") gives the value pi/6.
The same value is entry #19 ("pi/3", D = -3, a_0 = 2/3) of the table of S(H,a,P;1/2) in
H. Cohen, J. Guillera, arXiv:2101.12592, Section 6.  Both sources state it WITHOUT proof.

It is also the point (a,b,c) = (-1/2, 0, 1/2) of the three-parameter summand S(a,b,c;n) of
K. C. Au, arXiv:2602.08721v4, Example 8 (WZ-seed Seed10(a+3b-n+1/2, 2b+1, c+k+n)).

Checked here (exact rational arithmetic unless said otherwise):
  (i)   Sun's k-th term at k = n+1  ==  (1/18) (5+8n) [1, 3/4, 5/4; 3/2, 3/2, 3/2]_n 9^-n;
  (ii)  Guillera's shifted series  sum_n [1,3/4,5/4;3/2,3/2,3/2]_n (a + b/2 + b n) z^(n+1/2)  equals
        (1/(6 sqrt3)) sum_n (5+8n)[...]_n 9^-n,  hence "G_2(1/2,q) = pi/6"  <=>  T2   (numerically: 60 digits);
  (iii) Cohen-Guillera's S(H,a,P;1/2) = sum_n P(n+1/2) H_{n+1/2} / 9^(n+1/2), P(x) = 1+8x, equals
        a_0 |D| sqrt(-D) (pi/3)/pi^2 = 2 sqrt3/(3 pi)                      (numerically: 60 digits);
  (iv)  Au's S(-1/2, 0, 1/2; n) == -(5+8n) [1,3/4,5/4;3/2,3/2,3/2]_n 9^-n.
"""
from fractions import Fraction as Fr
from math import comb
from mpmath import mp, mpf, sqrt, pi, gamma, rf

def poch(a, m):
    r = Fr(1)
    for i in range(m): r *= a + i
    return r
def kern(n): return poch(Fr(1), n)*poch(Fr(3, 4), n)*poch(Fr(5, 4), n)/poch(Fr(3, 2), n)**3
def sun(k):  return Fr((8*k-3)*comb(4*k, 2*k), k*(4*k-1)*9**k*comb(2*k, k)**2)

ok1 = all(sun(n+1) == Fr(1, 18)*(5+8*n)*kern(n)/9**n for n in range(300))
print('(i)   Sun term(k=n+1) == (1/18)(5+8n)[1,3/4,5/4;3/2,3/2,3/2]_n 9^-n  for n < 300:', ok1)

def au(n, a=Fr(-1, 2), b=Fr(0), c=Fr(1, 2)):
    return (Fr(1, 9)**n*(2*a-6*b-6*c-8*n-1)*poch(-a-3*b+Fr(1, 2), n)*poch(-a+3*b+Fr(1, 2), n)*poch(a+3*b+3*c+Fr(1, 2), 2*n)
            / (poch(1-2*a, 2*n)*poch(c+1, n)*poch(2*b+c+1, n)))
ok4 = all(au(n) == -(5+8*n)*kern(n)/9**n for n in range(200))
print('(iv)  Au 2602.08721 Ex. 8: S(-1/2,0,1/2;n) == -(5+8n)[...]_n 9^-n  for n < 200:', ok4)
ok4b = all(au(n, Fr(0), Fr(0), Fr(0)) == -(1+8*n)*poch(Fr(1, 2), n)*poch(Fr(1, 4), n)*poch(Fr(3, 4), n)/(poch(Fr(1), n)**3*9**n) for n in range(100))
print('      (sanity) S(0,0,0;n) == -(1+8n)(1/2)_n(1/4)_n(3/4)_n/(n!^3 9^n):', ok4b)

mp.dps = 70
T = sum(mpf((5+8*n)*kern(n).numerator)/(kern(n).denominator*9**n) for n in range(80))
print('      sum_n (5+8n)[...]_n 9^-n - pi sqrt3 =', mp.nstr(T - pi*sqrt(3), 5))
a_, b_, z = 1/(2*sqrt(mpf(3))), 8/(2*sqrt(mpf(3))), mpf(1)/9
G2 = sum(mpf(kern(n).numerator)/kern(n).denominator*(a_ + b_/2 + b_*n)*z**(n+mpf(1)/2) for n in range(80))
ok2 = abs(G2 - pi/6) < mpf(10)**(-60)
print('(ii)  Guillera 1802.02023, Table 2:  G_2(1/2, e^{-pi sqrt6}) - pi/6 =', mp.nstr(G2 - pi/6, 5), '->', ok2)
def Hx(x): return rf(mpf(1)/2, x)*rf(mpf(1)/4, x)*rf(mpf(3)/4, x)/rf(mpf(1), x)**3
S = sum((1+8*(n+mpf(1)/2))*Hx(n+mpf(1)/2)/mpf(9)**(n+mpf(1)/2) for n in range(80))
ok3 = abs(S - 2*sqrt(mpf(3))/(3*pi)) < mpf(10)**(-60)
print('(iii) Cohen-Guillera 2101.12592, #19:  S(H,a,P;1/2) - (2/3)*3*sqrt3*(pi/3)/pi^2 =', mp.nstr(S - 2*sqrt(mpf(3))/(3*pi), 5), '->', ok3)
print('CONCLUSION: T2 is Guillera\'s (2018) conjectured half-shift of Ramanujan\'s (8n+1)/9^n series:', ok1 and ok2 and ok3 and ok4 and ok4b)
