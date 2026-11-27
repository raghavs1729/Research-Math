#!/usr/bin/env python3
"""
THEOREM 6 (Z.-W. Sun, arXiv:2506.01870v6, Conjecture 3.2, second identity):

   sum_{k>=0} (3k+2)^2 16^k H_k / ((2k+1)(6k+1)(6k+5) binom(6k,3k))
        =  (5/2) K - pi log(432) / (6 sqrt 3),        K = L(2,(-3/.)),  H_k = 1 + 1/2 + ... + 1/k.

Proof by differentiating a two-parameter extension of the WZ pair of Theorem 1; verification script.

Two-parameter family (n,k >= 0 integers, (a,b) real near (0,0)):
  F(n,k;a,b) = (-1)^(n+k) 3^-(k+n+b) (3/4)^(a+n+1/2)
               Gamma(k+3n+2a+b+1) Gamma(k+b-a+1/2) Gamma(2n+2a+3/2) Gamma(n+a+1)^2
               / ( Gamma(2n+2a+1) Gamma(1/2-a-n) Gamma(k+3n+2a+b+3/2) Gamma(k+n+b+1) Gamma(n+a+3/4) Gamma(n+a+5/4) )
  R(n,k;a,b) = 3 (k+3n+2a+b+2)^2 / ( 4 (k+3n+2a+b+2)^2 - 1 ),   G = R F
  (seed 2F1(2A,1-A;2A+1/2;-1/3) at A = 1/2+a+n, summation index k+b+n; at a=b=0: (sqrt6/pi) x the pair of Thm 1).

Facts used (RESULTS.md, Theorem 6):
  (i)   WZ equation for all (a,b)                                            [exact, below]
  (ii)  sum_n G(n,0;a,b) = sum_k F(0,k;a,b) for (a,b) near 0                   [numeric check below]
  (iii) sum_k F(0,k;a,0) = 1/sqrt2 for all small a  (Ebisu, Mem. AMS 248 (2017), formula (1,3,3-1)(xiv):
        2F1(2A,1-A;2A+1/2;-1/3) = sqrt(pi) Gamma(2A+1/2)/(3^A Gamma(A+1/2)^2))  [numeric check below]
  (iv)  d/db sum_k F(0,k;0,b) at b=0 equals -(sqrt6/(2 pi)) (2 log3 S1 + 4 S2),
        S1 = sum (-1/3)^k/(2k+1) = pi/(2 sqrt3),  S2 = sum (-1/3)^k/(2k+1)^2 = sqrt3 Ti_2(1/sqrt3)
        = -(sqrt3 pi/12) log 3 + (15/8) K            (inverse tangent integral; Lewin)
  (v)   (d/da - 2 d/db) log F(n,0;a,b) at 0 = 3 H_n + log 432   [digamma identities, checked numerically]
        and R(n,0;a,b) depends on (a,b) only through 2a+b.
  Hence  sum_n G(n,0) (3H_n + log 432) = (15/2)(sqrt6/pi) K, i.e. the theorem (G(n,0) = (sqrt6/pi) t(n)).
"""
import sympy as sp
from sympy import Rational as Q
n, k, a, b = sp.symbols('n k a b')
u = k + 3*n + 2*a + b
sigma = -Q(1, 3)*(u + 1)*(k + b - a + Q(1, 2))/((u + Q(3, 2))*(k + n + b + 1))
rho = (-Q(1, 3)*Q(3, 4)*(u + 1)*(u + 2)*(u + 3)*(2*n + 2*a + Q(3, 2))*(2*n + 2*a + Q(5, 2))*(n + a + 1)**2*(-(a + n + Q(1, 2)))
       / ((2*n + 2*a + 1)*(2*n + 2*a + 2)*(u + Q(3, 2))*(u + Q(5, 2))*(u + Q(7, 2))*(k + n + b + 1)*(n + a + Q(3, 4))*(n + a + Q(5, 4))))
R = 3*(u + 2)**2/(4*(u + 2)**2 - 1)
assert sp.cancel(sp.together(rho - 1 - (R.subs(k, k + 1)*sigma - R))) == 0
print("(i) WZ equation with both parameters: identically 0: OK")

import mpmath as mp
mp.mp.dps = 50
def F(nn, kk, aa, bb):
    aa, bb = mp.mpf(aa), mp.mpf(bb)
    return ((-1)**(nn + kk)*mp.mpf(3)**(-(kk + nn + bb))*(mp.mpf(3)/4)**(aa + nn + mp.mpf(1)/2)
            * mp.gamma(kk + 3*nn + 2*aa + bb + 1)*mp.gamma(kk + bb - aa + mp.mpf(1)/2)*mp.gamma(2*nn + 2*aa + mp.mpf(3)/2)*mp.gamma(nn + aa + 1)**2
            / (mp.gamma(2*nn + 2*aa + 1)*mp.gamma(mp.mpf(1)/2 - aa - nn)*mp.gamma(kk + 3*nn + 2*aa + bb + mp.mpf(3)/2)*mp.gamma(kk + nn + bb + 1)
               * mp.gamma(nn + aa + mp.mpf(3)/4)*mp.gamma(nn + aa + mp.mpf(5)/4)))
def Rv(nn, kk, aa, bb):
    uu = kk + 3*nn + 2*mp.mpf(aa) + mp.mpf(bb)
    return 3*(uu + 2)**2/(4*(uu + 2)**2 - 1)
# closed form has the ratios sigma, rho
for (nn, kk, aa, bb) in [(0, 0, '0.1', '0.2'), (2, 3, '-0.07', '0.11'), (1, 5, '0.03', '-0.2')]:
    sv = sigma.subs({n: nn, k: kk, a: Q(aa), b: Q(bb)}); rv = rho.subs({n: nn, k: kk, a: Q(aa), b: Q(bb)})
    assert abs(F(nn, kk + 1, aa, bb)/F(nn, kk, aa, bb) - sv) < mp.mpf(10)**-40
    assert abs(F(nn + 1, kk, aa, bb)/F(nn, kk, aa, bb) - rv) < mp.mpf(10)**-40
SG = lambda aa, bb: mp.fsum(Rv(nn, 0, aa, bb)*F(nn, 0, aa, bb) for nn in range(0, 80))
SF = lambda aa, bb: mp.fsum(F(0, kk, aa, bb) for kk in range(0, 120))
for (aa, bb) in [('0', '0'), ('0.05', '0'), ('0', '0.07'), ('0.03', '-0.04'), ('-0.02', '0.05')]:
    assert abs(SG(aa, bb) - SF(aa, bb)) < mp.mpf(10)**-40
print("(ii) sum_n G(n,0;a,b) = sum_k F(0,k;a,b) at 5 parameter points (40 digits): OK")
for aa in ('0', '0.05', '-0.08', '0.11'):
    assert abs(SF(aa, 0) - 1/mp.sqrt(2)) < mp.mpf(10)**-40
A = mp.mpf('0.37')
assert abs(mp.hyp2f1(2*A, 1 - A, 2*A + mp.mpf(1)/2, -mp.mpf(1)/3) - mp.sqrt(mp.pi)*mp.gamma(2*A + mp.mpf(1)/2)/(3**A*mp.gamma(A + mp.mpf(1)/2)**2)) < mp.mpf(10)**-40
print("(iii) sum_k F(0,k;a,0) = 1/sqrt2 (Ebisu's evaluation): OK")
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
S1 = mp.nsum(lambda j: (-mp.mpf(1)/3)**j/(2*j + 1), [0, mp.inf]); S2 = mp.nsum(lambda j: (-mp.mpf(1)/3)**j/(2*j + 1)**2, [0, mp.inf])
assert abs(S1 - mp.pi/(2*mp.sqrt(3))) < mp.mpf(10)**-40
assert abs(S2 - (-(mp.sqrt(3)*mp.pi/12)*mp.log(3) + mp.mpf(15)/8*K)) < mp.mpf(10)**-40
dFb = mp.diff(lambda bb: SF(0, bb), 0)
assert abs(dFb + (mp.sqrt(6)/(2*mp.pi))*(2*mp.log(3)*S1 + 4*S2)) < mp.mpf(10)**-30
assert abs(dFb + mp.mpf(15)/4*mp.sqrt(6)/mp.pi*K) < mp.mpf(10)**-30
print("(iv) S1, S2 = sqrt3 Ti_2(1/sqrt3) and d/db sum_k F(0,k;0,b)|_0 = -(15/4)(sqrt6/pi) K: OK")
H = lambda m: mp.harmonic(m)
for nn in (0, 1, 2, 7):
    la = mp.diff(lambda aa: mp.log(abs(F(nn, 0, aa, 0))), 0); lb = mp.diff(lambda bb: mp.log(abs(F(nn, 0, 0, bb))), 0)
    assert abs(la - 2*lb - (3*H(nn) + mp.log(432))) < mp.mpf(10)**-30
print("(v) (d/da - 2 d/db) log F(n,0) = 3 H_n + log 432: OK")
# final
mp.mp.dps = 120
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
t = lambda m: mp.mpf((3*m + 2)**2*16**m)/((2*m + 1)*(6*m + 1)*(6*m + 5)*mp.binomial(6*m, 3*m))
S = mp.fsum(t(m)*mp.harmonic(m) for m in range(0, 260))
val = mp.mpf(5)/2*K - mp.pi*mp.log(432)/(6*mp.sqrt(3))
print("final: |lhs - (5K/2 - pi log(432)/(6 sqrt3))| =", mp.nstr(abs(S - val), 5))
assert abs(S - val) < mp.mpf(10)**-110
print("ALL CHECKS PASSED")
