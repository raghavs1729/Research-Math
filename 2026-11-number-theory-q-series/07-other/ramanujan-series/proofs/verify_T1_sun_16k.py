#!/usr/bin/env python3
"""
THEOREM 1 (Z.-W. Sun, arXiv:2506.01870v6, Conjecture 3.2, first identity):

      sum_{k>=0} (3k+2)^2 16^k / ((2k+1)(6k+1)(6k+5) binom(6k,3k))  =  pi/(2 sqrt 3).

WZ proof, verification script.  Everything below is exact (sympy rationals) except section 4.

WZ pair (n,k >= 0 integers; (a)_m = rising factorial):

  f(n,k) = (-1/3)^k (1/4)^n * (1/2)_n (3/2)_{2n} (1)_n^2 (1)_{k+3n} (1/2)_k
                              / ( (1)_{2n} (3/4)_n (5/4)_n (3/2)_{k+3n} (1)_{k+n} )
  R(n,k) = 3 (k+3n+2)^2 / ((2k+6n+3)(2k+6n+5)),        g(n,k) = R(n,k) f(n,k).

Origin: the strange evaluation 2F1(2a, 1-a; 2a+1/2; -1/3) (A. Ebisu, Special values of the
hypergeometric series, Mem. AMS 2017, family (1,3,3-1)) used as a WZ seed in the sense of Au,
with a -> 1/2 + n, k -> k + n.

Checked here:
  1. f(n,k+1)/f(n,k) = sigma(n,k), f(n+1,k)/f(n,k) = rho(n,k) (explicit rational functions), by exact
     evaluation of the closed form on a grid, and the compatibility identity as polynomials;
  2. the WZ equation  rho - 1 = R(n,k+1) sigma - R(n,k)  as an identity of rational functions
     (this is  f(n+1,k) - f(n,k) = g(n,k+1) - g(n,k));
  3. g(n,0) = (3n+2)^2 16^n/((2n+1)(6n+1)(6n+5) binom(6n,3n)) exactly for n < 40, and as a ratio identity;
     f(0,k) = (-1/3)^k/(2k+1);  the bounds used in the limit argument;
  4. numerics to 100 digits.
"""
import sympy as sp
from sympy import Rational as Q, rf, binomial

n, k = sp.symbols('n k')

def f(nn, kk):
    return (Q(-1, 3)**kk * Q(1, 4)**nn * rf(Q(1, 2), nn) * rf(Q(3, 2), 2*nn) * rf(1, nn)**2 * rf(1, kk + 3*nn) * rf(Q(1, 2), kk)
            / (rf(1, 2*nn) * rf(Q(3, 4), nn) * rf(Q(5, 4), nn) * rf(Q(3, 2), kk + 3*nn) * rf(1, kk + nn)))

R = 3*(k + 3*n + 2)**2 / ((2*k + 6*n + 3)*(2*k + 6*n + 5))
sigma = -Q(1, 3) * (k + 3*n + 1) * (k + Q(1, 2)) / ((k + 3*n + Q(3, 2)) * (k + n + 1))
rho = (Q(1, 4) * (n + Q(1, 2)) * (2*n + Q(3, 2)) * (2*n + Q(5, 2)) * (n + 1)**2
       * (k + 3*n + 1) * (k + 3*n + 2) * (k + 3*n + 3)
       / ((2*n + 1) * (2*n + 2) * (n + Q(3, 4)) * (n + Q(5, 4))
          * (k + 3*n + Q(3, 2)) * (k + 3*n + Q(5, 2)) * (k + 3*n + Q(7, 2)) * (k + n + 1)))

# 1. ratios of the closed form
for nn in range(0, 6):
    for kk in range(0, 6):
        assert f(nn, kk + 1) == f(nn, kk) * sigma.subs({n: nn, k: kk})
        assert f(nn + 1, kk) == f(nn, kk) * rho.subs({n: nn, k: kk})
print("1. ratio formulas sigma, rho agree with the closed form on 0<=n,k<6: OK")

# 2. WZ equation as a rational-function identity
wz = sp.cancel(sp.together(rho - 1 - (R.subs(k, k + 1) * sigma - R)))
assert wz == 0
print("2. WZ equation rho - 1 = R(n,k+1) sigma - R(n,k): identically 0: OK")

# 3. boundary data
def t(m):
    return Q((3*m + 2)**2 * 16**m, (2*m + 1)*(6*m + 1)*(6*m + 5)*binomial(6*m, 3*m))
for nn in range(40):
    assert R.subs({n: nn, k: 0}) * f(nn, 0) == t(nn)
# ratio identity proving g(n,0) = t(n) for all n (both sides equal 4/15... at n=0 and have the same ratio)
tr = ((3*n + 5)**2 * 16 * (2*n + 1) * (6*n + 1) * (6*n + 5)
      / ((3*n + 2)**2 * (2*n + 3) * (6*n + 7) * (6*n + 11))
      * (3*n + 1) * (3*n + 2) * (3*n + 3) * (3*n + 1) * (3*n + 2) * (3*n + 3)
      / ((6*n + 1) * (6*n + 2) * (6*n + 3) * (6*n + 4) * (6*n + 5) * (6*n + 6)))
gr = sp.cancel(R.subs({n: n + 1, k: 0}) / R.subs(k, 0) * rho.subs(k, 0))
assert sp.cancel(gr - tr) == 0
print("3a. g(n,0) = (3n+2)^2 16^n/((2n+1)(6n+1)(6n+5) C(6n,3n)): n<40 exactly, and ratio identity: OK")
for kk in range(40):
    assert f(0, kk) == Q(-1, 3)**kk / (2*kk + 1)
print("3b. f(0,k) = (-1/3)^k/(2k+1) for k<40: OK")
# bounds: |sigma| <= 1/3 for n,k>=0 ;  rho(n,0) <= 1/3 for n >= 0 (used: sum_k |f(N,k)| <= (3/2) f(N,0) -> 0)
num, den = sp.fraction(sp.cancel(sp.Rational(1, 3) - (-sigma)))
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n, k).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n, k).coeffs())
num, den = sp.fraction(sp.cancel(sp.Rational(1, 3) - rho.subs(k, 0)))
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n).coeffs())
print("3c. 0 <= -sigma(n,k) <= 1/3 and 0 < rho(n,0) <= 1/3 (numerator/denominator of the differences have"
      " nonnegative coefficients): OK")

# 4. numerics
import mpmath as mp
mp.mp.dps = 120
S = mp.nsum(lambda m: mp.mpf((3*int(m) + 2)**2 * 16**int(m)) / ((2*int(m) + 1)*(6*int(m) + 1)*(6*int(m) + 5)*mp.binomial(6*int(m), 3*int(m))), [0, mp.inf])
print("4. |lhs - pi/(2 sqrt 3)| =", mp.nstr(abs(S - mp.pi/(2*mp.sqrt(3))), 5), " (120-digit arithmetic)")
assert abs(S - mp.pi/(2*mp.sqrt(3))) < mp.mpf(10)**(-110)
print("ALL CHECKS PASSED")
