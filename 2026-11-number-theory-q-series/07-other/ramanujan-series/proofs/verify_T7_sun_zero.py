#!/usr/bin/env python3
"""
THEOREM 7 (Z.-W. Sun, MathOverflow 506695, 6 Jan 2026):

   sum_{k>=1} binom(2k,k) 16^k (1656k^4-2604k^3+1158k^2-229k+15)
              / ( k^2 (2k-1)^2 (6k-1)(6k-5) binom(6k,3k) binom(3k,k) )  =  0.

WZ proof, verification script (exact sympy arithmetic except section 4).

WZ pair, supported on the triangle 0 <= k < n  (Phi(n,k) = 0 for k >= n):
  Phi(n,k) = (-1)^(k+1) (4k+4n+3) 2^(1-6n) k! (4n)! (2n)!
             / ( pi (n-k-1)! (2k+2n+1)(2k+2n+2) Gamma(k+2n+3/2) Gamma(k+3n+3/2) ),
  Gam(n,k) = R(n,k) Phi(n,k) for k < n,   Gam(n,n) = -Phi(n+1,n),   Gam(n,k) = 0 for k > n,
  R(n,k)   = P(n,k) / ( 3(6n+1)(6n+5)(k-n)(2k+4n+3)(2k+6n+3)(2k+6n+5)(4k+4n+3) ),  P below.
Origin: Au's seed "Seed3" (Dougall 5F4 with b -> b/2, c -> b/2+1/2):
  F(nu,k) = Seed3(1/2-2nu, -2nu, 1/2-3nu, k+1/2+2nu); Phi, Gam are the (finite) values at integer nu = n
  (F has the factor Gamma(k-nu+1)/Gamma(-2nu); the WZ equation for real nu passes to the limit).

Proof outline (RESULTS.md, Theorem 7): summing the WZ equation over 0<=k<=N, 0<=n<N gives
     sum_{n<N} Gam(n,0) = - sum_{k<N} Phi(N,k),
|Phi(N,k+1)/Phi(N,k)| <= 1/24, |Phi(N,0)| -> 0, and Gam(n,0) = -(1/(6 pi^2)) t(n+1), t = Sun's summand.
"""
import sympy as sp
from sympy import Rational as Q, binomial
n, k = sp.symbols('n k')
P = ((576*n + 480)*k**6 + (6912*n**2 + 10080*n + 3600)*k**5 + (31104*n**3 + 66720*n**2 + 46576*n + 10480)*k**4
     + (59904*n**4 + 180000*n**3 + 195472*n**2 + 90464*n + 14920)*k**3
     + (25920*n**5 + 145152*n**4 + 263360*n**3 + 211992*n**2 + 77796*n + 10454)*k**2
     + (-62208*n**6 - 165120*n**5 - 137152*n**4 - 7880*n**3 + 44364*n**2 + 21138*n + 2861)*k
     + (-52992*n**7 - 208128*n**6 - 334416*n**5 - 279712*n**4 - 127288*n**3 - 29432*n**2 - 2609*n + 12))
R = P/(3*(6*n+1)*(6*n+5)*(k-n)*(2*k+4*n+3)*(2*k+6*n+3)*(2*k+6*n+5)*(4*k+4*n+3))
# ratios of Phi (valid on the support):
sigma = -(4*k+4*n+7)*(k+1)*(n-k-1)*(2*k+2*n+1)*(2*k+2*n+2)/((4*k+4*n+3)*(2*k+2*n+3)*(2*k+2*n+4)*(k+2*n+Q(3,2))*(k+3*n+Q(3,2)))
rho = ((4*k+4*n+7)*(4*n+1)*(4*n+2)*(4*n+3)*(4*n+4)*(2*n+1)*(2*n+2)*(2*k+2*n+1)*(2*k+2*n+2)
       / (64*(4*k+4*n+3)*(n-k)*(2*k+2*n+3)*(2*k+2*n+4)*(k+2*n+Q(3,2))*(k+2*n+Q(5,2))*(k+3*n+Q(3,2))*(k+3*n+Q(5,2))*(k+3*n+Q(7,2))))

def Phi(nn, kk):
    if kk >= nn: return sp.Integer(0)
    return ((-1)**(kk+1)*(4*kk+4*nn+3)*Q(2)**(1-6*nn)*sp.factorial(kk)*sp.factorial(4*nn)*sp.factorial(2*nn)
            / (sp.pi*sp.factorial(nn-kk-1)*(2*kk+2*nn+1)*(2*kk+2*nn+2)*sp.gamma(kk+2*nn+Q(3,2))*sp.gamma(kk+3*nn+Q(3,2))))
def Gam(nn, kk):
    if kk > nn: return sp.Integer(0)
    if kk == nn: return -Phi(nn+1, nn)
    return R.subs({n: nn, k: kk})*Phi(nn, kk)

# 1. ratios agree with the closed form
for nn in range(1, 6):
    for kk in range(0, nn):
        if kk+1 < nn: assert sp.simplify(Phi(nn, kk+1)/Phi(nn, kk) - sigma.subs({n: nn, k: kk})) == 0
        assert sp.simplify(Phi(nn+1, kk)/Phi(nn, kk) - rho.subs({n: nn, k: kk})) == 0
print("1. ratio formulas agree with the closed form of Phi: OK")
# 2. WZ equation: generic rational identity (k <= n-2), and exact check of ALL cases on a grid incl. the boundary k = n-1, n
assert sp.cancel(sp.together(rho - 1 - (R.subs(k, k+1)*sigma - R))) == 0
for nn in range(0, 7):
    for kk in range(0, 9):
        assert sp.simplify(Phi(nn+1, kk) - Phi(nn, kk) - (Gam(nn, kk+1) - Gam(nn, kk))) == 0
# boundary case k = n-1 for all n >= 1 as a rational identity:
#   Phi(n+1,n-1) - Phi(n,n-1) = -Phi(n+1,n) - R(n,n-1) Phi(n,n-1)
# divide by Phi(n,n-1):  rho(n,n-1) - 1 = -rho(n,n-1) * [Phi(n+1,n)/Phi(n+1,n-1)] - R(n,n-1),
# and Phi(n+1,n)/Phi(n+1,n-1) = sigma(n+1,n-1):
lhs = rho.subs(k, n-1) - 1
rhs = -rho.subs(k, n-1)*sigma.subs({n: n+1, k: n-1}, simultaneous=True) - R.subs(k, n-1)
assert sp.cancel(sp.together(lhs - rhs)) == 0
print("2. WZ equation: generic identity, boundary identity (k = n-1), and all cases on the grid n<7, k<9: OK")
# 3. boundary data
t = lambda m: Q(binomial(2*m, m)*16**m*(1656*m**4-2604*m**3+1158*m**2-229*m+15), m**2*(2*m-1)**2*(6*m-1)*(6*m-5)*binomial(6*m, 3*m)*binomial(3*m, m))
for nn in range(0, 12):
    assert sp.simplify(Gam(nn, 0) + t(nn+1)/(6*sp.pi**2)) == 0
# ratio identity for Gam(n,0) = R(n,0)Phi(n,0), n >= 1:
m = n + 1
poly = lambda mm: 1656*mm**4-2604*mm**3+1158*mm**2-229*mm+15
tr = ((2*m+1)*(2*m+2)/(m+1)**2*16*poly(m+1)*m**2*(2*m-1)**2*(6*m-1)*(6*m-5)
      / (poly(m)*(m+1)**2*(2*m+1)**2*(6*m+5)*(6*m+1))
      / ((6*m+1)*(6*m+2)*(6*m+3)*(6*m+4)*(6*m+5)*(6*m+6)/((3*m+1)*(3*m+2)*(3*m+3))**2)
      / ((3*m+1)*(3*m+2)*(3*m+3)/((m+1)*(2*m+1)*(2*m+2))))
gr = sp.cancel(R.subs({n: n+1, k: 0})/R.subs(k, 0)*rho.subs(k, 0))
assert sp.cancel(gr - tr) == 0
print("3a. Gam(n,0) = -(1/(6 pi^2)) t(n+1): n<12 exactly, and ratio identity: OK")
# bounds: |sigma(N,k)| <= 1/24 for 0<=k<=N-2 (write N = k+2+j, j>=0);  |Phi(N,0)| ratio rho(N,0) <= 1/6
j = sp.symbols('j')
d = sp.cancel(sp.together(Q(1, 24) + sigma.subs(n, k+2+j))); num, den = sp.fraction(d)
pn, pd = sp.Poly(sp.expand(num), k, j), sp.Poly(sp.expand(den), k, j)
if all(c <= 0 for c in pd.coeffs()): pn, pd = -pn, -pd
assert all(c >= 0 for c in pn.coeffs()) and all(c > 0 for c in pd.coeffs())
d = sp.cancel(sp.together(Q(1, 6) - rho.subs({k: 0, n: n+1}))); num, den = sp.fraction(d)
pn, pd = sp.Poly(sp.expand(num), n), sp.Poly(sp.expand(den), n)
if all(c <= 0 for c in pd.coeffs()): pn, pd = -pn, -pd
assert all(c >= 0 for c in pn.coeffs()) and all(c > 0 for c in pd.coeffs())
print("3b. |Phi(N,k+1)/Phi(N,k)| <= 1/24 on the support; Phi(N+1,0)/Phi(N,0) <= 1/6 (N>=1): OK")
# 4. finite identity and numerics
for N in (1, 2, 5, 9):
    lhs = sum(Gam(nn, 0) for nn in range(N)); rhs = -sum(Phi(N, kk) for kk in range(N))
    assert sp.simplify(lhs - rhs) == 0
print("4a. finite identity sum_{n<N} Gam(n,0) = -sum_{k<N} Phi(N,k) for N = 1,2,5,9 (exact): OK")
import mpmath as mp
mp.mp.dps = 130
tm = lambda m: mp.binomial(2*m, m)*mp.mpf(16)**m*(1656*m**4-2604*m**3+1158*m**2-229*m+15)/(mp.mpf(m)**2*(2*m-1)**2*(6*m-1)*(6*m-5)*mp.binomial(6*m, 3*m)*mp.binomial(3*m, m))
S = mp.fsum(tm(m) for m in range(1, 200))
print("4b. |sum| =", mp.nstr(abs(S), 5), "(first term %s)" % mp.nstr(tm(1), 10))
assert abs(S) < mp.mpf(10)**-110
print("ALL CHECKS PASSED")
