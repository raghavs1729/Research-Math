#!/usr/bin/env python3
"""
Exact checks behind PiSeries/SunZero.lean (Theorem 7 of ramanujan-series/RESULTS.md, Sun's zero-sum series,
MathOverflow 506695), in the normalisation used in Lean:  a(n,k) = 6 pi^2 Phi(n,k)  (rational), a(1,0) = t(1) = -32/75.

  sigma = sN/sD, rho = rN/rD, R = PN/RD   (half-integers cleared: (k+2n+3/2) = (2k+4n+3)/2, ...)

Checks: (1) the cleared forms agree with certificates/T7.gp; (2) the three polynomial identities that Lean
proves by `ring` (compatibility, WZ, boundary row k = n-1); (3) closed form of a(n,0) and -R(n,0) a(n,0) = t(n+1);
(4) the certificates of the inequalities (polynomials with nonnegative coefficients), printed in Lean syntax.
"""
import sympy as sp
from sympy import Rational as Q, binomial
from lean_poly import lean
n, k, m, v, i = sp.symbols('n k m v i')
PN = ((576*n + 480)*k**6 + (6912*n**2 + 10080*n + 3600)*k**5 + (31104*n**3 + 66720*n**2 + 46576*n + 10480)*k**4
      + (59904*n**4 + 180000*n**3 + 195472*n**2 + 90464*n + 14920)*k**3
      + (25920*n**5 + 145152*n**4 + 263360*n**3 + 211992*n**2 + 77796*n + 10454)*k**2
      + (-62208*n**6 - 165120*n**5 - 137152*n**4 - 7880*n**3 + 44364*n**2 + 21138*n + 2861)*k
      + (-52992*n**7 - 208128*n**6 - 334416*n**5 - 279712*n**4 - 127288*n**3 - 29432*n**2 - 2609*n + 12))
RD = 3*(6*n+1)*(6*n+5)*(k-n)*(2*k+4*n+3)*(2*k+6*n+3)*(2*k+6*n+5)*(4*k+4*n+3)
# T7.gp
sigma = -(4*k+4*n+7)*(k+1)*(n-k-1)*(2*k+2*n+1)*(2*k+2*n+2)/((4*k+4*n+3)*(2*k+2*n+3)*(2*k+2*n+4)*(k+2*n+Q(3,2))*(k+3*n+Q(3,2)))
rho = ((4*k+4*n+7)*(4*n+1)*(4*n+2)*(4*n+3)*(4*n+4)*(2*n+1)*(2*n+2)*(2*k+2*n+1)*(2*k+2*n+2)
       / (64*(4*k+4*n+3)*(n-k)*(2*k+2*n+3)*(2*k+2*n+4)*(k+2*n+Q(3,2))*(k+2*n+Q(5,2))*(k+3*n+Q(3,2))*(k+3*n+Q(5,2))*(k+3*n+Q(7,2))))
# Lean forms
sM = 4*(4*k+4*n+7)*(k+1)*(n-k-1)*(2*k+2*n+1)*(2*k+2*n+2)          # sigma = -sM/sD
sN = -sM
sD = (4*k+4*n+3)*(2*k+2*n+3)*(2*k+2*n+4)*(2*k+4*n+3)*(2*k+6*n+3)
rN = (4*k+4*n+7)*(4*n+1)*(4*n+2)*(4*n+3)*(4*n+4)*(2*n+1)*(2*n+2)*(2*k+2*n+1)*(2*k+2*n+2)
rD = 2*(4*k+4*n+3)*(n-k)*(2*k+2*n+3)*(2*k+2*n+4)*(2*k+4*n+3)*(2*k+4*n+5)*(2*k+6*n+3)*(2*k+6*n+5)*(2*k+6*n+7)
S = lambda e, nn, kk: e.subs({n: nn, k: kk}, simultaneous=True)
assert sp.cancel(sN/sD - sigma) == 0 and sp.cancel(rN/rD - rho) == 0
print("1. cleared forms of sigma, rho agree with T7.gp: OK")
# 2. polynomial identities (as stated in Lean)
assert sp.expand(rN*S(sN, n+1, k)*sD*S(rD, n, k+1) - sN*S(rN, n, k+1)*rD*S(sD, n+1, k)) == 0
assert sp.expand((rN - rD)*(S(RD, n, k+1)*sD*RD) - rD*(S(PN, n, k+1)*sN*RD - PN*(S(RD, n, k+1)*sD))) == 0
b = lambda e: S(e, m+1, m)
assert sp.expand((b(rN) - b(rD))*(S(sD, m+2, m)*b(RD)) - (-(b(rN)*S(sN, m+2, m)*b(RD)) - b(PN)*(b(rD)*S(sD, m+2, m)))) == 0
print("2. compat_poly, wz_poly, wz_boundary_poly: OK")
# 3. closed form of a(n,0), n >= 1
poly = lambda x: 1656*x**4 - 2604*x**3 + 1158*x**2 - 229*x + 15
t = lambda mm: Q(binomial(2*mm, mm)*16**mm*poly(mm), mm**2*(2*mm-1)**2*(6*mm-1)*(6*mm-5)*binomial(6*mm, 3*mm)*binomial(3*mm, mm))
cf = lambda nn: Q(-9*nn*(4*nn+3)*(6*nn+5)*16**(nn+1)*binomial(2*nn+2, nn+1),
                  (nn+1)**2*(2*nn+1)**2*(4*nn+1)*binomial(6*nn+6, 3*nn+3)*binomial(3*nn+3, nn+1))
A = t(1); assert A == Q(-32, 75)
for nn in range(1, 30):
    assert cf(nn) == A and -S(PN/RD, nn, 0)*A == t(nn+1)
    A = A*S(rN/rD, nn, 0)
assert sp.expand(S(PN, n, 0) + (2*n+1)*(4*n+1)*(4*n+3)*poly(n+1)) == 0
# ratio of cf(n+1)/cf(n) is rho(n,0) (binomial recurrences)
c2 = 2*(2*(n+1)+1)/(n+2)                                   # C(2n+4,n+2)/C(2n+2,n+1)
c3 = (3*(n+1)+1)*(3*(n+1)+2)*(3*(n+1)+3)/((n+2)*(2*(n+1)+1)*(2*(n+1)+2))
c6 = 8*(6*(n+1)+1)*(6*(n+1)+3)*(6*(n+1)+5)/((3*(n+1)+1)*(3*(n+1)+2)*(3*(n+1)+3))
rat = lambda x: -9*x*(4*x+3)*(6*x+5)/((x+1)**2*(2*x+1)**2*(4*x+1))
assert sp.cancel(rat(n+1)/rat(n)*16*c2/(c3*c6) - S(rN/rD, n, 0)) == 0
print("3. a(n,0) = -9n(4n+3)(6n+5) 16^(n+1) C(2n+2,n+1) / ((n+1)^2 (2n+1)^2 (4n+1) C(6n+6,3n+3) C(3n+3,n+1)),"
      " -R(n,0) a(n,0) = t(n+1): OK")
# 4. inequalities
def nonneg(e, syms):
    p = sp.Poly(sp.expand(e), *syms)
    assert all(c >= 0 for c in p.coeffs()), p
    return p.as_expr()
names = {k: '(k : ℝ)', v: '(v : ℝ)', m: '(m : ℝ)', i: '(i : ℝ)'}
e1 = nonneg((sD - 24*sM).subs(n, k+2+v), (k, v))
print("4a. sD - 24 sM at n = k+2+v:\n", lean(e1, {k: names[k], v: names[v]}))
e2 = nonneg(S(rD - 6*rN, m+1, 0), (m,))
print("4b. rD - 6 rN at (n,k) = (m+1,0):\n", lean(e2, {m: names[m]}))
e3 = nonneg(poly(i+2), (i,))
print("4c. poly(i+2):\n", lean(e3, {i: names[i]}))
print("ALL CHECKS PASSED")
