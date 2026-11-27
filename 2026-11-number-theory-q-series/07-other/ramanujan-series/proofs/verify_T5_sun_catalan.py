#!/usr/bin/env python3
"""
THEOREM 5 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 4.33(i), first identity):

   sum_{k>=1} (77k^3-86k^2+29k-3) (-32)^k / ( k^3 (3k-1)(3k-2) binom(3k,k) binom(4k,k)^2 )  =  -6 G,
   G = Catalan's constant.

WZ proof, verification script (exact sympy arithmetic except sections 4, 5).

Data (integers n,k >= 0):  psi(0,0) = 1,
  psi(n,k+1)/psi(n,k) = sigma(n,k) = (2k+n+1)(2k+n+2) / ((2k+2n+3)(2k+4n+3)),
  psi(n+1,k)/psi(n,k) = rho(n,k)   = -(n+1)(3n+1)(3n+2)(3n+3)(2k+n+1)
                                      / ((4n+1)(4n+3)(2k+2n+3)(2k+4n+3)(2k+4n+5)),
  R(n,k) = ((16n+12)k^2+(66n^2+86n+28)k+77n^3+145n^2+88n+17) / ((4n+1)(4n+3)(2k+4n+3)),   gamma = R psi.
Closed form: psi(n,k) = (-1)^n 2 pi^(3/2) n! (3n)! (2k+n)! 2^(-2k-5n-3)
                         / ( Gamma(k+n+3/2) Gamma(2n+1/2) Gamma(k+2n+3/2) ).
Origin: Au's seed "Seed1":  F(n,k) = Seed1(a = -n, c = 1+n, k + 1/2 + n) vanishes identically for integer n
because of a factor 1/Gamma(-n); psi is its first-order coefficient in n (same shift ratios, same certificate).

Proof outline (RESULTS.md, Theorem 5):  psi(0,k) = 4^k/((2k+1)^2 C(2k,k)), sum_k psi(0,k) = 2G;
gamma(n,K) -> 0; sum_k |psi(N,k)| <= 6 |psi(N,0)| -> 0;  gamma(n,0) = -(1/3) t(n+1).
"""
import sympy as sp
from sympy import Rational as Q, binomial
n, k = sp.symbols('n k')
sigma = (2*k+n+1)*(2*k+n+2)/((2*k+2*n+3)*(2*k+4*n+3))
rho = -(n+1)*(3*n+1)*(3*n+2)*(3*n+3)*(2*k+n+1)/((4*n+1)*(4*n+3)*(2*k+2*n+3)*(2*k+4*n+3)*(2*k+4*n+5))
R = ((16*n+12)*k**2+(66*n**2+86*n+28)*k+77*n**3+145*n**2+88*n+17)/((4*n+1)*(4*n+3)*(2*k+4*n+3))

assert sp.cancel(rho*sigma.subs(n, n+1) - sigma*rho.subs(k, k+1)) == 0
print("1. compatibility of the two ratios: OK")
assert sp.cancel(sp.together(rho - 1 - (R.subs(k, k+1)*sigma - R))) == 0
print("2. WZ equation rho - 1 = R(n,k+1) sigma - R(n,k): identically 0: OK")

# 3. boundary data
t = lambda m: Q((77*m**3-86*m**2+29*m-3)*(-32)**m, m**3*(3*m-1)*(3*m-2)*binomial(3*m, m)*binomial(4*m, m)**2)
pn0 = Q(1)
for nn in range(30):
    assert R.subs({n: nn, k: 0})*pn0 == -Q(1, 3)*t(nn+1)
    pn0 = pn0*rho.subs({n: nn, k: 0})
p0k = Q(1)
for kk in range(30):
    assert p0k == Q(4**kk, (2*kk+1)**2*binomial(2*kk, kk))
    p0k = p0k*sigma.subs({n: 0, k: kk})
# ratio identities valid for all n, k
assert sp.cancel(sigma.subs(n, 0) - 4*(2*k+1)**2*(k+1)**2/((2*k+3)**2*(2*k+1)*(2*k+2))) == 0
m = n + 1
tr = ((77*(m+1)**3-86*(m+1)**2+29*(m+1)-3)*(-32)*m**3*(3*m-1)*(3*m-2)
      / ((77*m**3-86*m**2+29*m-3)*(m+1)**3*(3*m+2)*(3*m+1))
      / ((3*m+1)*(3*m+2)*(3*m+3)/((m+1)*(2*m+1)*(2*m+2)))
      / ((4*m+1)*(4*m+2)*(4*m+3)*(4*m+4)/((m+1)*(3*m+1)*(3*m+2)*(3*m+3)))**2)
gr = sp.cancel(R.subs({n: n+1, k: 0})/R.subs(k, 0)*rho.subs(k, 0))
assert sp.cancel(gr - tr) == 0
print("3a. psi(0,k) = 4^k/((2k+1)^2 C(2k,k));  gamma(n,0) = -(1/3) t(n+1) (n<30 exactly + ratio identity): OK")

def nonneg_rational(expr, *vars_):
    """True if expr = P/Q with P having nonnegative and Q positive coefficients (after a common sign)"""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    pn, pd = sp.Poly(sp.expand(num), *vars_), sp.Poly(sp.expand(den), *vars_)
    if all(c <= 0 for c in pd.coeffs()):
        pn, pd = -pn, -pd
    return all(c >= 0 for c in pn.coeffs()) and all(c > 0 for c in pd.coeffs())
# bounds:  1 - sigma >= (n+3/4)/(k+2n+3/2);  |rho(n,0)| <= 1/8;  0 <= R <= 2k + 2n + 2
assert nonneg_rational(1 - sigma - (n+Q(3, 4))/(k+2*n+Q(3, 2)), n, k)
assert nonneg_rational(Q(1, 8) + rho.subs(k, 0), n)
assert nonneg_rational(-rho.subs(k, 0), n)
assert nonneg_rational(2*k + 2*n + 2 - R, n, k)
print("3b. 1-sigma >= (n+3/4)/(k+2n+3/2), |rho(n,0)| <= 1/8, R <= 2k+2n+2: OK")

# 4. Gamma closed form; sum_k psi(0,k) = 2G
import mpmath as mp
mp.mp.dps = 50
def Psi(nn, kk):
    return ((-1)**nn*2*mp.pi**mp.mpf(1.5)*mp.factorial(nn)*mp.factorial(3*nn)*mp.factorial(2*kk+nn)*mp.mpf(2)**(-2*kk-5*nn-3)
            / (mp.gamma(kk+nn+mp.mpf(3)/2)*mp.gamma(2*nn+mp.mpf(1)/2)*mp.gamma(kk+2*nn+mp.mpf(3)/2)))
assert abs(Psi(0, 0) - 1) < mp.mpf(10)**-45
for nn in range(4):
    for kk in range(4):
        assert abs(Psi(nn, kk+1)/Psi(nn, kk) - sigma.subs({n: nn, k: kk})) < mp.mpf(10)**-40
        assert abs(Psi(nn+1, kk)/Psi(nn, kk) - rho.subs({n: nn, k: kk})) < mp.mpf(10)**-40
# sum_{k>=0} 4^k/((2k+1)^2 C(2k,k)) = int_0^{pi/2} x/sin x dx = 2G
I = mp.quad(lambda th: th/mp.sin(th), [0, mp.pi/2])
assert abs(I - 2*mp.catalan) < mp.mpf(10)**-40
S0 = mp.nsum(lambda j: mp.mpf(4)**j/((2*j+1)**2*mp.binomial(2*j, j)), [0, mp.inf])
print("4. closed form consistent; int_0^{pi/2} x/sin x dx = 2G; sum_k psi(0,k) =", mp.nstr(S0, 12), "(2G =", mp.nstr(2*mp.catalan, 12), ", terms ~ k^(-3/2))")

# 5. numerics, 100+ digits
mp.mp.dps = 120
tm = lambda j: mp.mpf(77*j**3-86*j**2+29*j-3)*mp.mpf(-32)**j/(mp.mpf(j)**3*(3*j-1)*(3*j-2)*mp.binomial(3*j, j)*mp.binomial(4*j, j)**2)
S = mp.fsum(tm(j) for j in range(1, 120))
print("5. |lhs + 6G| =", mp.nstr(abs(S + 6*mp.catalan), 5))
assert abs(S + 6*mp.catalan) < mp.mpf(10)**-110
print("ALL CHECKS PASSED")
