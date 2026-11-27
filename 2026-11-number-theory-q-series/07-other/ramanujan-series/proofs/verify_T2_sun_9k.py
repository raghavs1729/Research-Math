#!/usr/bin/env python3
"""
THEOREM 2 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.9(i), first identity; MathOverflow 456779):

      sum_{k>=1} (8k-3) binom(4k,2k) / ( k (4k-1) 9^k binom(2k,k)^2 )  =  sqrt(3) pi / 18.

WZ proof, verification script (exact sympy arithmetic except section 4).

WZ pair (integers n,k >= 0), defined through its ratios:
  f(0,0) = 1/2,
  f(n,k+1)/f(n,k) = sigma(n,k) = (2k+1)(2k+n+3)(3k+n+1)(3k+n+2)(3k+n+3)
                                  / ((2k+n+1)(2k+2n+3)(3k+2n+3)(3k+2n+4)(3k+2n+5)),
  f(n+1,k)/f(n,k) = rho(n,k)   = (n+1)(4n+3)(4n+5)(2k+n+2)(3k+n+1)
                                  / (9(2n+1)(2k+n+1)(2k+2n+3)(3k+2n+3)(3k+2n+4)),
  R(n,k) = (3k+2n+2)(6k+8n+5) / (9(2n+1)(2k+n+1)),      g = R f.
Closed form:  f(n,k) = -(pi/sqrt 3) * F(n,k),
  F(n,k) = (2k+n+1) 3^(1/2-2n) n! Gamma(k+1/2) Gamma(3k+n+1)
           / ( Gamma(-2n-1/2) Gamma(k+n+3/2) Gamma(n+1/2) Gamma(3k+2n+3) ).
Origin: Au's seed "Seed4" (Dougall's 5F4 with b -> b/3, c -> b/3+1/3, d -> b/3+2/3; arXiv:2312.14051),
  F(n,k) = Seed4(a = -n, b = -1/2 - 2n, k + 1/2 + n).

Proof outline (see RESULTS.md):  f > 0;  sum_k f(0,k) = sum 1/((3k+1)(3k+2)) = pi/(3 sqrt 3);
g(n,K) -> 0 (K -> oo);  sum_k f(N,k) <= (3/2) f(N,0) -> 0;  g(n,0) = 2 t(n+1) with t the summand.
"""
import sympy as sp
from sympy import Rational as Q, binomial
n, k = sp.symbols('n k')
sigma = ((2*k+1)*(2*k+n+3)*(3*k+n+1)*(3*k+n+2)*(3*k+n+3)
         / ((2*k+n+1)*(2*k+2*n+3)*(3*k+2*n+3)*(3*k+2*n+4)*(3*k+2*n+5)))
rho = ((n+1)*(4*n+3)*(4*n+5)*(2*k+n+2)*(3*k+n+1)
       / (9*(2*n+1)*(2*k+n+1)*(2*k+2*n+3)*(3*k+2*n+3)*(3*k+2*n+4)))
R = (3*k+2*n+2)*(6*k+8*n+5) / (9*(2*n+1)*(2*k+n+1))

# 1. the two ratios are compatible (so f is well defined), and agree with the Gamma closed form
assert sp.cancel(rho*sigma.subs(n, n+1) - sigma*rho.subs(k, k+1)) == 0
import mpmath as mp
mp.mp.dps = 40
def Fg(nn, kk):
    return ((2*kk+nn+1)*mp.mpf(3)**(mp.mpf(1)/2-2*nn)*mp.factorial(nn)*mp.gamma(kk+mp.mpf(1)/2)*mp.gamma(3*kk+nn+1)
            / (mp.gamma(-2*nn-mp.mpf(1)/2)*mp.gamma(kk+nn+mp.mpf(3)/2)*mp.gamma(nn+mp.mpf(1)/2)*mp.gamma(3*kk+2*nn+3)))
for nn in range(4):
    for kk in range(4):
        assert abs(Fg(nn, kk+1)/Fg(nn, kk) - sigma.subs({n: nn, k: kk})) < mp.mpf(10)**-30
        assert abs(Fg(nn+1, kk)/Fg(nn, kk) - rho.subs({n: nn, k: kk})) < mp.mpf(10)**-30
assert abs(-(mp.pi/mp.sqrt(3))*Fg(0, 0) - mp.mpf(1)/2) < mp.mpf(10)**-30
print("1. compatibility rho(n,k) sigma(n+1,k) = sigma(n,k) rho(n,k+1): OK; Gamma closed form consistent: OK")

# 2. WZ equation
assert sp.cancel(sp.together(rho - 1 - (R.subs(k, k+1)*sigma - R))) == 0
print("2. WZ equation rho - 1 = R(n,k+1) sigma - R(n,k): identically 0: OK")

# 3. boundary data
# f(0,k) = 1/((3k+1)(3k+2)):
assert sp.cancel(sigma.subs(n, 0) - ((3*k+1)*(3*k+2))/((3*k+4)*(3*k+5))) == 0
# g(n,0) = 2 t(n+1): equality at n = 0 and equality of ratios
t = lambda m: Q((8*m-3)*binomial(4*m, 2*m), m*(4*m-1)*9**m*binomial(2*m, m)**2)
assert R.subs({n: 0, k: 0})*Q(1, 2) == 2*t(1)
m = n + 1
tr = ((8*m+5)*m*(4*m-1)/((8*m-3)*(m+1)*(4*m+3)*9)
      * ((4*m+1)*(4*m+2)*(4*m+3)*(4*m+4)/((2*m+1)*(2*m+2))**2)      # binom(4m+4,2m+2)/binom(4m,2m)
      / ((2*m+1)*(2*m+2)/(m+1)**2)**2)                              # (binom(2m+2,m+1)/binom(2m,m))^2
gr = sp.cancel(R.subs({n: n+1, k: 0})/R.subs(k, 0)*rho.subs(k, 0))
assert sp.cancel(gr - tr) == 0
fn0 = Q(1, 2)
for nn in range(30):
    assert R.subs({n: nn, k: 0})*fn0 == 2*t(nn+1)
    fn0 = fn0*rho.subs({n: nn, k: 0})
print("3a. f(0,k) = 1/((3k+1)(3k+2)); g(n,0) = 2 (8m-3)C(4m,2m)/(m(4m-1)9^m C(2m,m)^2), m=n+1: OK")
# bounds used in the limit argument
#   sigma(N,k) <= (2k+1)/(2k+5) for N >= 1  (so sum_k f(N,k) <= f(N,0) * 3 * sum 1/((2k+1)(2k+3)) = (3/2) f(N,0))
N1 = sp.symbols('N1')  # N = N1 + 1 >= 1
d = sp.cancel(sp.together((2*k+1)/(2*k+5) - sigma.subs(n, N1+1)))
num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), N1, k).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), N1, k).coeffs())
#   rho(n,0) <= 1/8
d = sp.cancel(sp.together(Q(1, 8) - rho.subs(k, 0))); num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n).coeffs())
#   R(n,k) <= k + 2 (crude) and f(n,K) K -> 0: f(n,K) <= f(n,0) 3/((2K+1)(2K+3)) for n>=1; for n = 0 f(0,K) K -> 0 directly.
d = sp.cancel(sp.together(k + 2 - R)); num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n, k).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n, k).coeffs())
print("3b. sigma(N,k) <= (2k+1)/(2k+5) (N>=1), rho(n,0) <= 1/8, R(n,k) <= k+2: OK")

# 4. numerics (100+ digits)
mp.mp.dps = 120
S = mp.nsum(lambda j: mp.mpf(8*int(j)-3)*mp.binomial(4*int(j), 2*int(j))/(int(j)*(4*int(j)-1)*mp.mpf(9)**int(j)*mp.binomial(2*int(j), int(j))**2), [1, mp.inf])
print("4. |lhs - sqrt(3) pi/18| =", mp.nstr(abs(S - mp.sqrt(3)*mp.pi/18), 5))
assert abs(S - mp.sqrt(3)*mp.pi/18) < mp.mpf(10)**(-110)
S0 = mp.nsum(lambda j: 1/((3*j+1)*(3*j+2)), [0, mp.inf])
assert abs(S0 - mp.pi/(3*mp.sqrt(3))) < mp.mpf(10)**(-100)
print("   sum 1/((3k+1)(3k+2)) = pi/(3 sqrt 3): OK")
print("ALL CHECKS PASSED")
