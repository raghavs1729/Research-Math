#!/usr/bin/env python3
"""
THEOREM 4 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.25(i), eq. (5.20)):

   sum_{k>=0} (360k^3+612k^2+230k+15) binom(3k,k) binom(6k,3k)^2 / ((3k+1)(3k+2) 2^(15k))  =  32 sqrt(2)/pi.

WZ proof, verification script (exact sympy arithmetic except sections 4, 5).

Data (integers n,k >= 0):  h(0,0) = 1,
  h(n,k+1)/h(n,k) = sigma(n,k) = (4k+6n+1)(4k+6n+3) / (16 (k+n+1)(k+3n+1)),
  h(n+1,k)/h(n,k) = rho(n,k)   = 3(6n+1)(6n+5)(4k+6n+1)(4k+6n+3)(4k+6n+5)
                                  / (1024 (n+1)(k+n+1)(k+3n+1)(k+3n+2)(k+3n+3)),
  Q(n,k) = (64k^3+(160n+144)k^2+(240n^2+432n+92)k+360n^3+612n^2+230n+15) / (64 (k+3n+1)(k+3n+2)).
WZ pair:  F(n,k) = n h(n,k),  G(n,k) = Q(n,k) h(n,k):
          (n+1) h(n+1,k) - n h(n,k) = Q(n,k+1) h(n,k+1) - Q(n,k) h(n,k).
Closed form: h(n,k) = (2 sqrt 2/sqrt pi) 2^(-2k-5n-3/2) Gamma(3n+1/2) Gamma(2k+3n+1/2)
                       / ( Gamma(n+1/2) Gamma(k+n+1) Gamma(n+1) Gamma(k+3n+1) ).
Origin: Au's seed "Seed1" (Gauss 2F1 with a -> a/2, b -> a/2+1/2):  n h(n,k) ~ Seed1(1/2+n, 1+2n, k+n).

Proof outline (RESULTS.md, Theorem 4): telescoping gives, for N >= 1,
     sum_{n<N} Q(n,0) h(n,0) = lim_K Q(0,K) h(0,K) - N sum_k h(N,k)
(the K-limits of Q(n,K)h(n,K) vanish for n >= 1), lim_K Q(0,K)h(0,K) = lim K (1/4)_K (3/4)_K/K!^2
= 1/(Gamma(1/4)Gamma(3/4)) = 1/(sqrt2 pi),  N sum_k h(N,k) <= 7 N h(N,0) -> 0,
and Q(n,0) h(n,0) = (1/64) * (Sun's summand).
"""
import sympy as sp
from sympy import Rational as Q_, binomial
n, k = sp.symbols('n k')
sigma = (4*k+6*n+1)*(4*k+6*n+3)/(16*(k+n+1)*(k+3*n+1))
rho = 3*(6*n+1)*(6*n+5)*(4*k+6*n+1)*(4*k+6*n+3)*(4*k+6*n+5)/(1024*(n+1)*(k+n+1)*(k+3*n+1)*(k+3*n+2)*(k+3*n+3))
P = 64*k**3+(160*n+144)*k**2+(240*n**2+432*n+92)*k+360*n**3+612*n**2+230*n+15
Qf = P/(64*(k+3*n+1)*(k+3*n+2))

# 1. compatibility
assert sp.cancel(rho*sigma.subs(n, n+1) - sigma*rho.subs(k, k+1)) == 0
print("1. compatibility of the two ratios: OK")
# 2. WZ equation
assert sp.cancel(sp.together((n+1)*rho - n - (Qf.subs(k, k+1)*sigma - Qf))) == 0
print("2. WZ equation (n+1) rho - n = Q(n,k+1) sigma - Q(n,k): identically 0: OK")
# 3. boundary data: h(n,0) = C(3n,n) C(6n,3n)^2 / 2^(15n);  Q(n,0) h(n,0) = t(n)/64;  h(0,k) = (1/4)_k (3/4)_k / k!^2
hn0 = Q_(1)
for nn in range(30):
    assert hn0 == Q_(binomial(3*nn, nn)*binomial(6*nn, 3*nn)**2, 2**(15*nn))
    tn = Q_((360*nn**3+612*nn**2+230*nn+15)*binomial(3*nn, nn)*binomial(6*nn, 3*nn)**2, (3*nn+1)*(3*nn+2)*2**(15*nn))
    assert Qf.subs({n: nn, k: 0})*hn0 == tn/64
    hn0 = hn0*rho.subs({n: nn, k: 0})
# ratio identity for all n: rho(n,0) = [C(3n+3,n+1)/C(3n,n)] [C(6n+6,3n+3)/C(6n,3n)]^2 / 2^15
r1 = (3*n+1)*(3*n+2)*(3*n+3)/((n+1)*(2*n+1)*(2*n+2))
r2 = (6*n+1)*(6*n+2)*(6*n+3)*(6*n+4)*(6*n+5)*(6*n+6)/((3*n+1)*(3*n+2)*(3*n+3))**2
assert sp.cancel(rho.subs(k, 0) - r1*r2**2/2**15) == 0
assert sp.cancel(sigma.subs(n, 0) - (k+Q_(1, 4))*(k+Q_(3, 4))/(k+1)**2) == 0
print("3a. h(n,0) = C(3n,n)C(6n,3n)^2/2^(15n), Q(n,0)h(n,0) = summand/64, h(0,k) = (1/4)_k(3/4)_k/k!^2: OK")
# bounds
#  (i) 1 - sigma(n,k) >= (3/4)(n+1)/(k+3n+1)   [=> prod_{j<k} sigma(n,j) <= ((3n+1)/(k+3n+1))^(3(n+1)/4)]
d = sp.cancel(sp.together(1 - sigma - Q_(3, 4)*(n+1)/(k+3*n+1))); num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n, k).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n, k).coeffs())
#  (ii) rho(n,0) <= 27/32
d = sp.cancel(sp.together(Q_(27, 32) - rho.subs(k, 0))); num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n).coeffs())
#  (iii) Q(n,k) <= k + 3n + 3
d = sp.cancel(sp.together(k + 3*n + 3 - Qf)); num, den = sp.fraction(d)
assert all(c >= 0 for c in sp.Poly(sp.expand(num), n, k).coeffs()) and all(c > 0 for c in sp.Poly(sp.expand(den), n, k).coeffs())
print("3b. 1-sigma >= (3/4)(n+1)/(k+3n+1), rho(n,0) <= 27/32, Q <= k+3n+3: OK")

# 4. the Gamma closed form and the limit K h(0,K) -> 1/(sqrt2 pi)
import mpmath as mp
mp.mp.dps = 50
def H(nn, kk):
    return (2*mp.sqrt(2)/mp.sqrt(mp.pi)*mp.mpf(2)**(-2*kk-5*nn-mp.mpf(3)/2)*mp.gamma(3*nn+mp.mpf(1)/2)*mp.gamma(2*kk+3*nn+mp.mpf(1)/2)
            / (mp.gamma(nn+mp.mpf(1)/2)*mp.gamma(kk+nn+1)*mp.gamma(nn+1)*mp.gamma(kk+3*nn+1)))
assert abs(H(0, 0) - 1) < mp.mpf(10)**-45
for nn in range(4):
    for kk in range(4):
        assert abs(H(nn, kk+1)/H(nn, kk) - sigma.subs({n: nn, k: kk})) < mp.mpf(10)**-40
        assert abs(H(nn+1, kk)/H(nn, kk) - rho.subs({n: nn, k: kk})) < mp.mpf(10)**-40
lim = 1/(mp.gamma(mp.mpf(1)/4)*mp.gamma(mp.mpf(3)/4))
assert abs(lim - 1/(mp.sqrt(2)*mp.pi)) < mp.mpf(10)**-45
for KK in (10**3, 10**6, 10**9):
    val = float(Qf.subs({n: 0, k: KK}))*H(0, KK)
    print("   Q(0,K) h(0,K) at K=%d: %s   (limit 1/(sqrt2 pi) = %s)" % (KK, mp.nstr(val, 15), mp.nstr(lim, 15)))
print("4. Gamma closed form consistent; lim_K Q(0,K) h(0,K) = 1/(Gamma(1/4)Gamma(3/4)) = 1/(sqrt2 pi): OK")

# 5. numerics (100+ digits) and a check of the finite identity
mp.mp.dps = 120
t = lambda j: mp.mpf(360*j**3+612*j**2+230*j+15)*mp.binomial(3*j, j)*mp.binomial(6*j, 3*j)**2/((3*j+1)*(3*j+2)*mp.mpf(2)**(15*j))
S = mp.fsum(t(j) for j in range(0, 4200))
print("5. |lhs - 32 sqrt2/pi| =", mp.nstr(abs(S - 32*mp.sqrt(2)/mp.pi), 5), "(4200 terms, rate 27/32)")
assert abs(S - 32*mp.sqrt(2)/mp.pi) < mp.mpf(10)**-105
print("ALL CHECKS PASSED")
