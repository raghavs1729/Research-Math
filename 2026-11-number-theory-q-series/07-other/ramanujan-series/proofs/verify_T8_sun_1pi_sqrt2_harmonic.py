#!/usr/bin/env python3
"""
THEOREM 8 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.25(i), second identity):

   sum_{k>=0} binom(3k,k) binom(6k,3k)^2 [ (360k^3+612k^2+230k+15)(H_{2k}-H_k) - 180k^2 + 36k + 23 ]
              / ((3k+1)(3k+2) 2^(15k))   =   48 sqrt(2) log(2) / pi.

Proof by differentiating a one-parameter deformation of the WZ pair of Theorem 4; verification script.

Deformation (s real, |s| <= 1/10):
  h_s(n,k) = 2^(-2k-5n+s-3/2) Gamma(3n+s+1/2) Gamma(2k+3n-s+1/2) / ( Gamma(n+s+1/2) Gamma(k+n-s+1) n! Gamma(k+3n+1) ),
  sigma_s(n,k) = (4k+6n-2s+1)(4k+6n-2s+3) / (16 (k+3n+1)(k+n-s+1)),
  rho_s(n,k)   = (6n+2s+1)(6n+2s+3)(6n+2s+5)(4k+6n-2s+1)(4k+6n-2s+3)(4k+6n-2s+5)
                 / (1024 (n+1)(2n+2s+1)(k+3n+1)(k+3n+2)(k+3n+3)(k+n-s+1)),
  Q_s(n,k)     = N_s(n,k) / (64 (2n+2s+1)(k+3n+1)(k+3n+2)),   N_s below,
  WZ pair F = n h_s, G = Q_s h_s:   (n+1) rho_s - n = Q_s(n,k+1) sigma_s - Q_s(n,k).
  (n h_s(n,k) = Seed1(1/2+s+n, 1+s+2n, k-s+n).)

Facts (RESULTS.md, Theorem 8):
  (i)   WZ equation, identically in s                                              [exact]
  (ii)  sum_n Q_s(n,0) h_s(n,0) = lim_K Q_s(0,K) h_s(0,K) = 1/(4 sqrt(pi)) for |s| <= 1/10
        (same argument as Theorem 4; the limit is K 2^(-2K+s-3/2) Gamma(2K-s+1/2)/(Gamma(K-s+1) K!) -> 1/(4 sqrt pi))  [numeric]
  (iii) d/ds log h_s(n,0) at 0 = log 2 - psi(n+1/2) + psi(n+1) = 3 log 2 - 2 (H_{2n} - H_n)        [numeric]
  (iv)  d/ds log Q_s(n,0) at 0 = (360n^2-72n-46)/(360n^3+612n^2+230n+15)                         [exact]
  Hence sum_n t(n) [ 3 log2 - 2(H_2n - H_n) + (360n^2-72n-46)/p(n) ] = 0 with t(n) = 64 Q_0(n,0) h_0(n,0)/h_0(0,0)
  the summand of Theorem 4; with Theorem 4 (sum t = 32 sqrt2/pi) this is the statement.
"""
import sympy as sp
from sympy import Rational as Q_
n, k, s = sp.symbols('n k s')
Ns = ((128*n + (128*s + 64))*k**3 + (320*n**2 + (384*s + 448)*n + (-192*s**2 + 192*s + 144))*k**2
      + (480*n**3 + (1056*s + 1104)*n**2 + (-352*s**2 + 736*s + 616)*n + (96*s**3 - 240*s**2 + 40*s + 92))*k
      + (720*n**4 + (1440*s + 1584)*n**3 + (-192*s**2 + 1440*s + 1072)*n**2 + (96*s**3 - 272*s**2 + 296*s + 260)*n
         + (-16*s**4 + 64*s**3 - 56*s**2 - 16*s + 15)))
Qs = Ns/(64*(2*n + 2*s + 1)*(k + 3*n + 1)*(k + 3*n + 2))
sigma = (4*k + 6*n - 2*s + 1)*(4*k + 6*n - 2*s + 3)/(16*(k + 3*n + 1)*(k + n - s + 1))
rho = ((6*n + 2*s + 1)*(6*n + 2*s + 3)*(6*n + 2*s + 5)*(4*k + 6*n - 2*s + 1)*(4*k + 6*n - 2*s + 3)*(4*k + 6*n - 2*s + 5)
       / (1024*(n + 1)*(2*n + 2*s + 1)*(k + 3*n + 1)*(k + 3*n + 2)*(k + 3*n + 3)*(k + n - s + 1)))
assert sp.cancel(rho*sigma.subs(n, n + 1) - sigma*rho.subs(k, k + 1)) == 0
assert sp.cancel(sp.together((n + 1)*rho - n - (Qs.subs(k, k + 1)*sigma - Qs))) == 0
print("(i) compatibility and WZ equation with the parameter s: identically 0: OK")
# at s = 0 this is the pair of Theorem 4
p = 360*n**3 + 612*n**2 + 230*n + 15
assert sp.cancel(Qs.subs({s: 0, k: 0}) - p/(64*(3*n + 1)*(3*n + 2))) == 0
# (iv)
dQ = sp.cancel((sp.diff(Qs.subs(k, 0), s)/Qs.subs(k, 0)).subs(s, 0))
assert sp.cancel(dQ - (360*n**2 - 72*n - 46)/p) == 0
print("(iv) d/ds log Q_s(n,0)|_0 = (360n^2-72n-46)/p(n): OK")

import mpmath as mp
mp.mp.dps = 50
def h(ss, nn, kk):
    ss = mp.mpf(ss)
    return (mp.mpf(2)**(-2*kk - 5*nn + ss - mp.mpf(3)/2)*mp.gamma(3*nn + ss + mp.mpf(1)/2)*mp.gamma(2*kk + 3*nn - ss + mp.mpf(1)/2)
            / (mp.gamma(nn + ss + mp.mpf(1)/2)*mp.gamma(kk + nn - ss + 1)*mp.factorial(nn)*mp.gamma(kk + 3*nn + 1)))
Qf = sp.lambdify((s, n, k), Qs, 'mpmath')
for ss in ('0.07', '-0.05'):
    for (nn, kk) in [(0, 0), (1, 2), (3, 1)]:
        assert abs(h(ss, nn, kk + 1)/h(ss, nn, kk) - sigma.subs({s: Q_(ss), n: nn, k: kk})) < mp.mpf(10)**-40
        assert abs(h(ss, nn + 1, kk)/h(ss, nn, kk) - rho.subs({s: Q_(ss), n: nn, k: kk})) < mp.mpf(10)**-40
Phi = lambda ss: mp.fsum(Qf(mp.mpf(ss), nn, 0)*h(ss, nn, 0) for nn in range(0, 1400))
c = 1/(4*mp.sqrt(mp.pi))
for ss in ('0', '0.1', '-0.1', '0.033'):
    assert abs(Phi(ss) - c) < mp.mpf(10)**-40
    lim = [mp.mpf(K)*mp.mpf(2)**(-2*K + mp.mpf(ss) - mp.mpf(3)/2)*mp.gamma(2*K - mp.mpf(ss) + mp.mpf(1)/2)/(mp.gamma(K - mp.mpf(ss) + 1)*mp.factorial(K)) for K in (10**4, 10**8)]
    assert abs(lim[1] - c) < mp.mpf(10)**-8
print("(ii) sum_n Q_s(n,0) h_s(n,0) = 1/(4 sqrt pi) at s = 0, +-0.1, 0.033 (40 digits); K-limit consistent: OK")
H = mp.harmonic
for nn in (0, 1, 5, 12):
    dl = mp.diff(lambda ss: mp.log(h(ss, nn, 0)), 0)
    assert abs(dl - (3*mp.log(2) - 2*(H(2*nn) - H(nn)))) < mp.mpf(10)**-30
print("(iii) d/ds log h_s(n,0)|_0 = 3 log 2 - 2 (H_2n - H_n): OK")
mp.mp.dps = 120
t = lambda j: mp.binomial(3*j, j)*mp.binomial(6*j, 3*j)**2/((3*j + 1)*(3*j + 2)*mp.mpf(2)**(15*j))
S = mp.fsum(t(j)*((360*j**3 + 612*j**2 + 230*j + 15)*(mp.harmonic(2*j) - mp.harmonic(j)) - 180*j**2 + 36*j + 23) for j in range(0, 4300))
val = 48*mp.sqrt(2)*mp.log(2)/mp.pi
print("final: |lhs - 48 sqrt2 log2/pi| =", mp.nstr(abs(S - val), 5))
assert abs(S - val) < mp.mpf(10)**-105
print("ALL CHECKS PASSED")
