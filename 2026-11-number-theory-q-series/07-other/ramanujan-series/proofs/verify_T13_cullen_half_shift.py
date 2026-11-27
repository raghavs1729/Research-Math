#!/usr/bin/env python3
"""
THEOREM 13 (H. Cohen and J. Guillera, arXiv:2101.12592, Section 6, item 13: the value at x = 1/2 of the shifted
Cullen series for 1/pi^4, stated there as an observation):

   sum_{n>=0} P(n+1/2) n!^16 (4n+2)! / (2n+1)!^10  =  2048 pi^4 / 15,      P(x) = 43680x^4+20632x^3+4340x^2+466x+21,

equivalently  sum_{n>=0} P(n+1/2) (1)_n^7 (3/4)_n (5/4)_n / ((3/2)_n^9 4096^n) = 1024 pi^4/15,
equivalently  S(H,a,P;1/2) = (2048/15)/pi^4  in the notation of Cohen-Guillera (Cullen's series is S(H,a,P;0) = 2048/pi^4).

Verification script for the proof in RESULTS.md.  With
   E(n,k) = (2k+3n+1/2) Gamma(n+1/2)^10 Gamma(k+n+1/2)^5 Gamma(k+4n) / ( Gamma(2n)^5 Gamma(k-n+3/2) Gamma(k+2n+1)^5 )
(Au's Dougall7F6 seed in the direction of his Example I, arXiv:2312.14051, stripped of an antiperiodic factor) and
its WZ certificate R:   -E(n+1,k) - E(n,k) = R(n,k+1) E(n,k+1) - R(n,k) E(n,k),   the identity
   T1(a) + T2(a) = T3(a),   T1 = sum_k E(a,k),  T2 = sum_m (-1)^m v(a+m),  T3 = sum_m (-1)^m R(a+m,0) E(a+m,0),
   v(n) = lim_k R(n,k) E(n,k),
holds for a near 1/2; at a = 1/2: T3 = 0, T1 = 2 zeta(3); and d/da at a = 1/2 gives the theorem.
"""
import sys, mpmath as mp
sys.path.insert(0, '../code')
from wz import *
from numerics import pari_to_func

Fs = seed_instance('Dougall7F6', dict(a=-1, b=-1, c=-1, d=-1, e=-1), 2, consts=dict(a='1/2', b='1/2', c='1/2', d='1/2', e='1/2'), k0=0)
R = wz_mate(Fs)
n, k = pari('n'), pari('k')
# ratios of E (exact rational functions)
sig = ((2*k+3*n+pari('5/2'))/(2*k+3*n+pari('1/2')) * (k+n+pari('1/2'))**5 * (k+4*n) / ((k-n+pari('3/2')) * (k+2*n+1)**5))
En = ((2*k+3*n+pari('7/2'))/(2*k+3*n+pari('1/2')) * (n+pari('1/2'))**10 * (k+n+pari('1/2'))**5
      * (k+4*n)*(k+4*n+1)*(k+4*n+2)*(k+4*n+3) * (k-n+pari('1/2'))
      / ((2*n)**5 * (2*n+1)**5 * (k+2*n+1)**5 * (k+2*n+2)**5))          # E(n+1,k)/E(n,k)
rho = -En                                                                   # U(m,k) = (-1)^m E(a+m,k)
assert rho - 1 - (pari.subst(R, 'k', k+1)*sig - R) == 0
assert rho*pari.subst(sig, 'n', n+1) - sig*pari.subst(rho, 'k', k+1) == 0
print("WZ equation  rho - 1 = R(n,k+1) sigma - R(n,k)  with rho = -E(n+1,k)/E(n,k), and compatibility: exact OK")
P = 43680*n**4 + 20632*n**3 + 4340*n**2 + 466*n + 21
assert pari.subst(R, 'k', 0) - (2*n-1)*P/(16384*n**4*(6*n+1)) == 0
print("R(n,0) = (2n-1) P(n) / (16384 n^4 (6n+1)): exact OK")
Rn, Rd = pari.numerator(R), pari.denominator(R)
assert pari.poldegree(Rn, 'k') == 9 and pari.poldegree(Rd, 'k') == 6
r = pari.pollead(Rn, 'k')/pari.pollead(Rd, 'k')
assert r + (820*n**2 + 180*n + 13)/(8192*n**5) == 0
print("lim_k R(n,k)/k^3 = -(820n^2+180n+13)/(8192 n^5): exact OK;  degree of R in k: 9/6;  monomials in numerator:",
      len(str(Rn).split('+')) + len(str(Rn).split('-')) - 1)
open('../certificates/T13.gp', 'w').write(
    "\\\\ Theorem 13 (Cohen-Guillera, half-shift of Cullen's 1/pi^4 series): sum P(n+1/2) n!^16 (4n+2)!/(2n+1)!^10 = 2048 pi^4/15\n"
    "\\\\ pair U(m,k) = (-1)^m E(a+m,k), V = R*U;  sg = E(n,k+1)/E(n,k), rho = -E(n+1,k)/E(n,k)\n"
    "sg = " + str(sig) + ";\nrho = " + str(rho) + ";\nR = " + str(R) + ";\n"
    "print(\"T13 WZ: \", rho - 1 - (subst(R,k,k+1)*sg - R) == 0);\n"
    "print(\"T13 compat: \", rho*subst(sg,n,n+1) - sg*subst(rho,k,k+1) == 0);\n"
    "print(\"T13 R(n,0): \", subst(R,k,0) == (2*n-1)*(43680*n^4+20632*n^3+4340*n^2+466*n+21)/(16384*n^4*(6*n+1)));\n")

# ---------------------------------------------------------------- numerics
mp.mp.dps = 45
Rf = pari_to_func(R, ['n', 'k'])
G, half = mp.gamma, mp.mpf(1)/2
def E(nn, kk):
    return (2*kk+3*nn+half)*G(nn+half)**10*G(kk+nn+half)**5*G(kk+4*nn)*mp.rgamma(kk-nn+1.5)/(G(2*nn)**5*G(kk+2*nn+1)**5)
sf = pari_to_func(sig, ['n', 'k']); ef = pari_to_func(En, ['n', 'k'])
for (a_, b_) in [(mp.mpf('0.37'), 0), (mp.mpf('0.5'), 3), (mp.mpf('1.61'), 2)]:
    assert abs(E(a_, b_+1)/E(a_, b_) - sf(a_, mp.mpf(b_))) < mp.mpf(10)**-35
    assert abs(E(a_+1, b_)/E(a_, b_) - ef(a_, mp.mpf(b_))) < mp.mpf(10)**-35
print("Gamma form of E has the ratios used: OK")
v = lambda nn: -(820*nn*nn + 180*nn + 13)/(8192*nn**5)*2*G(nn+half)**10/G(2*nn)**5
assert abs(v(mp.mpf('0.8')) - Rf(mp.mpf('0.8'), mp.mpf(10)**7)*E(mp.mpf('0.8'), mp.mpf(10)**7)) < mp.mpf(10)**-6
T1 = lambda a: mp.nsum(lambda kk: E(a, kk), [0, mp.inf])
T2 = lambda a: mp.nsum(lambda m: (-1)**int(m)*v(a+m), [0, mp.inf])
T3 = lambda a: mp.nsum(lambda m: (-1)**int(m)*Rf(a+m, mp.mpf(0))*E(a+m, 0), [0, mp.inf])
for a_ in ('0.3', '0.45', '0.55', '0.7'):
    a_ = mp.mpf(a_)
    assert abs(T1(a_) + T2(a_) - T3(a_)) < mp.mpf(10)**-30
print("T1(a) + T2(a) = T3(a) at a = 0.3, 0.45, 0.55, 0.7 (30 digits): OK")
# values at a = 1/2
for kk in range(6):
    assert abs(E(half, kk) - mp.mpf(2)/(kk+1)**3) < mp.mpf(10)**-40
    assert abs(mp.diff(lambda x: mp.log(E(x, kk)), half) + mp.mpf(9)/(2*(kk+1))) < mp.mpf(10)**-25
print("E(1/2,k) = 2/(k+1)^3 and d/dn log E(n,k) at n=1/2 equals -9/(2(k+1)): OK   =>  T1(1/2) = 2 zeta(3), T1'(1/2) = -9 zeta(4)")
H = mp.harmonic
az = lambda m: (mp.mpf(-1)/1024)**m*mp.rf(1, m)**5/mp.rf(mp.mpf(3)/2, m)**5
for m in range(6):
    assert abs((-1)**m*v(m+half) + az(m)*(77+250*m+205*m*m)/32) < mp.mpf(10)**-38
    dl = mp.diff(lambda x: mp.log(abs(v(x))), m+half)
    assert abs(dl - ((410*m+250)/mp.mpf(205*m*m+250*m+77) - 10*(H(2*m+1) - H(m)))) < mp.mpf(10)**-22
print("(-1)^m v(m+1/2) = -(1/32)(-1/1024)^m (1)_m^5/(3/2)_m^5 (205m^2+250m+77); log-derivative = (410m+250)/(205m^2+250m+77) - 10(H_{2m+1}-H_m): OK")
W1 = mp.nsum(lambda m: az(m)*(77+250*m+205*m*m), [0, mp.inf])
W2 = mp.nsum(lambda m: az(int(m))*((77+250*m+205*m*m)*(H(2*int(m)+1) - H(int(m))) - 41*m - 25), [0, mp.inf])
assert abs(W1 - 64*mp.zeta(3)) < mp.mpf(10)**-35 and abs(W2 - 48*mp.zeta(4)) < mp.mpf(10)**-35
print("cited identities (Amdeberhan-Zeilberger; Wei arXiv:2303.07887 / Au arXiv:2212.02986): 64 zeta(3) and 48 zeta(4): numerically OK")
tsum = lambda m: P_(m+half)*mp.factorial(m)**16*mp.factorial(4*m+2)/mp.factorial(2*m+1)**10
P_ = lambda x: 43680*x**4 + 20632*x**3 + 4340*x*x + 466*x + 21
for m in range(5):
    dw = mp.diff(lambda x: Rf(x, mp.mpf(0))*E(x, 0), m+half)
    assert abs((-1)**m*dw - tsum(m)/2048) < mp.mpf(10)**-18*max(1, abs(tsum(m)))
print("(-1)^m d/dn [R(n,0)E(n,0)] at n = m+1/2 equals P(m+1/2) m!^16 (4m+2)!/(2m+1)!^10 / 2048: OK")
print("T1' + T2' = -9 zeta(4) + 15 zeta(4) = pi^4/15 ; numerically T3'(1/2) - pi^4/15 =", mp.nstr(mp.diff(T3, half) - mp.pi**4/15, 5))
mp.mp.dps = 130
S = mp.nsum(lambda m: P_(m+half)*mp.rf(1, m)**7*mp.rf(mp.mpf(3)/4, m)*mp.rf(mp.mpf(5)/4, m)/(mp.rf(mp.mpf(3)/2, m)**9*mp.mpf(4096)**m), [0, mp.inf])
print("final: |sum P(n+1/2)(1)_n^7(3/4)_n(5/4)_n/((3/2)_n^9 4096^n) - 1024 pi^4/15| =", mp.nstr(abs(S - 1024*mp.pi**4/15), 5))
assert abs(S - 1024*mp.pi**4/15) < mp.mpf(10)**-120
print("ALL CHECKS PASSED")
