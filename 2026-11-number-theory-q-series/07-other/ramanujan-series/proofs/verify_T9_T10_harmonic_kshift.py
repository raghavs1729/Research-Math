#!/usr/bin/env python3
"""
THEOREM 9 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.9(i), second identity):
   sum_{k>=1} binom(4k,2k) [ (8k-3)(5 H_{2k-1} - 4 H_{k-1}) - 6 ] / ( k (4k-1) 9^k binom(2k,k)^2 )  =  3K/2.

THEOREM 10 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.25(ii), second identity):
   sum_{k>=0} binom(3k,k) binom(6k,3k)^2 [ P(k)(4H_{6k} - 3H_{3k} - H_k) + f(k) ] / ((3k+1)(3k+2) 2^(15k)) = 192 sqrt2 log2/pi,
   P(k) = 360k^3+612k^2+230k+15,  f(k) = (1296k^3+1980k^2+912k+139)/((3k+1)(3k+2)).

Both follow from the WZ pairs of Theorems 2 and 4 used with a REAL second variable x = k + c and
differentiation in c at c = 0 ("k-shift derivative").  Verification script.

T9.  Pair of Theorem 2 for real x > -1/2:
  f(n,x) = -(pi/sqrt3) (2x+n+1) 3^(1/2-2n) n! Gamma(x+1/2) Gamma(3x+n+1) / (Gamma(-2n-1/2) Gamma(x+n+3/2) Gamma(n+1/2) Gamma(3x+2n+3)),
  R(n,x) = (3x+2n+2)(6x+8n+5)/(9(2n+1)(2x+n+1)),  g = R f;  WZ holds for real x (same rational identity).
  Phi(c) := sum_n g(n,c) = sum_k f(0,k+c) = sum_k 1/((3k+3c+1)(3k+3c+2))   (|c| <= 1/10; proof as in Theorem 2)
  Phi'(0) = -3K;   d/dx log f(n,x)|_0 = 4H_n - 5H_{2n+1} + 1/(2n+2);   d/dx log R(n,x)|_0 = 6/(8n+5) - 1/(2n+2);
  g(n,0) = 2 t(n+1).  Hence sum_m t(m) [5H_{2m-1} - 4H_{m-1} - 6/(8m-3)] = 3K/2.

T10. Pair of Theorem 4 for real x > -1/4:
  h(n,x) = (2 sqrt2/sqrt pi) 2^(-2x-5n-3/2) Gamma(3n+1/2) Gamma(2x+3n+1/2) / (Gamma(n+1/2) Gamma(x+n+1) n! Gamma(x+3n+1)),
  Q(n,x) as in Theorem 4.   Phi(b) := sum_n Q(n,b) h(n,b) = lim_K Q(0,K+b) h(0,K+b) = 1/(sqrt2 pi)  (|b| <= 1/10)
  0 = Phi'(0) = sum_n Q(n,0)h(n,0) [ d/dx log h + d/dx log Q ]_{x=0},
  d/dx log h(n,x)|_0 = (4H_{6n} - 3H_{3n} - H_n) - 6 log 2,   d/dx log Q(n,x)|_0 = f(n)/P(n).
"""
import sympy as sp
from sympy import Rational as Q_
n, x = sp.symbols('n x')
import mpmath as mp
mp.mp.dps = 50
H = mp.harmonic
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9

# ---------------- T9
sigma = ((2*x+1)*(2*x+n+3)*(3*x+n+1)*(3*x+n+2)*(3*x+n+3)/((2*x+n+1)*(2*x+2*n+3)*(3*x+2*n+3)*(3*x+2*n+4)*(3*x+2*n+5)))
rho = ((n+1)*(4*n+3)*(4*n+5)*(2*x+n+2)*(3*x+n+1)/(9*(2*n+1)*(2*x+n+1)*(2*x+2*n+3)*(3*x+2*n+3)*(3*x+2*n+4)))
R = (3*x+2*n+2)*(6*x+8*n+5)/(9*(2*n+1)*(2*x+n+1))
assert sp.cancel(sp.together(rho - 1 - (R.subs(x, x+1)*sigma - R))) == 0
dR = sp.cancel((sp.diff(R, x)/R).subs(x, 0))
assert sp.cancel(dR - (6/(8*n+5) - 1/(2*n+2))) == 0
def f(nn, xx):
    xx = mp.mpf(xx)
    return (-(mp.pi/mp.sqrt(3))*(2*xx+nn+1)*mp.mpf(3)**(mp.mpf(1)/2-2*nn)*mp.factorial(nn)*mp.gamma(xx+mp.mpf(1)/2)*mp.gamma(3*xx+nn+1)
            / (mp.gamma(-2*nn-mp.mpf(1)/2)*mp.gamma(xx+nn+mp.mpf(3)/2)*mp.gamma(nn+mp.mpf(1)/2)*mp.gamma(3*xx+2*nn+3)))
Rf = sp.lambdify((n, x), R, 'mpmath')
for (nn, xx) in [(0, '0.3'), (2, '1.7'), (1, '-0.2')]:
    assert abs(f(nn, mp.mpf(xx)+1)/f(nn, xx) - sigma.subs({n: nn, x: Q_(xx)})) < mp.mpf(10)**-40
    assert abs(f(nn+1, xx)/f(nn, xx) - rho.subs({n: nn, x: Q_(xx)})) < mp.mpf(10)**-40
    assert abs(f(0, xx) - 1/((3*mp.mpf(xx)+1)*(3*mp.mpf(xx)+2))) < mp.mpf(10)**-40
for c in ('0', '0.1', '-0.07'):
    lhs = mp.fsum(Rf(nn, mp.mpf(c))*f(nn, c) for nn in range(0, 70))
    rhs = (mp.digamma(mp.mpf(c)+mp.mpf(2)/3) - mp.digamma(mp.mpf(c)+mp.mpf(1)/3))/3     # = sum_k 1/((3k+3c+1)(3k+3c+2))
    assert abs(lhs - rhs) < mp.mpf(10)**-40
assert abs((mp.polygamma(1, mp.mpf(2)/3) - mp.polygamma(1, mp.mpf(1)/3))/3 + 3*K) < mp.mpf(10)**-40
for nn in (0, 1, 4, 9):
    dl = mp.diff(lambda xx: mp.log(f(nn, xx)), 0)
    assert abs(dl - (4*H(nn) - 5*H(2*nn+1) + mp.mpf(1)/(2*nn+2))) < mp.mpf(10)**-30
print("T9: WZ for real x, Phi(c) = (psi(c+2/3)-psi(c+1/3))/3 at c = 0, 0.1, -0.07 (40 digits), Phi'(0) = -3K,")
print("    d/dx log R and d/dx log f at x=0: OK")
mp.mp.dps = 120
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
S9 = mp.fsum(mp.binomial(4*j, 2*j)*((8*j-3)*(5*mp.harmonic(2*j-1) - 4*mp.harmonic(j-1)) - 6)/(j*(4*j-1)*mp.mpf(9)**j*mp.binomial(2*j, j)**2) for j in range(1, 160))
print("T9 final: |lhs - 3K/2| =", mp.nstr(abs(S9 - 3*K/2), 5)); assert abs(S9 - 3*K/2) < mp.mpf(10)**-110

# ---------------- T10
mp.mp.dps = 50
k = x
sigma4 = (4*k+6*n+1)*(4*k+6*n+3)/(16*(k+n+1)*(k+3*n+1))
rho4 = 3*(6*n+1)*(6*n+5)*(4*k+6*n+1)*(4*k+6*n+3)*(4*k+6*n+5)/(1024*(n+1)*(k+n+1)*(k+3*n+1)*(k+3*n+2)*(k+3*n+3))
P4 = 64*k**3+(160*n+144)*k**2+(240*n**2+432*n+92)*k+360*n**3+612*n**2+230*n+15
Qf = P4/(64*(k+3*n+1)*(k+3*n+2))
assert sp.cancel(sp.together((n+1)*rho4 - n - (Qf.subs(k, k+1)*sigma4 - Qf))) == 0
P = 360*n**3+612*n**2+230*n+15
fS = (1296*n**3+1980*n**2+912*n+139)/((3*n+1)*(3*n+2))
dQ = sp.cancel((sp.diff(Qf, k)/Qf).subs(k, 0))
assert sp.cancel(dQ - fS/P) == 0
def h(nn, xx):
    xx = mp.mpf(xx)
    return (2*mp.sqrt(2)/mp.sqrt(mp.pi)*mp.mpf(2)**(-2*xx-5*nn-mp.mpf(3)/2)*mp.gamma(3*nn+mp.mpf(1)/2)*mp.gamma(2*xx+3*nn+mp.mpf(1)/2)
            / (mp.gamma(nn+mp.mpf(1)/2)*mp.gamma(xx+nn+1)*mp.factorial(nn)*mp.gamma(xx+3*nn+1)))
Ql = sp.lambdify((n, x), Qf, 'mpmath')
for b in ('0', '0.1', '-0.08'):
    val = mp.fsum(Ql(nn, mp.mpf(b))*h(nn, b) for nn in range(0, 1400))
    assert abs(val - 1/(mp.sqrt(2)*mp.pi)) < mp.mpf(10)**-40
for nn in (0, 1, 3, 8):
    dl = mp.diff(lambda xx: mp.log(h(nn, xx)), 0)
    assert abs(dl - (4*H(6*nn) - 3*H(3*nn) - H(nn) - 6*mp.log(2))) < mp.mpf(10)**-30
print("T10: WZ, Phi(b) = 1/(sqrt2 pi) at b = 0, 0.1, -0.08 (40 digits), d/dx log Q|_0 = f/P (exact), d/dx log h|_0: OK")
mp.mp.dps = 120
t = lambda j: mp.binomial(3*j, j)*mp.binomial(6*j, 3*j)**2/((3*j+1)*(3*j+2)*mp.mpf(2)**(15*j))
S10 = mp.fsum(t(j)*((360*j**3+612*j**2+230*j+15)*(4*mp.harmonic(6*j) - 3*mp.harmonic(3*j) - mp.harmonic(j))
                    + mp.mpf(1296*j**3+1980*j**2+912*j+139)/((3*j+1)*(3*j+2))) for j in range(0, 4300))
val = 192*mp.sqrt(2)*mp.log(2)/mp.pi
print("T10 final: |lhs - 192 sqrt2 log2/pi| =", mp.nstr(abs(S10 - val), 5)); assert abs(S10 - val) < mp.mpf(10)**-104
print("ALL CHECKS PASSED")
