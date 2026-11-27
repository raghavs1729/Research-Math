#!/usr/bin/env python3
"""
THEOREM 3 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.1(i), second identity; MathOverflow 456444):

      sum_{k>=1} (4k-1) 3^k / ( (2k-1) k^2 binom(2k,k) binom(3k,k) )  =  2 K,
      K = L(2, (-3/.)) = sum_{m>=0} ( 1/(3m+1)^2 - 1/(3m+2)^2 ).

WZ proof with a continuous second variable; verification script.

WZ pair (n >= 0 integer, x > -1 real):
   F(n,x) = c_n (3/4)^x (-x)_{2n} Gamma(x+n+1)/Gamma(x+n+3/2),
   c_n    = 9^(-n) Gamma(3n+3/2) n!^2 / ( 4 (2n)! (3n)! Gamma(n+5/6) Gamma(n+7/6) ),
   R(n,x) = ( 8n^2 + 6n - (2n+3) x - x^2 ) / ((3n+1)(3n+2)),         G = R F.
Origin: Ebisu's evaluation 2F1(2a, 3a-1/2; 3a; 3/4) (Mem. AMS 2017, family (1,3,3-1)) as a WZ seed,
F(n,x) = seed(a = 1/2 + n, k = x - 2n).

Steps of the proof (RESULTS.md, Theorem 3) and what is checked here:
  (a) WZ equation F(n+1,x)-F(n,x) = G(n,x+1)-G(n,x): exact rational-function identity   [exact]
  (b) d/dx G(n,x) at x=0 equals -(3/(4 pi)) * (4m-1)3^m/((2m-1)m^2 C(2m,m)C(3m,m)), m=n+1   [exact]
  (c) D_n := sum_{k>=0} d/dx F(n,x)|_{x=k};  D_{n+1}-D_n = -dG/dx(n,0)                      [numeric]
  (d) D_0 = (3/(8 pi)) I,  I = int_0^1 log(3t/4)/((1-3t/4) sqrt(1-t)) dt = -4K              [numeric, 60 digits]
  (e) D_N = (c_N/sqrt(pi)) int_0^1 t^N (1-t)^(-1/2) Psi_N(3t/4) dt,
      Psi_N(y) = -(2N)! y^(2N) (1-y)^(-2N-1) int_0^{(1-y)/y} s^(2N)/(1+s) ds,
      |Psi_N(y)| <= (2N)!/((2N+1) y),  hence |D_N| <= c_N (2N)! 4 B(N,1/2)/(3 sqrt(pi) (2N+1)) -> 0   [numeric]
  (f) final identity to 100 digits.
"""
import sympy as sp
from sympy import Rational as Q
n, x = sp.symbols('n x')

# ratios of F:  sigma = F(n,x+1)/F(n,x),  rho = F(n+1,x)/F(n,x)
sigma = Q(3, 4)*(x + 1)*(x + n + 1)/((x + 1 - 2*n)*(x + n + Q(3, 2)))
# c_{n+1}/c_n:
cr = (Q(1, 9)*(3*n + Q(3, 2))*(3*n + Q(5, 2))*(3*n + Q(7, 2))*(n + 1)**2
      / ((2*n + 1)*(2*n + 2)*(3*n + 1)*(3*n + 2)*(3*n + 3)*(n + Q(5, 6))*(n + Q(7, 6))))
rho = cr*(2*n - x)*(2*n + 1 - x)*(x + n + 1)/(x + n + Q(3, 2))
R = (8*n**2 + 6*n - (2*n + 3)*x - x**2)/((3*n + 1)*(3*n + 2))

# (a)
assert sp.cancel(sp.together(rho - 1 - (R.subs(x, x + 1)*sigma - R))) == 0
print("(a) WZ equation rho - 1 = R(n,x+1) sigma - R(n,x): identically 0: OK")

import mpmath as mp
mp.mp.dps = 70
def c(N): return (mp.mpf(1)/9)**N*mp.gamma(3*N + mp.mpf(3)/2)*mp.factorial(N)**2/(4*mp.factorial(2*N)*mp.factorial(3*N)*mp.gamma(N + mp.mpf(5)/6)*mp.gamma(N + mp.mpf(7)/6))
def F(N, X): return c(N)*(mp.mpf(3)/4)**X*mp.rf(-X, 2*N)*mp.gamma(X + N + 1)/mp.gamma(X + N + mp.mpf(3)/2)
def Rf(N, X): return (8*N**2 + 6*N - (2*N + 3)*X - X**2)/mp.mpf((3*N + 1)*(3*N + 2))
def G(N, X): return Rf(N, X)*F(N, X)
# closed-form ratios agree with the Gamma expressions, WZ equation numerically at real x
for (N, X) in [(0, mp.mpf('0.37')), (1, mp.mpf('2.6')), (3, mp.mpf('9.21')), (2, mp.mpf('1.5'))]:
    assert abs(F(N, X + 1)/F(N, X) - sigma.subs({n: N, x: X})) < mp.mpf(10)**-50
    assert abs(F(N + 1, X)/F(N, X) - rho.subs({n: N, x: X})) < mp.mpf(10)**-50
    assert abs(F(N + 1, X) - F(N, X) - (G(N, X + 1) - G(N, X))) < mp.mpf(10)**-50
print("    ratios sigma, rho agree with the Gamma closed form; WZ holds numerically at non-integer x: OK")

# (b)  exact: G'(n,0) = -R(n,0) c_n (2n-1)! n!/Gamma(n+3/2)  (n>=1),  G'(0,0) = -(3/2) c_0/Gamma(3/2)
#      claim: G'(n,0) = -(3/(4 pi)) t(n+1).  Check n = 0 exactly and the ratio identity for n >= 1 -> n+1.
c0 = sp.gamma(Q(3, 2))/(4*sp.gamma(Q(5, 6))*sp.gamma(Q(7, 6)))
t = lambda m: Q((4*m - 1)*3**m, (2*m - 1)*m**2*sp.binomial(2*m, m)*sp.binomial(3*m, m))
lhs0 = sp.simplify(-Q(3, 2)*c0/sp.gamma(Q(3, 2)))            # = -3/(8 Gamma(5/6)Gamma(7/6)) = -9/(8 pi)
assert abs(sp.N(lhs0 + 3/(4*sp.pi)*t(1), 50)) < 1e-45
assert abs(sp.N(sp.gamma(Q(5, 6))*sp.gamma(Q(7, 6)) - sp.pi/3, 60)) < 1e-55   # = (1/6) Gamma(1/6)Gamma(5/6) = pi/3 (reflection)
# ratio for n >= 1:  G'(n+1,0)/G'(n,0) = [R(n+1,0)/R(n,0)] * cr * (2n+1)(2n)(n+1)/(n+3/2)
gr = sp.cancel(R.subs({n: n + 1, x: 0})/R.subs(x, 0)*cr*(2*n + 1)*(2*n)*(n + 1)/(n + Q(3, 2)))
m = n + 1
tr = sp.cancel((4*m + 3)*3*(2*m - 1)*m**2/((4*m - 1)*(2*m + 1)*(m + 1)**2)
               * (m + 1)**2/((2*m + 1)*(2*m + 2)) * (m + 1)*(2*m + 1)*(2*m + 2)/((3*m + 1)*(3*m + 2)*(3*m + 3)))
assert sp.cancel(gr - tr) == 0
dG = lambda N: mp.diff(lambda X: G(N, X), 0)
tm = lambda mm: mp.mpf((4*mm - 1)*3**mm)/((2*mm - 1)*mm**2*mp.binomial(2*mm, mm)*mp.binomial(3*mm, mm))
for N in range(0, 5):
    assert abs(dG(N) + 3/(4*mp.pi)*tm(N + 1)) < mp.mpf(10)**-40
print("(b) dG/dx(n,0) = -(3/(4 pi)) t(n+1): n=0 value, exact ratio identity, numeric n<5: OK")

# (c),(d),(e)
mp.mp.dps = 40
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
def D(N, kmax=500): return mp.fsum(mp.diff(lambda X: F(N, X), kk) for kk in range(kmax))
D0 = D(0)
I = mp.quad(lambda u: 8*mp.log(3*(1 - u*u)/4)/(1 + 3*u*u), mp.linspace(0, 1, 11))   # after t = 1-u^2
assert abs(I + 4*K) < mp.mpf(10)**-30 and abs(D0 - 3/(8*mp.pi)*I) < mp.mpf(10)**-30
print("(d) I = -4K and D_0 = (3/(8 pi)) I = -3K/(2 pi):", mp.nstr(D0, 25), mp.nstr(-3*K/(2*mp.pi), 25))
def Psi(N, y): return -mp.factorial(2*N)*y**(2*N)/(1 - y)**(2*N + 1)*mp.quad(lambda s: s**(2*N)/(1 + s), [0, (1 - y)/y])
def D_int(N): return c(N)/mp.sqrt(mp.pi)*mp.quad(lambda tt: tt**N*(1 - tt)**(-mp.mpf(1)/2)*Psi(N, 3*tt/4), mp.linspace(0, 1, 9))
Ds = [D0]
for N in (1, 2, 3):
    DN = D(N); Ds.append(DN)
    bound = c(N)*mp.factorial(2*N)*4*mp.beta(N, mp.mpf(1)/2)/(3*mp.sqrt(mp.pi)*(2*N + 1))
    assert abs(DN - D_int(N)) < mp.mpf(10)**-20 and abs(DN) <= bound
    assert abs(Ds[N] - Ds[N - 1] + dG(N - 1)) < mp.mpf(10)**-28          # (c)
    print("(c),(e) N=%d: D_N = %s = integral formula; |D_N| <= bound %s; D_N - D_{N-1} = -dG(N-1,0): OK" % (N, mp.nstr(DN, 12), mp.nstr(bound, 6)))
for y in (mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.74')):
    for N in (1, 4):
        # Psi_N(y) = y^{2N} d^{2N}/dy^{2N} [ log y/(1-y) ]
        direct = y**(2*N)*mp.diff(lambda z: mp.log(z)/(1 - z), y, 2*N)
        assert abs(direct - Psi(N, y)) < mp.mpf(10)**-12*abs(direct) and abs(Psi(N, y)) <= mp.factorial(2*N)/((2*N + 1)*y)
print("(e) closed form and bound for Psi_N: OK")

# (f)
mp.mp.dps = 120
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
S = mp.nsum(lambda j: mp.mpf((4*int(j) - 1)*3**int(j))/((2*int(j) - 1)*int(j)**2*mp.binomial(2*int(j), int(j))*mp.binomial(3*int(j), int(j))), [1, mp.inf])
print("(f) |lhs - 2K| =", mp.nstr(abs(S - 2*K), 5))
assert abs(S - 2*K) < mp.mpf(10)**-110
print("ALL CHECKS PASSED")
