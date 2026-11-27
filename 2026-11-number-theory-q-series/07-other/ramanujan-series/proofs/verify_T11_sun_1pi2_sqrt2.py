#!/usr/bin/env python3
"""
THEOREM 11 (Z.-W. Sun, arXiv:2307.03086v10, Conjecture 5.33(i), eq. (5.31), dated 2025-06-26):

   sum_{k>=0} (1344k^3+944k^2+156k+9) binom(2k,k)^2 binom(4k,2k)^3 / ((2k+1) 2^(22k))  =  64 sqrt(2)/pi^2.

Structure of the proof (RESULTS.md, Theorem 11):
  PART A (algebraic, exact).  A WZ pair (a, gam) supported on the triangle 0 <= k < n:
     a(n,k) = C_n (2k+1) Gamma(4k+2n+1/2) / ( (n-k-1)!^2 (n+k)!^2 Gamma(4k-2n+9/2) )   (0 <= k < n),   a(n,k) = 0 (k >= n),
     C_n    = 2^(6-8n) Gamma(n+3/4)^2 Gamma(n+1/4) / ( sqrt2 pi Gamma(n-1/4) Gamma(2n-1) ),
     gam(n,k) = R(n,k) a(n,k) (k < n),  gam(n,n) = -a(n+1,n),  gam(n,k) = 0 (k > n),
     R = NUM / (131072 n (2n-1)(2n+1)(4n-1)(k-n)^2 (2k+1)),
  with a(n+1,k) - a(n,k) = gam(n,k+1) - gam(n,k) for all n >= 1, k >= 0, and gam(n,0) = -t(n)/256 (n >= 1),
  a(1,0) = t(0)/256.  Telescoping (n from 1 to N-1, k from 0 to N):   sum_{n<N} t(n) = 256 T_N,  T_N := sum_{k<N} a(N,k).
  (Origin: eps^-2 coefficient of Seed9(1-2nu, 1/4-nu, 3/4-nu, k+nu), nu = n + eps; Au's Seed9.)
  PART B (asymptotic).  T_N -> sqrt2/(4 pi^2): Laplace's method for the log-concave sequence k -> a(N,k).
This script checks PART A exactly and the ingredients of PART B exactly (algebra) and numerically.
"""
import sympy as sp
from sympy import Rational as Q, binomial, gamma, sqrt, pi, factorial
n, k = sp.symbols('n k')
NUM = ((3145728*n + 1310720)*k**6 + (-6291456*n**2 + 3670016*n + 2621440)*k**5 + (-786432*n**3 - 11337728*n**2 - 344064*n + 1724416)*k**4
       + (9437184*n**4 + 786432*n**3 - 8716288*n**2 - 1884160*n + 376832)*k**3
       + (-8060928*n**5 + 10010624*n**4 + 3973120*n**3 - 3434496*n**2 - 892672*n - 15936)*k**2
       + (2752512*n**6 - 6324224*n**5 + 1916928*n**4 + 2363392*n**3 - 525824*n**2 - 170112*n - 12672)*k
       + (-344064*n**7 + 1134592*n**6 - 922624*n**5 - 195328*n**4 + 318016*n**3 - 1680*n**2 - 10044*n - 945))
R = NUM/(131072*n*(2*n-1)*(2*n+1)*(4*n-1)*(k-n)**2*(2*k+1))
rf4 = lambda z: z*(z+1)*(z+2)*(z+3)
sigma = (2*k+3)*(n-k-1)**2*rf4(4*k+2*n+Q(1,2))/((2*k+1)*(n+k+1)**2*rf4(4*k-2*n+Q(9,2)))
Cratio = (n+Q(3,4))**2*(n+Q(1,4))/(256*(n-Q(1,4))*(2*n-1)*(2*n))
rho = Cratio*(4*k+2*n+Q(1,2))*(4*k+2*n+Q(3,2))*(4*k-2*n+Q(5,2))*(4*k-2*n+Q(7,2))/((n-k)**2*(n+k+1)**2)

def a(nn, kk):
    if kk >= nn: return sp.Integer(0)
    C = Q(2)**(6-8*nn)*gamma(nn+Q(3,4))**2*gamma(nn+Q(1,4))/(sqrt(2)*pi*gamma(nn-Q(1,4))*gamma(2*nn-1))
    return C*(2*kk+1)*gamma(4*kk+2*nn+Q(1,2))/(factorial(nn-kk-1)**2*factorial(nn+kk)**2*gamma(4*kk-2*nn+Q(9,2)))
def gam(nn, kk):
    if kk > nn: return sp.Integer(0)
    if kk == nn: return -a(nn+1, nn)
    return R.subs({n: nn, k: kk})*a(nn, kk)
z = lambda e: sp.simplify(sp.gammasimp(sp.expand_func(e)))
# A1. ratios
for nn in range(1, 5):
    for kk in range(0, nn):
        if kk+1 < nn: assert z(a(nn, kk+1)/a(nn, kk) - sigma.subs({n: nn, k: kk})) == 0
        assert z(a(nn+1, kk)/a(nn, kk) - rho.subs({n: nn, k: kk})) == 0
print("A1. sigma, rho agree with the closed form of a(n,k): OK")
# A2. WZ: generic identity, boundary k = n-1, all cases on a grid
assert sp.cancel(sp.together(rho - 1 - (R.subs(k, k+1)*sigma - R))) == 0
lhs = rho.subs(k, n-1) - 1
rhs = -rho.subs(k, n-1)*sigma.subs({n: n+1, k: n-1}, simultaneous=True) - R.subs(k, n-1)
assert sp.cancel(sp.together(lhs - rhs)) == 0
for nn in range(1, 6):
    for kk in range(0, 8):
        assert z(a(nn+1, kk) - a(nn, kk) - (gam(nn, kk+1) - gam(nn, kk))) == 0
print("A2. WZ equation: generic identity, boundary identity (k = n-1), all cases on the grid 1<=n<6, k<8: OK")
# A3. gam(n,0) = -t(n)/256 (n>=1), a(1,0) = t(0)/256
t = lambda m: Q((1344*m**3+944*m**2+156*m+9)*binomial(2*m, m)**2*binomial(4*m, 2*m)**3, (2*m+1)*2**(22*m))
assert z(a(1, 0) - t(0)/256) == 0
for nn in range(1, 9):
    assert z(gam(nn, 0) + t(nn)/256) == 0
tr = ((1344*(n+1)**3+944*(n+1)**2+156*(n+1)+9)*(2*n+1)/((1344*n**3+944*n**2+156*n+9)*(2*n+3)*2**22)
      * ((2*n+1)*(2*n+2)/(n+1)**2)**2 * ((4*n+1)*(4*n+2)*(4*n+3)*(4*n+4)/((2*n+1)*(2*n+2))**2)**3)
gr = sp.cancel(R.subs({n: n+1, k: 0})/R.subs(k, 0)*rho.subs(k, 0))
assert sp.cancel(gr - tr) == 0
print("A3. gam(n,0) = -t(n)/256 (n = 1..8 exactly + ratio identity), a(1,0) = t(0)/256: OK")
for N in (1, 2, 3, 5):
    assert z(sum(t(m) for m in range(N)) - 256*sum(a(N, kk) for kk in range(N))) == 0
print("A4. finite identity sum_{n<N} t(n) = 256 T_N for N = 1,2,3,5 (exact): OK")

# ---------------- PART B: ingredients
x = sp.symbols('x', positive=True)
phi = (-6*sp.log(2) + 2*(2*x+1)*sp.log(2*x+1) - 2*(1-x)*sp.log(1-x) - 2*(1+x)*sp.log(1+x) - 2*(2*x-1)*sp.log(2*x-1))
xs = sp.sqrt(3)/2
assert sp.simplify(sp.exp(sp.diff(phi, x)) - ((1-x)**2*(2*x+1)**4/((1+x)**2*(2*x-1)**4))) == 0      # phi' = log rho_inf
assert sp.simplify(((1-x)*(2*x+1)**2 - (1+x)*(2*x-1)**2).subs(x, xs)) == 0                             # phi'(x*) = 0
assert sp.simplify(sp.diff(phi, x, 2).subs(x, xs) + 24) == 0                                           # phi''(x*) = -24
assert abs(sp.N(phi.subs(x, xs), 60)) < 1e-55                                                          # phi(x*) = 0
A = 2**Q(5,2)*pi**Q(-5,2)*x*(1-x)/((1+x)*(2*x-1)**4)
assert sp.simplify(A.subs(x, xs) - sp.sqrt(3)/(sp.sqrt(2)*pi**Q(5,2))) == 0
assert sp.simplify(A.subs(x, xs)*sp.sqrt(pi/12) - sp.sqrt(2)/(4*pi**2)) == 0
print("B1. phi' = log[(1-x)^2(2x+1)^4/((1+x)^2(2x-1)^4)], phi'(x*) = 0 at x* = sqrt3/2, phi''(x*) = -24, phi(x*) = 0,")
print("    psi0(x*) = sqrt3/(sqrt2 pi^(5/2)), psi0(x*) sqrt(pi/12) = sqrt2/(4 pi^2): OK")
# B1 (exact). shape of k -> a(N,k) for every N from 4 to 400, in rational arithmetic, from r(k) = sigma(N,k):
#   a > 0 (r > 0), r strictly decreasing on I = {N/2-1 <= k <= N-1}, r >= 1 for every k < min I,
#   and r(k) >= ((N-k-1)(N+2k)^2/((N+k+1)(N-2k)^2))^2 for 1 <= k <= N/2-2.
from fractions import Fraction as Fr
def rr(N, kk):
    num = Fr(2*kk+3)*(N-kk-1)**2; den = Fr(2*kk+1)*(N+kk+1)**2
    for i in range(4):
        num *= Fr(8*kk+4*N+1+2*i, 2); den *= Fr(8*kk-4*N+9+2*i, 2)
    return num/den
for N in range(4, 401):
    rs = [rr(N, kk) for kk in range(0, N-1)]
    assert all(v > 0 for v in rs)
    kmin = -(-(N-2)//2)                      # ceil(N/2 - 1)
    assert all(rs[i+1] < rs[i] for i in range(kmin, N-2))
    assert all(rs[i] >= 1 for i in range(0, kmin))
    for kk in range(1, (N-4)//2 + 1):        # k <= N/2 - 2
        assert rs[kk] >= (Fr((N-kk-1)*(N+2*kk)**2, (N+kk+1)*(N-2*kk)**2))**2 >= 1
print("B1. (exact, N = 4..400) a > 0, log-concave on I = [N/2-1, N-1], r >= 1 below I, lower bound for r: OK")
# B2. polynomial inequality used for monotonicity below N/2:  (N-k-1)(N+2k)^2 - (N+k+1)(N-2k)^2 = 6N^2 k - 8k^3 - 2N^2 - 8k^2
Nn = sp.symbols('N')
assert sp.expand((Nn-k-1)*(Nn+2*k)**2 - (Nn+k+1)*(Nn-2*k)**2 - (6*Nn**2*k - 8*k**3 - 2*Nn**2 - 8*k**2)) == 0
print("B2. algebraic identity for the region k < N/2: OK")
import mpmath as mp
mp.mp.dps = 30
def la(N, kk):   # log a(N,k)
    N = mp.mpf(N); kk = mp.mpf(kk); lg = mp.loggamma
    logC = (6-8*N)*mp.log(2) - mp.log(mp.sqrt(2)*mp.pi) + 2*lg(N+mp.mpf(3)/4) + lg(N+mp.mpf(1)/4) - lg(N-mp.mpf(1)/4) - lg(2*N-1)
    return logC + mp.log(2*kk+1) + lg(4*kk+2*N+mp.mpf(1)/2) - 2*lg(N-kk) - 2*lg(N+kk+1) - mp.re(lg(4*kk-2*N+mp.mpf(9)/2))
target = mp.sqrt(2)/(4*mp.pi**2)
for N in (40, 200, 1000, 5000, 20000):
    vals = [la(N, kk) for kk in range(N)]
    T = mp.fsum(mp.exp(v) for v in vals)
    k0 = int(mp.floor(mp.sqrt(3)/2*N))
    # floating-point re-check away from the junction (the exact check is B1 above)
    kI = int(N//2 + 2)
    d = [vals[i+1] - vals[i] for i in range(N-1)]
    assert all(d[i+1] < d[i] for i in range(kI, N-2))
    assert all(d[i] > 0 for i in range(0, N//2 - 2))
    lap = mp.exp(vals[k0])*mp.sqrt(mp.pi*N/12)
    print("B3. N=%6d: T_N - sqrt2/(4pi^2) = %s;  a(k0) sqrt(N)/psi0(x*) - 1 = %s;  T_N/(a(k0) sqrt(pi N/12)) - 1 = %s"
          % (N, mp.nstr(T - target, 3), mp.nstr(mp.exp(vals[k0])*mp.sqrt(N)/(mp.sqrt(3)/(mp.sqrt(2)*mp.pi**2.5)) - 1, 3), mp.nstr(T/lap - 1, 3)))
mp.mp.dps = 120
tm = lambda m: mp.mpf(1344*m**3+944*m**2+156*m+9)*mp.binomial(2*m, m)**2*mp.binomial(4*m, 2*m)**3/((2*m+1)*mp.mpf(2)**(22*m))
S = mp.fsum(tm(m) for m in range(0, 80))
print("final: |lhs - 64 sqrt2/pi^2| =", mp.nstr(abs(S - 64*mp.sqrt(2)/mp.pi**2), 5))
assert abs(S - 64*mp.sqrt(2)/mp.pi**2) < mp.mpf(10)**-110
print("ALL CHECKS PASSED")
