"""WP-A sub-claims (exact rational arithmetic / sympy), used by PiSeries/SunPi2/PartA.lean:
   1. NUM(n,0) = -(4n-7)(4n-5)(4n-3)(4n-1)(1344n^3+944n^2+156n+9)
   2. A n = a(n,0) = 512 n^3 (2n-1) C(2n,n)^2 C(4n,2n)^3 / (2^(22n) (4n-7)(4n-5)(4n-3)),  -R(n,0) A n = t(n)/256
   3. compatibility rho(n,k) sigma(n+1,k) = sigma(n,k) rho(n,k+1)  (k+1 < n)
   4. diagonal cf(k+1,k) (cancelled form), cf(k+2,k+1) = cf(k+1,k) rho(k+1,k) sigma(k+2,k)
   5. column step cf(N+1,k) = cf(N,k) rho(N,k) for k < N, N+1 <= 2k+2 (N = 2d+e+1, k = d+e)"""
import sympy as sp
from fractions import Fraction as Fr
from math import comb, factorial as f
n, k = sp.symbols('n k')
NUM = ((3145728*n + 1310720)*k**6 + (-6291456*n**2 + 3670016*n + 2621440)*k**5 + (-786432*n**3 - 11337728*n**2 - 344064*n + 1724416)*k**4
       + (9437184*n**4 + 786432*n**3 - 8716288*n**2 - 1884160*n + 376832)*k**3
       + (-8060928*n**5 + 10010624*n**4 + 3973120*n**3 - 3434496*n**2 - 892672*n - 15936)*k**2
       + (2752512*n**6 - 6324224*n**5 + 1916928*n**4 + 2363392*n**3 - 525824*n**2 - 170112*n - 12672)*k
       + (-344064*n**7 + 1134592*n**6 - 922624*n**5 - 195328*n**4 + 318016*n**3 - 1680*n**2 - 10044*n - 945))
assert sp.expand(NUM.subs(k, 0) + (4*n-7)*(4*n-5)*(4*n-3)*(4*n-1)*(1344*n**3+944*n**2+156*n+9)) == 0
print("1. factorisation of NUM(n,0): OK")
def sigma(n,k): return Fr((2*k+3)*(n-k-1)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k+4*n+5)*(8*k+4*n+7), (2*k+1)*(n+k+1)**2*(8*k-4*n+9)*(8*k-4*n+11)*(8*k-4*n+13)*(8*k-4*n+15))
def rho(n,k): return Fr((4*n+1)*(4*n+3)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k-4*n+5)*(8*k-4*n+7), 131072*n*(2*n-1)*(4*n-1)*(n-k)**2*(n+k+1)**2)
def t(m): return Fr((1344*m**3+944*m**2+156*m+9)*comb(2*m,m)**2*comb(4*m,2*m)**3, (2*m+1)*2**(22*m))
A = {1: Fr(9,256)}
for m in range(1, 60): A[m+1] = A[m]*rho(m,0)
for m in range(1, 60):
    N0 = int(NUM.subs({n: m, k: 0})); RD = 131072*m*(2*m-1)*(2*m+1)*(4*m-1)*m**2
    assert A[m] == Fr(512*m**3*(2*m-1)*comb(2*m,m)**2*comb(4*m,2*m)**3, 2**(22*m)*(4*m-7)*(4*m-5)*(4*m-3))
    assert -Fr(N0, RD)*A[m] == t(m)/256
print("2. closed form of A n, -R(n,0) A n = t(n)/256 (n < 60): OK")
for nn in range(2, 40):
    for kk in range(0, nn-1):
        assert rho(nn,kk)*sigma(nn+1,kk) == sigma(nn,kk)*rho(nn,kk+1)
print("3. compatibility (k+1 < n < 40): OK")
def cf(N,k): return Fr(2**12*(4*N-1)*(2*k+1)*f(4*N)*f(8*k+4*N)*f(4*k-2*N+4), 2**(22*N)*f(2*N)*f(2*N-2)*f(4*k+2*N)*f(8*k-4*N+8)*f(N-k-1)**2*f(N+k)**2)
for kk in range(0, 30):
    assert cf(kk+1,kk) == Fr(2**12*(4*kk+3)*(2*kk+1)*f(12*kk+4), 2**(22*(kk+1))*f(2*kk)*f(6*kk+2)*f(2*kk+1)**2)
    assert cf(kk+2,kk+1) == cf(kk+1,kk)*rho(kk+1,kk)*sigma(kk+2,kk)
assert cf(1,0) == Fr(9,256)
print("4. diagonal of cf (k < 30): OK")
for d in range(0, 12):
    for e in range(0, 12):
        N, kk = 2*d+e+1, d+e
        assert cf(N+1,kk) == cf(N,kk)*rho(N,kk)
print("5. column step of cf (d, e < 12): OK")
