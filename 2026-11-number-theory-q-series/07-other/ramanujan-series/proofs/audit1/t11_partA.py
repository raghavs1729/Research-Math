#!/usr/bin/env python3
"""Independent exact (Fraction) re-check of Theorem 11, Part A.
a(n,k) is rational:  C_n = 2^(6-8n) (3/4)_n^2 (1/4)_n / ( -4 (-1/4)_n (2n-2)! )   [Gamma(1/4)Gamma(3/4) = pi sqrt2, Gamma(-1/4) = -4 Gamma(3/4)]
Gamma(4k+2n+1/2)/Gamma(4k-2n+9/2) = prod_{i=0}^{4n-5} (4k-2n+9/2+i)  (no poles: half-integers)."""
from fractions import Fraction as Fr
from math import factorial, comb

def rf(a, m):
    p = Fr(1)
    for i in range(m): p *= (a + i)
    return p

def C(n):
    return Fr(2)**(6 - 8*n)*rf(Fr(3, 4), n)**2*rf(Fr(1, 4), n)/(-4*rf(Fr(-1, 4), n)*factorial(2*n - 2))

def a(n, k):
    if k >= n: return Fr(0)
    return C(n)*(2*k + 1)*rf(Fr(8*k - 4*n + 9, 2), 4*n - 4)/(factorial(n - k - 1)**2*factorial(n + k)**2)

def NUM(n, k):
    return ((3145728*n + 1310720)*k**6 + (-6291456*n**2 + 3670016*n + 2621440)*k**5 + (-786432*n**3 - 11337728*n**2 - 344064*n + 1724416)*k**4
       + (9437184*n**4 + 786432*n**3 - 8716288*n**2 - 1884160*n + 376832)*k**3
       + (-8060928*n**5 + 10010624*n**4 + 3973120*n**3 - 3434496*n**2 - 892672*n - 15936)*k**2
       + (2752512*n**6 - 6324224*n**5 + 1916928*n**4 + 2363392*n**3 - 525824*n**2 - 170112*n - 12672)*k
       + (-344064*n**7 + 1134592*n**6 - 922624*n**5 - 195328*n**4 + 318016*n**3 - 1680*n**2 - 10044*n - 945))

def R(n, k):
    return Fr(NUM(n, k), 131072*n*(2*n - 1)*(2*n + 1)*(4*n - 1)*(k - n)**2*(2*k + 1))

def gam(n, k):
    if k > n: return Fr(0)
    if k == n: return -a(n + 1, n)
    return R(n, k)*a(n, k)

def t(m):
    return Fr((1344*m**3 + 944*m**2 + 156*m + 9)*comb(2*m, m)**2*comb(4*m, 2*m)**3, (2*m + 1)*2**(22*m))

# positivity
for n in range(1, 61):
    for k in range(n):
        assert a(n, k) > 0, (n, k)
print("a(n,k) > 0 for 1 <= n <= 60, 0 <= k < n: OK")

# ratio formulas of the text
def sigma(n, k):
    return Fr((2*k + 3)*(n - k - 1)**2)*rf(Fr(8*k + 4*n + 1, 2), 4)/((2*k + 1)*(n + k + 1)**2*rf(Fr(8*k - 4*n + 9, 2), 4))
def rho(n, k):
    return ((n + Fr(3, 4))**2*(n + Fr(1, 4))*(4*k + 2*n + Fr(1, 2))*(4*k + 2*n + Fr(3, 2))*(4*k - 2*n + Fr(5, 2))*(4*k - 2*n + Fr(7, 2))
            / (512*n*(n - Fr(1, 4))*(2*n - 1)*(n - k)**2*(n + k + 1)**2))
for n in range(1, 41):
    for k in range(n):
        assert a(n, k + 1) == sigma(n, k)*a(n, k), (n, k)     # includes k = n-1: sigma = 0
        assert a(n + 1, k) == rho(n, k)*a(n, k), (n, k)
print("sigma, rho of the text agree with a(n,k) for n <= 40 (all k < n, incl. k = n-1): OK")

# WZ equation, all four cases, on a big grid
bad = []
for n in range(1, 41):
    for k in range(0, 46):
        if a(n + 1, k) - a(n, k) != gam(n, k + 1) - gam(n, k): bad.append((n, k))
print("WZ equation a(n+1,k)-a(n,k) = gam(n,k+1)-gam(n,k) on 1<=n<=40, 0<=k<=45: failures =", bad)

# boundary
assert a(1, 0) == t(0)/256
for n in range(1, 41):
    assert gam(n, 0) == -t(n)/256, n
print("a(1,0) = t(0)/256 and gam(n,0) = -t(n)/256 for n = 1..40: OK")

# finite identity (A), N = 1..40; print N = 1..6
for N in range(1, 41):
    S = sum(t(m) for m in range(N)); T = sum(a(N, k) for k in range(N))
    assert S == 256*T, N
    if N <= 6: print("N=%d: sum_{n<N} t(n) = 256 T_N = %s  (~%.15f)" % (N, S, float(S)))
import math
print("64 sqrt2/pi^2 = %.15f" % (64*math.sqrt(2)/math.pi**2))
# telescoping bookkeeping exactly as in the text: n = 1..N-1, k = 0..N
for N in range(1, 12):
    lhs = sum(a(n + 1, k) - a(n, k) for n in range(1, N) for k in range(0, N + 1))
    rhs = sum(gam(n, k + 1) - gam(n, k) for n in range(1, N) for k in range(0, N + 1))
    assert lhs == rhs
    assert lhs == sum(a(N, k) for k in range(N)) - a(1, 0)
    assert rhs == -sum(gam(n, 0) for n in range(1, N))
    assert all(gam(n, N + 1) == 0 for n in range(1, N))
print("telescoping bookkeeping (n=1..N-1, k=0..N) for N <= 11: OK")
