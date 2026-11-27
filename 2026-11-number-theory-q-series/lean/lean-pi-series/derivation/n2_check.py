"""Numerical/exact checks of the reformulated ingredients for N2 (Sun's 64 sqrt2/pi^2 series)."""
from fractions import Fraction as Fr
from math import comb, factorial as f
import mpmath as mp
mp.mp.dps = 30
def sigma(n,k): return Fr((2*k+3)*(n-k-1)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k+4*n+5)*(8*k+4*n+7), (2*k+1)*(n+k+1)**2*(8*k-4*n+9)*(8*k-4*n+11)*(8*k-4*n+13)*(8*k-4*n+15))
def rho(n,k): return Fr((4*n+1)*(4*n+3)**2*(8*k+4*n+1)*(8*k+4*n+3)*(8*k-4*n+5)*(8*k-4*n+7), 131072*n*(2*n-1)*(4*n-1)*(n-k)**2*(n+k+1)**2)
def t(m): return Fr((1344*m**3+944*m**2+156*m+9)*comb(2*m,m)**2*comb(4*m,2*m)**3, (2*m+1)*2**(22*m))
A = {1: Fr(9,256)}
for n in range(1, 400): A[n+1] = A[n]*rho(n,0)
def arow(n):
    r = [A[n]]
    for k in range(n-1): r.append(r[-1]*sigma(n,k))
    return r
# 1. column relation, T_N, closed form
rows = {n: arow(n) for n in range(1, 41)}
for n in range(1, 40):
    for k in range(n): assert rows[n+1][k] == rows[n][k]*rho(n,k)
    assert sum(rows[n]) == sum(t(m) for m in range(n))/256
    assert all(x > 0 for x in rows[n])
def cf(N,k): return Fr(2**12*(4*N-1)*(2*k+1)*f(4*N)*f(8*k+4*N)*f(4*k-2*N+4), 2**(22*N)*f(2*N)*f(2*N-2)*f(4*k+2*N)*f(8*k-4*N+8)*f(N-k-1)**2*f(N+k)**2)
for n in range(1, 41):
    for k in range(n):
        if 4*k-2*n+4 >= 0: assert rows[n][k] == cf(n,k), (n,k)
print("1. a(n+1,k)=rho a(n,k); T_N = sum_{m<N} t(m)/256; a>0; factorial closed form for 4k-2N+4>=0: OK")
# 2. limits
xs = mp.sqrt(3)/2
for N in (50, 100, 200, 399):
    r = arow(N); T = sum(r); k0 = int(mp.floor(xs*N))
    print(N, 'T_N', mp.nstr(mp.mpf(T.numerator)/T.denominator, 12), 'target', mp.nstr(mp.sqrt(2)/(4*mp.pi**2), 12),
          ' a(k0)sqrtN/psi0', mp.nstr(mp.mpf(r[k0].numerator)/r[k0].denominator*mp.sqrt(N)/(mp.sqrt(3)/(mp.sqrt(2)*mp.pi**2.5)), 8),
          ' -N log sigma(k0+5)/(24*(k0+5-x*N))', mp.nstr(-N*mp.log(mp.mpf(sigma(N,k0+5).numerator)/sigma(N,k0+5).denominator)/(24*(k0+5-xs*N)), 6))
# 3. rho bounds: below (x*-d)n rho<1, above >1
for n in (30, 100, 300):
    lo = max(float(rho(n,k)) for k in range(0, int((float(xs)-0.02)*n)+1)); hi = min(float(rho(n,k)) for k in range(int((float(xs)+0.02)*n)+1, n))
    print(n, 'max rho below window', lo, ' min rho above window', hi)
