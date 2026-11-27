import math, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 30
def logr(N, k):
    s = math.log1p(2.0/(2*k + 1)) + 2*math.log((N - k - 1)/(N + k + 1))
    z = 4*k - 2*N + 4.5
    for i in range(4): s += math.log1p((4*N - 4)/(z + i))
    return s
bad_window = bad_up = bad_dn = None
for N in range(4, 20001):
    k0 = math.floor(math.sqrt(3)/2*N); J = math.floor(N**0.6)
    if k0 + J > N - 2: bad_window = N; continue
    if not logr(N, k0 + J) <= -12*J/N: bad_up = N
    if not logr(N, k0 - J - 1) >= 12*J/N: bad_dn = N
print("last N<=20000 with k0+J > N-2:", bad_window, "; last N where r(k0+J) <= e^{-12J/N} fails (window inside):", bad_up, "; last N where r(k0-J-1) >= e^{12J/N} fails:", bad_dn)
# 1/(1-e^{-u}) <= 2/u, u = 12J/N  (needs u <= 1.59)
print("u = 12J/N <= 1.59 from N =", next(N for N in range(4, 10**6) if all(12*math.floor(M**0.6)/M <= 1.59 for M in range(N, N + 2000))))
# Theorem 4: direct summation of sum_k h(n,k)/h(n,0)
for n in (1, 2, 40, 200):
    p = mp.mpf(1); S = mp.mpf(0)
    for k in range(0, 200000):
        S += p; p *= mp.mpf((4*k + 6*n + 1)*(4*k + 6*n + 3))/(16*(k + n + 1)*(k + 3*n + 1))
        if p < mp.mpf(10)**-28*S and k > 50: break
    print("T4: n=%d  sum_k h(n,k)/h(n,0) ~ %s (k up to %d)" % (n, mp.nstr(S, 10), k))
