#!/usr/bin/env python3
"""Exact Ramanujan tau(n), n <= N, from Delta = q * prod (1-q^n)^24 = q * (eta3)^8,
eta3 = prod (1-q^n)^3 = sum_{k>=0} (-1)^k (2k+1) q^{k(k+1)/2}  (Jacobi).
Exact integer arithmetic (FLINT fmpz_poly).  Returns list tau[0..N], tau[0]=0."""
import sys, time
from flint import fmpz_poly

def tau_list(N):
    c = [0]*N                      # coefficients of eta3 mod q^N
    k = 0
    while k*(k+1)//2 < N:
        c[k*(k+1)//2] = (2*k+1) if k % 2 == 0 else -(2*k+1)
        k += 1
    e3 = fmpz_poly(c)
    e6 = e3.mul_low(e3, N)
    e12 = e6.mul_low(e6, N)
    e24 = e12.mul_low(e12, N)      # coefficient of q^(n-1) is tau(n)
    co = e24.coeffs()
    co += [0]*(N-len(co))
    return [0] + [int(v) for v in co]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**5
    t0 = time.time()
    t = tau_list(N)
    print("N =", N, "time", round(time.time()-t0, 2), "s")
    print("tau(1..12) =", t[1:13])
    assert t[1:13] == [1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944]
    # Ramanujan congruence tau(n) = sigma_11(n) mod 691 at a few n; multiplicativity checks
    def sigma11(n):
        return sum(d**11 for d in range(1, n+1) if n % d == 0)
    for n in [2,3,5,7,30,97,360,1001,4999]:
        if n <= N: assert (t[n]-sigma11(n)) % 691 == 0, n
    for (a,b) in [(2,3),(3,5),(7,11),(16,81),(101,103)]:
        if a*b <= N: assert t[a*b] == t[a]*t[b]
    for p in [2,3,5,7,11,13,101]:
        if p*p <= N: assert t[p*p] == t[p]**2 - p**11
    print("checks passed")
