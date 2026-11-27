#!/usr/bin/env python3
"""Theorem 11 (B1), exact rational arithmetic, r(k) = sigma(N,k) from the theorem."""
from fractions import Fraction as Fr
import sys

def rf4(z): return z*(z + 1)*(z + 2)*(z + 3)
def r(N, k):
    return Fr((2*k + 3)*(N - k - 1)**2)*rf4(Fr(8*k + 4*N + 1, 2))/((2*k + 1)*(N + k + 1)**2*rf4(Fr(8*k - 4*N + 9, 2)))

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
fail = {"pos": [], "logconc": [], "below": [], "ineq": [], "poly": [], "k0case": [], "odd": [], "stated_prod": []}
minmargin_below = None
for N in range(1, NMAX + 1):
    # I = {k integer : N/2 - 1 <= k <= N-1}
    minI = max(0, -((-(N - 2))//2))          # ceil(N/2 - 1), clipped at 0
    rr = [r(N, k) for k in range(N)]
    # positivity of a(N,k), 0 <= k < N  <=>  r(k) > 0 for 0 <= k <= N-2 (a(N,0) > 0 checked in t11_partA.py)
    if not all(rr[k] > 0 for k in range(N - 1)): fail["pos"].append(N)
    if rr[N - 1] != 0: fail["pos"].append(("last", N))
    # log-concavity on I: r nonincreasing for consecutive k, k+1 in I (r(k) defined for k <= N-1)
    for k in range(minI, N - 1):
        if not rr[k + 1] < rr[k]: fail["logconc"].append((N, k))
    # form of r as stated: (1+2/(2k+1)) ((N-k-1)/(N+k+1))^2 prod (1 + (4N-4)/(4k-2N+9/2+i)); on I: 4k-2N+9/2 >= 1/2
    for k in range(minI, N):
        z = Fr(8*k - 4*N + 9, 2)
        if z < Fr(1, 2): fail["stated_prod"].append((N, k, "z"))
        p = (1 + Fr(2, 2*k + 1))*Fr(N - k - 1, N + k + 1)**2
        for i in range(4): p *= 1 + (4*N - 4)/(z + i)
        if p != rr[k]: fail["stated_prod"].append((N, k))
    # r(k) >= 1 for ALL k < min I
    for k in range(0, minI):
        if rr[k] < 1: fail["below"].append((N, k))
        m = rr[k] - 1
        if N >= 4 and (minmargin_below is None or m < minmargin_below[0]): minmargin_below = (m, N, k)
    # the text's route: 0 <= k <= N/2 - 2
    for k in range(0, N//2 + 1):
        if 2*k <= N - 4:
            z = Fr(8*k - 4*N + 9, 2)
            if not all(z + i < 0 for i in range(4)): fail["ineq"].append((N, k, "sign"))
            bound = (Fr((N - k - 1)*(N + 2*k)**2, (N + k + 1)*(N - 2*k)**2))**2
            if rr[k] < bound: fail["ineq"].append((N, k))
            if k == 0:
                b0 = 3*Fr(N - 1, N + 1)**2
                if not (rr[0] >= b0 >= 1): fail["k0case"].append(N)
        if 1 <= k and 2*k <= N:
            pol = 6*N*N*k - 8*k**3 - 2*N*N - 8*k*k
            if pol != (N - k - 1)*(N + 2*k)**2 - (N + k + 1)*(N - 2*k)**2: fail["poly"].append((N, k, "id"))
            if pol < 0: fail["poly"].append((N, k))
    if N % 2 == 1 and N >= 3:
        k = (N - 3)//2
        z = Fr(8*k - 4*N + 9, 2)
        if [z + i for i in range(4)] != [Fr(-3, 2), Fr(-1, 2), Fr(1, 2), Fr(3, 2)] or rr[k] < 1: fail["odd"].append(N)
        # is this index covered by neither I nor k <= N/2-2 ?
        assert not (2*k <= N - 4) and k < minI
for key, v in fail.items():
    print("%-12s failures: %s" % (key, v[:20]), "(total %d)" % len(v))
print("smallest r(k)-1 below I for 4 <= N <= %d: %.6g at N=%d, k=%d" % (NMAX, float(minmargin_below[0]), minmargin_below[1], minmargin_below[2]))
