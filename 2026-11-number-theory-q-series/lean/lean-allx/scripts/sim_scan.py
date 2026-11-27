#!/usr/bin/env python3
"""Simulate the integer scan: V(k) = sum_{p<=k} lo(p) (k-p) with
lo(p) = a(p) for p = 3 mod 4, -(a(p)+1) for p = 1 mod 4, a(p) = floor(M log2 p).
Reports min of V(k) / k^{3/2} (in units log2/M) and whether V(k) > 0 for 4 <= k <= N."""
import sys, math
import numpy as np
N = int(float(sys.argv[1])); M = int(sys.argv[2])
sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(N ** 0.5) + 1):
    if sieve[i]: sieve[i * i::i] = False
up = um = vp = vm = 0
bad = 0; worst = None
for n in range(N + 1):
    vp += up; vm += um
    if sieve[n]:
        a = (n ** M).bit_length() - 1
        r = n % 4
        if r == 3: up += a
        elif r == 1: um += a + 1
    if n >= 4:
        if not vm < vp:
            bad += 1
            if bad < 5: print("FAIL at", n, vp - vm)
        g = (vp - vm) * math.log(2) / M / n ** 1.5
        if worst is None or g < worst[0]: worst = (g, n)
        if n in (10**4, 10**5, 10**6, 10**7): print(n, g, flush=True)
print("N", N, "M", M, "fails", bad, "min G_lower", worst)
