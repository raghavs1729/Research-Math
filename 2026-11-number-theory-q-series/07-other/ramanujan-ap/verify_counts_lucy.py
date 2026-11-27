#!/usr/bin/env python3
"""Independent recomputation of pi(x;q,a), q = 3, 4, at given x, by a combinatorial method
(no sieve to x): pi(x) from primecount (Deleglise-Rivat/Gourdon), and
S_chi(x) = sum_{p <= x} chi(p) by the "Lucy" dynamic programme for a completely multiplicative
function (O(x^{3/4}) integer operations, exact int64 arithmetic in numpy):
   S(v, p) = S(v, p-1) - chi(p) [ S(v/p, p-1) - S(p-1, p-1) ]   for v >= p^2,
over v in { floor(x/k) }, starting from S(v, 1) = sum_{2 <= k <= v} chi(k).
Then pi(x;4,1) = (pi(x) - 1 + S_chi4)/2, pi(x;4,3) = (pi(x) - 1 - S_chi4)/2 (the prime 2 excluded),
     pi(x;3,1) = (pi(x) - 1 + S_chi3)/2, pi(x;3,2) = (pi(x) - 1 - S_chi3)/2 (the prime 3 excluded).
Usage: verify_counts_lucy.py x [x ...]"""
import sys, subprocess, math
import numpy as np
def chi_tables(q):
    if q == 4: f = {0: 0, 1: 1, 2: 0, 3: -1}
    else: f = {0: 0, 1: 1, 2: -1}
    pre = [0]
    for r in range(1, q + 1): pre.append(pre[-1] + f[r % q])          # sum_{k=1}^{r} chi(k), r <= q
    return f, pre
def char_prime_sum(x, q):
    f, pre = chi_tables(q)
    F = lambda v: (v // q) * pre[q] + np.take(np.array(pre), v % q)   # sum_{k<=v} chi(k) (pre[q] = 0)
    r = math.isqrt(x)
    small = F(np.arange(0, r + 1, dtype=np.int64)) - 1                # S(v), v = 0..r   (index v); v=0 unused
    small[0] = 0
    ks = np.arange(1, r + 1, dtype=np.int64)
    large = F(x // ks) - 1                                            # S(x//k), k = 1..r (index k-1)
    sieve = np.ones(r + 1, dtype=bool); sieve[:2] = False
    for p in range(2, r + 1):
        if not sieve[p]: continue
        sieve[p * p::p] = False
        fp = f[p % q]
        if fp == 0: continue
        sp = int(small[p - 1])
        p2 = p * p
        # large part: k with x//k >= p^2  <=>  k <= x // p^2
        kmax = min(r, x // p2)
        k = ks[:kmax]
        kp = k * p
        inl = kp <= r
        nl = int(np.count_nonzero(inl))                                # kp increasing: first nl are "large" lookups
        vals = np.empty(kmax, dtype=np.int64)
        vals[:nl] = large[kp[:nl] - 1]
        vals[nl:] = small[x // kp[nl:]]
        large[:kmax] -= fp * (vals - sp)
        # small part: v from r down to p^2
        if p2 <= r:
            v = np.arange(p2, r + 1, dtype=np.int64)
            small[p2:] -= fp * (small[v // p] - sp)
    return int(large[0])
if __name__ == "__main__":
    for xs in sys.argv[1:]:
        x = int(xs)
        pix = int(subprocess.run(["primecount", str(x)], capture_output=True, text=True, check=True).stdout)
        s4, s3 = char_prime_sum(x, 4), char_prime_sum(x, 3)
        assert (pix - 1 + s4) % 2 == 0 and (pix - 1 + s3) % 2 == 0
        print(f"x = {x}: pi = {pix}; pi(x;4,1) = {(pix - 1 + s4) // 2}, pi(x;4,3) = {(pix - 1 - s4) // 2}, "
              f"pi(x;3,1) = {(pix - 1 + s3) // 2}, pi(x;3,2) = {(pix - 1 - s3) // 2}")
