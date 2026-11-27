#!/usr/bin/env python3
"""explicit17.py q logC [kmax] -- evaluate both sides of AF inequality (3), a = 1, for the first
kmax primorials in S_{q,1}:   LHS_k = N/(phi(N) (log(phi(q) log N))^{1/phi(q)}),  RHS = 1/C(q,1).
(3) holds iff LHS > RHS.  logC is passed in (from mertens_C.py) with its error bar."""
import sys
import mpmath as mp
def isprime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))

mp.mp.dps = 40
q = int(sys.argv[1]); logC = mp.mpf(sys.argv[2]); kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 8
phi = q - 1
RHS = mp.exp(-logC)
ps = []; n = 1
while len(ps) < kmax:
    n += q
    if isprime(n):
        ps.append(n)
N = 1; ratio = mp.mpf(1)
for k, p in enumerate(ps, 1):
    N *= p; ratio *= mp.mpf(p) / (p - 1)
    LHS = ratio / mp.log(phi * mp.log(N)) ** (mp.mpf(1) / phi)
    print("k=%d p=%d N=%d LHS=%s RHS=1/C=%s  log(LHS/RHS)=%s  %s" % (
        k, p, N, mp.nstr(LHS, 12), mp.nstr(RHS, 12), mp.nstr(mp.log(LHS / RHS), 6),
        "(3) holds" if LHS > RHS else "(3) FAILS"))
