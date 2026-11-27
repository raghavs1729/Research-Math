#!/usr/bin/env python3
"""Check apscan_star against an exact computation (Fractions + Arb) for x <= X.
Usage: validate_star.py star_output X"""
import sys, math, collections
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
import smallx
ctx.prec = 128
X = int(float(sys.argv[2])); E = arb(1).exp()
CQ = [1, 3, 3, 4, 4]; CA = [0, 1, 2, 1, 3]; PHI = [1, 2, 2, 2, 2]
P = smallx.primes_upto(X)
bmax = int(math.floor(100 * math.log10(X))) - 1
scan = collections.Counter()
for line in open(sys.argv[1]):
    w = line.split()
    if w and w[0] == 'B' and int(w[4]) <= bmax: scan[(int(w[3]), int(w[4]))] += int(w[6])
bad = 0
for c in range(5):
    q, a, phi = CQ[c], CA[c], PHI[c]
    ev = []                                   # (position, weight) of prime powers in the class
    for p in P:
        if q % p == 0 and q > 1: continue
        v, k = p, 1
        while v <= X:
            if q == 1 or v % q == a: ev.append((v, Fr(1, k)))
            v *= p; k += 1
    ev.sort()
    lagpos = [(E * v, w) for v, w in ev]
    N = Fr(0); M = Fr(0); j = 0
    nf = collections.Counter(); tot = 0
    for v, w in ev:
        N += w
        while j < len(lagpos) and lagpos[j][0] < v:
            M += lagpos[j][1]; j += 1
        if v < 7: continue
        if M == 0: f = True
        else:
            val = arb(fmpq(N.numerator, N.denominator)) ** 2 * phi * arb(v).log() - E * v * arb(fmpq(M.numerator, M.denominator))
            assert val > 0 or val < 0
            f = bool(val > 0)
        if f:
            nf[int(math.floor(100 * math.log10(v) + 1e-12))] += 1; tot += 1
    for b in range(bmax + 1):
        if nf[b] != scan[(c, b)]:
            bad += 1; print("MISMATCH", c, b, nf[b], scan[(c, b)])
    print(f"class ({q},{a}): failing events (primes and prime powers >= 7) in complete bins: {sum(nf[b] for b in range(bmax+1))}")
print("mismatches:", bad)
