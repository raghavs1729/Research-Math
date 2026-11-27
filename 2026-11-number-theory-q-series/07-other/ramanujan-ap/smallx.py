#!/usr/bin/env python3
"""Exact structure of the failure set F = {x > 1 : pi(x;q,a)^2 >= (e x/(phi log x)) pi(x/e;q,a)}
for small x, independently of the C code (plain sieve + Arb ball arithmetic).

On (1, e): x/log x decreases, on [e, oo) it increases.  For x < e*p1 (p1 the least prime of
the class) pi(x/e) = 0 and the inequality fails trivially: F contains (1, e*p1).
For x >= e the function G(x) = n^2 - (e x/(phi log x)) m is non-increasing between
consecutive primes of the class, jumps up at primes p (n -> n+1) and down at x = e*l
(m -> m+1), so F is a union of intervals each starting at a prime of the class.

Usage: smallx.py X [q a]   -> components of F in (1, X], counts of failing primes / starts
                             per bin floor(100 log10 p) for comparison with apscan.
"""
import sys, collections
from flint import arb, ctx
ctx.prec = 128
E = arb(1).exp()
PHI = {1: 1, 3: 2, 4: 2, 5: 4, 8: 4, 12: 4}
def primes_upto(N):
    s = bytearray([1]) * (N + 1); s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]
def sign_G(n, m, phi, x):          # sign of n^2 phi log x - e x m  (x arb or int)
    x = arb(x)
    v = arb(n) ** 2 * phi * x.log() - E * x * m
    if v > 0: return 1
    if v < 0: return -1
    raise ArithmeticError("undecided sign")
def root(n, m, phi, lo):            # x >= lo >= e with x/log x = phi n^2/(e m); float is enough for display
    import math
    t = phi * n * n / (math.e * m)
    a, b = max(lo, math.e), max(lo, math.e) * 2 + 10
    while b / math.log(b) < t: b *= 2
    for _ in range(200):
        mid = (a + b) / 2
        if mid / math.log(mid) < t: a = mid
        else: b = mid
    return a
def analyse(q, a, X, P):
    phi = PHI[q]
    cls = [p for p in P if (q == 1 or (p % q == a and q % p != 0))]
    # events: (position as arb, kind, prime)
    ev = [(arb(p), 0, p) for p in cls if p <= X] + [(E * l, 1, l) for l in cls if float((E * l).mid()) <= X]
    ev.sort(key=lambda t: float(t[0].mid()))
    # check ordering is certified
    for i in range(len(ev) - 1):
        assert ev[i][0] < ev[i + 1][0]
    n = m = 0
    comps = []          # (start, end, endkind)
    active = True; start = 1.0      # F contains (1, e*p1): G = n^2 >= 0 with m = 0
    nfail = collections.Counter(); nstart = collections.Counter(); failing = []
    import math
    for i, (pos, kind, pr) in enumerate(ev):
        # state on [previous event, pos): (n, m).  Did an active failure end before pos?
        if active and m > 0:
            # failing just before pos?  G(pos-) with current (n, m)
            if sign_G(n, m, phi, pos) < 0:
                comps.append((start, root(n, m, phi, start), "root")); active = False
        if kind == 0:
            n += 1
            f = (m == 0) or sign_G(n, m, phi, pr) >= 0
            if f:
                b = int(math.floor(100 * math.log10(pr) + 1e-12))
                nfail[b] += 1; failing.append(pr)
                if not active:
                    active = True; start = float(pr); nstart[b] += 1
        else:
            m += 1
            if active and sign_G(n, m, phi, pos) < 0:
                comps.append((start, float(pos.mid()), "e*%d" % pr)); active = False
    return comps, nfail, nstart, failing, active
if __name__ == "__main__":
    X = int(float(sys.argv[1]))
    P = primes_upto(X)
    sel = [(int(sys.argv[2]), int(sys.argv[3]))] if len(sys.argv) > 3 else [(1, 0), (3, 1), (3, 2), (4, 1), (4, 3)]
    for q, a in sel:
        comps, nfail, nstart, failing, active = analyse(q, a, X, P)
        print(f"class ({q},{a}), x <= {X}: {len(comps)} completed components, failing primes {len(failing)}, still failing at X: {active}")
        print("  first components:", "; ".join(f"[{s:.6g}, {e:.6g}{']' if k=='root' else ')'}" for s, e, k in comps[:8]))
        for b in sorted(set(nfail) | set(nstart)):
            print(f"SBIN {q} {a} {b} {nfail[b]} {nstart[b]}")
