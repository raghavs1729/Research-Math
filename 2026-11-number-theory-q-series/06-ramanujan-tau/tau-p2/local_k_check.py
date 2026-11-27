#!/usr/bin/env python3
"""
Independent check (Python, exact integers; no code shared with local_k.gp) of the local obstruction
for tau(p^2) = c, c in {463, 601}:   t^2 - c = p^11,   t + sqrt(c) = eps^k * g^11  (k mod 11).

* eps: the fundamental unit of the maximal order of Q(sqrt c), computed here from the continued
  fraction of sqrt(c) (least solution of x^2 - c y^2 = +-1), and, for c = 1 mod 4, from the least
  solution of x^2 - c y^2 = +-4 (brute force over y), whichever is smaller.
* For a prime q and k in 0..10:  A(q, k) = { t mod q : (t + sqrt c) eps^(-k) is an 11th power in
  (O_K/q)^* }.  Arithmetic in F_q[y]/(y^2 - c) with pairs (works for split and inert q alike;
  "11th power" is decided by brute-force enumeration of all 11th powers of units).
* At q = 23 only t = 2 is allowed (Ramanujan's congruence, see TAU2.md).
Output: for each k the primes q <= QMAX with A(q, k) empty.
"""
from math import isqrt
import sys

def pell_unit(c):
    # least (x, y), y > 0, with x^2 - c y^2 = +-1 via continued fraction of sqrt(c)
    a0 = isqrt(c); m, d, a = 0, 1, a0
    h0, h1, k0, k1 = 1, a0, 0, 1
    while h1 * h1 - c * k1 * k1 not in (1, -1):
        m = d * a - m; d = (c - m * m) // d; a = (a0 + m) // d
        h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0
    return h1, k1

def fundamental_unit(c):
    x, y = pell_unit(c)                       # in Z[sqrt c]
    if c % 4 == 1:
        # units (X + Y sqrt c)/2 with X, Y odd: X^2 - c Y^2 = +-4; search Y up to 2y
        for Y in range(1, 2 * y + 1, 2):
            for s in (4, -4):
                X2 = c * Y * Y + s
                X = isqrt(X2)
                if X * X == X2:
                    return ("half", X, Y)     # smaller unit exists
            if Y > 10**7: break               # bounded search; reported below
    return ("int", x, y)

def mul(a, b, c, q): return ((a[0]*b[0] + c*a[1]*b[1]) % q, (a[0]*b[1] + a[1]*b[0]) % q)
def power(a, n, c, q):
    r = (1, 0)
    while n:
        if n & 1: r = mul(r, a, c, q)
        a = mul(a, a, c, q); n >>= 1
    return r
def norm(a, c, q): return (a[0]*a[0] - c*a[1]*a[1]) % q
def inv(a, c, q):
    n = pow(norm(a, c, q), -1, q)
    return (a[0]*n % q, (-a[1])*n % q)

def run(c, QMAX):
    kind, X, Y = fundamental_unit(c)
    print(f"c = {c}: fundamental unit type {kind}: ({X} + {Y} sqrt({c})){'/2' if kind=='half' else ''}, norm {(X*X - c*Y*Y)//(4 if kind=='half' else 1)}")
    killers = {k: [] for k in range(11)}
    q = 2
    while q < QMAX:
        q += 1
        if any(q % r == 0 for r in range(2, isqrt(q) + 1)) or (q*q - 1) % 11 or q == c or q == 2: continue
        s = isqrt(q**11 + c); assert s*s != q**11 + c
        half = pow(2, -1, q)
        e = ((X * half) % q, (Y * half) % q) if kind == "half" else (X % q, Y % q)
        units = [(x, y) for x in range(q) for y in range(q) if norm((x, y), c, q) != 0]
        P11 = {power(g, 11, c, q) for g in units}
        if len(P11) == len(units): continue        # every unit is an 11th power
        T = [2] if q == 23 else range(q)
        einv = inv(e, c, q)
        for k in range(11):
            ek = power(einv, k, c, q)
            if not any(norm((t, 1), c, q) != 0 and mul((t % q, 1), ek, c, q) in P11 for t in T):
                killers[k].append(q)
    for k in range(11): print(f"  k = {k:2d}: killed by q in {killers[k]}")
    print("  surviving k:", [k for k in range(11) if not killers[k]])

if __name__ == "__main__":
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 240
    run(463, QMAX); run(601, QMAX)
