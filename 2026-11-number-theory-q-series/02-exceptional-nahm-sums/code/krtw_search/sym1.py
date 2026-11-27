# Test the "symmetrized" (termwise in the multiset of squares) form of  K_j * C_1 = B_{j+1} * C_j
from fractions import Fraction
from itertools import permutations, combinations_with_replacement
import sys
q = Fraction(3, 7)
def poch(n, base=1):
    if n < 0: return None
    r = Fraction(1)
    for i in range(1, n+1): r *= (1 - q**(base*i))
    return r
def inv(n, base=1):
    p = poch(n, base)
    return Fraction(0) if p is None else 1/p
def wL(M, s):   # M has j+1 entries (last = t free), s has j entries
    j = len(s)
    w = Fraction(1)
    for i in range(j):
        ri = M[i]-s[i]
        rn = (M[i+1]-s[i+1]) if i+1 < j else 0
        w *= inv(ri-rn)
        w *= inv(s[i]-(s[i-1] if i > 0 else 0), 2)
        if w == 0: return w
    if M[j-1]-s[j-1] < 0: return Fraction(0)
    return w*inv(M[j], 2)
def wR(M, s):
    j = len(s)
    w = Fraction(1)
    for i in range(j):
        w *= inv(M[i]-M[i+1])
        w *= inv(s[i]-(s[i-1] if i > 0 else 0), 2)
        if w == 0: return w
    return w*inv(M[j], 2)
def test(j, maxpart):
    bad = 0; tot = 0
    for ms in combinations_with_replacement(range(maxpart+1), 2*j+1):
        L = Fraction(0); R = Fraction(0)
        for p in set(permutations(ms)):
            M, s = p[:j+1], p[j+1:]
            L += wL(M, s); R += wR(M, s)
        tot += 1
        if L != R:
            bad += 1
            if bad <= 5: print("  differs for multiset", ms, float(L), float(R))
    print("j=%d maxpart=%d: %d multisets, %d differ"%(j, maxpart, tot, bad))
test(1, 5)
test(2, 4)
