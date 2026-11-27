"""AUDIT-4 part 6 (own code): brute-force test of the FIRST (verbal) form of conj:496, independent of every theorem of the source.
For given n and EVERY obstruction B=(a,n-a), 0<=a<=n, compute f_B(A) for ALL grid points A by dynamic programming over
B-avoiding paths (no binomial formulas) and determine the set of maximisers.
fwd[x][y] = # B-avoiding monotone paths (0,0)->(x,y);  # B-avoiding paths (x,y)->(n,n) = fwd[n-y][n-x]
(reverse the path and reflect in the antidiagonal, which fixes B);  f_B(x,y) = fwd[x][y]*fwd[n-y][n-x].
usage: a6_bruteforce_form1.py n1 n2 ...
"""
import sys
from itertools import accumulate
from multiprocessing import Pool
from math import comb
from fractions import Fraction as F

def rho(n):
    return max(F(n-2*a+1, n-a)*comb(n, a)*comb(n-2, a-1)/F(comb(2*n-2, n-1), n) for a in range(1, (n-1)//2+1))

def run(n):
    corner = {(1, 1), (n-1, n-1)}
    bad = []          # obstructions for which the maximiser set does not contain both (1,1),(n-1,n-1)
    notexact = []     # obstructions for which the maximiser set is not exactly {(1,1),(n-1,n-1)}
    for a in range(0, n+1):
        bx, by = a, n-a
        rows = []
        prev = None
        for x in range(n+1):
            if x == 0:
                row = [1]*(n+1)
                if bx == 0:
                    for y in range(by, n+1): row[y] = 0
            else:
                if x == bx:
                    left = list(accumulate(prev[:by])); right = list(accumulate(prev[by+1:]))
                    row = left+[0]+right
                else:
                    row = list(accumulate(prev))
            rows.append(row); prev = row
        assert rows[bx][by] == 0
        if a == 1: assert rows[n][n] == comb(2*n, n)-comb(n, 1)**2      # sanity: total number of B-avoiding paths
        mx = -1; arg = []
        for x in range(n+1):
            rx = rows[x]
            for y in range(n+1):
                if (x == 0 and y == 0) or (x == n and y == n): continue
                v = rx[y]*rows[n-y][n-x]
                if v > mx: mx = v; arg = [(x, y)]
                elif v == mx: arg.append((x, y))
        s = set(arg)
        if not corner <= s: bad.append((a, sorted(s)))
        if s != corner: notexact.append((a, sorted(s)))
    return n, bad, notexact

if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]]
    with Pool(min(4, len(ns))) as p:
        for n, bad, notexact in p.imap(run, ns):
            r = rho(n)
            print("n=%d  rho(n)=%.10f  (rho<1: %s)" % (n, float(r), r < 1))
            print("   conjecture form 1 [for every B on x+y=n, 0<=a<=n, max of f_B attained at (1,1) and (n-1,n-1)]: %s" % (not bad))
            print("   maximiser set is exactly {(1,1),(n-1,n-1)} for every B: %s" % (not notexact))
            if bad: print("   obstructions violating form 1 (a, maximisers):", bad[:6], "... total", len(bad))
            elif notexact: print("   ties:", notexact[:6])
            sys.stdout.flush()
