"""AUDIT-4, part 1 (own code). Exact checks of
 (i)  source Prop. antidiag-criterion: f_B(1,0)-f_B(1,1) = G(n,a)-D(n), by brute-force path counting (DP on the grid);
 (ii) the two printed forms of G, D (statement form and bracket form) agree;
 (iii) R10 formula (1);
 (iv) source Lemma lem:monotone quotient;
 (v) the first form of the conjecture vs rho(n)<1 by brute force over ALL points A, for 5<=n<=NB.
"""
import sys
from fractions import Fraction as F
from math import comb

def paths_dp(n, B):
    """number of monotone paths (0,0)->(n,n) avoiding B through each point A: fwd[A]*bwd[A]."""
    fwd = [[0]*(n+1) for _ in range(n+1)]
    for x in range(n+1):
        for y in range(n+1):
            if (x, y) == B: continue
            if x == 0 and y == 0: fwd[x][y] = 1; continue
            fwd[x][y] = (fwd[x-1][y] if x else 0) + (fwd[x][y-1] if y else 0)
    bwd = [[0]*(n+1) for _ in range(n+1)]
    for x in range(n, -1, -1):
        for y in range(n, -1, -1):
            if (x, y) == B: continue
            if x == n and y == n: bwd[x][y] = 1; continue
            bwd[x][y] = (bwd[x+1][y] if x < n else 0) + (bwd[x][y+1] if y < n else 0)
    return [[fwd[x][y]*bwd[x][y] for y in range(n+1)] for x in range(n+1)]

def G_stmt(n, a): return F(n-2*a+1, n-a)*comb(n, a)*comb(n-2, a-1)
def G_brk(n, a):  return (2 - F(n-1, n-a))*comb(n, a)*comb(n-2, a-1)
def D_stmt(n):    return F(comb(2*n-2, n-1), n)
def D_brk(n):     return 2*comb(2*n-2, n-1) - comb(2*n-1, n-1)
def R(n, a):      return G_stmt(n, a)/D_stmt(n)
def R1(n, a):
    u = n-2*a
    return F((u+1)*a*2*(2*n-1)*comb(n, a)**2, (n-1)*n*comb(2*n, n))

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 60
# (ii),(iii),(iv) for 2<=n<=400, ALL 1<=a<=n-1 (formula (1) only needs 2a<n but holds for all a)
for n in range(2, 401):
    assert D_stmt(n) == D_brk(n)
    for a in range(1, n):
        assert G_stmt(n, a) == G_brk(n, a)
        assert R(n, a) == R1(n, a), (n, a)
    for a in range(2, (n+1)//2):
        assert 2*a < n
        assert R(n, a-1)/R(n, a) == F(a*(a-1)*(n-2*a+3), (n-a+1)**2*(n-2*a+1))
print("(ii),(iii),(iv): G,D forms agree; formula (1) == G/D; lem:monotone quotient -- exact, 2<=n<=400: OK")

# (i),(v)
mism = []
for n in range(2, NB+1):
    C = {(1,0),(0,1),(1,1),(n-1,n),(n,n-1),(n-1,n-1)}
    rho = max([R(n, a) for a in range(1, (n+1)//2) if 2*a < n], default=None)
    allcorner = True; onC = True
    for a in range(0, n+1):
        B = (a, n-a)
        f = paths_dp(n, B)
        if 1 <= a <= n-1:
            assert f[1][0]-f[1][1] == G_stmt(n, a)-D_stmt(n), (n, a)
            assert f[1][1] == f[n-1][n-1] and f[1][0] == f[n][n-1] and f[0][1] == f[n-1][n]
        pts = [(x, y) for x in range(n+1) for y in range(n+1) if (x, y) not in ((0,0),(n,n),B)]
        mx = max(f[x][y] for x, y in pts)
        arg = {(x, y) for x, y in pts if f[x][y] == mx}
        if not arg <= C: onC = False
        # "maximum is attained at (1,1) and (n-1,n-1)"
        if not ({(1,1),(n-1,n-1)} <= arg): allcorner = False
    if n >= 5:
        # source: for n>=9, [exists antidiagonal B whose max is at boundary points] <=> rho>1
        if n >= 9:
            assert onC, n
            assert allcorner == (rho <= 1), (n, float(rho))
        mism.append((n, onC, allcorner, rho > 1, rho == 1))
print("(i) f_B(1,0)-f_B(1,1) = G-D by brute-force DP, 2<=n<=%d, all 1<=a<=n-1: OK" % NB)
print("(v) for 9<=n<=%d: max over ALL A (all antidiagonal B incl. (0,n),(n,0), a=n/2, 2a>n) lies in C, and" % NB)
print("    [max attained at (1,1)&(n-1,n-1) for every antidiagonal B] <=> rho(n)<=1 : OK")
print("    n: (max on C?, conj-form-1 holds?, rho>1?, rho==1?) for small n:", [m for m in mism if m[0] <= 12])
print("    n<=%d with conj form 1 true:" % NB, [m[0] for m in mism if m[2]])
