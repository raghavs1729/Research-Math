import sympy as sp, sys
from spec1 import Gam, poch, x, q
def coeffs(n, m):
    d = n-m
    ds = list(range(d, d+2*m+1))
    cs = {}
    # triangular solve: R = d' successively
    for R in ds:
        val = Gam(n, m, R)
        for dp in ds:
            if dp < R:
                val -= cs[dp]*poch(q, R)/(poch(q, R-dp)*poch(x*q, R+dp))
        cs[R] = sp.factor(sp.cancel(val*poch(x*q, 2*R)/poch(q, R)))
    return cs
for (n, m) in [(1, 1), (2, 1), (3, 1), (2, 2)]:
    cs = coeffs(n, m)
    print("(n,m)=(%d,%d)" % (n, m))
    for dp, c in cs.items():
        print("   c_%d =" % dp, c)
    # check sum rule at one more R
    R = n+m+1
    tot = sum(c*poch(q, R)/(poch(q, R-dp)*poch(x*q, R+dp)) for dp, c in cs.items())
    print("   check R=%d:" % R, sp.simplify(tot - Gam(n, m, R)) == 0)
    sys.stdout.flush()
