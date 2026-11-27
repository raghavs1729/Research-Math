# spectral coefficients Gamma_{n,m}(R;x) (times (x^2q^2;q^2)_inf), exact rational functions
import sympy as sp
from functools import lru_cache
x, q = sp.symbols('x q')
def poch(a, n, b=None):
    b = q if b is None else b
    r = sp.Integer(1)
    for i in range(n): r *= (1 - a*b**i)
    return r
def qbin(n, k):
    if k < 0 or k > n: return sp.Integer(0)
    return sp.cancel(poch(q, n)/(poch(q, k)*poch(q, n-k)))
def g(d, X):
    # conjectural g_d(x)*(x^2q^2;q^2)_inf ; d=0 -> 1
    if d == 0: return sp.Integer(1)
    return (-1)**d*X**d*q**(2*d*d-d)*poch(X**2, d, q**2)*(1+q**d)*(1+X*q**d)*(1-X*q**(2*d))/(poch(q**2, d, q**2)*(1-X**2))
@lru_cache(None)
def Gam(n, m, R):
    """Gamma~_{n,m}(R;x) as expression in x,q  (n>=m>=0)"""
    d = n-m
    if R < d: return sp.Integer(0)
    if m == 0:
        return sp.factor(g(d, x)*poch(q, R)/(poch(q, R-d)*poch(x*q, R+d)))
    lam = x**(n+m)*q**(n*n+m*m)
    rhs = 0
    for r in range(R):
        rhs += qbin(R, r)*x**r*q**(r*r)*Gam(n, m, r)
    for b in range(1, m+1):
        for r in range(R+1):
            inner = Gam(n-b, m-b, r).subs(x, x*q**(2*b))*poch(x**2*q**2, 2*b, q**2)
            rhs += qbin(R, r)*x**(r+2*b)*q**((r+b)**2+b*b)*inner/poch(q**2, b, q**2)
    return sp.factor(sp.cancel(rhs/(lam - x**R*q**(R*R))))
if __name__ == "__main__":
    import sys
    for (n, m) in [(1,1),(2,1),(2,2),(3,1)]:
        for R in range(n-m, n-m+4):
            print("Gamma~_{%d,%d}(R=%d) ="%(n,m,R), Gam(n,m,R))
        sys.stdout.flush()
