# arXiv:2502.21050 (Y. Wang, Y. Zhang): Hankel determinants H_n(M(x)^r), M = Motzkin g.f.  Conjectures 1, 2.
# H_n(A) = det [a_{i+j}]_{0<=i,j<n}, H_0 = 1.   Authors verified r <= 27.
import flint, sys
from fractions import Fraction
def motzkin_power(r, D):
    m = [1, 1]
    for n in range(2, D): m.append(((2*n+1)*m[n-1] + (3*n-3)*m[n-2]) // (n+2))
    M = flint.fmpz_poly(m); P = flint.fmpz_poly([1])
    for _ in range(r): P = flint.fmpz_poly((P*M).coeffs()[:D])
    c = [int(x) for x in P.coeffs()]
    return c + [0]*(D-len(c))
def hankels(r, nmax):
    a = motzkin_power(r, 2*nmax + 2)
    H = {0: 1}
    return a, H
def H(a, n):
    if n == 0: return 1
    if n < 0: return None
    return int(flint.fmpz_mat([[a[i+j] for j in range(n)] for i in range(n)]).det())
def test(r, nn):
    N = 3*r*nn + 2*r + 4
    a = motzkin_power(r, 2*N + 2)
    cache = {}
    def h(n):
        if n not in cache: cache[n] = H(a, n)
        return cache[n]
    bad = []
    if r % 3 == 0:
        for n in range(0, 3*nn):
            if h(r*n) != (-1)**n*(n+1)**(r-1) or h(r*n+1) != (-1)**n*(n+1)**(r-1): bad.append(('C1a', r, n, h(r*n), h(r*n+1)))
            if n >= 1:
                s = h(r*n+2) + h(r*n-1)
                if abs(s) != (r*(r-3)//2)*(n+1)**(r-1): bad.append(('C1b', r, n, s, Fraction(s, (n+1)**(r-1))))
    else:
        for n in range(0, nn):
            v = [h(3*r*n), h(3*r*n+1), h(3*r*n+r), h(3*r*n+r+1)]
            if len(set(v)) != 1 or abs(v[0]) != 1: bad.append(('C2a', r, n, v))
            if h(3*r*n+2*r) != 0 or h(3*r*n+2*r+1) != 0: bad.append(('C2b', r, n, h(3*r*n+2*r), h(3*r*n+2*r+1)))
            s = h(3*r*n+2*r-1) + h(3*r*n+2*r+2)
            if abs(s) != (2*r*(n+1))**(r-2): bad.append(('C2c', r, n, s, Fraction(s, (2*r*(n+1))**(r-2))))
            if n >= 1:
                s = h(3*r*n+2) + h(3*r*n-1)
                if 2*abs(s) != r*(r-3): bad.append(('C2d', r, n, s))
    return bad
if __name__ == '__main__':
    for r in range(int(sys.argv[1]), int(sys.argv[2])+1):
        b = test(r, int(sys.argv[3]))
        print(r, 'OK' if not b else b[:6], flush=True)
