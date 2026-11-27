# Repaired forms of Conjectures 1, 2 of arXiv:2502.21050, tested.
from motzkin_hankel import motzkin_power, H
import sys
def test(r, nn):
    N = 3*r*nn + 2*r + 4
    a = motzkin_power(r, 2*N + 2); cache = {}
    def h(n):
        if n not in cache: cache[n] = H(a, n)
        return cache[n]
    bad = []
    if r % 3 == 0:
        sg = -1 if (r*(r-1)//2) % 2 else 1
        for n in range(0, 3*nn):
            if h(r*n) != sg**n*(n+1)**(r-1) or h(r*n+1) != sg**n*(n+1)**(r-1): bad.append(('R1a', r, n))
            if n >= 1 and abs(h(r*n+2) + h(r*n-1)) != (r*(r-3)//2)*(n+1)**(r-1): bad.append(('R1b', r, n))
    else:
        for n in range(0, nn):
            v = [h(3*r*n), h(3*r*n+1), h(3*r*n+r), h(3*r*n+r+1)]
            if v[0] != v[1] or v[2] != v[3] or abs(v[0]) != 1 or abs(v[2]) != 1: bad.append(('R2a', r, n, v))
            if h(3*r*n+2*r) != 0 or h(3*r*n+2*r+1) != 0: bad.append(('R2b', r, n))
            x, y = h(3*r*n+2*r-1), h(3*r*n+2*r+2)
            if x != -y or abs(x) != (2*r*(n+1))**(r-2): bad.append(('R2c', r, n, x, y))
            if n >= 1 and 2*abs(h(3*r*n+2) + h(3*r*n-1)) != r*(r-3): bad.append(('R2d', r, n))
    return bad
for r in range(int(sys.argv[1]), int(sys.argv[2])+1):
    b = test(r, int(sys.argv[3])); print(r, 'repaired statements OK' if not b else b[:6], flush=True)
