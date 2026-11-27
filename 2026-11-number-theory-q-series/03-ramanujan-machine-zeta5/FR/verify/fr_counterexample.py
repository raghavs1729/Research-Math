"""Exact check of the counterexample to Conjecture 1.3 of arXiv:2111.04468 as printed (degenerate range):
d = 2, B = -3/4, m = 1, x = (-2, 1, -2, 2)  =>  p_n = q_n = 0 for all n >= 1 (also for the integer inflation B=-12, m=4)."""
from fractions import Fraction as Fr
def run(B, m, xs, d, N=200):
    X1, X2 = xs[:d], xs[d:]
    def prod(X, n):
        r = Fr(1)
        for x in X: r *= (n - x)
        return r
    a = lambda n: (B / m) * prod(X1, n + 1) - m * prod(X2, n)
    b = lambda n: B * prod(xs, n)
    p = [Fr(1), a(0)]; q = [Fr(0), Fr(1)]
    for n in range(1, N + 1):
        p.append(a(n) * p[-1] + b(n) * p[-2]); q.append(a(n) * q[-1] + b(n) * q[-2])
    return p[1:], q[1:]
for (B, m) in ((Fr(-3, 4), Fr(1)), (Fr(-12), Fr(4))):
    p, q = run(B, m, [Fr(-2), Fr(1), Fr(-2), Fr(2)], 2)
    assert p[0] != 0 and q[0] == 1 and all(x == 0 for x in p[1:]) and all(x == 0 for x in q[1:])
    print("B = %s, m = %s: p_0 = %s, q_0 = 1, and p_n = q_n = 0 for 1 <= n <= 200" % (B, m, p[0]))
print("PASS")
