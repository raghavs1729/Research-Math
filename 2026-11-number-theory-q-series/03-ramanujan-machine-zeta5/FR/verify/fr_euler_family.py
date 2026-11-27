"""Exact tests for Theorem A (Conjecture 1.3 of arXiv:2111.04468: the Euler family has factorial reduction).

For random rational B, m, x_1..x_2d we inflate to an integer PCF, compute p_n, q_n by the recurrence and check
  (1) q_n = S_n := sum_{k=0}^n Pi2(k) prod_{i=k+2}^{n+1} H1(i),   p_n = Pi1(n+1) + H2(0) S_n          (exactly)
  (2) GCD_n := gcd(p_n, q_n) = |Pi1(n+1)| / D_n,  D_n = reduced denominator of sigma_n = sum_{k<=n} Pi2(k)/Pi1(k+1)
  (3) D_n divides the explicit bound Delta_n of Lemma 1 (so D_n <= C^n), and (1/n) log(GCD_n / n!^d) stays bounded.
All arithmetic exact (Fractions)."""
from fractions import Fraction as Fr
from math import gcd, factorial, log, lcm
import random, sys

def poly_prod(xs, n):
    r = Fr(1)
    for x in xs:
        r *= (n - x)
    return r

def run(d, B, m, xs, N=60, verbose=False):
    X1, X2 = xs[:d], xs[d:]
    # inflate: c with c*H1, c*H2 integer valued on integers: H1(n) = (B/m) prod (n-x_i), H2(n) = -m prod (n-x_i)
    den1 = 1
    for x in X1: den1 *= x.denominator
    den2 = 1
    for x in X2: den2 *= x.denominator
    c = lcm((B / m).denominator * den1, m.denominator * den2)
    H1 = lambda n: c * (B / m) * poly_prod(X1, n)
    H2 = lambda n: -c * m * poly_prod(X2, n)
    a = lambda n: H1(n + 1) + H2(n)
    b = lambda n: -H1(n) * H2(n)
    for n in range(0, 5):
        assert a(n).denominator == 1 and b(n).denominator == 1
        assert b(n) == c * c * B * poly_prod(xs, n) and a(n) == c * ((B / m) * poly_prod(X1, n + 1) - m * poly_prod(X2, n))
    p = [Fr(1), a(0)]; q = [Fr(0), Fr(1)]
    Pi1 = [Fr(1)]; Pi2 = [Fr(1)]
    for i in range(1, N + 3):
        Pi1.append(Pi1[-1] * H1(i)); Pi2.append(Pi2[-1] * H2(i))
    out = []
    sigma = Fr(0)
    for n in range(0, N + 1):
        if n >= 1:
            p.append(a(n) * p[-1] + b(n) * p[-2]); q.append(a(n) * q[-1] + b(n) * q[-2])
        pn, qn = p[n + 1], q[n + 1]
        S = sum(Pi2[k] * Pi1[n + 1] / Pi1[k + 1] for k in range(n + 1))
        assert qn == S and pn == Pi1[n + 1] + H2(0) * S, n
        sigma += Pi2[n] / Pi1[n + 1]
        assert S == Pi1[n + 1] * sigma
        G = gcd(int(pn), int(qn))
        assert pn.denominator == 1 and qn.denominator == 1
        assert Fr(G) == abs(Pi1[n + 1]) / sigma.denominator, (n, G)
        out.append(G)
    # growth
    vals = [(log(out[n]) - d * sum(log(i) for i in range(1, n + 1))) / n for n in (N // 2, N)]
    return vals

random.seed(1)
def rr(lo=-5, hi=5, dens=(1, 1, 2, 3, 4)):
    return Fr(random.randint(lo, hi), random.choice(dens))
cnt = 0
worst = 0.0
while cnt < 60:
    d = random.choice([1, 1, 2, 2, 3])
    B = rr(-6, 6, (1, 1, 2, 3)); m = rr(-4, 4, (1, 1, 2, 3))
    if B == 0 or m == 0: continue
    xs = [rr() for _ in range(2 * d)]
    if any(x.denominator == 1 and x >= 1 for x in xs[:d]):      # H1 must not vanish at positive integers
        continue
    v = run(d, B, m, xs, N=50 if d < 3 else 36)
    worst = min(worst, min(v))
    cnt += 1
print("Theorem A checks: %d random Euler-family PCFs, identities (1),(2) exact for all n; min (1/n) log(GCD_n/n!^d) = %.3f" % (cnt, worst))
# the paper's own example: PCF[2n^2+2n+1, -n^4] (B=-1, m=-1, x=0)
print("example b=-n^4, a=n^2+(n+1)^2:", run(2, Fr(-1), Fr(-1), [Fr(0)] * 4, N=60))
print("PASS")
