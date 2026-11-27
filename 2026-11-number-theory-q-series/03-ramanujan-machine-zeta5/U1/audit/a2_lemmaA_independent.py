"""AUDIT a2. Independent check of Lemma A (no sympy, no differential-ring reduction; different method from
verify/v1_operator_identity.py).

Method: truncated Taylor series arithmetic (own implementation, 100-digit mpmath numbers) around a random point x0.
  * u  := the solution of x(1-x)u'' + (1-3x/2)u' - u/16 = 0 with RANDOM u(x0), u'(x0)   (so: an arbitrary solution,
          not only 2F1(1/4,1/4;1;x));  built by the Taylor recursion of the ODE;
  * s  := sqrt(1-x)  (binomial series);
  * T  := any primitive of 1/(x s u^2), with RANDOM T(x0);
  * h  := exp(c T) + T^7 + 1/(T-d)   (c, d random; a test function that is not a polynomial of degree <= 4),
          h^(5) = c^5 exp(cT) + 2520 T^2 - 120/(T-d)^6;
  * L applied literally as theta^5 - 2x(theta+1/2)((theta+1/2)^4 + (5/8)(theta+1/2)^2 + 1/32) + x^2 (theta+1)^5,
    theta = x d/dx acting on Taylor series in e = x - x0.
We compare the Taylor coefficients of L[u^4 h(T)] and of h^(5)(T)/(s u^6).
Controls (non-vacuity): the same test must FAIL if 5/8 is replaced by 3/5, if 1/32 is replaced by 1/16, if the
exponent u^6 is replaced by u^5, or if sqrt(1-x) is dropped.
Also: exact check (Fractions) of the action of L on power series, and of Sym^4 at x = 0 via log-free solution.
"""
from mpmath import mp, mpf, mpc, sqrt, rand, exp as mexp, nstr
from fractions import Fraction as Fr
mp.dps = 100
M = 16     # series are kept modulo e^(M+1)

class S:
    def __init__(self, c): self.c = list(c) + [mpf(0)] * (M + 1 - len(c))
    @staticmethod
    def const(a): return S([a])
    def __add__(self, o):
        o = o if isinstance(o, S) else S.const(o); return S([x + y for x, y in zip(self.c, o.c)])
    __radd__ = __add__
    def __neg__(self): return S([-x for x in self.c])
    def __sub__(self, o): return self + (-(o if isinstance(o, S) else S.const(o)))
    def __rsub__(self, o): return (-self) + o
    def __mul__(self, o):
        if not isinstance(o, S): return S([x * o for x in self.c])
        r = [mpf(0)] * (M + 1)
        for i, x in enumerate(self.c):
            for j in range(M + 1 - i): r[i + j] += x * o.c[j]
        return S(r)
    __rmul__ = __mul__
    def inv(self):
        r = [1 / self.c[0]] + [mpf(0)] * M
        for n in range(1, M + 1):
            r[n] = -sum(self.c[k] * r[n - k] for k in range(1, n + 1)) / self.c[0]
        return S(r)
    def __truediv__(self, o): return self * (o.inv() if isinstance(o, S) else 1 / o)
    def __pow__(self, k):
        r = S.const(mpf(1))
        for _ in range(k): r = r * self
        return r
    def d(self):  return S([(k + 1) * self.c[k + 1] for k in range(M)] + [mpf(0)])     # loses the top order
    def integ(self, c0): return S([c0] + [self.c[k] / (k + 1) for k in range(M)])
    def exp(self):       # exp via f' = g' f
        r = [mexp(self.c[0])] + [mpf(0)] * M
        for n in range(1, M + 1):
            r[n] = sum(k * self.c[k] * r[n - k] for k in range(1, n + 1)) / n
        return S(r)
    def sqrt(self):
        r = [sqrt(self.c[0])] + [mpf(0)] * M
        for n in range(1, M + 1):
            r[n] = (self.c[n] - sum(r[k] * r[n - k] for k in range(1, n))) / (2 * r[0])
        return S(r)

def run(x0, al=mpf(5) / 8, be=mpf(1) / 32, upow=6, use_s=True):
    X = S([x0, 1])                                  # x = x0 + e
    # arbitrary solution u of the hypergeometric equation: coefficient recursion via fixed-point passes
    u0, u1 = rand() + mpf(1) / 2, rand() - mpf(1) / 2
    u = S([u0, u1])
    lead = (X * (1 - X)).inv()
    for _ in range(M):
        upp = (u * (mpf(1) / 16) - (1 - mpf(3) / 2 * X) * u.d()) * lead       # u'' from the ODE
        u = S([u0, u1] + [upp.c[k] / ((k + 1) * (k + 2)) for k in range(M - 1)])
    resid = X * (1 - X) * u.d().d() + (1 - mpf(3) / 2 * X) * u.d() - u * (mpf(1) / 16)
    assert max(abs(c) for c in resid.c[:M - 2]) < mpf(10) ** -50
    s = (1 - X).sqrt()
    T = (X * s * u * u).inv().integ(rand() * 3 - 1)
    c, dd = rand() + mpf(1) / 3, 7 + rand()
    h = (T * c).exp() + T ** 7 + (T - dd).inv()
    h5 = (T * c).exp() * c ** 5 + 2520 * T ** 2 - 120 * ((T - dd) ** 6).inv()
    f = u ** 4 * h
    th = lambda g: X * g.d()
    def shift_pow(g, a, k):                          # (theta + a)^k g
        for _ in range(k): g = th(g) + g * a
        return g
    half = mpf(1) / 2
    g1 = shift_pow(f, half, 1)
    Lf = shift_pow(f, 0, 5) - 2 * X * (shift_pow(g1, half, 4) + al * shift_pow(g1, half, 2) + be * g1) \
         + X * X * shift_pow(f, 1, 5)
    rhs = h5 / (u ** upow) / (s if use_s else 1)
    K = M - 6                                        # 5 derivatives lose 5 orders
    scale = max(abs(x) for x in rhs.c[:K])
    return max(abs(x - y) for x, y in zip(Lf.c[:K], rhs.c[:K])) / scale

ok = True
for x0 in (mpf('0.137'), mpf('0.5'), mpf('0.83'), mpf('0.971'), mpc('0.3', '0.4'), mpf('-1.7')):
    e = run(x0)
    print("x0 = %-14s  max rel. deviation of first %d Taylor coefficients of L[u^4 h(T)] - h^(5)(T)/(sqrt(1-x) u^6): %s"
          % (nstr(x0, 5), M - 6, nstr(e, 3)))
    ok &= e < mpf(10) ** -40
assert ok
print("controls (must be large):")
x0 = mpf('0.37')
for name, kw in (("5/8 -> 3/5", dict(al=mpf(3) / 5)), ("1/32 -> 1/16", dict(be=mpf(1) / 16)),
                 ("u^6 -> u^5", dict(upow=5)), ("drop sqrt(1-x)", dict(use_s=False))):
    e = run(x0, **kw); print("   %-16s deviation %s" % (name, nstr(e, 3))); assert e > mpf(10) ** -6

# exact: action of L on power series, and a_n/8 = 2T^5 + (5/4)T^3 + T/16
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
for n in range(-3, 40):
    Tn = Fr(2 * n + 1, 2)
    assert Fr(a(n), 8) == 2 * Tn**5 + Fr(5, 4) * Tn**3 + Fr(1, 16) * Tn
    assert Fr(a(n), 8) == 2 * Tn * (Tn**4 + Fr(5, 8) * Tn**2 + Fr(1, 32))
print("a_n/8 = 2(n+1/2)((n+1/2)^4 + (5/8)(n+1/2)^2 + 1/32): exact PASS")
print("a2 PASS")
