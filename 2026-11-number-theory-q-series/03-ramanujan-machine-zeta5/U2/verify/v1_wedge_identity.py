"""V1 (exact, sympy).  Let L4 = theta^4 - x (theta+1/2)^4, theta = x d/dx.  For ANY two solutions f, h of L4 y = 0,
       m := f*theta(h) - h*theta(f)
satisfies L5 m = 0, where
       L5 = theta^5 - 2x (theta+1/2)((theta+1/2)^4 + (theta+1/2)^2 + 1/16) + x^2 (theta+1)^5 .
Computation in the differential ring Q(x)[f0..f3, h0..h3], fk = theta^k f, with theta^4 f = x (theta+1/2)^4 f solved
for f4."""
import sympy as sp
x = sp.symbols('x')
f = sp.symbols('f0:4'); h = sp.symbols('h0:4')
half = sp.Rational(1, 2)
def top(y):      # theta^4 y from the equation: (1-x) y4 = x (2 y3 + 3/2 y2 + 1/2 y1 + y0/16)
    return x * (4 * half * y[3] + 6 * half**2 * y[2] + 4 * half**3 * y[1] + half**4 * y[0]) / (1 - x)
def theta(e):
    r = x * sp.diff(e, x)
    for y in (f, h):
        for k in range(3):
            r += sp.diff(e, y[k]) * y[k + 1]
        r += sp.diff(e, y[3]) * top(y)
    return sp.together(r)
m = f[0] * h[1] - h[0] * f[1]
th = [m]
for _ in range(5):
    th.append(sp.cancel(theta(th[-1])))
def shifted(c, k):
    return sum(sp.binomial(k, j) * c**(k - j) * th[j] for j in range(k + 1))
L5m = th[5] - 2 * x * (shifted(half, 5) + shifted(half, 3) + sp.Rational(1, 16) * shifted(half, 1)) + x**2 * shifted(1, 5)
res = sp.cancel(sp.together(L5m))
print("L5[f theta h - h theta f] =", res)
assert res == 0
# a_n/8 for U2 equals G(n) := 2T(T^4 + T^2 + 1/16), T = n + 1/2
n = sp.symbols('n'); T = n + half
a = 8 * (n**5 + (n + 1)**5) - 12 * (n**3 + (n + 1)**3) + 7 * (2 * n + 1)
assert sp.expand(a / 8 - 2 * T * (T**4 + T**2 + sp.Rational(1, 16))) == 0
print("a_n/8 = 2T(T^4+T^2+1/16), T = n+1/2: OK")
print("V1 PASS")
