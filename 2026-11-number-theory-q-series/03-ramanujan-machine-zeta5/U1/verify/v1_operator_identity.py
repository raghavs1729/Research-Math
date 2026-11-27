"""V1. Differential-algebra identity behind Theorem U1 (exact, sympy).

Let u be ANY solution of the hypergeometric equation of 2F1(1/4,1/4;1;x),
      x(1-x) u'' + (1 - 3x/2) u' - u/16 = 0,
let T be a function with T' = 1/(x sqrt(1-x) u^2), and let h be an arbitrary (5 times differentiable) function.
With  theta = x d/dx,  al = 5/8, be = 1/32  and
      L = theta^5 - 2x (theta+1/2)((theta+1/2)^4 + al (theta+1/2)^2 + be) + x^2 (theta+1)^5
we verify
      L[ u^4 * h(T) ] = h^(5)(T) / ( sqrt(1-x) * u^6 ).
In particular L[u^4 T^j] = 0 for j = 0..4 (L is the symmetric 4th power of the hypergeometric operator)."""
import sympy as sp
x, s, u, v = sp.symbols('x s u v')          # s = sqrt(1-x), v = u'
h = sp.symbols('h0:6')                       # h_k = h^(k)(T)
al, be = sp.Rational(5, 8), sp.Rational(1, 32)
Tp = 1 / (x * s * u**2)
def D(e):
    r = sp.diff(e, x) + sp.diff(e, s) * (-1 / (2 * s)) + sp.diff(e, u) * v \
        + sp.diff(e, v) * (-((1 - sp.Rational(3, 2) * x) * v - u / 16) / (x * (1 - x)))
    for k in range(5):
        r += sp.diff(e, h[k]) * h[k + 1] * Tp
    return r
def red(e):
    e = sp.together(e)
    num, den = sp.fraction(e)
    num = sp.rem(sp.expand(num), s**2 - (1 - x), s)
    den = sp.rem(sp.expand(den), s**2 - (1 - x), s)
    return sp.cancel(num / den)
theta = lambda e: red(x * D(e))
f = u**4 * h[0]
th = [f]
for k in range(5):
    th.append(theta(th[-1]))
half = sp.Rational(1, 2)
# (theta+1/2)^k f  and (theta+1)^5 f via binomial expansion in theta
def shifted(c, k):
    return sum(sp.binomial(k, j) * c**(k - j) * th[j] for j in range(k + 1))
Lf = th[5] - 2 * x * (shifted(half, 5) + al * shifted(half, 3) + be * shifted(half, 1)) + x**2 * shifted(1, 5)
diff = red(Lf - h[5] / (s * u**6))
num = sp.fraction(sp.together(diff))[0]
num = sp.rem(sp.expand(num), s**2 - (1 - x), s)
print("L[u^4 h(T)] - h^(5)(T)/(sqrt(1-x) u^6) =", sp.simplify(num))
assert sp.simplify(num) == 0
print("V1 PASS")
