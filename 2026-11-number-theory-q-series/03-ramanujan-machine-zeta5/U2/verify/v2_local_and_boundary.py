"""V2 (exact).  Local theory of L4 = theta^4 - x(theta+1/2)^4 at x = 1 and the boundary terms of Mellin transforms.

(a) In s = 1 - x (theta = -(1-s) d/ds) the coefficient of s^p in L4(sum c_m s^m) is
        (p+1)(p+2)^2(p+3) c_{p+3} + (terms in c_0..c_{p+2}),
    so c_0, c_1, c_2 are free and determine a formal (hence convergent) holomorphic solution: indicial exponents 0,1,1,2
    with a 3-dimensional space of holomorphic solutions.
(b) nu := s + (3/4)s^2 + ..., r2 := s^2 + ... (first coefficients printed), and theta^k values at x=1.
(c) Lemma M: for F holomorphic at x=1 with F_i = theta^i F(1), integration by parts gives
        int_0^1 x^(u-1) L4 F dx = u^4 Fhat(u) - (u+1/2)^4 Fhat(u+1) - (F_2 + F_1/2 - u(u+1/2) F_0),
    Fhat(u) = int_0^1 F x^(u-1) dx.  For nu the bracket is 0, for r2 it is 2.
(d) For mu := nu*theta(r2) - r2*theta(nu) and g(u) = int_0^1 mu x^u dx:
        int_0^1 x^(u-1) L5 mu dx = -u^5 g(u-1) + atil(u) g(u) - (u+1)^5 g(u+1) + B5(u),
    and B5(u) == 0 identically (so g satisfies the U2 recurrence in the continuous variable u).
All arithmetic exact (Fractions / sympy rationals)."""
from fractions import Fraction as Fr
from math import comb
import sympy as sp

M = 14
half = Fr(1, 2)
def theta(c):
    return [(m * c[m]) - ((m + 1) * c[m + 1] if m + 1 < len(c) else 0) for m in range(len(c))]
def L4(c):
    th = [c]
    for _ in range(4): th.append(theta(th[-1]))
    sh = [sum(comb(4, j) * half**(4 - j) * th[j][m] for j in range(5)) for m in range(len(c))]
    xs = [sh[m] - (sh[m - 1] if m >= 1 else 0) for m in range(len(c))]       # multiply by x = 1 - s
    return [th[4][m] - xs[m] for m in range(len(c))]
# (a) pivot
for p in range(0, 8):
    e = [Fr(0)] * (p + 6); e[p + 3] = Fr(1)
    d = L4(e)[p]
    assert d == (p + 1) * (p + 2)**2 * (p + 3), (p, d)
    # no dependence of power p on c_m for m > p+3
    for mm in range(p + 4, p + 6):
        e = [Fr(0)] * (p + 8); e[mm] = Fr(1)
        assert L4(e)[p] == 0
print("(a) PASS: pivot of c_{p+3} at power s^p is (p+1)(p+2)^2(p+3) != 0")
def solve(c0, c1, c2):
    c = [Fr(c0), Fr(c1), Fr(c2)] + [Fr(0)] * (M - 3)
    for p in range(0, M - 5):
        w = c[:p + 6]; w[p + 3] = Fr(0); w[p + 4] = Fr(0); w[p + 5] = Fr(0)
        c[p + 3] = -L4(w)[p] / ((p + 1) * (p + 2)**2 * (p + 3))
    return c
nu = solve(0, 1, Fr(3, 4)); r2 = solve(0, 0, 1)
res = L4(nu); assert all(v == 0 for v in res[:M - 6]); res = L4(r2); assert all(v == 0 for v in res[:M - 6])
print("(b) nu =", " + ".join("%s s^%d" % (nu[i], i) for i in range(1, 6)), "+ ...")
print("    r2 =", " + ".join("%s s^%d" % (r2[i], i) for i in range(2, 6)), "+ ...")
def thvals(c, k):
    out = []; cur = c
    for _ in range(k):
        out.append(cur[0]); cur = theta(cur)
    return out
print("    theta^i nu(1), i=0..3:", thvals(nu, 4), "   theta^i r2(1):", thvals(r2, 4))
# (c),(d): symbolic boundary polynomials
u = sp.symbols('u')
def mellin_boundary(polys_by_degree, nF):
    """operator sum_d x^d P_d(theta); returns (list of multipliers of Fhat-type integrals per d, boundary polynomial in F_i)
       using  int x^c theta^j F = sum_{i<j} (-(c+1))^(j-1-i) F_i + (-(c+1))^j int x^c F,   c = u - 1 + d."""
    F = sp.symbols('F0:%d' % nF)
    th = sp.symbols('th')
    bdry = 0; mult = []
    for d, P in enumerate(polys_by_degree):
        P = sp.Poly(sp.expand(P), th)
        c1 = u + d          # c + 1
        mult.append(sp.expand(P.as_expr().subs(th, -c1)))
        for (j,), coef in P.terms():
            for i in range(j):
                bdry += coef * (-c1)**(j - 1 - i) * F[i]
    return mult, sp.expand(bdry), F
th = sp.symbols('th'); h = sp.Rational(1, 2)
mult4, B4, F = mellin_boundary([th**4, -(th + h)**4], 4)
assert sp.expand(mult4[0] - u**4) == 0 and sp.expand(mult4[1] + (u + h)**4) == 0
assert sp.expand(B4 + (F[2] + F[1] / 2 - u * (u + h) * F[0])) == 0
print("(c) PASS: boundary term of L4 is -(F_2 + F_1/2 - u(u+1/2) F_0)")
tn = thvals(nu, 4); tr = thvals(r2, 4)
assert tn[2] + tn[1] / 2 == 0 and tn[0] == 0 and tr[2] + tr[1] / 2 == 2 and tr[0] == 0
print("    nu: F_2 + F_1/2 = 0, F_0 = 0;   r2: F_2 + F_1/2 = 2, F_0 = 0")
G = 2 * (th + h) * ((th + h)**4 + (th + h)**2 + sp.Rational(1, 16))
mult5, B5, F5 = mellin_boundary([th**5, -G, (th + 1)**5], 5)
atil = 2 * (u + h) * ((u + h)**4 + (u + h)**2 + sp.Rational(1, 16))
# with g(u) = int mu x^u: int x^(u-1) mu = g(u-1), int x^u mu = g(u), int x^(u+1) mu = g(u+1)
assert sp.expand(mult5[0] + u**5) == 0 and sp.expand(mult5[1] - atil) == 0 and sp.expand(mult5[2] + (u + 1)**5) == 0
# mu = nu*theta(r2) - r2*theta(nu) as a series in s
def mul(a, b):
    r = [Fr(0)] * M
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < M: r[i + j] += x * y
    return r
mu = [p - q for p, q in zip(mul(nu, theta(r2)), mul(r2, theta(nu)))]
print("(d) mu =", " + ".join("%s s^%d" % (mu[i], i) for i in range(0, 6)), "+ ...")
tm = thvals(mu, 5)
print("    theta^i mu(1), i=0..4:", tm)
B5val = sp.expand(B5.subs({F5[i]: sp.Rational(tm[i].numerator, tm[i].denominator) for i in range(5)}))
print("    B5(u) =", B5val)
assert B5val == 0
print("V2 PASS")
