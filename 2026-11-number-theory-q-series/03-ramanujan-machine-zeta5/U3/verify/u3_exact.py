"""U3, exact checks (Fractions / sympy).  See ../THEOREM.md.
(a) local solutions at x = 1 (s = 1 - x):  nu = s + (3/4)s^2 + ...,  r = s^2 + ...,  r* = 1 - s^2/16 + ...;
    theta^i r*(1) = (1, 0, -1/8);  Lemma 3 right side for r*:  -1/8 - u(u+1/2).
(b) R*(u) := 1/u - (2Y(u) - Z(u))/(16 u^4) satisfies u^4 R*(u) - (u+1/2)^4 R*(u+1) = -1/8 - u(u+1/2), given
    Y(u) - k(u) Y(u+1) = 1 and Z(u) - k(u) Z(u+1) = 1/(u+1),  k(u) = ((u+1/2)/(u+1))^4.
(c) mu_rho := nu theta(rho) - rho theta(nu):  the boundary polynomial B5(u) of L5 vanishes identically for EVERY rho
    holomorphic at x = 1 (checked on the basis rho = e0, e1, e2 of holomorphic solutions).
(d) gauge map G(X)_n = P1(n) X_n + 32 (n+1)^5 X_{n+1} applied to beta_j(n) := int_0^1 s^j x^n dx = j! n!/(n+j+1)!:
    mu* = m0 + m1 s + m2 s^2 + m3 s^3 + O(s^4)  ==>  sum_{j<=3} m_j G(beta_j)_n = O(1)   (the n^2 and n^1 terms cancel);
    mu  = -s^2 - (3/2)s^3 + O(s^4)              ==>  sum_{j=2,3} m_j G(beta_j)_n = 128 n + O(1).
(e) telescoping: H(k) := a_k^2 k^3/(k - 1/2), a_k = ((1/2)_k/k!)^2:  H(k+1) - H(k) = -(1/8) a_k^2 (1/(k+1) + 1/(k-1/2)).
(f) U3 sequences are increasing: at3(n) > (n+1)^5 + n^5.
"""
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
    xs = [sh[m] - (sh[m - 1] if m >= 1 else 0) for m in range(len(c))]
    return [th[4][m] - xs[m] for m in range(len(c))]
def solve(c0, c1, c2):
    c = [Fr(c0), Fr(c1), Fr(c2)] + [Fr(0)] * (M - 3)
    for p in range(0, M - 5):
        w = c[:p + 6]; w[p + 3] = Fr(0); w[p + 4] = Fr(0); w[p + 5] = Fr(0)
        c[p + 3] = -L4(w)[p] / ((p + 1) * (p + 2)**2 * (p + 3))
    return c
def thvals(c, k):
    out = []; cur = c
    for _ in range(k):
        out.append(cur[0]); cur = theta(cur)
    return out
def mul(a, b):
    r = [Fr(0)] * M
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < M: r[i + j] += x * y
    return r
nu = solve(0, 1, Fr(3, 4)); r2 = solve(0, 0, 1); rs = solve(1, 0, Fr(-1, 16))
for c in (nu, r2, rs): assert all(v == 0 for v in L4(c)[:M - 6])
print("(a) r* =", " + ".join("%s s^%d" % (rs[i], i) for i in range(0, 6)), "+ ...")
ts = thvals(rs, 4); print("    theta^i r*(1), i=0..3:", ts)
assert ts[:3] == [1, 0, Fr(-1, 8)]
u = sp.symbols('u'); h = sp.Rational(1, 2)
# (b)
Y0, Y1, Z0s, Z1s = sp.symbols('Y0 Y1 Z0 Z1')      # Y(u), Y(u+1), Z(u), Z(u+1)
kk = ((u + h) / (u + 1))**4
Rst = lambda uu, Yv, Zv: 1 / uu - (2 * Yv - Zv) / (16 * uu**4)
expr = u**4 * Rst(u, Y0, Z0s) - (u + h)**4 * Rst(u + 1, Y1, Z1s)
expr = expr.subs({Y1: (Y0 - 1) / kk, Z1s: (Z0s - 1 / (u + 1)) / kk})
assert sp.simplify(expr - (-sp.Rational(1, 8) - u * (u + h))) == 0
assert ts[2] + ts[1] / 2 == Fr(-1, 8) and ts[0] == 1
print("(b) PASS: R*(u) = 1/u - (2Y - Z)/(16u^4) solves the Lemma-3 recurrence of r*  (right side -1/8 - u(u+1/2))")
# (c)
def mellin_boundary(polys_by_degree, nF):
    F = sp.symbols('F0:%d' % nF); th = sp.symbols('th'); bdry = 0; mult = []
    for d, P in enumerate(polys_by_degree):
        P = sp.Poly(sp.expand(P), th); c1 = u + d
        mult.append(sp.expand(P.as_expr().subs(th, -c1)))
        for (j,), coef in P.terms():
            for i in range(j): bdry += coef * (-c1)**(j - 1 - i) * F[i]
    return mult, sp.expand(bdry), F
th = sp.symbols('th')
Gp = 2 * (th + h) * ((th + h)**4 + (th + h)**2 + sp.Rational(1, 16))
mult5, B5, F5 = mellin_boundary([th**5, -Gp, (th + 1)**5], 5)
def mu_of(rho): return [p - q for p, q in zip(mul(nu, theta(rho)), mul(rho, theta(nu)))]
for name, rho in (("e0", solve(1, 0, 0)), ("e1", solve(0, 1, 0)), ("e2", solve(0, 0, 1)), ("r*", rs)):
    tm = thvals(mu_of(rho), 5)
    val = sp.expand(B5.subs({F5[i]: sp.Rational(tm[i].numerator, tm[i].denominator) for i in range(5)}))
    assert val == 0, (name, val)
print("(c) PASS: B5(u) == 0 for mu_rho, rho = e0, e1, e2 (hence for all rho holomorphic at 1), and for r*")
mus = mu_of(rs); mu = mu_of(r2)
print("    mu* =", " + ".join("%s s^%d" % (mus[i], i) for i in range(0, 5)), "+ ...;  theta^i mu*(1):", thvals(mus, 5))
print("    mu  =", " + ".join("%s s^%d" % (mu[i], i) for i in range(0, 5)), "+ ...")
# (d)
n = sp.symbols('n', positive=True)
P1 = lambda m: -(2 * m + 1) * (16 * m**4 + 56 * m**3 + 76 * m**2 + 46 * m + 11)
def beta(j, m):
    r = sp.Integer(sp.factorial(j))
    for i in range(1, j + 2): r = r / (m + i)
    return r
Gb = lambda j: sp.together(P1(n) * beta(j, n) + 32 * (n + 1)**5 * beta(j, n + 1))
tot_s = sp.together(sum(sp.Rational(mus[j].numerator, mus[j].denominator) * Gb(j) for j in range(4)))
tot = sp.together(sum(sp.Rational(mu[j].numerator, mu[j].denominator) * Gb(j) for j in range(4)))
num, den = sp.fraction(sp.cancel(tot_s))
dn, dd = sp.degree(num, n), sp.degree(den, n)
print("(d) sum_{j<=3} m*_j G(beta_j)_n = rational function with deg num - deg den = %d;  limit = %s" % (dn - dd, sp.limit(tot_s, n, sp.oo)))
assert dn < dd            # strict: O(1/n), needed in Prop. 4.1
assert sp.simplify(tot_s - (6712*n**4 + 23412*n**3 + 31402*n**2 + 18297*n + 4015) / (8*(n+1)*(n+2)*(n+3)*(n+4)*(n+5))) == 0
assert sp.limit(tot_s * n, n, sp.oo) == 839
print("    individual orders: G(beta_0) ~ %s n^2, G(beta_1) ~ %s n^2, G(beta_2) ~ %s n, G(beta_3) -> %s"
      % (sp.limit(Gb(0) / n**2, n, sp.oo), sp.limit(Gb(1) / n**2, n, sp.oo), sp.limit(Gb(2) / n, n, sp.oo), sp.limit(Gb(3), n, sp.oo)))
l1 = sp.limit(tot / n, n, sp.oo); print("    for mu: sum_{j<=3} m_j G(beta_j)_n / n ->", l1); assert l1 == 128
# uniqueness of c2: with rho = r* + lam r the sum has slope 128 lam
# (e)
k = sp.symbols('k')
lhs = (k + h)**3 / (k + 1) - k**3 / (k - h)        # [H(k+1) - H(k)]/a_k^2, using a_{k+1}^2 = a_k^2 ((k+1/2)/(k+1))^4
assert sp.simplify(lhs + sp.Rational(1, 8) * (1 / (k + 1) + 1 / (k - h))) == 0
print("(e) PASS: H(k+1) - H(k) = -(1/8) a_k^2 (1/(k+1) + 1/(k-1/2))")
# (f)
T = n + h
at3 = 2 * T * (T**4 + 5 * T**2 + sp.Rational(25, 16))
d = sp.expand(at3 - (n + 1)**5 - n**5); print("(f) at3(n) - (n+1)^5 - n^5 =", sp.factor(d))
assert all(c > 0 for c in sp.Poly(d, n).all_coeffs())
# (g) certificate for the gauge lemma: R_3(G(X))_n = P_1(n-1) R_2(X)_n + 32 (n+1)^5 R_2(X)_{n+1} for arbitrary X
Xm1, X0, X1, X2 = sp.symbols('Xm1 X0 X1 X2')
Th = lambda m: m + h
at2 = lambda m: 2 * Th(m) * (Th(m)**4 + Th(m)**2 + sp.Rational(1, 16))
at3f = lambda m: 2 * Th(m) * (Th(m)**4 + 5 * Th(m)**2 + sp.Rational(25, 16))
Gm = lambda m, a_, b_: P1(m) * a_ + 32 * (m + 1)**5 * b_
R3 = (n + 1)**5 * Gm(n + 1, X1, X2) + n**5 * Gm(n - 1, Xm1, X0) - at3f(n) * Gm(n, X0, X1)
R2n = (n + 1)**5 * X1 + n**5 * Xm1 - at2(n) * X0
R2n1 = (n + 2)**5 * X2 + (n + 1)**5 * X0 - at2(n + 1) * X1
assert sp.expand(R3 - (P1(n - 1) * R2n + 32 * (n + 1)**5 * R2n1)) == 0
print("(g) PASS: R_3(G X)_n = P_1(n-1) R_2(X)_n + 32 (n+1)^5 R_2(X)_{n+1}")
print("U3 exact checks PASS")
