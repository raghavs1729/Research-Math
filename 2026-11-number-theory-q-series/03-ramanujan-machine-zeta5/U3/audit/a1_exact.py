"""Audit of U3 (independent exact checks).  Written by the referee; does not import the verify/ scripts.
 (A) a_n/8 = at3(n); at3 - (n+1)^5 - n^5; p_n, q_n of the CF vs 8^{n+1}(n+1)!^5 (w3, v3)_{n+1} for n < 40.
 (B) gauge identity (own elimination), G(w) = w3, G(v) = v3 + 4 w3 for n < 200.
 (C) local solutions at x = 1 via an explicit coefficient recurrence (derived here symbolically), r*, mu*,
     theta^i values, Lemma-3 boundary term for r*, B5 for mu (U2) and mu* by my own integration by parts.
 (D) Prop 3.1 functional relation for R*(u).  (E) Lemma 3.3 rational function + 128 lambda n.
 (F) telescoping H(k).  (G) Section 7 coefficient extraction redone, Legendre solve, the three forms.
"""
import sympy as sp
from fractions import Fraction as Fr
from math import factorial
n, u, k, X0, X1 = sp.symbols('n u k X0 X1')
h = sp.Rational(1, 2)
at2 = lambda m: 2*(m+h)*((m+h)**4 + (m+h)**2 + sp.Rational(1, 16))
at3 = lambda m: 2*(m+h)*((m+h)**4 + 5*(m+h)**2 + sp.Rational(25, 16))
an = 8*(n**5+(n+1)**5) + 20*(n**3+(n+1)**3) - 5*(2*n+1)
assert sp.expand(an/8 - at3(n)) == 0
print("(A) a_n/8 = at3(n): OK;  at3-(n+1)^5-n^5 =", sp.factor(sp.expand(at3(n)-(n+1)**5-n**5)))
toF = lambda e: Fr(int(sp.numer(e)), int(sp.denom(e)))
at2F = lambda m: toF(at2(sp.Integer(m))); at3F = lambda m: toF(at3(sp.Integer(m)))
def seqF(at, x0, x1, N):
    S = [Fr(x0), Fr(x1)]
    for mm in range(1, N):
        S.append((at(mm) * S[mm] - mm**5 * S[mm-1]) / (mm+1)**5)
    return S
N = 202
w = seqF(at2F, 1, at2F(0), N+1); v = seqF(at2F, 0, Fr(1, 8), N+1)
w3 = seqF(at3F, 1, at3F(0), N); v3 = seqF(at3F, 0, Fr(1, 8), N)
aI = lambda mm: 8*(mm**5+(mm+1)**5) + 20*(mm**3+(mm+1)**3) - 5*(2*mm+1)
p = [1, aI(0)]; q = [0, 1]          # p_-1, p_0 ; q_-1, q_0
for mm in range(1, 41):
    p.append(aI(mm)*p[-1] - 64*mm**10*p[-2]); q.append(aI(mm)*q[-1] - 64*mm**10*q[-2])
for mm in range(-1, 40):
    fac = 8**(mm+1) * factorial(mm+1)**5
    assert p[mm+1] == fac * w3[mm+1] and q[mm+1] == fac * v3[mm+1], mm
assert all(w3[i+1] > w3[i] for i in range(N-1)) and all(v3[i+1]*w3[i] > v3[i]*w3[i+1] for i in range(N-1))
assert all((w3[i]*v3[i+1]-w3[i+1]*v3[i]) == Fr(1, 8*(i+1)**5) for i in range(N-1))
print("    p_{n-1}=8^n n!^5 w3_n, q_{n-1}=8^n n!^5 v3_n (n<=40); w3, v3/w3 increasing; Casoratian 1/(8(n+1)^5): OK")
# (B)
P1 = lambda mm: -(2*mm+1)*(16*mm**4+56*mm**3+76*mm**2+46*mm+11)
Xm1 = (at2(n)*X0 - (n+1)**5*X1)/n**5
X2 = (at2(n+1)*X1 - (n+1)**5*X0)/(n+2)**5
G = lambda mm, a, b: P1(mm)*a + 32*(mm+1)**5*b
res = sp.cancel((n+1)**5*G(n+1, X1, X2) + n**5*G(n-1, Xm1, X0) - at3(n)*G(n, X0, X1))
assert res == 0
GI = lambda S, mm: P1(mm)*S[mm] + 32*(mm+1)**5*S[mm+1]
for mm in range(N):
    assert GI(w, mm) == w3[mm] and GI(v, mm) == v3[mm] + 4*w3[mm], mm
print("(B) gauge identity: OK; G(w)=w3, G(v)=v3+4w3 for n<%d: OK; G(w)_0,1 =" % N, GI(w,0), GI(w,1),
      "; G(v)_0,1 =", GI(v,0), GI(v,1))
# (C) L4 acting on s^m, theta = -(1-s) d/ds, computed symbolically in m
m = sp.symbols('m')
def th_apply(terms):
    out = {}
    for j, c in terms.items():
        e = m + j
        out[j-1] = out.get(j-1, 0) - c*e
        out[j] = out.get(j, 0) + c*e
    return out
def add(a, b, cb=1):
    o = dict(a)
    for j, c in b.items(): o[j] = o.get(j, 0) + cb*c
    return o
def times_x(a):
    o = {}
    for j, c in a.items():
        o[j] = o.get(j, 0) + c; o[j+1] = o.get(j+1, 0) - c
    return o
Pw = [{0: sp.Integer(1)}]
for _ in range(4): Pw.append(th_apply(Pw[-1]))
sh = {}
for i in range(5): sh = add(sh, Pw[i], sp.binomial(4, i)*h**(4-i))
Lop = add(Pw[4], times_x(sh), -1)
Lop = {j: sp.factor(sp.expand(c)) for j, c in Lop.items() if sp.expand(c) != 0}
print("(C) L4 s^m = sum_j P_j(m) s^(m+j), j in", sorted(Lop), "; P_-3(m) =", Lop[-3])
def local(c0, c1, c2, M=30):
    c = [Fr(c0), Fr(c1), Fr(c2)]
    for pp in range(M):
        tot = Fr(0)
        for j, Pj in Lop.items():
            idx = pp - j
            if idx != pp + 3 and 0 <= idx < len(c): tot += toF(Pj.subs(m, idx)) * c[idx]
        c.append(-tot / toF(Lop[-3].subs(m, pp+3)))
    return c
def theta(c): return [mm*c[mm] - ((mm+1)*c[mm+1] if mm+1 < len(c) else 0) for mm in range(len(c))]
def thv(c, kk):
    out = []; cur = c
    for _ in range(kk): out.append(cur[0]); cur = theta(cur)
    return out
def mul(a, b, M=25):
    r = [Fr(0)]*M
    for i in range(min(M, len(a))):
        for j in range(min(M-i, len(b))): r[i+j] += a[i]*b[j]
    return r
nu = local(0, 1, Fr(3, 4)); r = local(0, 0, 1); rs = local(1, 0, Fr(-1, 16))
print("    nu =", nu[:4], " r =", r[:5], " r* =", rs[:5])
mus = [a-b for a, b in zip(mul(nu, theta(rs)), mul(rs, theta(nu)))]
mu = [a-b for a, b in zip(mul(nu, theta(r)), mul(r, theta(nu)))]
print("    mu* =", mus[:5], " theta^i mu*(1) =", thv(mus[:20], 5))
print("    mu  =", mu[:5], " theta^i mu(1) =", thv(mu[:20], 5))
print("    theta^i r*(1) =", thv(rs[:20], 3), " theta^i nu(1) =", thv(nu[:20], 3))
# Own integration by parts: for operator sum_d x^d Q_d(theta), int_0^1 x^(u-1) x^d theta^j F dx with c = u-1+d:
#   int x^c theta F = F(1) - (c+1) int x^c F  (boundary at 0 vanishes for c > -1, F = O(log^k))
#   => int x^c theta^j F = sum_{i<j} (-(c+1))^(j-1-i) F_i + (-(c+1))^j int x^c F
def bdry(Qs, Fvals):
    B = 0
    for d, Q in enumerate(Qs):
        Qp = sp.Poly(sp.expand(Q), sp.Symbol('t'))
        c1 = u + d
        for (j,), co in Qp.terms():
            for i in range(j): B += co * (-c1)**(j-1-i) * Fvals[i]
    return sp.expand(B)
t = sp.Symbol('t')
B4 = bdry([t**4, -(t+h)**4], [sp.Rational(x.numerator, x.denominator) for x in thv(rs[:20], 4)])
print("    Lemma-3 boundary for r*: int x^(u-1) L4 r* = u^4 rhat(u) - (u+1/2)^4 rhat(u+1) + B, B =", sp.factor(B4),
      "  => u^4 rhat(u) - (u+1/2)^4 rhat(u+1) =", sp.expand(-B4))
assert sp.expand(-B4 - (-sp.Rational(1, 8) - u*(u+h))) == 0
Q5 = [t**5, -2*(t+h)*((t+h)**4 + (t+h)**2 + sp.Rational(1, 16)), (t+1)**5]
for name, MM in (("mu", mu), ("mu*", mus)):
    fv = [sp.Rational(x.numerator, x.denominator) for x in thv(MM[:22], 5)]
    print("    B5(u) for %s =" % name, bdry(Q5, fv))
    assert bdry(Q5, fv) == 0
# general rho = a e0 + b e1 + c e2
a_, b_, c_ = sp.symbols('a b c')
def localsym(c0, c1, c2, M=12):
    c = [c0, c1, c2]
    for pp in range(M):
        tot = 0
        for j, Pj in Lop.items():
            idx = pp - j
            if idx != pp + 3 and 0 <= idx < len(c): tot += Pj.subs(m, idx) * c[idx]
        c.append(sp.expand(-tot / Lop[-3].subs(m, pp+3)))
    return c
rg = localsym(a_, b_, c_)
nuS = [sp.Rational(x.numerator, x.denominator) for x in nu[:15]]
def thS(c): return [sp.expand(mm*c[mm] - ((mm+1)*c[mm+1] if mm+1 < len(c) else 0)) for mm in range(len(c))]
def mulS(a, b, M=10):
    rr = [0]*M
    for i in range(M):
        for j in range(M-i): rr[i+j] += a[i]*b[j]
    return [sp.expand(x) for x in rr]
mug = [sp.expand(x-y) for x, y in zip(mulS(nuS, thS(rg)), mulS(rg, thS(nuS)))]
cur = mug; fv = []
for _ in range(5): fv.append(cur[0]); cur = thS(cur)
print("    B5 for nu theta(rho) - rho theta(nu), rho general holomorphic:", bdry(Q5, fv))
# (D)
Yu, Yu1, Zu, Zu1 = sp.symbols('Yu Yu1 Zu Zu1')
kap = ((u+h)/(u+1))**4
Rs = lambda uu, Y_, Z_: (uu**3 - Y_/8 + Z_/16)/uu**4        # 2 Y*(u)/u^4
ex = u**4*Rs(u, Yu, Zu) - (u+h)**4*Rs(u+1, Yu1, Zu1)
ex = ex.subs({Yu1: (Yu-1)/kap, Zu1: (Zu-1/(u+1))/kap})
assert sp.simplify(ex + sp.Rational(1, 8) + u*(u+h)) == 0
# check (1.2) for Z termwise: A_{k+1}(u)/(u+k+2) = kappa(u) A_k(u+1)/((u+1)+k+1)  -> Z(u) - kappa Z(u+1) = A_0/(u+1)
print("(D) R*(u) satisfies u^4R*(u) - (u+1/2)^4 R*(u+1) = -1/8 - u(u+1/2): OK")
# (E)
beta = lambda j, mm: sp.factorial(j)*sp.prod([1/(mm+i) for i in range(1, j+2)])
Gb = lambda j: P1(n)*beta(j, n) + 32*(n+1)**5*beta(j, n+1)
cn = sp.expand(P1(n) + 32*(n+1)**5); print("(E) c(n) =", cn)
ms = [sp.Rational(x.numerator, x.denominator) for x in mus[:4]]
tot = sp.cancel(sum(ms[j]*Gb(j) for j in range(4)))
nu_, de_ = sp.fraction(tot)
print("    sum m*_j G(beta_j)_n: deg num =", sp.degree(nu_, n), ", deg den =", sp.degree(de_, n),
      "; n*value ->", sp.limit(n*tot, n, sp.oo))
for j in range(4): print("      G(beta_%d)_n ~" % j, sp.limit(Gb(j)/n**(2-j if j < 3 else 0), n, sp.oo), "* n^%d" % (2-j if j < 3 else 0))
lam = sp.Symbol('lam')
mr = [sp.Rational(x.numerator, x.denominator) for x in mu[:4]]
totr = sp.cancel(sum(mr[j]*Gb(j) for j in range(4)))
print("    for mu (r): G(...)_n/n ->", sp.limit(totr/n, n, sp.oo), "; constant term:", sp.limit(totr-128*n, n, sp.oo))
# (F)
lhs = (k+h)**3/(k+1) - k**3/(k-h)
assert sp.simplify(lhs + sp.Rational(1, 8)*(1/(k+1) + 1/(k-h))) == 0
print("(F) telescoping identity: OK")
# (G) Section 7, own bookkeeping
e = sp.Symbol('e'); pi = sp.pi
t0, t1, t2, t3, z0, z1, z2, z3, zeta3, zeta5 = sp.symbols('tau0 tau1 tau2 tau3 zt0 zt1 zt2 zt3 zeta3 zeta5')
N0, N1, N2, N3 = sp.symbols('N0:4')
Nser = N0 + N1*e + N2*e**2 + N3*e**3
NNm = (e*sp.cot(pi*e))**4
NNm = sp.series(NNm, e, 0, 4).removeO()
YN = lambda x: t0 + t1*x + t2*x**2 + t3*x**3
ZN = lambda x: z0 + z1*x + z2*x**2 + z3*x**3
# Y(-e) = N(-e) (Y/N)(-e);  N(e) Y(-e) = N(e)N(-e) (Y/N)(-e)
NYs = sp.expand(-Nser*e**3/2 + NNm*(-YN(-e)/16 + ZN(-e)/32))
hz = sum(sp.binomial(-3, j)*(2**(3+j)-1)*sp.zeta(3+j)*e**j for j in range(4))
hz = hz.subs(sp.zeta(3), zeta3).subs(sp.zeta(5), zeta5)
c1 = sp.expand(NYs).coeff(e, 1)
c3Xi = sp.expand(2*hz*NYs).coeff(e, 3)
rel = {t2: 2*pi**2*t0, t3: sp.Rational(4, 3)*pi**2*t1, N0: 1/pi**2}
S0, S1, Z0, e1 = sp.symbols('S0 S1 Z0 e1')
e0v = pi**2*Z0; e2v = e0v + 8
zrel = {z0: e0v, z1: pi*e1, z2: pi**2*(e0v+e2v), z3: sp.Rational(4, 3)*pi**3*e1 - 16*pi**2}
sub2 = {t0: pi**2*S0, t1: 8*pi**2*S1}
mu4s = sp.simplify((-c1/6).subs(rel).subs(zrel).subs(sub2))
gs0 = sp.simplify((-2*c3Xi).subs(rel).subs(zrel).subs(sub2))
print("(G) zeta(3)-coefficient of g*(0):", sp.simplify(sp.diff(gs0, zeta3)), "; [e^3] N Y*(-e) =",
      sp.simplify(sp.expand(NYs).coeff(e, 3).subs(rel).subs(zrel)))
mu4 = 4*S1/(3*pi**2); X2U = sp.Rational(93, 4)*zeta5 + pi**6*S0/(64*S1); g0 = 192*mu4*X2U
sol = sp.solve(sp.Eq(g0*mu4s - mu4*gs0, sp.Rational(1, 6)), e1); assert len(sol) == 1
Kk = 4 + pi**2*(Z0 - 2*S0)
ell3 = sp.simplify((gs0/(192*mu4s) - 4).subs(e1, sol[0]))
final = sp.Rational(93, 4)*zeta5 - 4 + pi**6/64*S0*Kk/(1 + S1*Kk)
print("    e1 =", sp.factor(sol[0]), ";  ell3 - stated formula =", sp.simplify(ell3 - final))
inv2 = X2U
third_stated = inv2 - 4 + pi**6*S0/(64*S1*(1 + S1*Kk))
third_fixed = inv2 - 4 - pi**6*S0/(64*S1*(1 + S1*Kk))
print("    third form as stated  - formula =", sp.simplify(third_stated - final))
print("    third form, sign flipped - formula =", sp.simplify(third_fixed - final))
print("a1 done")
