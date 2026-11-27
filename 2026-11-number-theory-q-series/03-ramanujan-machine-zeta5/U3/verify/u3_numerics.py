"""U3: numerical confirmation of every analytic step of ../THEOREM.md (not part of the proof).
 (1) r* = 2 sum_j (-1)^j eta*_{3-j} y_j on (0,1), eta*_i = [u^i] Y*(u), Y*(u) = u^3/2 - Y(u)/16 + Z(u)/32 (vs local series at 1/2);
 (2) Mellin transform of r*: 2 Y*(u)/u^4;
 (3) g*(u) = int_0^1 mu* x^u dx: functional equation; g*(0) = -2[e^3] Xi(e)Y*(-e); mu4* = -(1/6)[e^1] N(e)Y*(-e);
     g*(1) = atil(0) g*(0) - 24 mu4*;  g*(n) = g*(0) w_n - 192 mu4* v_n;  G(g*)_n -> 0, G(g)_n/n -> 128;
 (4) bilateral sum T_Z(u) = (1+t^2)(e0 + e1 t + e2 t^2), t = tan(pi u); e0 = pi^2 Z0, e2 = e0 + 8; zeta_i relations;
 (5) Legendre relation g(0) mu4* - mu4 g*(0) = 1/6;  e1 = 8(1 + 4 S1 + pi^2 S1 Z0)/(pi S0);
 (6) final formula against the 48-digit value of the continued fraction (and its direct evaluation)."""
import os, pickle
from mpmath import mp, mpf, mpc, pi, gamma, zeta, log, quad, hyper, taylor, tan, cos, factorial
from fractions import Fraction as Fr
from math import comb
here = os.path.dirname(os.path.abspath(__file__)); code = os.path.join(here, '..', '..', 'code')
mp.dps = 50
half = mpf(1) / 2
K = 260
def pmul(a, b):
    r = [mpf(0)] * 4
    for i in range(4):
        for j in range(4 - i): r[i + j] += a[i] * b[j]
    return r
def pinv(a):
    r = [1 / a[0], mpf(0), mpf(0), mpf(0)]
    for n in range(1, 4): r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r
Acoef = [[mpf(1), mpf(0), mpf(0), mpf(0)]]
for k in range(K):
    f = pmul([k + half, mpf(1), mpf(0), mpf(0)], pinv([mpf(k + 1), mpf(1), mpf(0), mpf(0)]))
    f2 = pmul(f, f); Acoef.append(pmul(Acoef[-1], pmul(f2, f2)))
def frob(x, mmax=4):
    lx = log(x); xe = [lx**i / factorial(i) for i in range(4)]
    vals = [[mpf(0)] * 4 for _ in range(mmax)]; xp = mpf(1)
    for k in range(K):
        cur = pmul(Acoef[k], xe)
        for m in range(mmax):
            for j in range(4): vals[m][j] += cur[j] * xp
            cur = pmul(cur, [mpf(k), mpf(1), mpf(0), mpf(0)])
        xp *= x
        if xp < mpf(10)**(-mp.dps - 5): break
    return vals
M = K
def theta(c): return [(m * c[m]) - ((m + 1) * c[m + 1] if m + 1 < len(c) else 0) for m in range(len(c))]
def L4w(c, p):
    th = [c]
    for _ in range(4): th.append(theta(th[-1]))
    sh = [sum(comb(4, j) * half**(4 - j) * th[j][m] for j in range(5)) for m in range(len(c))]
    return th[4][p] - (sh[p] - (sh[p - 1] if p >= 1 else 0))
def solve(c0, c1, c2):
    c = [mpf(c0), mpf(c1), mpf(c2)] + [mpf(0)] * (M - 3)
    for p in range(0, M - 5):
        w = c[max(0, p - 1):p + 6]; off = max(0, p - 1)
        w[p + 3 - off] = mpf(0); w[p + 4 - off] = mpf(0); w[p + 5 - off] = mpf(0)
        full = [mpf(0)] * off + w
        c[p + 3] = -L4w(full, p) / ((p + 1) * (p + 2)**2 * (p + 3))
    return c
nu_c = solve(0, 1, mpf(3) / 4); r2_c = solve(0, 0, 1); rs_c = solve(1, 0, -mpf(1) / 16)
def loc(c, x, mmax=4):
    s = 1 - x; out = []; cur = c
    for m in range(mmax):
        out.append(sum(cur[i] * s**i for i in range(len(cur) - 8))); cur = theta(cur)
    return out
Nf = lambda e: (gamma(1 + e) / gamma(half + e))**4
Yf = lambda e: hyper([1, half + e, half + e, half + e, half + e], [1 + e] * 4, 1)
Zf = lambda e: hyper([1, half + e, half + e, half + e, half + e], [1 + e, 1 + e, 1 + e, 2 + e], 1) / (1 + e)
Ztf = lambda v: hyper([1, half + v, half + v, half + v, v - half], [1 + v] * 4, 1) / (v - half)    # sum_j A_j(v)/(v+j-1/2)
Ysf = lambda e: e**3 / 2 - Yf(e) / 16 + Zf(e) / 32
Xif = lambda e: 2 * zeta(3, half + e) * Nf(e)
mp.dps = 90
Nt = taylor(Nf, 0, 3); Xt = taylor(Xif, 0, 3); Yt = taylor(Yf, 0, 3); Yst = taylor(Ysf, 0, 3)
tau = taylor(lambda e: Yf(e) / Nf(e), 0, 3); zt = taylor(lambda e: Zf(e) / Nf(e), 0, 3)
mp.dps = 50
Nt = [+t for t in Nt]; Xt = [+t for t in Xt]; eta = [+t for t in Yt]; etas = [+t for t in Yst]; tau = [+t for t in tau]; zt = [+t for t in zt]
c_nu = [(-1)**j * Nt[3 - j] for j in range(4)]
c_r2 = [2 * (-1)**j * eta[3 - j] for j in range(4)]
c_rs = [2 * (-1)**j * etas[3 - j] for j in range(4)]
tol = mpf(10)**-22
x0 = half; FV = frob(x0); lr = loc(rs_c, x0)
d1 = max(abs(sum(c_rs[j] * FV[m][j] for j in range(4)) - lr[m]) for m in range(4))
print("(1) max |Frobenius form - local series| for r* at x=1/2 (values, theta^1..3): %s" % mp.nstr(d1, 3)); assert d1 < tol
def fn(c_fro, c_loc, x, m=0):
    if x < half:
        FVx = frob(x, m + 1); return sum(c_fro[j] * FVx[m][j] for j in range(4))
    return loc(c_loc, x, m + 1)[m]
nu = lambda x, m=0: fn(c_nu, nu_c, x, m); r2 = lambda x, m=0: fn(c_r2, r2_c, x, m); rs = lambda x, m=0: fn(c_rs, rs_c, x, m)
mu = lambda x: nu(x) * r2(x, 1) - r2(x) * nu(x, 1)
mus = lambda x: nu(x) * rs(x, 1) - rs(x) * nu(x, 1)
mp.dps = 30
pts = [0, mpf(1) / 1000, mpf(1) / 50, mpf(1) / 5, half, mpf(4) / 5, 1]
for u in (mpf(13) / 10, mpf(2)):
    c_ = quad(lambda x: rs(x) * x**(u - 1), pts); d_ = 2 * Ysf(u) / u**4
    print("(2) u=%s: r*^ - 2Y*(u)/u^4 = %s" % (mp.nstr(u, 3), mp.nstr(c_ - d_, 3))); assert abs(c_ - d_) < mpf(10)**-15
g = lambda u: quad(lambda x: mu(x) * x**u, pts); gs = lambda u: quad(lambda x: mus(x) * x**u, pts)
atil = lambda u: 2 * (u + half) * ((u + half)**4 + (u + half)**2 + mpf(1) / 16)
for u in (mpf(7) / 10, mpf(19) / 10):
    r = u**5 * gs(u - 1) + (u + 1)**5 * gs(u + 1) - atil(u) * gs(u)
    print("(3) functional equation residual of g* at u=%s: %s" % (mp.nstr(u, 3), mp.nstr(r, 3))); assert abs(r) < mpf(10)**-13
NG = 31
gv = [g(n) for n in range(NG + 1)]; gsv = [gs(n) for n in range(NG + 1)]
alt3 = lambda X, Yc: sum(X[3 - i] * Yc[i] * (-1)**i for i in range(4))          # [e^3] X(e) Yc(-e)
mu4 = (Nt[0] * eta[1] - Nt[1] * eta[0]) / 6; mu4s = (Nt[0] * etas[1] - Nt[1] * etas[0]) / 6
print("    g*(0) + 2[e^3]Xi(e)Y*(-e) = %s ;  g*(1) - atil(0)g*(0) + 24 mu4* = %s ;  mu*(x)/log^4 x at 1e-12: %s vs mu4* = %s"
      % (mp.nstr(gsv[0] + 2 * alt3(Xt, etas), 3), mp.nstr(gsv[1] - atil(0) * gsv[0] + 24 * mu4s, 3), mp.nstr(mus(mpf(10)**-12) / log(mpf(10)**-12)**4, 8), mp.nstr(mu4s, 8)))
assert abs(gsv[0] + 2 * alt3(Xt, etas)) < mpf(10)**-15 and abs(gsv[1] - atil(0) * gsv[0] + 24 * mu4s) < mpf(10)**-13
f2 = lambda m: Fr(8 * (m**5 + (m + 1)**5) - 12 * (m**3 + (m + 1)**3) + 7 * (2 * m + 1), 8)
def seq(at, x0, x1, N):
    s = [Fr(x0), Fr(x1)]
    for m in range(1, N): s.append((at(m) * s[m] - m**5 * s[m - 1]) / (m + 1)**5)
    return s
w = seq(f2, 1, f2(0), NG + 1); v = seq(f2, 0, Fr(1, 8), NG + 1)
toM = lambda q: mpf(q.numerator) / q.denominator
for n in (2, 5, 9):
    d = gsv[n] - (gsv[0] * toM(w[n]) - 192 * mu4s * toM(v[n]))
    print("    g*(%d) - (g*(0) w_n - 192 mu4* v_n) = %s" % (n, mp.nstr(d, 3))); assert abs(d) < mpf(10)**-12
P1 = lambda m: -(2 * m + 1) * (16 * m**4 + 56 * m**3 + 76 * m**2 + 46 * m + 11)
G = lambda s, m: P1(m) * s[m] + 32 * (m + 1)**5 * s[m + 1]
print("    G(g*)_n, n = 5, 10, 20, 30: %s   (theory: O(1/n); minimal U3 solution, ~ n^-5)" % [mp.nstr(G(gsv, m), 4) for m in (5, 10, 20, 30)])
print("    G(g)_n / n, n = 5, 10, 20, 30: %s   (-> 128)" % [mp.nstr(G(gv, m) / m, 8) for m in (5, 10, 20, 30)])
assert abs(G(gsv, 30)) < mpf(10)**-3 and abs(G(gv, 30) / 30 - 128) < 20
# (4)
mp.dps = 40
S0 = eta[0]; D = eta[1] - 8 * log(2) * eta[0]; S1 = D / 8; Z0 = Zf(mpf(0))
e0 = zt[0]; e1 = zt[1] / pi; e2 = zt[2] / pi**2 - e0
print("(4) e0 - pi^2 Z0 = %s ; e2 - e0 - 8 = %s ; zeta_3 - (4 pi^3/3) e1 + 16 pi^2 = %s"
      % (mp.nstr(e0 - pi**2 * Z0, 3), mp.nstr(e2 - e0 - 8, 3), mp.nstr(zt[3] - 4 * pi**3 / 3 * e1 + 16 * pi**2, 3)))
assert abs(e0 - pi**2 * Z0) < mpf(10)**-25 and abs(e2 - e0 - 8) < mpf(10)**-18 and abs(zt[3] - 4 * pi**3 / 3 * e1 + 16 * pi**2) < mpf(10)**-15
X1 = hyper([half, half, half, -half], [1, 1, 1], 1)
print("    Z0 - 2 sum a_k^2/(1-2k) + 8/pi^2 = %s" % mp.nstr(Z0 - 2 * X1 + 8 / pi**2, 3)); assert abs(Z0 - 2 * X1 + 8 / pi**2) < mpf(10)**-30
mp.dps = 25
for uu in (mpf(1) / 5, mpf(-3) / 10, mpc('0.15', '0.2')):
    vv = half - uu
    T = Zf(uu) / Nf(uu) - tan(pi * uu)**4 * Ztf(vv) / Nf(vv)
    t = tan(pi * uu); R = (1 + t * t) * (e0 + e1 * t + e2 * t * t)
    print("    T_Z(u) - (1+t^2)(e0+e1 t+e2 t^2) at u=%s: %s" % (mp.nstr(uu, 3), mp.nstr(abs(T - R), 3))); assert abs(T - R) < mpf(10)**-13
# (5)
mp.dps = 30
leg = gv[0] * mu4s - mu4 * gsv[0]
print("(5) g(0) mu4* - mu4 g*(0) = %s  (1/6 = %s)" % (mp.nstr(leg, 20), mp.nstr(mpf(1) / 6, 20))); assert abs(leg - mpf(1) / 6) < mpf(10)**-14
e1f = 8 * (1 + 4 * S1 + pi**2 * S1 * Z0) / (pi * S0)
print("    e1 = %s ; 8(1 + 4 S1 + pi^2 S1 Z0)/(pi S0) = %s" % (mp.nstr(e1, 25), mp.nstr(e1f, 25))); assert abs(e1 - e1f) < mpf(10)**-18
print("    192 mu4* = %s ; 8(1 + S1 K)/(pi^4 S0) = %s" % (mp.nstr(192 * mu4s, 20), mp.nstr(8 * (1 + S1 * (4 + pi**2 * (Z0 - 2 * S0))) / (pi**4 * S0), 20)))
# (6)
mp.dps = 48
C = {k: [mpf(x) for x in vv_] for k, vv_ in pickle.load(open(os.path.join(code, 'u_consts.pkl'), 'rb')).items()}
U3 = C['U3'][0]
mp.dps = 60
Ls = {k: mpf(v_) for k, v_ in pickle.load(open(os.path.join(code, 'f8_L.pkl'), 'rb')).items()}
S0h = 16 * Ls['L2'] / pi**2; S1h = -2 * Ls['L1']; Z0h = hyper([half] * 4 + [1], [1, 1, 1, 2], 1)
Kc = 4 + pi**2 * (Z0h - 2 * S0h)
val = mpf(93) / 4 * zeta(5) - 4 + pi**6 / 64 * S0h * Kc / (1 + S1h * Kc)
print("(6) K = %s ;  1/U3(formula) = %s ;  U3(formula) = %s" % (mp.nstr(Kc, 30), mp.nstr(val, 40), mp.nstr(1 / val, 45)))
print("    U3 (table value, 48 digits)                = %s ;  difference = %s" % (mp.nstr(U3, 45), mp.nstr(U3 - 1 / val, 3)))
assert abs(U3 - 1 / val) < mpf(10)**-40
val2 = gsv[0] / (192 * mu4s) - 4
mp.dps = 20
print("    Level A: g*(0)/(192 mu4*) - 4 = %s  vs 1/U3 = %s" % (mp.nstr(val2, 18), mp.nstr(1 / U3, 18)))
print("U3 numerics PASS")
