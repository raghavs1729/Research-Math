"""V4. Numerical confirmation of every analytic step of the proof of Theorem U2 (not part of the proof).
 (1) Frobenius coordinates: nu = sum_j (-1)^j N_{3-j} y_j,  r2 = 2 sum_j (-1)^j eta_{3-j} y_j  on (0,1)
     (N(e) = (Gamma(1+e)/Gamma(1/2+e))^4, Y(e) = sum_k ((1/2+e)_k/(1+e)_k)^4 = sum eta_i e^i), checked at x = 1/2
     against the local series at x = 1, values and theta-derivatives.
 (2) Mellin transforms: nu^(u) = Gamma(u)^4/Gamma(u+1/2)^4, r2^(u) = 2Y(u)/u^4 (quadrature), u = 1.3, 2.
 (3) g(u) = int_0^1 mu x^u dx satisfies u^5 g(u-1) + (u+1)^5 g(u+1) = atil(u) g(u) (u = 0.7, 1.9);
     g(0) = -2 [e^3] Xi(e) Y(-e);  mu ~ mu4 log^4 x at 0 with mu4 = (N_0 eta_1 - N_1 eta_0)/6;  g(1) = atil(0) g(0) - 24 mu4.
 (4) Bilateral sum: T(u) := Y(u)/N(u) + tan^4(pi u) Y(1/2-u)/N(1/2-u) = sec^2(pi u)[tau_0 sec^2(pi u) + (tau_1/pi) tan(pi u)],
     tau_i = [u^i] Y/N;  consequences tau_2 = 2 pi^2 tau_0, tau_3 = (4 pi^2/3) tau_1.
 (5) Final: 1/U2 = (93/4) zeta(5) + (pi^6/8) Y(0)/(Y'(0) - 8 log2 Y(0)) against the 300-digit value of the fraction;
     Y(0) = 16 L(f,2)/pi^2,  Y'(0) - 8 log 2 Y(0) = -16 L(f,1)."""
import sys, os, pickle
from mpmath import mp, mpf, mpc, pi, gamma, zeta, log, quad, hyper, taylor, tan, cos, factorial, matrix, lu_solve
here = os.path.dirname(os.path.abspath(__file__))
code = os.path.join(here, '..', '..', 'code')
mp.dps = 50
half = mpf(1) / 2
K = 260
# ---- Frobenius basis at 0 (truncated eps-series), values of theta^m y_j at x
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
    return vals            # vals[m][j] = theta^m y_j (x)
# ---- local series at x = 1
M = K
def theta(c): return [(m * c[m]) - ((m + 1) * c[m + 1] if m + 1 < len(c) else 0) for m in range(len(c))]
from math import comb
def L4w(c, p):       # coefficient of s^p in L4(c), c a short window
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
nu_c = solve(0, 1, mpf(3) / 4); r2_c = solve(0, 0, 1)
def loc(c, x, mmax=4):
    s = 1 - x; out = []; cur = c
    for m in range(mmax):
        out.append(sum(cur[i] * s**i for i in range(len(cur) - 8))); cur = theta(cur)
    return out
# ---- generating functions
Nf = lambda e: (gamma(1 + e) / gamma(half + e))**4
Yf = lambda e: hyper([1, half + e, half + e, half + e, half + e], [1 + e] * 4, 1)
Xif = lambda e: 2 * zeta(3, half + e) * Nf(e)
mp.dps = 90
Nt = taylor(Nf, 0, 3); Xt = taylor(Xif, 0, 3); Yt = taylor(Yf, 0, 3)
mp.dps = 50
Nt = [+t for t in Nt]; Xt = [+t for t in Xt]; eta = [+t for t in Yt]
c_nu = [(-1)**j * Nt[3 - j] for j in range(4)]
c_r2 = [2 * (-1)**j * eta[3 - j] for j in range(4)]
tol = mpf(10)**-22
# (1)
x0 = half
FV = frob(x0); ln = loc(nu_c, x0); lr = loc(r2_c, x0)
d1 = max(abs(sum(c_nu[j] * FV[m][j] for j in range(4)) - ln[m]) for m in range(4))
d2 = max(abs(sum(c_r2[j] * FV[m][j] for j in range(4)) - lr[m]) for m in range(4))
print("(1) max |Frobenius form - local series| at x=1/2 (values and theta^1..3):  nu: %s   r2: %s" % (mp.nstr(d1, 3), mp.nstr(d2, 3)))
assert d1 < tol and d2 < tol
def fn(c_fro, c_loc, x, m=0):
    if x < half: 
        FVx = frob(x, m + 1); return sum(c_fro[j] * FVx[m][j] for j in range(4))
    return loc(c_loc, x, m + 1)[m]
nu = lambda x, m=0: fn(c_nu, nu_c, x, m); r2 = lambda x, m=0: fn(c_r2, r2_c, x, m)
mu = lambda x: nu(x) * r2(x, 1) - r2(x) * nu(x, 1)
mp.dps = 30
pts = [0, mpf(1) / 1000, mpf(1) / 50, mpf(1) / 5, half, mpf(4) / 5, 1]
# (2)
for u in (mpf(13) / 10, mpf(2)):
    a_ = quad(lambda x: nu(x) * x**(u - 1), pts); b_ = gamma(u)**4 / gamma(u + half)**4
    c_ = quad(lambda x: r2(x) * x**(u - 1), pts); d_ = 2 * Yf(u) / u**4
    print("(2) u=%s: nu^ - Gamma ratio = %s ;  r2^ - 2Y(u)/u^4 = %s" % (mp.nstr(u, 3), mp.nstr(a_ - b_, 3), mp.nstr(c_ - d_, 3)))
    assert abs(a_ - b_) < mpf(10)**-15 and abs(c_ - d_) < mpf(10)**-15
# (3)
g = lambda u: quad(lambda x: mu(x) * x**u, pts)
atil = lambda u: 2 * (u + half) * ((u + half)**4 + (u + half)**2 + mpf(1) / 16)
for u in (mpf(7) / 10, mpf(19) / 10):
    r = u**5 * g(u - 1) + (u + 1)**5 * g(u + 1) - atil(u) * g(u)
    print("(3) functional equation residual at u=%s: %s" % (mp.nstr(u, 3), mp.nstr(r, 3)))
    assert abs(r) < mpf(10)**-13
g0 = g(0); g1 = g(1)
alt3 = sum(Xt[3 - i] * eta[i] * (-1)**i for i in range(4))
mu4 = (Nt[0] * eta[1] - Nt[1] * eta[0]) / 6
xs = mpf(10)**-12
print("    g(0) + 2[e^3]Xi(e)Y(-e) = %s ;  mu(x)/log^4 x at x=1e-12: %s vs mu4 = %s ;  g(1) - atil(0)g(0) + 24 mu4 = %s"
      % (mp.nstr(g0 + 2 * alt3, 3), mp.nstr(mu(xs) / log(xs)**4, 8), mp.nstr(mu4, 8), mp.nstr(g1 - atil(0) * g0 + 24 * mu4, 3)))
assert abs(g0 + 2 * alt3) < mpf(10)**-15 and abs(g1 - atil(0) * g0 + 24 * mu4) < mpf(10)**-13
print("    PCF = 192 mu4/g(0) = %s" % mp.nstr(192 * mu4 / g0, 20))
# (4)
mp.dps = 50
mp.dps = 90
tau = taylor(lambda e: Yf(e) / Nf(e), 0, 3)
mp.dps = 40
tau = [+t for t in tau]
print("(4) tau_2 - 2 pi^2 tau_0 = %s ;  tau_3 - (4 pi^2/3) tau_1 = %s" % (mp.nstr(tau[2] - 2 * pi**2 * tau[0], 3), mp.nstr(tau[3] - 4 * pi**2 / 3 * tau[1], 3)))
assert abs(tau[2] - 2 * pi**2 * tau[0]) < mpf(10)**-18 and abs(tau[3] - 4 * pi**2 / 3 * tau[1]) < mpf(10)**-16
mp.dps = 25
for uu in (mpf(1) / 5, mpf(-3) / 10, mpc('0.15', '0.2')):
    T = Yf(uu) / Nf(uu) + tan(pi * uu)**4 * Yf(half - uu) / Nf(half - uu)
    R = (tau[0] / cos(pi * uu)**2 + tau[1] / pi * tan(pi * uu)) / cos(pi * uu)**2
    print("    T(u) - closed form at u=%s: %s" % (mp.nstr(uu, 3), mp.nstr(abs(T - R), 3)))
    assert abs(T - R) < mpf(10)**-14
# (5)
mp.dps = 300
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open(os.path.join(code, 'u_consts.pkl'), 'rb')).items()}
Ls = {k: mpf(v) for k, v in pickle.load(open(os.path.join(code, 'f8_L.pkl'), 'rb')).items()}
U2 = C['U2'][0]; L1, L2 = Ls['L1'], Ls['L2']
Y0 = eta[0]; D = eta[1] - 8 * log(2) * eta[0]
val = mpf(93) / 4 * zeta(5) + pi**6 / 8 * Y0 / D
print("(5) 1/U2 - [(93/4) zeta(5) + (pi^6/8) Y(0)/(Y'(0) - 8 log2 Y(0))] = %s   (Y, Y' from mpmath, ~25 digits)" % mp.nstr(1 / U2 - val, 3))
print("    Y(0) - 16 L(f,2)/pi^2 = %s ;  Y'(0) - 8 log2 Y(0) + 16 L(f,1) = %s" % (mp.nstr(Y0 - 16 * L2 / pi**2, 3), mp.nstr(D + 16 * L1, 3)))
print("    U2 - 8 L(f,1)/(186 zeta(5) L(f,1) - pi^4 L(f,2)) = %s   (300-digit data)" % mp.nstr(U2 - 8 * L1 / (186 * zeta(5) * L1 - pi**4 * L2), 3))
assert abs(1 / U2 - val) < mpf(10)**-18 and abs(U2 - 8 * L1 / (186 * zeta(5) * L1 - pi**4 * L2)) < mpf(10)**-290
print("V4 PASS")
