"""Numerical confirmation of every step of ../THEOREM-2.md (Zudilin's relations (6) for (r,t) = (1/2,1/4) and (1/2,1/3)).
Not part of the proof.  Sections [A] t = 1/4,  [B] t = 1/3;  see the statements printed."""
from mpmath import (mp, mpf, mpc, pi, gamma, taylor, hyper, jtheta, exp, sqrt, log, quad, inf, nsum, nprod, diff, nstr, gammainc, psi)
mp.dps = 30
half = mpf(1) / 2
TOL = mpf(10)**-15
def Fj(r, t):
    r = mpf(r); t = mpf(t); pref = 1 / (gamma(r) * gamma(1 - r) * gamma(t) * gamma(1 - t))
    def G(e):
        c = gamma(r + e) * gamma(1 - r + e) * gamma(t + e) * gamma(1 - t + e) / gamma(1 + e)**4
        return pref * c * hyper([r + e, 1 - r + e, t + e, 1 - t + e, 1], [1 + e] * 4, 1)
    old = mp.dps; mp.dps = old + 40; v = taylor(G, 0, 3); mp.dps = old
    return [+x for x in v]
def Useries(t, x):
    """U_0, U_1, U_2 at x: Taylor coefficients in eps of sum_n u_t(n+eps) x^(n+eps)."""
    t = mpf(t)
    def G(e):
        u0 = gamma(half + e) * gamma(t + e) * gamma(1 - t + e) / (gamma(half) * gamma(t) * gamma(1 - t) * gamma(1 + e)**3)
        return u0 * x**e * hyper([half + e, t + e, 1 - t + e, 1], [1 + e] * 3, x)
    old = mp.dps; mp.dps = old + 30; v = taylor(G, 0, 2); mp.dps = old
    return [+c for c in v]
NQ = 500
def eta_coeffs(expo):
    lead = sum(d * e for d, e in expo.items()); assert lead % 24 == 0
    co = [0] * (NQ + 1); co[0] = 1
    for d, e in expo.items():
        for m in range(1, NQ // d + 1):
            step = d * m
            for _ in range(abs(e)):
                if e > 0:
                    for i in range(NQ, step - 1, -1): co[i] -= co[i - step]
                else:
                    for i in range(step, NQ + 1): co[i] += co[i - step]
    sh = lead // 24
    return {n + sh: co[n] for n in range(NQ + 1 - sh) if co[n]}
def comb(*terms):
    out = {}
    for mult, expo in terms:
        for n, v in eta_coeffs(expo).items(): out[n] = out.get(n, 0) + mult * v
    return {n: v for n, v in out.items() if v and n <= NQ - 20}
def qval(coef, tau):                 # sum a_n e^{2 pi i n tau}
    q = exp(2j * pi * tau); return sum(c * q**n for n, c in coef.items())
def Lfe(coef, lev, sign, s, split=1.0, k=4):
    """L(f,s) assuming f(i/(lev t)) = sign * lev^(k/2) t^k f(it); 'split' moves the cut point (result must not depend on it)."""
    t0 = mpf(split) / sqrt(lev)
    tot = mpf(0)
    for n, c in sorted(coef.items()):
        a = 2 * pi * n
        if a * t0 > 90: break
        tot += c * (gammainc(s, a * t0) / a**s + sign * lev**(k / mpf(2) - s) * gammainc(k - s, a / (lev * t0)) / a**(k - s))
    return tot * (2 * pi)**s / gamma(s)
def th(y): 
    q = exp(-pi * y); return jtheta(2, 0, q), jtheta(3, 0, q), jtheta(4, 0, q)
# ------------------------------------------------------------------ [A]  t = 1/4
print("[A] (r,t) = (1/2,1/4)")
f2 = comb((1, {2: 4, 4: 4}))
f1 = comb((1, {4: 16, 2: -4, 8: -4}))
chi4 = lambda n: [0, 1, 0, -1][n % 4]
assert all(f1.get(n, 0) == f2.get(n, 0) * chi4(n) for n in range(1, 400))
print("  f1 := eta4^16/(eta2^4 eta8^4) = f2 (x) chi_{-4}, f2 = eta2^4 eta4^4: coefficients agree for n < 400")
F = Fj(half, mpf(1) / 4)
for y in (mpf('1.6'), mpf('2.3')):
    t2, t3, t4 = th(y); P = t3**4 + t2**4; X = 16 * t2**4 * t3**4 * t4**8 / P**4
    U = Useries(mpf(1) / 4, X)
    d0 = abs(U[0] - P); d1 = abs(U[1] + pi * y * P); d2 = abs(U[2] - (pi**2 * y**2 / 2 + pi**2) * P)
    Xf = lambda yy: (lambda a, b, c: 16 * a**4 * b**4 * c**8 / (b**4 + a**4)**4)(*th(yy))
    dX = diff(Xf, y); d3 = abs(dX + 4 * pi * t2**2 * t3**2 * t4**4 / P * sqrt(X * (1 - X)))
    d4 = abs(t2**2 * t3**2 * t4**4 - 4 * qval(f2, 1j * y / 4))
    print("  y = %s: X = %s; |U0 - P| = %s, |U1 + pi y P| = %s, |U2 - (pi^2 y^2/2 + pi^2) P| = %s, dX/dy: %s, G1 = 4 f2(iy/4): %s" % (nstr(y, 3), nstr(X, 8), nstr(d0, 2), nstr(d1, 2), nstr(d2, 2), nstr(d3, 2), nstr(d4, 2)))
    assert max(d0, d1, d2, d4) < TOL and d3 < mpf(10)**-12
t2, t3, t4 = th(sqrt(2)); assert abs(16 * t2**4 * t3**4 * t4**8 / (t3**4 + t2**4)**4 - 1) < TOL
G1 = lambda y: (lambda a, b, c: a**2 * b**2 * c**4)(*th(y))
pA = [lambda y: mpf(1), lambda y: -pi * y, lambda y: pi**2 * y**2 / 2 + pi**2]
for j in range(3):
    v = 4 * quad(lambda y: pA[j](y) * G1(y), [sqrt(2), 2, 4, 8, 16, 40])
    print("  F_%d(1) = %s ; 4 int_{sqrt2}^oo p_%d(y) G1(iy) dy = %s" % (j, nstr(F[j], 20), j, nstr(v, 20))); assert abs(v - F[j]) < TOL
f2t = lambda t: qval(f2, 1j * t).real
r8 = 1 / sqrt(mpf(8))
rho0 = quad(lambda t: f2t(t) * (1 + 8 * t * t), [r8, 1, 2, 4, 8]); rho1 = 2 * quad(lambda t: f2t(t) * t, [r8, 1, 2, 4, 8])
L21 = 2 * pi * rho0; L22 = 4 * pi**2 * rho1
print("  L(f2,1) = %s, L(f2,2) = %s (Mellin integrals folded with the Fricke relation)" % (nstr(L21, 20), nstr(L22, 20)))
print("    diffs:", nstr(F[1] + 32 * L22 / pi, 3), nstr(F[2] - 32 * pi * L21, 3), nstr(F[3] - 7 * pi**2 / 3 * F[1], 3))
assert abs(F[1] + 32 * L22 / pi) < TOL and abs(F[2] - 32 * pi * L21) < TOL and abs(F[3] - 7 * pi**2 / 3 * F[1]) < mpf(10)**-10
print("  F_1 = -(32/pi) L(f2,2), F_2 = 32 pi L(f2,1), F_3 = (7 pi^2/3) F_1 : OK")
tt = mpc('0.11', '0.63'); f2c = lambda z: qval(f2, z)
assert abs(f2c(tt / (4 * tt + 1)) / ((4 * tt + 1)**4 * f2c(tt)) + 1) < mpf(10)**-12 and abs(f2c(-1 / (8 * tt)) / (64 * tt**4 * f2c(tt)) - 1) < mpf(10)**-12
print("  f2(tau/(4tau+1)) = -(4tau+1)^4 f2(tau), f2(-1/(8tau)) = 64 tau^4 f2(tau): OK at a sample point")
L1 = [Lfe(f1, 16, 1, s) for s in (1, 2, 3)]; L1b = [Lfe(f1, 16, 1, s, split=1.25) for s in (1, 2, 3)]
assert max(abs(a - b) for a, b in zip(L1, L1b)) < mpf(10)**-15
assert abs(L1[0] - 4 / pi * L22) < TOL and abs(L1[1] - pi * L21) < TOL and abs(L1[2] - pi / 2 * L22) < TOL
print("  L(f1,1) = (4/pi) L(f2,2), L(f1,2) = pi L(f2,1), L(f1,3) = (pi/2) L(f2,2): OK   (L(f1,k) = %s)" % [nstr(v, 15) for v in L1])
rat = [L1[i] / F[i + 1] for i in range(3)]
print("  L(f1,k)/F_k(1) = %s   (expected -1/8, 1/32, -3/448)" % [nstr(v, 18) for v in rat])
assert abs(rat[0] + mpf(1) / 8) < TOL and abs(rat[1] - mpf(1) / 32) < TOL and abs(rat[2] + mpf(3) / 448) < mpf(10)**-10
# ------------------------------------------------------------------ [B]  t = 1/3
print("[B] (r,t) = (1/2,1/3)")
gt = comb((1, {2: 6, 6: 2}), (9, {2: 3, 6: 2, 18: 3}))
fz = comb((1, {6: 14, 2: -3, 18: -3}), (-3, {2: 3, 6: 2, 18: 3}))
chi3 = lambda n: [0, 1, -1][n % 3]
assert all(fz.get(n, 0) == gt.get(n, 0) * chi3(n) for n in range(1, 400))
print("  Zudilin's f = eta6^14/(eta2^3 eta18^3) - 3 eta2^3 eta6^2 eta18^3 = g~ (x) chi_{-3},  g~ = eta2^6 eta6^2 + 9 eta2^3 eta6^2 eta18^3: n < 400")
def eta(tau):
    q = exp(2j * pi * tau); return exp(2j * pi * tau / 24) * nprod(lambda n: 1 - q**n, [1, inf])
def acub(tau):
    q = exp(2j * pi * tau); return nsum(lambda m, n: q**(m * m + m * n + n * n), [-30, 30], [-30, 30])
g0 = lambda tau: acub(tau) * eta(tau)**3 * eta(3 * tau)**3
tt = mpc('0.13', '0.71')
assert abs(g0(tt) - qval(gt, tt / 2)) < mpf(10)**-15
assert abs(g0(tt / (3 * tt + 1)) / ((3 * tt + 1)**4 * g0(tt)) + 1) < mpf(10)**-12 and abs(g0(tt + 1) / g0(tt) + 1) < mpf(10)**-12 and abs(g0(-1 / (3 * tt)) / (9 * tt**4 * g0(tt)) - 1) < mpf(10)**-12
print("  g0 := a(q) eta(tau)^3 eta(3tau)^3 = g~(tau/2); g0(tau/(3tau+1)) = -(3tau+1)^4 g0(tau), g0(tau+1) = -g0(tau), g0(-1/(3tau)) = 9 tau^4 g0(tau): OK")
F = Fj(half, mpf(1) / 3)
def cub(y):
    tau = 1j * y; a = acub(tau).real; e1 = eta(tau).real; e3 = eta(3 * tau).real
    return a, e1**3 / e3, 3 * e3**3 / e1, e1, e3
for y in (mpf('0.8'), mpf('1.3')):
    a, b, c, e1, e3 = cub(y); hh = c**3 / a**3; x = 4 * hh * (1 - hh)
    U = Useries(mpf(1) / 3, x)
    d0 = abs(U[0] - a * a); d1 = abs(U[1] + 2 * pi * y * a * a); d2 = abs(U[2] - (2 * pi**2 * y**2 + 2 * pi**2 / 3) * a * a)
    xf = lambda yy: (lambda A, B, C, E1, E3: 4 * (C / A)**3 * (1 - (C / A)**3))(*cub(yy))
    dx = diff(xf, y); d3 = abs(dx + 4 * pi * 3 * sqrt(3) * e1**3 * e3**3 / a * sqrt(x * (1 - x)))
    print("  y = %s: h = %s, x = %s; |a^3 - b^3 - c^3| = %s, |U0 - a^2| = %s, |U1 + 2 pi y a^2| = %s, |U2 - (2 pi^2 y^2 + 2 pi^2/3) a^2| = %s, dx/dy: %s" % (nstr(y, 3), nstr(hh, 6), nstr(x, 6), nstr(abs(a**3 - b**3 - c**3), 2), nstr(d0, 2), nstr(d1, 2), nstr(d2, 2), nstr(d3, 2)))
    assert max(d0, d1, d2) < TOL and d3 < mpf(10)**-12
a, b, c, e1, e3 = cub(1 / sqrt(mpf(3))); assert abs(c**3 / a**3 - half) < TOL
g0y = lambda y: qval(gt, 1j * y / 2).real
pB = [lambda y: mpf(1), lambda y: -2 * pi * y, lambda y: 2 * pi**2 * (y * y + mpf(1) / 3)]
r3 = 1 / sqrt(mpf(3))
for j in range(3):
    v = 12 * sqrt(3) * quad(lambda y: pB[j](y) * g0y(y), [r3, 1, 2, 4, 8, 16, 40])
    print("  F_%d(1) = %s ; 12 sqrt3 int_{1/sqrt3}^oo p_%d(y) g0(iy) dy = %s" % (j, nstr(F[j], 20), j, nstr(v, 20))); assert abs(v - F[j]) < TOL
rr0 = quad(lambda y: g0y(y) * (1 + 3 * y * y), [r3, 1, 2, 4, 8, 16, 40]); rr1 = 2 * quad(lambda y: g0y(y) * y, [r3, 1, 2, 4, 8, 16, 40])
Lg1 = pi * rr0; Lg2 = pi**2 * rr1
print("  L(g~,1) = %s, L(g~,2) = %s" % (nstr(Lg1, 20), nstr(Lg2, 20)))
print("    diffs:", nstr(F[1] + 12 * sqrt(3) / pi * Lg2, 3), nstr(F[2] - 8 * sqrt(3) * pi * Lg1, 3), nstr(F[3] - 5 * pi**2 / 3 * F[1], 3))
assert abs(F[1] + 12 * sqrt(3) / pi * Lg2) < TOL and abs(F[2] - 8 * sqrt(3) * pi * Lg1) < TOL and abs(F[3] - 5 * pi**2 / 3 * F[1]) < mpf(10)**-10
print("  F_1 = -(12 sqrt3/pi) L(g~,2), F_2 = 8 sqrt3 pi L(g~,1), F_3 = (5 pi^2/3) F_1 : OK")
Lf = [Lfe(fz, 36, 1, s) for s in (1, 2, 3)]; Lfb = [Lfe(fz, 36, 1, s, split=1.25) for s in (1, 2, 3)]
assert max(abs(a_ - b_) for a_, b_ in zip(Lf, Lfb)) < mpf(10)**-15
assert abs(Lf[0] - 4 * sqrt(3) / pi * Lg2) < TOL and abs(Lf[1] - 4 * sqrt(3) * pi / 9 * Lg1) < TOL and abs(Lf[2] - 2 * sqrt(3) * pi / 9 * Lg2) < TOL
print("  L(f,1) = (4 sqrt3/pi) L(g~,2), L(f,2) = (4 sqrt3 pi/9) L(g~,1), L(f,3) = (2 sqrt3 pi/9) L(g~,2): OK   (L(f,k) = %s)" % [nstr(v, 15) for v in Lf])
rat = [Lf[i] / F[i + 1] for i in range(3)]
print("  L(f,k)/F_k(1) = %s   (expected -1/3, 1/18, -1/90)" % [nstr(v, 18) for v in rat])
assert abs(rat[0] + mpf(1) / 3) < TOL and abs(rat[1] - mpf(1) / 18) < TOL and abs(rat[2] + mpf(1) / 90) < mpf(10)**-10
print("PASS")
