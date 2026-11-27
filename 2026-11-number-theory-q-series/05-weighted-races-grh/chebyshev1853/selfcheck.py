#!/usr/bin/env python3
"""
selfcheck.py -- numerical checks (mpmath quadrature, 25 digits) of the elementary integrals and inequalities
used by hand in CHEBYSHEV1853.md.  Prints "ok"/"FAIL" per item.
"""
import mpmath as mp
mp.mp.dps = 25
ok_all = True


def chk(name, a, b, tol=mp.mpf(10) ** -15):
    global ok_all
    good = abs(a - b) <= tol * max(1, abs(b))
    ok_all = ok_all and good
    print("%-88s %s   (%s vs %s)" % (name, "ok" if good else "FAIL", mp.nstr(a, 14), mp.nstr(b, 14)))


def chk_true(name, cond):
    global ok_all
    ok_all = ok_all and bool(cond)
    print("%-88s %s" % (name, "ok" if cond else "FAIL"))


# Theorem 4
chk("c_a: int_0^1 (1-(1+v)^(-1/2)) dv/v = 2 log((1+sqrt2)/2)",
    mp.quad(lambda v: (1 - (1 + v) ** -0.5) / v, [0, 1]), 2 * mp.log((1 + mp.sqrt(2)) / 2))
A = lambda V: -mp.log((mp.sqrt(1 + V) + 1) / (mp.sqrt(1 + V) - 1))
chk("int_1^64 dv/(v sqrt(1+v)) = A(64) - A(1)", mp.quad(lambda v: 1 / (v * mp.sqrt(1 + v)), [1, 4, 16, 64]), A(64) - A(1))
chk("A(1) = -2 log(1+sqrt2)", A(1), -2 * mp.log(1 + mp.sqrt(2)))
chk("A(inf) - A(1) - c_a = 2 log 2", -A(1) - 2 * mp.log((1 + mp.sqrt(2)) / 2), 2 * mp.log(2))
x = mp.mpf(1000)
chk("int_0^1 dv/(v log^2(x/v)) = 1/log x  (x = 1000)", mp.quad(lambda a: 1 / (mp.log(x) + a) ** 2, [0, mp.inf]), 1 / mp.log(x))
lam = mp.mpf('0.01')
chk("int_0^1 u lam^(u-1) du = (1 - lam + lam log lam)/(lam log^2 lam)", mp.quad(lambda u: u * lam ** (u - 1), [0, 1]),
    (1 - lam + lam * mp.log(lam)) / (lam * mp.log(lam) ** 2))
chk_true("1/Gamma(u) >= u on [0,1] (grid)", all(mp.rgamma(mp.mpf(i) / 200) >= mp.mpf(i) / 200 for i in range(0, 201)))
# section 2.5
y = mp.mpf(50)
for k in (2, 3, 5):
    g = lambda u: mp.exp(-u ** k / y)
    dg = lambda u: -k * u ** (k - 1) / y * mp.exp(-u ** k / y)
    chk("int_0^oo g_k = Gamma(1+1/k) y^(1/k), k=%d" % k, mp.quad(g, [0, 1, 3, 10, mp.inf]), mp.gamma(1 + mp.mpf(1) / k) * y ** (mp.mpf(1) / k))
    chk("int_0^oo sqrt(u)(-g_k') = Gamma(1+1/(2k)) y^(1/(2k)), k=%d" % k, mp.quad(lambda u: mp.sqrt(u) * -dg(u), [0, 1, 3, 10, mp.inf]),
        mp.gamma(1 + mp.mpf(1) / (2 * k)) * y ** (mp.mpf(1) / (2 * k)))
F = lambda v: v * mp.exp(-v * v) + mp.sqrt(mp.pi) / 2 * mp.erfc(v)
a, b = mp.mpf(3), mp.mpf(11)
chk("int_a^b u(-g_2') du = sqrt(y) [F(a/sqrt y) - F(b/sqrt y)]", mp.quad(lambda u: u * 2 * u / y * mp.exp(-u * u / y), [a, b]),
    mp.sqrt(y) * (F(a / mp.sqrt(y)) - F(b / mp.sqrt(y))))
chk("F'(v) = -2 v^2 e^{-v^2} (v = 0.7)", mp.diff(F, mp.mpf('0.7')), -2 * mp.mpf('0.49') * mp.exp(-mp.mpf('0.49')), tol=mp.mpf(10) ** -12)
U = mp.mpf(4)
chk_true("int_U^oo e^{-u^3/y} du <= y e^{-U^3/y}/(3 U^2)", mp.quad(lambda u: mp.exp(-u ** 3 / y), [U, 2 * U, mp.inf]) <= y * mp.exp(-U ** 3 / y) / (3 * U ** 2))
T = mp.mpf(30)
chk("int_T^oo t^2 e^{-t/y} dt = y e^{-T/y}(T^2 + 2Ty + 2y^2)", mp.quad(lambda t: t * t * mp.exp(-t / y), [T, 10 * T, mp.inf]),
    y * mp.exp(-T / y) * (T * T + 2 * T * y + 2 * y * y))
# Lemma 3
chk("sum 1/((2n-1)2n(2n+1)) = log 2 - 1/2", mp.nsum(lambda n: 1 / ((2 * n - 1) * 2 * n * (2 * n + 1)), [1, mp.inf]), mp.log(2) - mp.mpf(1) / 2)
chk("y^-2 int_1^oo t e^{-t/y} dt = (1+1/y) e^{-1/y}", mp.quad(lambda t: t * mp.exp(-t / y), [1, 10 * y, mp.inf]) / y ** 2, (1 + 1 / y) * mp.exp(-1 / y))
chk("y^-2 int_1^oo log t e^{-t/y} dt = E_1(1/y)/y", mp.quad(lambda t: mp.log(t) * mp.exp(-t / y), [1, 10 * y, mp.inf]) / y ** 2, mp.e1(1 / y) / y)
chk_true("0 < E_1(u) < e^{-u} log(1+1/u) (u = 1e-3 .. 10)", all(0 < mp.e1(u) < mp.exp(-u) * mp.log(1 + 1 / u) for u in [mp.mpf(10) ** (mp.mpf(j) / 10) for j in range(-30, 11)]))
chk_true("(1+u) e^{-u} >= 1 - u^2/2 (u in [0,1])", all((1 + u) * mp.exp(-u) >= 1 - u * u / 2 for u in [mp.mpf(j) / 100 for j in range(0, 101)]))
chk("Gamma(rho+2)/(rho(rho+1)) = Gamma(rho), rho = 1/2+6i", mp.gamma(mp.mpc(2.5, 6)) / (mp.mpc(0.5, 6) * mp.mpc(1.5, 6)), mp.gamma(mp.mpc(0.5, 6)))
chk("|Gamma(1/2+it)|^2 = pi/cosh(pi t), t = 6.02", abs(mp.gamma(mp.mpc(0.5, 6.02))) ** 2, mp.pi / mp.cosh(mp.pi * mp.mpf('6.02')))
h = lambda t: (mp.mpf(1) / 4 + t * t) * mp.sqrt(mp.pi / mp.cosh(mp.pi * t))
chk_true("(1/4+t^2) sqrt(pi/cosh(pi t)) decreasing on [2, 60] (grid)", all(h(mp.mpf(j) / 10) > h(mp.mpf(j + 1) / 10) for j in range(20, 600)))
# section 2.6
chk_true("e^{-v} <= 1/(e v)  (equality at v = 1)", all(mp.exp(-v) <= (1 + mp.mpf(10) ** -20) / (mp.e * v) for v in [mp.mpf(j) / 10 for j in range(1, 200)]))
# kernel F''
xx = mp.mpf(7)
Fk = lambda t: mp.exp(-t / xx) / mp.log(t)
t0 = mp.mpf(5)
chk("F''(t) formula (x = 7, t = 5)", mp.diff(Fk, t0, 2),
    mp.exp(-t0 / xx) * (1 / (xx ** 2 * mp.log(t0)) + 2 / (xx * t0 * mp.log(t0) ** 2) + 1 / (t0 ** 2 * mp.log(t0) ** 2) + 2 / (t0 ** 2 * mp.log(t0) ** 3)),
    tol=mp.mpf(10) ** -12)
# Mellin
p, s = mp.mpf(7), mp.mpf('1.7')
chk("int_0^oo e^{-p/x} x^{-s-1} dx = Gamma(s) p^{-s}", mp.quad(lambda t: mp.exp(-p * t) * t ** (s - 1), [0, 1, mp.inf]), mp.gamma(s) * p ** (-s))
# Proposition 8 (section 2.11): the three upper bounds for mu and the elementary constants
g0 = mp.euler
mu = lambda lam: mp.quad(lambda u: lam ** (u - 1) * mp.rgamma(u), [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128])
chk_true("(mu1) mu(lam) <= 1/(lam (log(1/lam) - gamma0)^2) + 3.1 for lam = 1e-8 .. 1/4",
         all(mu(l) <= 1 / (l * (mp.log(1 / l) - g0) ** 2) + mp.mpf('3.1') for l in [mp.mpf(10) ** -8, mp.mpf(10) ** -4, mp.mpf('0.01'), mp.mpf('0.1'), mp.mpf('0.25')]))
chk_true("(mu2) mu(lam) <= 8 on [1/4, 1] (grid)", all(mu(mp.mpf(j) / 20) <= 8 for j in range(5, 21)))
chk_true("(mu3) mu(lam) <= 0.75 + 1.13 lam e^lam for lam = 1, 2, 5, 10, 20",
         all(mu(mp.mpf(l)) <= mp.mpf('0.75') + mp.mpf('1.13') * l * mp.exp(l) for l in (1, 2, 5, 10, 20)))
chk_true("lam log^2(1/lam) mu(lam) -> 1: value at 1e-30 within [1, 1.02]", 1 <= mp.mpf(10) ** -30 * mp.log(mp.mpf(10) ** 30) ** 2 * mu(mp.mpf(10) ** -30) <= mp.mpf('1.02'))
chk_true("int_0^1 u e^{gamma0 u} du <= 0.75", mp.quad(lambda u: u * mp.exp(g0 * u), [0, 1]) <= mp.mpf('0.75'))
chk_true("1/Gamma(u) <= u e^{gamma0 u} on [0,1] (grid)", all(mp.rgamma(mp.mpf(j) / 200) <= mp.mpf(j) / 200 * mp.exp(g0 * mp.mpf(j) / 200) * (1 + mp.mpf(10) ** -20) for j in range(0, 201)))
chk_true("Gamma >= 0.8856 on [1,2] (grid)", all(mp.gamma(1 + mp.mpf(j) / 1000) >= mp.mpf('0.8856') for j in range(0, 1001)))
chk_true("-A(eta) <= log(6/eta) for eta in (0,1]", all(-A(mp.mpf(j) / 100) <= mp.log(600 / mp.mpf(j)) for j in range(1, 101)))
chk_true("1/(log 4 - gamma0) <= 1.24", 1 / (mp.log(4) - g0) <= mp.mpf('1.24'))
print("ALL OK:", ok_all)
