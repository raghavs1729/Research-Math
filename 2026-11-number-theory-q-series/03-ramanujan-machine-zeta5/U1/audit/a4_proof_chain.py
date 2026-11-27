"""AUDIT a4. Independent numerical/exact checks of the proof chain of THEOREM.md, sections 3-5.
(1) exact: recurrences, initial values, Proposition 1, Casoratian D_n = 1/(8(n+1)^5), w_n >= n^(-3/2)/16, w_n >= u_n.
(2) Q := sum v_n x^n (v_n from the recurrence, n <= 6000) equals theta3(it)^8 Phi(it) at x = 4 lam (1-lam), t > 1
    (this tests (3.2), analyticity/initial values and 'Q = Q~ on (0,1)' at points far from x = 0).
(3) (4.2)/(4.3): Q'(x)/P'(x) -> l as x -> 1-, both from the power series and from phi + g phi_t/g_t; l from the
    Casoratian series (accelerated) vs (4.3) vs (372 zeta(5) - pi^5)/96.
(4) section 5: integral of (t-1)^3(t+1)e^{-pi n t}; (5.1); M(s) vs the Dirichlet-series expression at s = 4, 5 and at
    non-integer / complex s (tests the analytic continuation, not only the two values); regularity at s = 4, 5.
(5) Remarks 1 and 4(iii)."""
from fractions import Fraction as Fr
from math import factorial
from mpmath import mp, mpf, mpc, exp, pi, zeta, gamma, quad, inf, nstr, ellipk, hyp2f1, sqrt, findroot
mp.dps = 50
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
at = lambda n: Fr(a(n), 8)
# ---------- (1) exact
N = 120
w = [Fr(1)]; v = [Fr(0), Fr(1, 8)]
w.append(at(0) * w[0])
for n in range(1, N):
    w.append((at(n) * w[n] - n**5 * w[n - 1]) / (n + 1)**5)
    v.append((at(n) * v[n] - n**5 * v[n - 1]) / (n + 1)**5)
u = [Fr(1)]
for k in range(N): u.append(u[-1] * Fr((4 * k + 1)**2, (4 * k + 4)**2))
def mul(A, B): return [sum(A[i] * B[n - i] for i in range(n + 1)) for n in range(N + 1)]
u2 = mul(u, u); u4 = mul(u2, u2)
assert u4 == w, "w_n = [x^n] 2F1(1/4,1/4;1;x)^4"
p = [1, a(0)]; q = [0, 1]
for n in range(1, N):
    p.append(a(n) * p[-1] - 64 * n**10 * p[-2]); q.append(a(n) * q[-1] - 64 * n**10 * q[-2])
for n in range(N + 1):
    assert p[n] == 8**n * factorial(n)**5 * w[n] and q[n] == 8**n * factorial(n)**5 * v[n]
for n in range(N):
    assert w[n] * v[n + 1] - w[n + 1] * v[n] == Fr(1, 8 * (n + 1)**5)
for n in range(1, N + 1):
    assert w[n] >= u[n] and u[n]**2 * 256 * n**3 >= 1 and v[n] > 0      # u_n >= n^(-3/2)/16  <=> 256 n^3 u_n^2 >= 1
print("(1) exact, n <= %d: w = [x^n]u^4; Prop. 1 (p_{n-1}, q_{n-1}); Casoratian = 1/(8(n+1)^5); w_n >= u_n >= n^(-3/2)/16; v_n > 0: PASS" % N)

# ---------- theta side
def th2(t): qq = exp(-pi * t); return 2 * sum(qq ** ((n + mpf(1) / 2) ** 2) for n in range(0, 40))
def th3(t): qq = exp(-pi * t); return 1 + 2 * sum(qq ** (n * n) for n in range(1, 40))
def th4(t): qq = exp(-pi * t); return 1 + 2 * sum((-1) ** n * qq ** (n * n) for n in range(1, 40))
def dth3(t): qq = exp(-pi * t); return 2 * sum(-pi * n * n * qq ** (n * n) for n in range(1, 40))
sig5 = lambda n: sum(d ** 5 for d in range(1, n + 1) if n % d == 0)
NA = 260
AN = [0] + [16 * (sig5(n) - (65 * sig5(n // 2) if n % 2 == 0 else 0) + (64 * sig5(n // 4) if n % 4 == 0 else 0)) for n in range(1, NA)]
F = lambda t: th2(t)**4 * th4(t)**4 * (th4(t)**4 - th2(t)**4)
phi = lambda t: sum(AN[n] * exp(-pi * n * t) / mpf(n)**5 for n in range(1, NA)) / 2
dphi = lambda t: sum(-pi * AN[n] * exp(-pi * n * t) / mpf(n)**4 for n in range(1, NA)) / 2
mp.dps = 140
ell_claim = (372 * zeta(5) - pi**5) / 96
mp.dps = 50

# ---------- (2) power series with coefficients from the recurrence vs modular expression
NN = 6000
mp.dps = 140      # the Richardson step below has weights ~ 1e45
wf = [mpf(1), mpf(1) / 4]; vf = [mpf(0), mpf(1) / 8]
for n in range(1, NN):
    an8 = mpf(a(n)) / 8
    wf.append((an8 * wf[n] - mpf(n)**5 * wf[n - 1]) / mpf(n + 1)**5)
    vf.append((an8 * vf[n] - mpf(n)**5 * vf[n - 1]) / mpf(n + 1)**5)
def ps(c, x, der=0):
    if der == 0: return sum(c[n] * x**n for n in range(len(c)))
    return sum(n * c[n] * x**(n - 1) for n in range(1, len(c)))
for t in (mpf(3), mpf(2), mpf('1.5'), mpf('1.25')):
    l = th2(t)**4 / th3(t)**4; x = 4 * l * (1 - l)
    dP = ps(wf, x) - th3(t)**8; dQ = ps(vf, x) - th3(t)**8 * phi(t)
    print("(2) t=%-5s x=%s   sum w_n x^n - th3^8 = %s    sum v_n x^n - th3^8 Phi = %s   (x^N = %s)"
          % (nstr(t, 4), nstr(x, 12), nstr(dP, 3), nstr(dQ, 3), nstr(x**NN, 3)))
    assert abs(dP) < mpf(10)**-25 and abs(dQ) < mpf(10)**-25
# ---------- (3)
# l from the Casoratian series, with Richardson acceleration of the partial sums (tail ~ c1/n + c2/n^2 + ...)
part = [mpf(0)]
for n in range(0, NN - 1):
    part.append(part[-1] + 1 / (8 * mpf(n + 1)**5 * wf[n] * wf[n + 1]))        # part[m] = v_m/w_m
assert abs(part[NN - 1] - vf[NN - 1] / wf[NN - 1]) < mpf(10)**-40
k = 14; Nt = NN - 1
ell_series = sum((-1)**(k - j) * mpf(Nt - k + j)**k / (factorial(j) * factorial(k - j)) * part[Nt - k + j] for j in range(k + 1))
mp.dps = 50
ell_43 = sum(AN[n] * exp(-pi * n) * (2 + pi * n) / (4 * mpf(n)**5) for n in range(1, NA))
print("(3) l (Casoratian series, Richardson k=14) - claim = %s ;  l (4.3) - claim = %s ; 8*l - (372z5-pi^5)/12 = %s"
      % (nstr(ell_series - ell_claim, 3), nstr(ell_43 - ell_claim, 3), nstr(8 * ell_series - (372 * zeta(5) - pi**5) / 12, 3)))
assert abs(ell_series - ell_claim) < mpf(10)**-30 and abs(ell_43 - ell_claim) < mpf(10)**-45
print("    Q'(x)/P'(x) as x -> 1-  [modular formula phi + g phi_t/g_t  |  power series, N=%d]  minus l:" % NN)
for t in (mpf('1.5'), mpf('1.2'), mpf('1.05'), mpf('1.01'), mpf('1.001'), mpf('1.00001')):
    l = th2(t)**4 / th3(t)**4; x = 4 * l * (1 - l)
    r_mod = phi(t) + th3(t)**8 * dphi(t) / (8 * th3(t)**7 * dth3(t))
    r_ps = ps(vf, x, 1) / ps(wf, x, 1) if x**NN < mpf(10)**-20 else None
    print("      t=%-8s 1-x=%s  modular: %s   series: %s" % (nstr(t, 7), nstr(1 - x, 3), nstr(r_mod - ell_claim, 5),
          "n/a" if r_ps is None else nstr(r_ps - ell_claim, 5) + "  (modular-series = %s)" % nstr(r_mod - r_ps, 3)))
# growth of P': P'(x) ~ C/sqrt(1-x) -> infinity
# ---------- (4)
for n in (1, 2, 5):
    I = quad(lambda t: (t - 1)**3 * (t + 1) * exp(-pi * n * t), [1, 2, 4, inf])
    assert abs(I - 12 * exp(-pi * n) * (2 + pi * n) / (pi * n)**5) < mpf(10)**-45
I = quad(lambda t: (t - 1)**3 * (t + 1) * F(t), [1, 1.5, 2, 3, 5, 9, inf])
print("(4) (pi^5/48) int_1^oo (t-1)^3(t+1)F(it)dt - claim =", nstr(pi**5 / 48 * I - ell_claim, 3))
assert abs(pi**5 / 48 * I - ell_claim) < mpf(10)**-40
Mint = lambda s: quad(lambda t: (t**(s - 1) - t**(5 - s)) * F(t), [1, 1.5, 2, 3, 5, 9, inf])
Mdir = lambda s: 16 * gamma(s) * pi**(-s) * zeta(s) * zeta(s - 5) * (1 - 65 * mpf(2)**(-s) + 64 * mpf(4)**(-s))
for s in (mpf(8), mpf(7), mpf(5), mpf(4), mpf('2.5'), mpf('0.5'), mpf('-2.5'), mpc(3, 2), mpc('5.5', -1)):
    d = Mint(s) - Mdir(s); print("    M(s) [integral over (1,oo)] - 16 Gamma(s) pi^-s zeta(s) zeta(s-5)(1-65 2^-s+64 4^-s) at s=%s: %s" % (nstr(s, 4), nstr(abs(d), 3)))
    assert abs(d) < mpf(10)**-38
# direct Mellin integral over (0,oo) at s = 8 (no use of the functional equation): theta series with 120 terms, lower
# limit 0.02 (the neglected piece is < 0.02^8 * max|F| on (0,0.02), and |F(it)| < 1e-40 there: checked at 0.02)
def Fbig(t):
    qq = exp(-pi * t)
    a2 = 2 * sum(qq ** ((n + mpf(1) / 2) ** 2) for n in range(0, 120)); a4 = 1 + 2 * sum((-1) ** n * qq ** (n * n) for n in range(1, 120))
    return a2**4 * a4**4 * (a4**4 - a2**4)
assert abs(Fbig(mpf('0.02'))) < mpf(10)**-40
s = mpf(8); M0 = quad(lambda t: t**(s - 1) * Fbig(t), [0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, inf])
dser = sum(AN[n] / mpf(n)**20 for n in range(1, NA))
assert abs(M0 - Mdir(s)) < mpf(10)**-38
print("    int_0^oo t^7 F(it) dt - Mdir(8) = %s ;  sum a_n n^-20 - 16 z(20)z(15)(1-65/2^20+64/4^20) = %s"
      % (nstr(M0 - Mdir(s), 3), nstr(dser - 16 * zeta(20) * zeta(15) * (1 - 65 * mpf(2)**-20 + 64 * mpf(4)**-20), 3)))
print("    M(5) - 186 zeta(5)/pi^5 = %s ;  M(4) - 1/4 = %s ;  symmetry M(s) + M(6-s) at s=4.3: %s"
      % (nstr(Mint(mpf(5)) - 186 * zeta(5) / pi**5, 3), nstr(Mint(mpf(4)) - mpf(1) / 4, 3), nstr(Mint(mpf('4.3')) + Mint(mpf('1.7')), 3)))
assert Fr(16 * 24) * Fr(-1, 2) * (1 - Fr(65, 32) + Fr(64, 1024)) == 186
assert Fr(16 * 6) * Fr(1, 90) * Fr(-1, 12) * (1 - Fr(65, 16) + Fr(64, 256)) == Fr(1, 4)
assert Fr(1, 48) * 186 == Fr(372, 96) and Fr(1, 48) * Fr(1, 2) == Fr(1, 96)
print("    rational arithmetic of M(5), M(4), final combination: exact PASS")
# ---------- (5) remarks
r1 = sum((AN[n] // 16) * (2 + pi * n) * exp(-pi * n) / mpf(n)**5 for n in range(1, NA))
print("(5) Remark 1 series - (372 zeta(5) - pi^5)/384 =", nstr(r1 - (372 * zeta(5) - pi**5) / 384, 3))
J = quad(lambda z: ellipk(z)**3 * (2 * ellipk(1 - z) - ellipk(z)) * (1 - 2 * z), [0, mpf(1) / 4, mpf(1) / 2, mpf(3) / 4, 1])
print("    Remark 4(iii): 3/int_0^1 K^3 (2K' - K)(1-2z) dz - 1/l =", nstr(3 / J - 1 / ell_claim, 3))
# convergence constant c = 1/(8 C^2 l^2), C = lim n^(3/2) w_n = 4 pi^(3/2)/(Gamma(3/4)^6 Gamma(1/4)^2)
C = 4 * pi**mpf(1.5) / (gamma(mpf(3) / 4)**6 * gamma(mpf(1) / 4)**2)
print("    n^(3/2) w_n at n=%d: %s  vs C = %s ;  predicted c = 1/(8 C^2 l^2) = %s" % (NN, nstr(mpf(NN)**1.5 * wf[NN], 10), nstr(C, 10), nstr(1 / (8 * C**2 * ell_claim**2), 10)))
print("a4 PASS")
