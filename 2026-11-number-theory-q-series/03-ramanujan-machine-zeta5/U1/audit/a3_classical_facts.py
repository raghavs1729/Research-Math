"""AUDIT a3. High-precision numerical test of the classical facts (C1)-(C7) exactly as stated in THEOREM.md
(q = e^{i pi tau}, tau = i t), of Lemma B, and of the parametrisation in section 3.  Theta functions are summed
directly from their definitions in THEOREM.md section 1 (mpmath.jtheta is used only as a cross-check)."""
from mpmath import mp, mpf, exp, pi, hyp2f1, diff, zeta, gamma, sqrt, jtheta, nstr, quad, inf, ellipk
mp.dps = 60
TOL = mpf(10) ** -50
def th2(t, c=1):   # c rescales tau: theta_2(c*i*t)
    q = exp(-pi * t * c); return 2 * sum(q ** ((n + mpf(1) / 2) ** 2) for n in range(0, 60))
def th3(t, c=1):
    q = exp(-pi * t * c); return 1 + 2 * sum(q ** (n * n) for n in range(1, 60))
def th4(t, c=1):
    q = exp(-pi * t * c); return 1 + 2 * sum((-1) ** n * q ** (n * n) for n in range(1, 60))
lam = lambda t: th2(t) ** 4 / th3(t) ** 4
sig5 = lambda n: sum(d ** 5 for d in range(1, n + 1) if n % d == 0)
def E6(q):   # E_6(q) = 1 - 504 sum sigma_5(n) q^n
    return 1 - 504 * sum(sig5(n) * q ** n for n in range(1, 400))
ts = [mpf('0.6'), mpf('0.8'), mpf(1), mpf('1.3'), mpf(2), mpf('3.1')]
def chk(name, vals):
    m = max(abs(v) for v in vals); print("%-62s max dev %s" % (name, nstr(m, 3))); assert m < TOL, name
# cross-check own theta sums against mpmath (nome convention q = e^{-pi t})
chk("own theta_j vs mpmath.jtheta(j,0,e^{-pi t})", [th2(t) - jtheta(2, 0, exp(-pi * t)) for t in ts] +
    [th3(t) - jtheta(3, 0, exp(-pi * t)) for t in ts] + [th4(t) - jtheta(4, 0, exp(-pi * t)) for t in ts])
chk("(C1) th3^4 = th2^4 + th4^4", [th3(t) ** 4 - th2(t) ** 4 - th4(t) ** 4 for t in ts])
chk("(C2) lambda(i) = 1/2", [lam(mpf(1)) - mpf(1) / 2])
grid = [mpf(k) / 20 for k in range(4, 120)]          # t in [0.2, 6): larger t would need more than 60 terms
assert all(lam(grid[i]) > lam(grid[i + 1]) for i in range(len(grid) - 1)) and 0 < lam(grid[-1]) < lam(grid[0]) < 1
print("(C2) lambda(it) decreasing on a grid t in [0.2,6), values in (0,1): ok   lam(0.2)=%s lam(5.95)=%s" % (nstr(lam(grid[0]), 8), nstr(lam(grid[-1]), 8)))
chk("(C2) 2F1(1/2,1/2;1;lambda(it)) = th3(it)^2  (all t tested, also t<1)", [hyp2f1(.5, .5, 1, lam(t)) - th3(t) ** 2 for t in ts])
chk("(C2) K(k) = (pi/2) th3^2, k^2 = lambda (mpmath ellipk(m=lambda))", [ellipk(lam(t)) - pi / 2 * th3(t) ** 2 for t in ts])
# exact termwise derivatives of the theta series in t
def dth2(t):
    q = exp(-pi * t); return 2 * sum(-pi * (n + mpf(1) / 2) ** 2 * q ** ((n + mpf(1) / 2) ** 2) for n in range(0, 60))
def dth3(t):
    q = exp(-pi * t); return 2 * sum(-pi * n * n * q ** (n * n) for n in range(1, 60))
dlam = lambda t: 4 * th2(t) ** 3 * dth2(t) / th3(t) ** 4 - 4 * th2(t) ** 4 * dth3(t) / th3(t) ** 5
chk("(C3) d lambda(it)/dt = -pi lambda(1-lambda) th3^4 (termwise derivative)", [dlam(t) + pi * lam(t) * (1 - lam(t)) * th3(t) ** 4 for t in ts])
chk("(C4) th3(i/t)^2 = t th3(it)^2", [th3(1 / t) ** 2 - t * th3(t) ** 2 for t in ts])
chk("(C4) th2(i/t)^2 = t th4(it)^2", [th2(1 / t) ** 2 - t * th4(t) ** 2 for t in ts])
chk("(C4) th4(i/t)^2 = t th2(it)^2", [th4(1 / t) ** 2 - t * th2(t) ** 2 for t in ts])
A = lambda t: th3(t) ** 4; B = lambda t: th4(t) ** 4; C = lambda t: th2(t) ** 4
chk("(C5) E6(q^2) = (1/2)(C+A)(A+B)(B-C)", [E6(exp(-2 * pi * t)) - (C(t) + A(t)) * (A(t) + B(t)) * (B(t) - C(t)) / 2 for t in ts])
chk("(C5) th3(2tau)^2 = (th3^2+th4^2)/2", [th3(t, 2) ** 2 - (th3(t) ** 2 + th4(t) ** 2) / 2 for t in ts])
chk("(C5) th4(2tau)^2 = th3 th4", [th4(t, 2) ** 2 - th3(t) * th4(t) for t in ts])
chk("(C5) th2(2tau)^2 = (th3^2-th4^2)/2", [th2(t, 2) ** 2 - (th3(t) ** 2 - th4(t) ** 2) / 2 for t in ts])
h_ = mpf(1) / 2
chk("(C5) th3(tau/2)^2 = th3^2+th2^2", [th3(t, h_) ** 2 - th3(t) ** 2 - th2(t) ** 2 for t in ts])
chk("(C5) th4(tau/2)^2 = th3^2-th2^2", [th4(t, h_) ** 2 - th3(t) ** 2 + th2(t) ** 2 for t in ts])
chk("(C5) th2(tau/2)^2 = 2 th2 th3", [th2(t, h_) ** 2 - 2 * th2(t) * th3(t) for t in ts])
zs = [mpf(0), mpf('0.01'), mpf('0.2'), mpf('0.4'), mpf('0.5')]
chk("(C6) 2F1(1/4,1/4;1;4z(1-z)) = 2F1(1/2,1/2;1;z), 0<=z<=1/2", [hyp2f1(.25, .25, 1, 4 * z * (1 - z)) - hyp2f1(.5, .5, 1, z) for z in zs])
z = mpf('0.7'); print("     (C6) is FALSE for z > 1/2 as expected (range restriction is needed): dev at z=0.7 = %s" % nstr(hyp2f1(.25, .25, 1, 4 * z * (1 - z)) - hyp2f1(.5, .5, 1, z), 5))
part = sum(mpf(sig5(n)) / mpf(n) ** 12 for n in range(1, 3000)); print("(C7) sum_{n<3000} sigma_5(n) n^-12 - zeta(12)zeta(7) = %s" % nstr(part - zeta(12) * zeta(7), 5))
assert abs(part - zeta(12) * zeta(7)) < mpf(10) ** -15
chk("(C7) zeta(0) = -1/2, zeta(-1) = -1/12, zeta(4) = pi^4/90", [zeta(0) + mpf(1) / 2, zeta(-1) + mpf(1) / 12, zeta(4) - pi ** 4 / 90])

# ---- Lemma B
F = lambda t: C(t) * B(t) * (B(t) - C(t))
an = lambda n: 16 * (sig5(n) - (65 * sig5(n // 2) if n % 2 == 0 else 0) + (64 * sig5(n // 4) if n % 4 == 0 else 0))
AN = [0] + [an(n) for n in range(1, 400)]
assert AN[1:4] == [16, -512, 3904]
Fser = lambda t: sum(AN[n] * exp(-pi * n * t) for n in range(1, 400))
chk("Lemma B: F = sum a_n q^n with a_n = 16(s5(n)-65 s5(n/2)+64 s5(n/4))", [F(t) - Fser(t) for t in ts])
chk("Lemma B: 63 F = -2E6(q)+130E6(q^2)-128E6(q^4)", [63 * F(t) + 2 * E6(exp(-pi * t)) - 130 * E6(exp(-2 * pi * t)) + 128 * E6(exp(-4 * pi * t)) for t in ts])
chk("Lemma B: F(i/t) = -t^6 F(it)", [F(1 / t) + t ** 6 * F(t) for t in ts])
chk("Lemma B (proof): E6(q^4) = (A+B)(36AB-(A+B)^2)/64", [E6(exp(-4 * pi * t)) - (A(t) + B(t)) * (36 * A(t) * B(t) - (A(t) + B(t)) ** 2) / 64 for t in ts])
chk("Lemma B (proof): E6(q)   = (A+C)((A+C)^2-36AC)", [E6(exp(-pi * t)) - (A(t) + C(t)) * ((A(t) + C(t)) ** 2 - 36 * A(t) * C(t)) for t in ts])
assert max(abs(AN[n]) - 2080 * sig5(n) for n in range(1, 400)) <= 0
print("|a_n| <= 2080 sigma_5(n) for n < 400 (trivially true in general: 16(1+65+64) = 2080): ok")
# exact polynomial identity, by hand-expanded integer arithmetic on random integers (no CAS)
import random
for _ in range(50):
    b, c = random.randint(-10**6, 10**6), random.randint(-10**6, 10**6); a = b + c
    e2 = (c + a) * (a + b) * (b - c) * 32          # 64*E6(tau)
    e4 = (a + b) * (36 * a * b - (a + b) ** 2)      # 64*E6(2tau)
    e1 = 64 * (a + c) * ((a + c) ** 2 - 36 * a * c) # 64*E6(tau/2)
    assert 64 * 63 * c * b * (b - c) == -2 * e1 + 130 * e2 - 128 * e4
print("polynomial identity 63 CB(B-C) = -2E6(tau/2)+130E6(tau)-128E6(2tau) on A=B+C (50 random integer points, degree 3): ok")

# ---- section 3 parametrisation: x = 4 lam (1-lam), sqrt(1-x) = 1-2 lam (t>1), u(x) = th3^2, dx/dt
for t in (mpf('1.05'), mpf('1.3'), mpf(2), mpf('3.1')):
    l = lam(t); x = 4 * l * (1 - l)
    assert abs(sqrt(1 - x) - (1 - 2 * l)) < TOL and l < mpf(1) / 2
    assert abs(hyp2f1(.25, .25, 1, x) - th3(t) ** 2) < TOL
for t in (mpf('1.05'), mpf('1.3'), mpf(2)):
    x = 4 * lam(t) * (1 - lam(t)); d = 4 * (1 - 2 * lam(t)) * dlam(t)
    assert abs(d + pi * x * sqrt(1 - x) * th3(t) ** 4) < TOL
t = mpf('0.8'); l = lam(t); x = 4 * l * (1 - l)
print("section 3: u(x)=th3(it)^2, sqrt(1-x)=1-2lam, dx/dt=-pi x sqrt(1-x) u^2 for t>1: ok;  (for t<1 they fail, e.g. t=0.8: u(x)-th3^2 = %s)" % nstr(hyp2f1(.25, .25, 1, x) - th3(t) ** 2, 5))
# q = lam/16 + O(lam^2):
t = mpf(12); print("nome: q/(lam/16) at t=12: %s ;  x/(64 q): %s" % (nstr(exp(-pi * t) / (lam(t) / 16), 15), nstr(4 * lam(t) * (1 - lam(t)) / (64 * exp(-pi * t)), 15)))
# g_t(1) = -2 g(1)
print("g_t(1)/g(1) =", nstr(8 * dth3(mpf(1)) / th3(mpf(1)), 40), " (claim: -2)")
assert abs(8 * dth3(mpf(1)) / th3(mpf(1)) + 2) < TOL
print("a3 PASS")
