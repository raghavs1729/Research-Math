"""V3. Numerical confirmation for Theorem U1 (not part of the proof).
(a) the continued fraction, evaluated independently of the proof (asymptotic matching, code/asym_match.py),
    agrees with 96/(372 zeta(5) - pi^5) to 300 digits;
(b) direct evaluation of the convergents p_n/q_n for n up to 20000 (plain recurrence) approaches the same value
    at the predicted rate ~ c/n;
(c) F(i/t) = -t^6 F(it); the Mellin moments M_s = int_0^oo t^(s-1) F(it) dt = int_1^oo (t^(s-1) - t^(5-s)) F(it) dt:
    M_4 = 1/4, M_5 = 186 zeta(5)/pi^5; and (pi^5/48) int_1^oo (t-1)^3 (t+1) F(it) dt = (372 zeta(5) - pi^5)/96;
(d) Phi(i) - (i/2) Phi'(i) computed from the q-series equals (372 zeta(5) - pi^5)/96."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # asym_match.py is copied next to this script
from fractions import Fraction as Fr
from mpmath import mp, mpf, zeta, pi, quad, exp, inf, jtheta
import asym_match as am

mp.dps = 320
target = 96 / (372 * zeta(5) - pi**5)
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
L, coef, errs = am.limit(5, 8, am.at5(Fr(-15, 8), Fr(9, 8)), Fr(-3, 2), Fr(-5, 2), a, 170, 220)
print("(a) |PCF - 96/(372 zeta(5) - pi^5)| =", mp.nstr(abs(L - target), 5))
assert abs(L - target) < mpf(10) ** (-295)

mp.dps = 40
p0, p1, q0, q1 = mpf(1), mpf(a(0)), mpf(0), mpf(1)
for n in range(1, 20001):
    p0, p1 = p1, a(n) * p1 - 64 * n**10 * p0
    q0, q1 = q1, a(n) * q1 - 64 * n**10 * q0
    s = abs(q1); p0 /= s; p1 /= s; q0 /= s; q1 /= s
    if n in (10, 100, 1000, 10000, 20000):
        print("(b) n = %5d   p_n/q_n - target = %s   n*(diff) = %s" % (n, mp.nstr(p1 / q1 - 96 / (372 * zeta(5) - pi**5), 8), mp.nstr(n * (p1 / q1 - 96 / (372 * zeta(5) - pi**5)), 8)))

mp.dps = 50
def F(t):
    q = exp(-pi * t)
    t2, t4 = jtheta(2, 0, q)**4, jtheta(4, 0, q)**4
    return t2 * t4 * (t4 - t2)
# symmetry F(i/t) = -t^6 F(it), checked at a few points; then M_s = int_1^oo (t^(s-1) - t^(5-s)) F(it) dt
for t in (mpf(1) / 2, mpf(7) / 10, mpf(9) / 10, mpf(13) / 10):
    assert abs(F(1 / t) + t**6 * F(t)) < mpf(10)**-40
M = {s: quad(lambda t: (t**(s - 1) - t**(5 - s)) * F(t), [1, 2, 4, 8, 16, inf]) for s in (4, 5)}
print("(c) M_4 - 1/4 =", mp.nstr(M[4] - mpf(1) / 4, 5), "  M_5 - 186 zeta(5)/pi^5 =", mp.nstr(M[5] - 186 * zeta(5) / pi**5, 5))
assert abs(M[4] - mpf(1) / 4) < mpf(10)**-40 and abs(M[5] - 186 * zeta(5) / pi**5) < mpf(10)**-40
I = quad(lambda t: (t - 1)**3 * (t + 1) * F(t), [1, 2, 4, 8, 16, inf])
print("    (pi^5/48) int_1^oo (t-1)^3 (t+1) F(it) dt - (372 zeta(5) - pi^5)/96 =", mp.nstr(pi**5 / 48 * I - (372 * zeta(5) - pi**5) / 96, 5))
assert abs(pi**5 / 48 * I - (372 * zeta(5) - pi**5) / 96) < mpf(10)**-40
sig5 = lambda n: sum(d**5 for d in range(1, n + 1) if n % d == 0)
an = lambda n: 16 * (sig5(n) - (65 * sig5(n // 2) if n % 2 == 0 else 0) + (64 * sig5(n // 4) if n % 4 == 0 else 0))
q = exp(-pi)
Phi = sum(mpf(an(n)) * q**n / (2 * mpf(n)**5) for n in range(1, 60))
dPhi = sum(mpf(an(n)) * (1j * pi * n) * q**n / (2 * mpf(n)**5) for n in range(1, 60))      # d/dtau at tau = i
val = Phi - (1j / 2) * dPhi
print("(d) Phi(i) - (i/2)Phi'(i) - (372 zeta(5) - pi^5)/96 =", mp.nstr(abs(val - (372 * zeta(5) - pi**5) / 96), 5))
assert abs(val - (372 * zeta(5) - pi**5) / 96) < mpf(10)**-40
print("V3 PASS")
