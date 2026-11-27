"""Checks for Theorems C, D of ../THEOREM.md ("only if" direction, degree (1,2)):  a_n = A n + k, b_n = B n^2 + C n + E.
D := A^2 + 4B,  omega := 2kB - A(B + C).
(a) algebra (sympy): with N = x + k z + B(1-x')z^2, Delta = 1 - A z - B z^2 = (1-l1 z)(1-l2 z), x + x' = -C/B:
       N/Delta = const + al1/(1-l1 z) + al2/(1-l2 z),   al1 - al2 = omega/(B (l1 - l2));
    and if A^2 + 4B = 0 (l = A/2):  N(1/l) = omega/(2 B l).
(b) Theorem C by brute force mod p: for random parameters and every prime p < 400 with p not dividing 2 B omega and
    D not a nonzero square mod p, there is NO n <= 3000 with p | GCD_n and p | GCD_{n+1}.
    (For comparison: number of primes p < 400 with D a nonzero square mod p for which such n exists.)
(c) the degenerate examples and counterexamples of section 5 (exact).
(d) Corollary: log min(GCD_n, GCD_{n+1}) - (1/2) log n! stays below C n in examples."""
import sympy as sp, random
from math import gcd, lgamma, log
from fractions import Fraction as Fr
# (a)
x, xp, k, A, B, z, l1, l2 = sp.symbols('x xp k A B z l1 l2')
N = x + k * z + B * (1 - xp) * z**2
# Delta = (1 - l1 z)(1 - l2 z): A = l1 + l2, B = -l1 l2
subs = {A: l1 + l2, B: -l1 * l2}
al1 = sp.simplify((N / (1 - l2 * z)).subs(z, 1 / l1).subs(subs))          # coefficient of 1/(1 - l1 z)
al2 = sp.simplify((N / (1 - l1 * z)).subs(z, 1 / l2).subs(subs))
omega_over_B = 2 * k - A * (1 - x - xp)                                     # omega/B, since (B + C)/B = 1 - x - x'
assert sp.simplify(al1 - al2 - (omega_over_B / (l1 - l2)).subs(subs)) == 0
lam = sp.symbols('lam')
assert sp.simplify(N.subs(z, 1 / lam).subs(B, -lam**2) - (2 * k - 2 * lam * (1 - x - xp)) / (2 * lam)) == 0
print("(a) PASS: al1 - al2 = omega/(B(l1 - l2));  double root: N(1/l) = omega/(2 B l)")
# (b)
def primes(M):
    s = [True] * M; out = []
    for i in range(2, M):
        if s[i]:
            out.append(i)
            for j in range(i * i, M, i): s[j] = False
    return out
PR = primes(400)
def persistent(Ai, ki, Bi, Ci, Ei, p, NMAX=3000):
    p0, p1, q0, q1 = 1 % p, ki % p, 0, 1 % p
    prev = (p1 == 0 and q1 == 0)
    for n in range(1, NMAX + 1):
        a = (Ai * n + ki) % p; b = ((Bi * n + Ci) * n + Ei) % p
        p0, p1 = p1, (a * p1 + b * p0) % p; q0, q1 = q1, (a * q1 + b * q0) % p
        cur = (p1 == 0 and q1 == 0)
        if cur and prev: return n
        prev = cur
    return None
random.seed(7)
tested = viol = splitpers = splittot = 0
cases = []
while len(cases) < 300:
    Ai = random.randint(1, 9); ki = random.randint(-9, 9); Bi = random.choice([i for i in range(-9, 10) if i]); Ci = random.randint(-9, 9); Ei = random.randint(-9, 9)
    om = 2 * ki * Bi - Ai * (Bi + Ci)
    if om == 0: continue
    cases.append((Ai, ki, Bi, Ci, Ei, om))
cases += [(2, 1, -1, 0, 3, None), (4, 3, -4, 1, 2, None), (6, 1, -9, 2, 5, None)]          # double root D = 0
for (Ai, ki, Bi, Ci, Ei, om) in cases:
    om = 2 * ki * Bi - Ai * (Bi + Ci); Dd = Ai * Ai + 4 * Bi
    if om == 0: continue
    for p in PR:
        if (2 * Bi * om) % p == 0: continue
        nonzero_square = (Dd % p != 0) and pow(Dd % p, (p - 1) // 2, p) == 1
        r = persistent(Ai, ki, Bi, Ci, Ei, p)
        if not nonzero_square:
            tested += 1
            if r is not None: viol += 1; print("VIOLATION", (Ai, ki, Bi, Ci, Ei), p, r)
        else:
            splittot += 1; splitpers += (r is not None)
print("(b) Theorem C: %d (case, prime) pairs with D not a nonzero square mod p tested, violations: %d" % (tested, viol))
print("    (for D a nonzero square mod p: two consecutive GCDs divisible by p occur in %d of %d pairs)" % (splitpers, splittot))
assert viol == 0
# (c)
def run(a, b, Nn):
    p = [1, a(0)]; q = [0, 1]
    for n in range(1, Nn + 1):
        p.append(a(n) * p[-1] + b(n) * p[-2]); q.append(a(n) * q[-1] + b(n) * q[-2])
    return p[1:], q[1:]
from math import factorial
p, q = run(lambda n: n + 1, lambda n: n * (n - 1), 60)
assert all(p[n] == q[n] and p[n] >= factorial(n + 1) for n in range(61))
print("(c) PCF[n+1, n(n-1)]: p_n = q_n >= (n+1)! for n <= 60  (terminating; FR; A^2+4B = 5, not in (Eq. 6), not in (Eq. 7))")
p, q = run(lambda n: n * n + 1, lambda n: (n - 1)**2 * (n * n - 2), 40)
assert all(p[n] == q[n] and p[n] >= factorial(n)**2 for n in range(41))
print("    PCF[n^2+1, (n-1)^2 (n^2-2)]: p_n = q_n >= n!^2 for n <= 40  (terminating; FR; b not decomposable over Q)")
p, q = run(lambda n: n - 3, lambda n: 2 * (n - 1) * (n - 3), 60)
assert all(p[n] == 0 and q[n] == 0 for n in range(2, 61))
print("    PCF[n-3, 2(n-1)(n-3)]  (family (Eq. 6): B = 2, m = 1, k = -3): p_n = q_n = 0 for 2 <= n <= 60")
p, q = run(lambda n: n - 1, lambda n: -(n - 1) * (n - 2), 60)
assert all(p[n] == 0 and q[n] == 0 for n in range(1, 61))
print("    PCF[n-1, -(n-1)(n-2)]  (family (Eq. 7): k = 1/2, x = (1,2)):       p_n = q_n = 0 for 1 <= n <= 60")
# (d)
for (Ai, ki, Bi, Ci, Ei) in ((1, 1, 1, 0, 2), (3, 1, 1, 1, 0), (2, 1, -1, 0, 3), (1, 0, 3, 1, 1), (2, 3, 2, 2, 1)):
    p, q = run(lambda n: Ai * n + ki, lambda n: (Bi * n + Ci) * n + Ei, 2001)
    row = []
    for n in (250, 500, 1000, 2000):
        g = min(gcd(p[n], q[n]), gcd(p[n + 1], q[n + 1]))
        row.append("%.3f" % ((log(g) - 0.5 * lgamma(n + 1)) / n))
    print("(d) (A,k;B,C,E) = %s, D = %d: (log min(GCD_n,GCD_{n+1}) - (1/2) log n!)/n at n = 250..2000: %s" % ((Ai, ki, Bi, Ci, Ei), Ai * Ai + 4 * Bi, row))
print("PASS")
