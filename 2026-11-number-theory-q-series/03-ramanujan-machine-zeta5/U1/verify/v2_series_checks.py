"""V2. Exact (rational-arithmetic) checks for Theorem U1.

(a) w_n := [x^n] 2F1(1/4,1/4;1;x)^4 satisfies (n+1)^5 w_{n+1} + n^5 w_{n-1} = atil(n) w_n, atil = a_n/8, and
    p_{n-1} = 8^n n!^5 w_n  (p_n the numerators of the continued fraction);  also w_n = sum_k A_k A_{n-k},
    A_k = binom(2k,k)^3/64^k (Clausen).
(b) theta-function identity (q = e^{i pi tau}):
       F := theta2^4 theta4^4 (theta4^4 - theta2^4) = (-2 E6(q) + 130 E6(q^2) - 128 E6(q^4))/63,
    i.e. F = sum a_n q^n with a_n = 16 (sigma5(n) - 65 sigma5(n/2) + 64 sigma5(n/4)),  checked to O(q^NQ),
    and the polynomial identity that proves it from Jacobi's and Landen's identities.
(c) Q(x) := theta3^8 * (1/2) sum a_n q^n/n^5, expanded in x = 4 lambda (1-lambda), has coefficients f_q(n) with
    q_{n-1} = 8^n n!^5 f_q(n)  (q_n the denominators of the continued fraction).
All arithmetic is exact (Fractions / integers)."""
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp

N = 40           # order in x
NQ = 300         # order in q for (b)
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
b = lambda n: -64 * n**10

def mul(A, B, n=N):
    C = [Fr(0)] * (n + 1)
    for i, x in enumerate(A[:n + 1]):
        if x:
            for j, y in enumerate(B[:n + 1 - i]):
                C[i + j] += x * y
    return C

# ---- (a)
u = [Fr(1)]
for k in range(N):
    u.append(u[-1] * Fr((4 * k + 1) ** 2, 16 * (k + 1) ** 2))     # (1/4)_k^2/k!^2
u2 = mul(u, u); w = mul(u2, u2)
A = [Fr(comb(2 * k, k) ** 3, 64 ** k) for k in range(N + 1)]
assert u2 == A, "Clausen: 2F1(1/4,1/4;1;x)^2 = 3F2(1/2,1/2,1/2;1,1;x)"
assert w == mul(A, A)
for n in range(0, N):
    lhs = (n + 1) ** 5 * w[n + 1] + (n ** 5 * w[n - 1] if n >= 1 else 0)
    assert lhs == Fr(a(n), 8) * w[n], n
p = [1, a(0)]; qd = [0, 1]                    # p_{-1},p_0 ; q_{-1},q_0
for n in range(1, N):
    p.append(a(n) * p[-1] + b(n) * p[-2]); qd.append(a(n) * qd[-1] + b(n) * qd[-2])
for n in range(0, N + 1):
    assert p[n] == 8 ** n * factorial(n) ** 5 * w[n], n       # p[n] = p_{n-1}
print("(a) PASS: recurrence for w_n (n < %d), p_{n-1} = 8^n n!^5 w_n (n <= %d), Clausen" % (N, N))

# ---- (b) q-series
def qmul(A, B, n=NQ):
    C = [0] * (n + 1)
    for i, x in enumerate(A):
        if x:
            for j in range(0, n + 1 - i):
                if B[j]:
                    C[i + j] += x * B[j]
    return C
def qpow(A, e):
    R = [1] + [0] * NQ
    for _ in range(e):
        R = qmul(R, A)
    return R
th3 = [0] * (NQ + 1); th4 = [0] * (NQ + 1); th2sq = None
k = 0
while k * k <= NQ:
    th3[k * k] += 1 if k == 0 else 2
    th4[k * k] += (1 if k == 0 else 2) * (-1) ** k
    k += 1
# theta2^4 = theta3^4 - theta4^4 (Jacobi); we use it as the definition here and cross-check with the product formula
th3_4 = qpow(th3, 4); th4_4 = qpow(th4, 4)
th2_4 = [x - y for x, y in zip(th3_4, th4_4)]
# cross-check: theta2^4 = 16 q prod (1-q^{2n})^4 (1+q^{2n})^8  (Jacobi triple product form), to O(q^NQ)
prod = [1] + [0] * NQ
for n in range(1, NQ // 2 + 1):
    for _ in range(8):
        for i in range(NQ, 2 * n - 1, -1):
            prod[i] += prod[i - 2 * n]
    for _ in range(4):
        for i in range(NQ, 2 * n - 1, -1):
            prod[i] -= prod[i - 2 * n]
assert th2_4 == [0] + [16 * c for c in prod[:NQ]], "Jacobi: theta3^4 = theta2^4 + theta4^4"
F = qmul(qmul(th2_4, th4_4), [x - y for x, y in zip(th4_4, th2_4)])
sig5 = lambda n: sum(d**5 for d in range(1, n + 1) if n % d == 0)
S5 = [0] + [sig5(n) for n in range(1, NQ + 1)]
acoef = [0] * (NQ + 1)
for n in range(1, NQ + 1):
    acoef[n] = 16 * (S5[n] - (65 * S5[n // 2] if n % 2 == 0 else 0) + (64 * S5[n // 4] if n % 4 == 0 else 0))
assert F == acoef
print("(b) PASS: theta2^4 theta4^4 (theta4^4-theta2^4) = sum a_n q^n, a_n = 16(s5(n) - 65 s5(n/2) + 64 s5(n/4)), n <= %d" % NQ)
# polynomial identity: with A=theta3^4, B=theta4^4, C=theta2^4 (at tau), A = B + C:
A_, B_, C_ = sp.symbols('A B C')
E6_tau = sp.Rational(1, 2) * (C_ + A_) * (A_ + B_) * (B_ - C_)
E6_2tau = sp.Rational(1, 64) * (A_ + B_) * (36 * A_ * B_ - (A_ + B_) ** 2)
E6_half = (A_ + C_) * ((A_ + C_) ** 2 - 36 * A_ * C_)
ident = sp.expand((63 * C_ * B_ * (B_ - C_) - (-2 * E6_half + 130 * E6_tau - 128 * E6_2tau)).subs(A_, B_ + C_))
assert ident == 0
print("    polynomial identity 63 F = -2 E6(tau/2) + 130 E6(tau) - 128 E6(2 tau) (mod A = B + C): PASS")
# cross-check the three theta expressions of E6 against the q-series of E6 (to O(q^60))
M = 60
def E6ser(step):
    r = [Fr(0)] * (M + 1); r[0] = Fr(1)
    for n in range(1, M // step + 1):
        r[n * step] -= 504 * S5[n]
    return r
ser = lambda L: [Fr(x) for x in L[:M + 1]]
A3, B3, C3 = ser(th3_4), ser(th4_4), ser(th2_4)
add = lambda X, Y, s=1: [x + s * y for x, y in zip(X, Y)]
sc = lambda c, X: [c * x for x in X]
m = lambda X, Y: mul(X, Y, M)
e_tau = sc(Fr(1, 2), m(m(add(C3, A3), add(A3, B3)), add(B3, C3, -1)))
assert e_tau == E6ser(2)
AB = add(A3, B3); AC = add(A3, C3)
e_2tau = sc(Fr(1, 64), m(AB, add(sc(36, m(A3, B3)), m(AB, AB), -1)))
assert e_2tau == E6ser(4)
e_half = m(AC, add(m(AC, AC), sc(36, m(A3, C3)), -1))
assert e_half == E6ser(1)
print("    theta expressions of E6(tau/2), E6(tau), E6(2tau) agree with their q-series to O(q^%d): PASS" % M)

# ---- (c) expansion of theta3^8 * Phi in x
# lambda = theta2^4/theta3^4 as series in q (order N+1), x = 4 lambda (1-lambda) ; invert to get q(x)
T = N + 2
def trunc(L): return [Fr(v) for v in L[:T + 1]]
def inv(A):   # 1/A, A[0] != 0
    R = [Fr(1) / A[0]] + [Fr(0)] * T
    for n in range(1, T + 1):
        R[n] = -sum(A[k] * R[n - k] for k in range(1, n + 1)) / A[0]
    return R
mT = lambda X, Y: mul(X, Y, T)
lam = mT(trunc(th2_4), inv(trunc(th3_4)))
xq = sc(4, mT(lam, add([Fr(1)] + [Fr(0)] * T, lam, -1)))           # x as series in q: 64 q + ...
assert xq[0] == 0 and xq[1] == 64
# power series in q -> coefficients in x: we need G(q) = sum g_n x(q)^n ; solve triangularly
def to_x(G):
    G = list(G); out = []
    xp = [Fr(1)] + [Fr(0)] * T
    for n in range(0, N + 1):
        c = G[n] / xp[n] if n > 0 else G[0]
        if n == 0:
            c = G[0]
        out.append(c)
        G = [g - c * y for g, y in zip(G, xp)]
        xp = mT(xp, xq)
    return out
th3_8 = mT(trunc(th3_4), trunc(th3_4))
Phi = [Fr(0)] + [Fr(acoef[n], 2 * n**5) for n in range(1, T + 1)]
Pser = to_x(th3_8)
Qser = to_x(mT(th3_8, Phi))
assert Pser == w[:N + 1], "theta3^8 = sum w_n x^n"
for n in range(0, N + 1):
    assert qd[n] == 8 ** n * factorial(n) ** 5 * Qser[n], n      # qd[n] = q_{n-1}
print("(c) PASS: theta3^8 = sum w_n x^n and theta3^8 * Phi = sum f_q(n) x^n with q_{n-1} = 8^n n!^5 f_q(n), n <= %d" % N)
print("V2 PASS")
