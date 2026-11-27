"""Referee symbolic checks (sympy), own code.
 S1  Prop 2.1: normal-ordered coefficients of H'E, E H~', G E, E G~ (symbols u=q^j, v=Q^b, P, Y)
 S2  (B-iv): Gamma^{(c)}_{d,e} formula reproduces phi_d(R;aq^c) for small d,c (exact, symbolic a,q at numeric R)
 S3  (B-v): tau phi_n = alpha phi_{n-1} + beta phi_n + gamma phi_{n+1}, symbolic in z=q^R, A=q^n
 S4  Prop 4.1: elimination and factorisation of the bracket
 S5  Cor 4.2: W_{N,M} = (xq;q)_{2M} V_{N,M} by iterating Prop 4.1 (N <= 10)
 S6  Prop 5.1: elimination, Fricke factorisation
 S7  Cor 5.2: iterating Prop 5.1 gives the closed form V# (N <= 10)
 S8  Lemma 7.1: V V# (x^2 q^{2L+2};Q)_{2M} (xq)_{2N} = t_{N,M}  (N <= 10)
 S9  (B4)-identity for E and the symmetrisation sum_{N>=M} lam^k t = sum_{n,m} lam_n^{k+1} lam_m^{k+1} bB_n bJ_m
"""
import sympy as sp

q, x, a, z, A, u, v, P, Y = sp.symbols('q x a z A u v P Y')
Q = q**2


def poch(c, base, n):
    r = sp.Integer(1)
    for i in range(n):
        r *= (1 - c * base**i)
    return r


def zero(e):
    return sp.simplify(sp.together(e)) == 0


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f"{name}: {'OK' if ok else 'FAIL'}")


# ---------------- S1 ----------------
# coefficient of S^j u^b times (q)_j (Q;Q)_b ; u = q^j, v = Q^b, with 1/(q)_{j-1} = (1-u)/(q)_j etc.
# H' E : H' = (Y^-1+x)P^-1 (1 - S) + xP(1+xY) - x^2 P Y U
HE = (1 / Y + x) / P * (1 - (1 - u)) + x * P * (1 + x * Y) - x**2 * P * Y * (1 - v)
# E H~' : move S^j U^b past f(P,Y): f(q^{-j}P, Q^b Y)
# H~' = (Y^-1+x)P^-1 - x q^-1 P(1+xY) S + x P(1+xY) - q^-2 P^-1 Y^-1 U
EH = (1 / (v * Y) + x) * u / P \
    - x / q * (P / u * q) * (1 + x * v * Y) * (1 - u) \
    + x * (P / u) * (1 + x * v * Y) \
    - q**-2 * (u / P) * (1 / (v * Y)) * Q * (1 - v)
# note: from S^{j-1}: P -> q^{-(j-1)}P = (q/u)P ; from U^{b-1}: Y -> Q^{b-1}Y = (v/Q)Y
target = u / (P * Y) + x * u / P + x * P + x**2 * v * P * Y
report("S1 H'E coefficient", zero(HE - target))
report("S1 E H~' coefficient", zero(EH - target))
# G
GE = (1 - (1 + q) * (1 - u) + q * (1 - u) * (1 - u / q)) / (P**2 * Y) + 1 / Y - x * (1 + 1 / q) * (1 - u) \
    + x**2 * (1 + P**2) * Y - x**2 * (1 + P**2) * Y * (1 - v)
# G~ = P^-2Y^-1 + Y^-1 - x(1+q^-1)S - x^2 q^-1(1+q^-1)P^2 Y S + x^2 q^-3 P^2 Y S^2 + x^2(1+P^2)Y - q^-2(1+P^-2)Y^-1 U
Pj = P / u          # q^{-j} P
Pj1 = P * q / u     # q^{-(j-1)} P
Pj2 = P * q**2 / u  # q^{-(j-2)} P
EG = (1 / (Pj**2 * v * Y)) + 1 / (v * Y) - x * (1 + 1 / q) * (1 - u) \
    - x**2 / q * (1 + 1 / q) * Pj1**2 * v * Y * (1 - u) \
    + x**2 * q**-3 * Pj2**2 * v * Y * (1 - u) * (1 - u / q) \
    + x**2 * (1 + Pj**2) * v * Y \
    - q**-2 * (1 + 1 / Pj**2) * (Q / (v * Y)) * (1 - v)
targetG = u**2 / (P**2 * Y) + 1 / Y - x * (1 + 1 / q) * (1 - u) + x**2 * v * (1 + P**2) * Y
report("S1 G E coefficient", zero(GE - targetG))
report("S1 E G~ coefficient", zero(EG - targetG))

# W H' = H~' W and W G = G~ W, checked as operators on a generic function at a generic point:
# use symbolic f(R,S) values via sympy Function and explicit w
R, S = sp.symbols('R S', integer=True)
f = sp.Function('f')
w = lambda r, s: x**(r + 2 * s) * q**((r + s)**2 + s**2)
Pq, Yq = q**R, q**(2 * S)
def Hp_on(F):
    return (1 / Yq + x) / Pq * (F(R, S) - F(R - 1, S)) + x * Pq * (1 + x * Yq) * F(R, S) - x**2 * Pq * Yq * F(R, S + 1)
def Ht_on(F):
    return (1 / Yq + x) / Pq * F(R, S) - x / q * Pq * (1 + x * Yq) * F(R - 1, S) + x * Pq * (1 + x * Yq) * F(R, S) \
        - q**-2 / (Pq * Yq) * F(R, S + 1)
WF = lambda r, s: w(r, s) * f(r, s)
lhs = w(R, S) * Hp_on(f)
rhs = Ht_on(WF)
report("S1 W H' = H~' W", sp.simplify(sp.expand(sp.powsimp(sp.expand(lhs - rhs), force=True))) == 0)
def G_on(F):
    return (F(R, S) - (1 + q) * F(R - 1, S) + q * F(R - 2, S)) / (Pq**2 * Yq) + F(R, S) / Yq \
        - x * (1 + 1 / q) * F(R - 1, S) + x**2 * (1 + Pq**2) * Yq * (F(R, S) - F(R, S + 1))
def Gt_on(F):
    return F(R, S) / (Pq**2 * Yq) + F(R, S) / Yq - x * (1 + 1 / q) * F(R - 1, S) \
        - x**2 / q * (1 + 1 / q) * Pq**2 * Yq * F(R - 1, S) + x**2 * q**-3 * Pq**2 * Yq * F(R - 2, S) \
        + x**2 * (1 + Pq**2) * Yq * F(R, S) - q**-2 * (1 + Pq**-2) / Yq * F(R, S + 1)
lhs = w(R, S) * G_on(f)
rhs = Gt_on(WF)
report("S1 W G = G~ W", sp.simplify(sp.expand(sp.powsimp(sp.expand(lhs - rhs), force=True))) == 0)


# ---------------- S2 ----------------
def phi(e, Rv, aa):
    if Rv < e:
        return sp.Integer(0)
    return 1 / (poch(q, q, Rv - e) * poch(aa * q, q, Rv + e))


def qbinom(n, k):
    if k < 0 or k > n:
        return 0
    return sp.cancel(poch(q, q, n) / (poch(q, q, k) * poch(q, q, n - k)))


def Gamma(c, d, e):
    return (1 - a * q**(2 * e)) / (1 - a) * poch(a, q, d + e) / poch(a * q**(c + 1), q, d + e) * (-1)**(e - d) \
        * q**((e - d) * (e - d - 1) // 2) * a**(e - d) * q**((e + d) * (e - d)) * qbinom(c, e - d)


okS2 = True
for c in range(0, 4):
    for d in range(0, 3):
        ssum = sum(Gamma(c, d, e) for e in range(d, d + c + 1))
        okS2 &= zero(ssum - poch(a * q, q, c))
        for Rv in range(0, d + c + 3):
            lhs = phi(d, Rv, a * q**c)
            rhs = sum(Gamma(c, d, e) * phi(e, Rv, a) for e in range(d, d + c + 1))
            okS2 &= zero(lhs - rhs)
report("S2 (B-iv) Gamma formula and coefficient sum (c<=3, d<=2)", okS2)
okS2b = all(zero(phi(d, Rv, a * q**2) - poch(a * q, q, 2) * phi(d + 1, Rv + 1, a)) for d in range(3) for Rv in range(6))
report("S2 phi_d(R;aq^2) = (aq)_2 phi_{d+1}(R+1;a)", okS2b)

# ---------------- S3 ----------------
n = sp.symbols('n')
al = A / q / ((1 - a * A**2 / q) * (1 - a * A**2))
be = -A / q * (1 + a * A**2) / ((1 - a * A**2 / q) * (1 - a * A**2 * q))
ga = a * A**3 / ((1 - a * A**2) * (1 - a * A**2 * q))
# with D = (q)_{R-n+1}(aq)_{R+n+1}: tau phi_n * D = z ; phi_{n-1} D = (1-a z A)(1-a z A q) ; phi_n D = (1-zq/A)(1-a z A q)
# phi_{n+1} D = (1 - z q/A)(1 - z/A)
pm, p0, pp = (1 - a * z * A) * (1 - a * z * A * q), (1 - z * q / A) * (1 - a * z * A * q), (1 - z * q / A) * (1 - z / A)
report("S3 (B-v) tau phi_n", zero(al * pm + be * p0 + ga * pp - z))
report("S3 alpha+beta+gamma = 0", zero(al + be + ga))
# phi_n(R+1) D = 1 ; q^{2R} phi_n(R+1) D = z^2 : solve and check coefficient sums
c1, c2, c3 = sp.symbols('c1 c2 c3')
for target, expect in ((1, 1), (z**2, 0)):
    eqs = sp.Poly(sp.expand(c1 * pm + c2 * p0 + c3 * pp - target), z).all_coeffs()
    sol = sp.solve(eqs, [c1, c2, c3], dict=True)[0]
    report(f"S3 coefficient sum for {target}", zero(sol[c1] + sol[c2] + sol[c3] - expect))
# numeric spot check of tau phi_n at actual R (n=1..3)
okS3 = True
for nn in range(1, 4):
    aln, ben, gan = [sp.simplify(t.subs(A, q**nn)) for t in (al, be, ga)]
    for Rv in range(0, 7):
        lhs = q**Rv * phi(nn, Rv + 1, a)
        rhs = aln * phi(nn - 1, Rv, a) + ben * phi(nn, Rv, a) + gan * phi(nn + 1, Rv, a)
        okS3 &= zero(lhs - rhs)
report("S3 tau phi_n at R=0..6, n=1..3", okS3)

# ---------------- S4 ----------------
eN, eM, mu0, mu1, mu2, v1 = sp.symbols('eN eM mu0 mu1 mu2 v1')
eps = eN * eM
g = eN**2 + eM**2 - 4 * x
sol = sp.solve([(1 + x) * mu1 - eps * mu0, (1 + x) * mu2 - eps * mu1 - x**2 * (1 - 1 / q) * v1,
                mu2 + (1 - x)**2 * mu0 - x**2 * v1 - g * mu0], [mu0, mu1, mu2], dict=True)[0]
pred = x**2 * (1 + x) * (x + 1 / q) * v1 / ((eN**2 - (1 + x)**2) * (eM**2 - (1 + x)**2))
report("S4 elimination", zero(sol[mu0] - pred))
NN, MM = sp.symbols('NN MM')
en = lambda t: 1 / t + x * t  # t = q^n
tN, tM = sp.symbols('tN tM')
report("S4 eps = e_N e_M", zero(tN**-1 * tM**-1 * (1 + x * tN**2) * (1 + x * tM**2) - en(tN) * en(tM)))
report("S4 g = e_N^2+e_M^2-4x", zero(tN**-2 + x**2 * tN**2 + tM**-2 + x**2 * tM**2 - (en(tN)**2 + en(tM)**2 - 4 * x)))
report("S4 e_n^2-(1+x)^2 = q^{-2n}(1-q^{2n})(1-x^2q^{2n})", zero(en(tN)**2 - (1 + x)**2 - tN**-2 * (1 - tN**2) * (1 - x**2 * tN**2)))


# ---------------- S5 ----------------
def Rfac(N, M, xx):
    return xx**2 * (1 + xx) * (xx + 1 / q) * q**(2 * N + 2 * M) * poch(xx * q, q, 2) / (
        (1 - q**(2 * N)) * (1 - q**(2 * M)) * (1 - xx**2 * q**(2 * N)) * (1 - xx**2 * q**(2 * M)))


def W_rec(N, M, xx):
    if M == 0:
        return sp.Integer(1)
    return Rfac(N, M, xx) * W_rec(N - 1, M - 1, xx * q**2)


def V(N, M, xx):
    return xx**(2 * M) * q**(2 * M * N + M * (2 * M - 1)) * poch(-xx, q, M) / (
        poch(Q, Q, M) * poch(Q**(N - M + 1), Q, M) * poch(xx * q**M, q, M) * poch(xx**2 * q**(2 * N), Q, M))


okS5 = all(zero(W_rec(N, M, x) - poch(x * q, q, 2 * M) * V(N, M, x)) for N in range(0, 11) for M in range(0, N + 1))
report("S5 Cor 4.2 (N<=10)", okS5)

# ---------------- S6 ----------------
rho_c = sp.symbols('c')  # c = q^rho + q^-rho
tr_ = sp.symbols('tr')   # q^rho
mu0, mu1, vn, eps_, g_ = sp.symbols('mu0 mu1 vn eps g')
aQ = x**2 * tr_**2 * q**2
sol = sp.solve([mu1 / tr_ + (x * (1 / tr_ + tr_) - eps_) * mu0 - x * tr_ * (1 - aQ) * vn,
                (1 / tr_**2 + 1) * mu1 - x * (1 + 1 / q) * (1 - aQ) * vn - g_ * mu0], [mu0, mu1], dict=True)[0]
cc = tr_ + 1 / tr_
pred = x * (1 - aQ) * (tr_**2 - 1 / q) * vn / (g_ + x * cc**2 - eps_ * cc)
report("S6 elimination", zero(sol[mu0] - pred))
# Fricke: with eps = e_N e_M, g = e_N^2+e_M^2-4x, e_n = q^-n + x q^n ; L = N-M
tL = sp.symbols('tL')
br = (en(tN)**2 + en(tM)**2 - 4 * x) + x * cc**2 - en(tN) * en(tM) * cc
fr = (tN * tM)**-1 * (tN / tM)**-1 * (1 - (tN / tM) / tr_) * (1 - (tN / tM) * tr_) * (1 - x * tN * tM / tr_) * (1 - x * tN * tM * tr_)
report("S6 Fricke factorisation of the bracket", zero(br - fr))


# ---------------- S7 ----------------
def Vs_rec(N, M):
    L = N - M
    vv = sp.Integer(1)
    for rho in range(L - 1, -1, -1):
        vv *= -x * q**(L + N + M - 1) * (1 - x**2 * q**(2 * rho + 2)) * (1 - q**(2 * rho + 1)) / (
            (1 - q**(L - rho)) * (1 - q**(L + rho)) * (1 - x * q**(N + M - rho)) * (1 - x * q**(N + M + rho)))
    return vv


def Vs(N, M):
    L = N - M
    return (-x)**L * q**(L * (2 * L - 1) + 2 * M * L) * poch(q, Q, L) / (poch(q, q, L) * poch(q**L, q, L)) * \
        poch(x**2 * Q, Q, L) / (poch(x * q**(2 * M + 1), q, L) * poch(x * q**(2 * M + L), q, L))


okS7 = all(zero(Vs_rec(N, M) - Vs(N, M)) for N in range(0, 11) for M in range(0, N + 1))
okS7b = all(zero(poch(q, Q, L) / (poch(q, q, L) * poch(q**L, q, L)) - (1 + q**L) / poch(Q, Q, L)) for L in range(1, 10))
report("S7 Cor 5.2 (N<=10)", okS7)
report("S7 (q;q^2)_L/((q)_L(q^L;q)_L) = (1+q^L)/(Q;Q)_L (1<=L<=9)", okS7b)


# ---------------- S8 ----------------
def t(N, M):
    lam = x**(N + M) * q**(N * N + M * M)
    return lam * (-1)**(N + M) * q**(N * N + M * M - N - M) * (q**N + q**M) * (1 + x * q**(N + M)) * \
        (1 - x * q**(2 * N)) * (1 - x * q**(2 * M)) * poch(x**2, Q, N) * poch(x**2, Q, M) / (
            poch(Q, Q, N) * poch(Q, Q, M) * (1 - x) * (1 - x**2) * (2 if N == M else 1))


def bB(n_):
    return (-1)**n_ * q**(n_ * n_) * poch(x**2, Q, n_) * (1 - x * q**(2 * n_)) / (poch(Q, Q, n_) * (1 - x))


def bJ(m_):
    return (-1)**m_ * q**(m_ * m_ - m_) * poch(x**2, Q, m_) * (1 - x**2 * q**(4 * m_)) / (poch(Q, Q, m_) * (1 - x**2))


okS8 = all(zero(V(N, M, x) * Vs(N, M) * poch(x**2 * q**(2 * (N - M) + 2), Q, 2 * M) * poch(x * q, q, 2 * N) - t(N, M))
           for N in range(0, 11) for M in range(0, N + 1))
report("S8 Lemma 7.1 (N<=10)", okS8)
okS9 = all(zero(t(N, M) - x**(N + M) * q**(N * N + M * M) * (bB(N) * bJ(M) + bB(M) * bJ(N)) / (2 if N == M else 1))
           for N in range(0, 9) for M in range(0, N + 1))
report("S9 t_{N,M} = lam (bB_N bJ_M + bB_M bJ_N)/(1+delta) (N<=8)", okS9)
print("ALL SYMBOLIC CHECKS OK" if ok_all else "SOME SYMBOLIC CHECK FAILED")
