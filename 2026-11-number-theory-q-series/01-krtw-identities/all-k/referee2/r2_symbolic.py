"""Referee 2, computation (6) and further algebra: sympy + exact rationals.

S1  Lemma 2.1 normal-ordered coefficients (symbols qj = q^j, Qb = Q^b, P, Y), and W h W^-1 = h~ on monomials.
S2  Lemma 4.1 for ALL N, M: reduction to rational identities in u = q^N, v = q^M, x:
      eps - mu_N - 2x = (u-v)(1+xu^2)(1-xuv)/(u^2 v), and rho_n := beta^B_n/beta^J_n = (1+x) u_n/(1+x u_n^2)
      (the (x^2;Q)_n/(Q;Q)_n parts cancel), so Lemma 4.1 <=> rho_N (eps-mu_N-2x) + rho_M (eps-mu_M-2x) = 0.
S3  Lemma 4.1 itself, symbolically in (q, x) for N, M <= 6, and exactly at random rational points for N, M <= 14.
S4  b^Omega_n: definition == factored form (n <= 10, symbolic), and value at x = 1 (from the factored form) == (-1)^n (q^{n^2-2n}+q^{n^2+2n}).
S5  Theorem 6 step 1: (1/(Q;Q)_n, (-1)^n q^{n^2} (a^2;Q)_n (1-aq^{2n})/((Q;Q)_n (1-a))) is a Bailey pair relative to a
    (general a, n <= 7, symbolic in a, q).
S6  Theorem 6 step 3: Gamma^(1) re-expansion phi_d(R;aq) = G_dd phi_d(R;a) + G_{d,d+1} phi_{d+1}(R;a) (symbolic, d <= 4, R <= 9),
    and alpha*_e computed from alpha''_e by that formula == lambda_e^2 b^Omega_e (e <= 10, symbolic).
S7  Theorem 6 steps 1-4 at FINITE n: beta^{(k)}_n (the multisum with 1/(q)_{n-N_1}) == sum_{r<=n} lambda_r^{k+1} b^Omega_r
    / ((q)_{n-r} (xq)_{n+r}), exact at rational points, k <= 6, n <= 7.  Positive control: the same for (B1) and B_{k+1}.
    Negative control: b^Omega with x-1 -> x+1.
"""
import random, itertools
from fractions import Fraction as Fr
import sympy as sp

random.seed(2026)
q, x, a, u, v, qj, Qb, P, Y, z = sp.symbols('q x a u v qj Qb P Y z')
lines = []


def pr(s):
    print(s, flush=True)
    lines.append(s)


def poch(c, base, n):
    r = sp.Integer(1)
    for j in range(n):
        r *= (1 - c * base ** j)
    return r


Qs = q ** 2
allok = True


def is_zero_rat(expr):
    """gcd-free zero test for a rational function: common denominator, expand numerator."""
    num, den = sp.fraction(sp.together(expr))
    return sp.expand(num) == 0


# ---------------------------------------------------------------- S1
hE = qj / P + x * P * Y
Eht = qj / P - x / qj * Qb * (1 - qj) * P * Y + x / qj * Qb * P * Y
xPYuE = x * P * Y * (1 - Qb)
ok = sp.simplify(hE - Eht - xPYuE) == 0
# W h W^-1 = h~ : check (W h F)(R,S) = (h~ W F)(R,S) pointwise, F an arbitrary function (symbols f_{R,S})
w = lambda r, s: x ** (r + 2 * s) * q ** ((r + s) ** 2 + s ** 2)
fs = {(r, s): sp.Symbol('f_%d_%d' % (r, s)) for r in range(-1, 6) for s in range(0, 5)}
Fv = lambda r, s: 0 if r < 0 else fs[(r, s)]
ok2 = True
for Rv in range(0, 6):
    for Sv in range(0, 5):
        WhF = w(Rv, Sv) * (q ** (-Rv) * (Fv(Rv, Sv) - Fv(Rv - 1, Sv)) + x * q ** (Rv + 2 * Sv) * Fv(Rv, Sv))
        htWF = (q ** (-Rv) * w(Rv, Sv) * Fv(Rv, Sv) - x / q * q ** Rv * q ** (2 * Sv) * (w(Rv - 1, Sv) * Fv(Rv - 1, Sv) if Rv >= 1 else 0)
                + x * q ** Rv * q ** (2 * Sv) * w(Rv, Sv) * Fv(Rv, Sv))
        if sp.expand(WhF - htWF) != 0:
            ok2 = False
pr("S1  Lemma 2.1 normal-ordered coefficients: hE - Eh~ = xPYuE: %s;  W h = h~ W (h~ = P^-1 - x q^-1 P Y Ssh + x P Y): %s"
   % ("OK" if ok else "FAIL", "OK" if ok2 else "FAIL"))
allok &= ok and ok2

# ---------------------------------------------------------------- S2
eps_uv = (1 + x * u ** 2) * (1 + x * v ** 2) / (u * v)
mu = lambda t: t ** -2 + x ** 2 * t ** 2
d1 = sp.cancel(eps_uv - mu(u) - 2 * x - (u - v) * (1 + x * u ** 2) * (1 - x * u * v) / (u ** 2 * v))
d2 = sp.cancel(eps_uv - mu(v) - 2 * x - (v - u) * (1 + x * v ** 2) * (1 - x * u * v) / (v ** 2 * u))
rho = lambda t: (1 + x) * t / (1 + x * t ** 2)
d3 = sp.cancel(rho(u) * (eps_uv - mu(u) - 2 * x) + rho(v) * (eps_uv - mu(v) - 2 * x))
# rho_n from the definitions, as a function of t = q^n: beta^B/beta^J = q^n (1-xq^{2n})(1-x^2)/((1-x)(1-x^2 q^{4n}))
t = sp.symbols('t')
d4 = sp.cancel(t * (1 - x * t ** 2) * (1 - x ** 2) / ((1 - x) * (1 - x ** 2 * t ** 4)) - rho(t))
ok = d1 == 0 and d2 == 0 and d3 == 0 and d4 == 0
pr("S2  Lemma 4.1 for all N,M (u=q^N, v=q^M): factorisations %s, %s; ratio beta^B/beta^J = (1+x)q^n/(1+xq^{2n}) %s; "
   "rho_N(eps-mu_N-2x)+rho_M(eps-mu_M-2x) = 0 %s" % (d1 == 0, d2 == 0, d4 == 0, d3 == 0))
allok &= ok


# ---------------------------------------------------------------- definitions
def betaB(n):
    return (-1) ** n * q ** (n * n) * poch(x ** 2, Qs, n) * (1 - x * q ** (2 * n)) / (poch(Qs, Qs, n) * (1 - x))


def betaJ(m):
    return (-1) ** m * q ** (m * m - m) * poch(x ** 2, Qs, m) * (1 - x ** 2 * q ** (4 * m)) / (poch(Qs, Qs, m) * (1 - x ** 2))


def lam(n):
    return x ** n * q ** (n * n)


def bOm(n):
    return betaB(n) * (q ** (-2 * n) + x - 1 + x ** 2 * q ** (2 * n)) / (x * (1 + x))


def bOm_fact(n):
    if n == 0:
        return sp.Integer(1)
    return ((-1) ** n * q ** (n * n - 2 * n) * poch(x ** 2 * Qs, Qs, n - 1) * (1 - x * q ** (2 * n))
            * (1 - (1 - x) * q ** (2 * n) + x ** 2 * q ** (4 * n)) / (x * poch(Qs, Qs, n)))


def eps(N, M):
    return q ** (-N - M) * (1 + x * q ** (2 * N)) * (1 + x * q ** (2 * M))


def tNM(N, M):
    return lam(N) * lam(M) * (betaB(N) * betaJ(M) + betaB(M) * betaJ(N)) / (2 if N == M else 1)


def lemma41_diff(N, M):
    return tNM(N, M) * (eps(N, M) - 1 - x) / (x * (1 + x)) - lam(N) * lam(M) * (bOm(N) * betaJ(M) + bOm(M) * betaJ(N)) / (2 if N == M else 1)


# ---------------------------------------------------------------- S3
bad = []
for N in range(0, 7):
    for M in range(0, N + 1):
        if not is_zero_rat(lemma41_diff(N, M)):
            bad.append((N, M))
pr("S3a Lemma 4.1 symbolic in (q,x), all 0 <= M <= N <= 6: %s" % ("OK" if not bad else "FAIL %s" % bad))
allok &= not bad


def lemma41_diff_exact(N, M, qq, xx):
    Qq = qq * qq

    def fp(c, base, n):
        r = Fr(1)
        for j in range(n):
            r *= (1 - c * base ** j)
        return r
    bB = lambda n: (-1) ** n * qq ** (n * n) * fp(xx * xx, Qq, n) * (1 - xx * qq ** (2 * n)) / (fp(Qq, Qq, n) * (1 - xx))
    bJ = lambda m: (-1) ** m * qq ** (m * m - m) * fp(xx * xx, Qq, m) * (1 - xx * xx * qq ** (4 * m)) / (fp(Qq, Qq, m) * (1 - xx * xx))
    lm = lambda n: xx ** n * qq ** (n * n)
    bO = lambda n: bB(n) * (qq ** (-2 * n) + xx - 1 + xx * xx * qq ** (2 * n)) / (xx * (1 + xx))
    ep = qq ** (-N - M) * (1 + xx * qq ** (2 * N)) * (1 + xx * qq ** (2 * M))
    dl = 2 if N == M else 1
    t = lm(N) * lm(M) * (bB(N) * bJ(M) + bB(M) * bJ(N)) / dl
    return t * (ep - 1 - xx) / (xx * (1 + xx)) - lm(N) * lm(M) * (bO(N) * bJ(M) + bO(M) * bJ(N)) / dl


bad = []
cnt = 0
for N in range(0, 21):
    for M in range(0, N + 1):
        for _ in range(4):
            qq = Fr(random.randint(1, 9), random.randint(10, 19)) * random.choice([1, -1])
            xx = Fr(random.randint(1, 9), random.randint(10, 23)) * random.choice([1, -1])
            cnt += 1
            if lemma41_diff_exact(N, M, qq, xx) != 0:
                bad.append((N, M, qq, xx))
pr("S3b Lemma 4.1 exact at random rational points, all 0 <= M <= N <= 20 (%d evaluations): %s" % (cnt, "OK" if not bad else "FAIL %s" % bad[:3]))
allok &= not bad

# ---------------------------------------------------------------- S4
bad = []
bad1 = []
for n in range(0, 11):
    if not is_zero_rat(bOm(n) - bOm_fact(n)):
        bad.append(n)
    val1 = sp.together(bOm_fact(n)).subs(x, 1)  # factored form has no pole at x = 1 (it equals the definition, checked above)
    target = 1 if n == 0 else (-1) ** n * (q ** (n * n - 2 * n) + q ** (n * n + 2 * n))
    if not is_zero_rat(val1 - target):
        bad1.append(n)
pr("S4  b^Omega_n definition == factored form (n <= 10): %s;  b^Omega_n(x=1) == (-1)^n(q^{n^2-2n}+q^{n^2+2n}) (n <= 10; "
   "from the factored form): %s" % ("OK" if not bad else "FAIL %s" % bad, "OK" if not bad1 else "FAIL %s" % bad1))
allok &= not bad and not bad1

# ---------------------------------------------------------------- S5
bad = []
for n in range(0, 8):
    beta_n = 1 / poch(Qs, Qs, n)
    rhs = sum((-1) ** r * q ** (r * r) * poch(a ** 2, Qs, r) * (1 - a * q ** (2 * r)) / (poch(Qs, Qs, r) * (1 - a))
              / (poch(q, q, n - r) * poch(a * q, q, n + r)) for r in range(n + 1))
    if not is_zero_rat(rhs - beta_n):
        bad.append(n)
pr("S5  Bailey pair of Theorem 6 step 1 for general a (n <= 7): %s" % ("OK" if not bad else "FAIL %s" % bad))
allok &= not bad


# ---------------------------------------------------------------- S6
def phiB(e, Rv, aa):
    if Rv < e:
        return sp.Integer(0)
    return 1 / (poch(q, q, Rv - e) * poch(aa * q, q, Rv + e))


G_dd = lambda d, aa: (1 - aa * q) / (1 - aa * q ** (2 * d + 1))
G_dd1 = lambda d, aa: -aa * q ** (2 * d + 1) * (1 - aa * q) / (1 - aa * q ** (2 * d + 1))
bad = []
for d in range(0, 5):
    for Rv in range(0, 10):
        diff = phiB(d, Rv, a * q) - G_dd(d, a) * phiB(d, Rv, a) - G_dd1(d, a) * phiB(d + 1, Rv, a)
        if not is_zero_rat(diff):
            bad.append((d, Rv))
ok6a = not bad


def alpha_pp(n):  # two (L_inf) steps relative to xq applied to the step-1 alpha with a = xq
    aa = x * q
    al = (-1) ** n * q ** (n * n) * poch(aa ** 2, Qs, n) * (1 - aa * q ** (2 * n)) / (poch(Qs, Qs, n) * (1 - aa))
    return (aa ** n * q ** (n * n)) ** 2 * al


bad2 = []
for e in range(0, 11):
    ast = alpha_pp(e) * G_dd(e, x) + (alpha_pp(e - 1) * G_dd1(e - 1, x) if e >= 1 else 0)
    if not is_zero_rat(ast - lam(e) ** 2 * bOm(e)):
        bad2.append(e)
pr("S6  Gamma^(1) re-expansion (d <= 4, R <= 9): %s;  alpha*_e (from alpha'' and Gamma^(1)) == lambda_e^2 b^Omega_e (e <= 10): %s"
   % ("OK" if ok6a else "FAIL %s" % bad, "OK" if not bad2 else "FAIL %s" % bad2))
allok &= ok6a and not bad2


# ---------------------------------------------------------------- S7 (exact Fractions)
def S7(k, nmax, qq, xx, which='omega', perturb=False):
    Qq = qq * qq

    def fp(c, base, n):
        r = Fr(1)
        for j in range(n):
            r *= (1 - c * base ** j)
        return r

    def betaB_f(n):
        return (-1) ** n * qq ** (n * n) * fp(xx * xx, Qq, n) * (1 - xx * qq ** (2 * n)) / (fp(Qq, Qq, n) * (1 - xx))

    def bO_f(n):
        c = 1 if (perturb and n >= 1) else -1
        return betaB_f(n) * (qq ** (-2 * n) + xx + c + xx * xx * qq ** (2 * n)) / (xx * (1 + xx))

    lin = (k, k + 1) if which == 'omega' else ()
    res = []
    for n in range(nmax + 1):
        # beta^{(k)}_n: sum over n >= N_1 >= ... >= N_{k+1} >= 0
        tot = Fr(0)
        for Ns in itertools.combinations_with_replacement(range(n, -1, -1), k + 1):
            # combinations_with_replacement of a decreasing range gives nonincreasing tuples
            term = xx ** sum(Ns) * qq ** (sum(t * t for t in Ns) + sum(Ns[j - 1] for j in lin))
            term /= fp(qq, qq, n - Ns[0])
            for i in range(k):
                term /= fp(qq, qq, Ns[i] - Ns[i + 1])
            term /= fp(Qq, Qq, Ns[k])
            tot += term
        if which == 'omega':
            al = lambda r: (xx ** r * qq ** (r * r)) ** (k + 1) * bO_f(r)
        else:
            al = lambda r: (xx ** r * qq ** (r * r)) ** (k + 1) * betaB_f(r)
        rhs = sum(al(r) / (fp(qq, qq, n - r) * fp(xx * qq, qq, n + r)) for r in range(n + 1))
        res.append(tot == rhs)
    return all(res)


bad = []
badB = []
negdet = 0
negtot = 0
for k in range(1, 7):
    for (qq, xx) in [(Fr(2, 7), Fr(3, 11)), (Fr(-3, 5), Fr(5, 13))]:
        if not S7(k, 7 if k <= 4 else 5, qq, xx, 'omega'):
            bad.append((k, qq, xx))
        if not S7(k, 7 if k <= 4 else 5, qq, xx, 'B'):
            badB.append((k, qq, xx))
        negtot += 1
        if not S7(k, 4, qq, xx, 'omega', perturb=True):
            negdet += 1
pr("S7  Theorem 6 Bailey pair at finite n (beta^{(k)}_n == sum_r lam_r^{k+1} b^Om_r/((q)_{n-r}(xq)_{n+r})), k <= 6, n <= 7 (k<=4) / 5, "
   "2 points: %s;  positive control (B1)-pair for B_{k+1}: %s;  negative control (x-1 -> x+1) detected %d/%d"
   % ("OK" if not bad else "FAIL %s" % bad, "OK" if not badB else "FAIL %s" % badB, negdet, negtot))
allok &= not bad and not badB and negdet == negtot

pr("RESULT: %s" % ("ALL OK" if allok else "SOME CHECK FAILED"))
