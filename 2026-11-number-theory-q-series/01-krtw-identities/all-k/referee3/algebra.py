#!/usr/bin/env python3
"""Referee 3, check (5): algebraic identities, exactly.

Part I (flint fmpz_mpoly in q, x; a rational function is num/den, equality by cross-multiplication):
  L7.1   V_{N,M} V#_{N,M} (x^2 q^{2L+2};Q)_{2M} (xq;q)_{2N} = lam_{N,M} (bB_N bJ_M + bB_M bJ_N)/(1+delta)
         (right side from the DEFINITIONS of beta^B, beta^J, not from the closed form of t)
  t      closed form of t_{N,M} (F-all-k Sec. 7, from (B4)) = definition
  C4.2   prod_{i<M} [Prop 4.1 factor at (N-i, M-i, x q^{2i})] = (xq;q)_{2M} V_{N,M}(x)
  C5.2   prod_{rho<L} [Prop 5.1 factor] = V#_{N,M}
  L4.1   t (eps-1-x)/(x(1+x)) = lam_N lam_M (bOm_N bJ_M + bOm_M bJ_N)/(1+delta)
  bOm    closed form (iv-all-k Sec. 0) = definition beta^B (q^{-2n}+x-1+x^2q^{2n})/(x(1+x));  value at x = 1
  a*     alpha*_e (iv-all-k Sec. 6 item 3, computed from alpha'' and Gamma^{(1)} of (B-iv)) = lam_e^2 bOm_e
  pair   Sec. 6 item 1: sum_{r<=n} alpha_r(a) phi_r(n;a) = 1/(Q;Q)_n for a = x q (and a = x)
Part II (sympy, generic symbols): (B-iv) Gamma formula and coefficient sum; (B-v) coefficients and the
  'likewise' sums; base-Q facts of Sec. 5; Prop 4.1 elimination and bracket factorisation; Prop 5.1
  elimination and Fricke factorisation; the factorisations in iv Lemma 4.1 and (B4).
"""
import sys
from fractions import Fraction as Fr
from flint import fmpz_mpoly_ctx
import sympy as sp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8

ctx = fmpz_mpoly_ctx.get(('q', 'x'), 'lex')
qg, xg = ctx.gens()
ONEP = ctx.from_dict({(0, 0): 1})


class MR:
    """rational function num/den in q, x with integer coefficients; monomials with negative exponents allowed"""

    def __init__(self, num, den=None):
        self.n = num
        self.d = ONEP if den is None else den

    def __mul__(self, o):
        if not isinstance(o, MR):
            o = MR(ctx.from_dict({(0, 0): o}))
        return MR(self.n * o.n, self.d * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, MR):
            o = MR(ctx.from_dict({(0, 0): o}))
        return MR(self.n * o.d, self.d * o.n)

    def __add__(self, o):
        return MR(self.n * o.d + o.n * self.d, self.d * o.d)

    def __sub__(self, o):
        return MR(self.n * o.d - o.n * self.d, self.d * o.d)

    def __eq__(self, o):
        return (self.n * o.d - o.n * self.d).is_zero()


def mono(eq, ex, c=1):
    """c q^eq x^ex (negative exponents -> denominator)"""
    num = {(max(eq, 0), max(ex, 0)): c}
    den = {(max(-eq, 0), max(-ex, 0)): 1}
    return MR(ctx.from_dict(num), ctx.from_dict(den))


def one_minus(eq, ex, s=1):
    """1 - s q^eq x^ex  (eq, ex >= 0)"""
    return MR(ctx.from_dict({(0, 0): 1}) - ctx.from_dict({(eq, ex): s}))


def one_plus(eq, ex):
    return one_minus(eq, ex, -1)


def poch(eq0, ex, step, n, s=1):
    """prod_{i<n} (1 - s x^ex q^{eq0 + i step})"""
    r = MR(ONEP)
    for i in range(n):
        r = r * one_minus(eq0 + i * step, ex, s)
    return r


ONE = MR(ONEP)


def betaB(n):
    # (-1)^n q^{n^2} (x^2;Q)_n (1-x q^{2n}) / ((Q;Q)_n (1-x))
    return mono(n * n, 0, (-1) ** n) * poch(0, 2, 2, n) * one_minus(2 * n, 1) / (poch(2, 0, 2, n) * one_minus(0, 1))


def betaJ(m):
    return mono(m * m - m, 0, (-1) ** m) * poch(0, 2, 2, m) * one_minus(4 * m, 2) / (poch(2, 0, 2, m) * one_minus(0, 2))


def lam(n):
    return mono(n * n, n)


def V(N, M):
    return mono(2 * M * N + M * (2 * M - 1), 2 * M) * poch(0, 1, 1, M, s=-1) / (
        poch(2, 0, 2, M) * poch(2 * (N - M + 1), 0, 2, M) * poch(M, 1, 1, M) * poch(2 * N, 2, 2, M))


def Vs(N, M):
    L = N - M
    r = mono(L * (2 * L - 1) + 2 * M * L, L, (-1) ** L)
    r = r * poch(1, 0, 2, L) / (poch(1, 0, 1, L) * poch(L, 0, 1, L))
    r = r * poch(2, 2, 2, L) / (poch(2 * M + 1, 1, 1, L) * poch(2 * M + L, 1, 1, L))
    return r


def t_def(N, M):
    s = betaB(N) * betaJ(M) + betaB(M) * betaJ(N)
    r = lam(N) * lam(M) * s
    return r / 2 if N == M else r


def t_closed(N, M):
    r = lam(N) * lam(M) * mono(N * N + M * M - N - M, 0, (-1) ** (N + M))
    r = r * (mono(N, 0) + mono(M, 0)) * one_plus(N + M, 1) * one_minus(2 * N, 1) * one_minus(2 * M, 1)
    r = r * poch(0, 2, 2, N) * poch(0, 2, 2, M) / (poch(2, 0, 2, N) * poch(2, 0, 2, M) * one_minus(0, 1) * one_minus(0, 2))
    return r / 2 if N == M else r


def eps(N, M):
    return mono(-N - M, 0) * one_plus(2 * N, 1) * one_plus(2 * M, 1)


def bOm_def(n):
    mu = mono(-2 * n, 0) + mono(1, 1, 0) * 0 + mono(0, 1) - ONE + mono(2 * n, 2)
    return betaB(n) * mu / (mono(0, 1) * one_plus(0, 1))


def bOm_closed(n):
    if n == 0:
        return ONE
    pol = ONE - (ONE - mono(0, 1)) * mono(2 * n, 0) + mono(4 * n, 2)
    return mono(n * n - 2 * n, -1, (-1) ** n) * poch(2, 2, 2, n - 1) * one_minus(2 * n, 1) * pol / poch(2, 0, 2, n)


def prop41_factor(N, M, s):
    """Prop 4.1 factor with x -> x q^{s}:  x^2(1+x)(x+q^{-1}) q^{2N+2M} (xq;q)_2 / ((1-q^{2N})(1-q^{2M})(1-x^2q^{2N})(1-x^2q^{2M}))"""
    X = mono(s, 1)                         # the substituted x
    num = X * X * (ONE + X) * (X + mono(-1, 0)) * mono(2 * N + 2 * M, 0) * one_minus(1 + s, 1) * one_minus(2 + s, 1)
    den = one_minus(2 * N, 0) * one_minus(2 * M, 0) * one_minus(2 * N + 2 * s, 2) * one_minus(2 * M + 2 * s, 2)
    return num / den


def prop51_factor(N, M, rho):
    L = N - M
    num = mono(L + N + M - 1, 1, -1) * one_minus(2 * rho + 2, 2) * one_minus(2 * rho + 1, 0)
    den = one_minus(L - rho, 0) * one_minus(L + rho, 0) * one_minus(N + M - rho, 1) * one_minus(N + M + rho, 1)
    return num / den


def alpha_pair(n, a_eq, a_ex):
    """alpha_n relative to a = q^{a_eq} x^{a_ex}: (-1)^n q^{n^2} (a^2;Q)_n (1-a q^{2n}) / ((Q;Q)_n (1-a))"""
    return mono(n * n, 0, (-1) ** n) * poch(2 * a_eq, 2 * a_ex, 2, n) * one_minus(a_eq + 2 * n, a_ex) / (
        poch(2, 0, 2, n) * one_minus(a_eq, a_ex))


def phi_rn(r, n, a_eq, a_ex):
    """phi_r(n; a) = 1/((q)_{n-r} (aq)_{n+r})"""
    return ONE / (poch(1, 0, 1, n - r) * poch(a_eq + 1, a_ex, 1, n + r))


def part1():
    ok = True
    res = {}
    for N in range(NMAX + 1):
        for M in range(N + 1):
            L = N - M
            lhs = V(N, M) * Vs(N, M) * poch(2 * L + 2, 2, 2, 2 * M) * poch(1, 1, 1, 2 * N)
            a = lhs == t_def(N, M)
            b = t_closed(N, M) == t_def(N, M)
            W = ONE
            for i in range(M):
                W = W * prop41_factor(N - i, M - i, 2 * i)
            c = W == poch(1, 1, 1, 2 * M) * V(N, M)
            P = ONE
            for rho in range(L):
                P = P * prop51_factor(N, M, rho)
            d = P == Vs(N, M)
            lhs41 = t_def(N, M) * (eps(N, M) - ONE - mono(0, 1)) / (mono(0, 1) * one_plus(0, 1))
            rhs41 = lam(N) * lam(M) * (bOm_def(N) * betaJ(M) + bOm_def(M) * betaJ(N))
            if N == M:
                rhs41 = rhs41 / 2
            e = lhs41 == rhs41
            res[(N, M)] = (a, b, c, d, e)
            ok &= a and b and c and d and e
    names = ["L7.1", "t closed", "C4.2", "C5.2", "L4.1"]
    for i, nm in enumerate(names):
        good = all(v[i] for v in res.values())
        print(f"{nm}: all N >= M >= 0, N <= {NMAX}: {'OK' if good else 'FAIL ' + str([k for k, v in res.items() if not v[i]])}")
    # negative control for L7.1: drop the factor (xq;q)_{2N} -> (xq;q)_{2N-1}
    nc = V(2, 1) * Vs(2, 1) * poch(2 * 1 + 2, 2, 2, 2) * poch(1, 1, 1, 3) == t_def(2, 1)
    print(f"   negative control (L7.1 with (xq;q)_3 instead of (xq;q)_4 at (2,1)): equal = {nc}")
    # bOm closed form, x = 1 values, alpha* = lam^2 bOm
    nb = 2 * NMAX + 4
    okb = all(bOm_def(n) == bOm_closed(n) for n in range(nb + 1))
    print(f"bOm closed form = definition, n <= {nb}: {'OK' if okb else 'FAIL'}")
    ok &= okb
    # alpha*: alpha''_e = (xq)^{2e} q^{2e^2} alpha_e(a = xq); Gamma^{(1)}_{d,d}(x) = (1-xq)/(1-xq^{2d+1}),
    # Gamma^{(1)}_{d,d+1}(x) = -x q^{2d+1}(1-xq)/(1-xq^{2d+1})  (re-derived in Part II from the general Gamma)
    def a2(e):
        return mono(2 * e + 2 * e * e, 2 * e) * alpha_pair(e, 1, 1)

    def G11(d, e):
        if e == d:
            return one_minus(1, 1) / one_minus(2 * d + 1, 1)
        return mono(2 * d + 1, 1, -1) * one_minus(1, 1) / one_minus(2 * d + 1, 1)
    oka = True
    for e in range(nb + 1):
        ast = G11(e, e) * a2(e)
        if e >= 1:
            ast = ast + G11(e - 1, e) * a2(e - 1)
        oka &= ast == lam(e) * lam(e) * bOm_closed(e)
    print(f"alpha*_e = lam_e^2 bOm_e, e <= {nb}: {'OK' if oka else 'FAIL'}")
    ok &= oka
    # Bailey pair of Sec. 6 item 1 (general a, tested at a = x q and a = x)
    okp = True
    for (aq_, ax_) in ((1, 1), (0, 1)):
        for n in range(7):
            s = MR(ctx.from_dict({}))
            for r in range(n + 1):
                s = s + alpha_pair(r, aq_, ax_) * phi_rn(r, n, aq_, ax_)
            okp &= s == ONE / poch(2, 0, 2, n)
    print(f"Bailey pair beta_n = 1/(Q;Q)_n, alpha_n(a), a in {{xq, x}}, n <= 6: {'OK' if okp else 'FAIL'}")
    ok &= okp
    return ok


def part2():
    ok = True
    q, x, a, z, t, u, v, E, Gs, v1 = sp.symbols('q x a z t u v E G v1')

    def P(c, base, n):
        return sp.Mul(*[(1 - c * base ** i) for i in range(n)])

    # (B-iv): sum_e Gamma^{(c)}_{d,e}(a) p_e(z) = (aq;q)_c with p_e(z) = (z q^{1-e};q)_{e-d} (a z q^{e+1};q)_{d+c-e}
    def qbin(n, k):
        return sp.Mul(*[(1 - q ** (n - i)) for i in range(k)]) / sp.Mul(*[(1 - q ** (i + 1)) for i in range(k)])

    def Gam(c, d, e):
        return ((1 - a * q ** (2 * e)) / (1 - a)) * (P(a, q, d + e) / P(a * q ** (c + 1), q, d + e)) * (-1) ** (e - d) \
            * q ** ((e - d) * (e - d - 1) // 2) * a ** (e - d) * q ** ((e + d) * (e - d)) * qbin(c, e - d)
    okg = True
    for c in range(0, 5):
        for d in range(0, 4):
            S = sum(Gam(c, d, e) * P(z * q ** (1 - e), q, e - d) * P(a * z * q ** (e + 1), q, d + c - e)
                    for e in range(d, d + c + 1))
            okg &= sp.simplify(sp.cancel(S - P(a * q, q, c))) == 0
            okg &= sp.cancel(sum(Gam(c, d, e) for e in range(d, d + c + 1)) - P(a * q, q, c)) == 0
    print(f"(B-iv) Gamma formula (identity of numerators over (q)_(R-d)(aq)_(R+d+c)) and sum = (aq)_c, c<=4, d<=3: "
          f"{'OK' if okg else 'FAIL'}")
    ok &= okg
    # (B-v): with t = q^n, z = q^R; numerators over D = (q)_{R-n+1}(aq)_{R+n+1}
    Nm1 = (1 - a * z * t) * (1 - a * z * t * q)
    N0 = (1 - z * q / t) * (1 - a * z * t * q)
    Np1 = (1 - z / t) * (1 - z * q / t)
    al = t / q / ((1 - a * t * t / q) * (1 - a * t * t))
    be = -t / q * (1 + a * t * t) / ((1 - a * t * t / q) * (1 - a * t * t * q))
    ga = a * t ** 3 / ((1 - a * t * t) * (1 - a * t * t * q))
    okv = sp.cancel(al * Nm1 + be * N0 + ga * Np1 - z) == 0 and sp.cancel(al + be + ga) == 0
    A1, A2, A3 = sp.symbols('A1 A2 A3')
    for target, ssum in ((sp.Integer(1), 1), (z * z, 0)):
        sol = sp.solve(sp.Poly(sp.expand(A1 * Nm1 + A2 * N0 + A3 * Np1 - target), z).all_coeffs(), [A1, A2, A3], dict=True)
        okv &= len(sol) == 1 and sp.cancel(sol[0][A1] + sol[0][A2] + sol[0][A3] - ssum) == 0
    print(f"(B-v) alpha_n, beta_n, gamma_n and the 'likewise' coefficient sums 1, 0: {'OK' if okv else 'FAIL'}")
    ok &= okv
    # base-Q facts (Sec. 5): Z = Q^S, T = Q^e; numerators over (Q;Q)_{S-e} (aQ;Q)_{S+e+1} resp. (aQ;Q)_{S+e+2}
    Qs, Z, T = sp.symbols('Qs Z T')
    n_e = 1 - a * T * Qs * Z           # phi^Q_e(S;a)
    n_e1 = 1 - Z / T                   # phi^Q_{e+1}(S;a)
    c1, c2 = sp.symbols('c1 c2')
    sol = sp.solve(sp.Poly(sp.expand(c1 * n_e + c2 * n_e1 - (1 - a * Qs)), Z).all_coeffs(), [c1, c2], dict=True)[0]
    okq = sp.cancel(sol[c1] + sol[c2] - (1 - a * Qs)) == 0
    be_ = (1 - a * Qs) / (1 / T - a * T * Qs)
    okq &= sp.cancel(be_ * (n_e - n_e1) - Z * (1 - a * Qs)) == 0
    # Q^S phi^Q_e(S;aQ^2): over (Q;Q)_{S-e}(aQ;Q)_{S+e+2}, numerator Z (aQ;Q)_2; basis numerators
    m0 = (1 - a * T * Qs * Z) * (1 - a * T * Qs ** 2 * Z)
    m1 = (1 - Z / T) * (1 - a * T * Qs ** 2 * Z)
    m2 = (1 - Z / T) * (1 - Z / (T * Qs))
    sol = sp.solve(sp.Poly(sp.expand(A1 * m0 + A2 * m1 + A3 * m2 - Z * (1 - a * Qs) * (1 - a * Qs ** 2)), Z).all_coeffs(),
                   [A1, A2, A3], dict=True)
    okq &= len(sol) == 1 and sp.cancel(sol[0][A1] + sol[0][A2] + sol[0][A3]) == 0
    print(f"base-Q facts of Sec. 5 (c1+c2 = 1-aQ, b_e, sum 0 for Q^S phi^Q_e(S;aQ^2)): {'OK' if okq else 'FAIL'}")
    ok &= okq
    # Prop 4.1 elimination
    m0, m1, m2 = sp.symbols('m0 m1 m2')
    sol = sp.solve([(1 + x) * m1 - E * m0, (1 + x) * m2 - E * m1 - x ** 2 * (1 - 1 / q) * v1,
                    m2 + (1 - x) ** 2 * m0 - x ** 2 * v1 - Gs * m0], [m0, m1, m2], dict=True)[0]
    claim = x ** 2 * (1 + x) * (x + 1 / q) * v1 / (E ** 2 + (1 + x) ** 2 * ((1 - x) ** 2 - Gs))
    oke = sp.cancel(sol[m0] - claim) == 0
    eN, eM = 1 / u + x * u, 1 / v + x * v
    br = (eN * eM) ** 2 + (1 + x) ** 2 * ((1 - x) ** 2 - (eN ** 2 + eM ** 2 - 4 * x))
    oke &= sp.expand(br - (eN ** 2 - (1 + x) ** 2) * (eM ** 2 - (1 + x) ** 2)) == 0
    oke &= sp.cancel(eN ** 2 - (1 + x) ** 2 - (1 - u * u) * (1 - x * x * u * u) / u ** 2) == 0
    oke &= sp.expand(eN ** 2 + eM ** 2 - 4 * x - (1 / u ** 2 + x * x * u * u + 1 / v ** 2 + x * x * v * v)) == 0
    print(f"Prop 4.1 elimination, bracket = (e_N^2-(1+x)^2)(e_M^2-(1+x)^2), e_n^2-(1+x)^2 = q^-2n(1-q^2n)(1-x^2q^2n): "
          f"{'OK' if oke else 'FAIL'}")
    ok &= oke
    # Prop 5.1 elimination (t = q^rho, c = t + 1/t, aQ = x^2 q^{2 rho + 2})
    cc = t + 1 / t
    aQ = x * x * t * t * q * q
    sol = sp.solve([m1 / t + (x * cc - E) * m0 - x * t * (1 - aQ) * v1,
                    (1 / t ** 2 + 1) * m1 - x * (1 + 1 / q) * (1 - aQ) * v1 - Gs * m0], [m0, m1], dict=True)[0]
    claim = x * (1 - aQ) * (t * t - 1 / q) * v1 / (Gs + x * cc ** 2 - E * cc)
    okf = sp.cancel(sol[m0] - claim) == 0
    Ev, Gv = eN * eM, eN ** 2 + eM ** 2 - 4 * x
    br = Gv + x * cc ** 2 - Ev * cc
    # q^L = u/v, q^{N+M} = u v
    fr = (v / u) * (1 / (u * v)) * (1 - u / v / t) * (1 - u / v * t) * (1 - x * u * v / t) * (1 - x * u * v * t)
    okf &= sp.cancel(br - fr) == 0
    print(f"Prop 5.1 elimination and Fricke factorisation of the bracket: {'OK' if okf else 'FAIL'}")
    ok &= okf
    # iv Lemma 4.1 factorisations and (B4) E-identity
    epsu = 1 / (u * v) + x * u / v + x * v / u + x * x * u * v
    ok41 = sp.cancel(epsu - (1 / u ** 2 + x * x * u * u) - 2 * x - (u - v) * (1 + x * u * u) * (1 - x * u * v) / (u * u * v)) == 0
    ok41 &= sp.cancel(epsu - (1 / v ** 2 + x * x * v * v) - 2 * x - (v - u) * (1 + x * v * v) * (1 - x * u * v) / (v * v * u)) == 0
    EE = (1 / v) * (1 - x * u * u) * (1 - x * x * v ** 4) + (1 / u) * (1 - x * v * v) * (1 - x * x * u ** 4)
    ok41 &= sp.cancel(EE - (u + v) * (1 + x * u * v) * (1 - x * u * u) * (1 - x * v * v) / (u * v)) == 0
    print(f"iv Lemma 4.1 factorisations of eps - mu - 2x, and the (B4) E-identity: {'OK' if ok41 else 'FAIL'}")
    ok &= ok41
    return ok


if __name__ == "__main__":
    o1 = part1()
    o2 = part2()
    print("ALL ALGEBRA CHECKS " + ("PASSED" if (o1 and o2) else "FAILED"))
