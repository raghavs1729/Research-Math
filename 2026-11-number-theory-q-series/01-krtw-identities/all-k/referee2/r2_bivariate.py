"""Referee 2, computations (2) and (5): bivariate checks in Z[[x,q]] mod (x^(D+1), q^L), exact integers.

Every series involved has q-exponent >= x-exponent termwise, so D = L-1 captures ALL powers of x mod q^L.

Built from definitions only:
  K_k, K^Omega_k (Omega literal), K^flat_k (extra q^{-r_1})   -- level DP and (k <= kbrute) naive enumeration
  C_k, J_k := C_k / C_1, B_{k+1}, B^{(k)}_{k+1}                 -- multisum DPs
  beta^B_n, beta^J_m, lambda_n, t_{N,M}, eps_{N,M} as in F-all-k.md section 7 / Lemma 3.2
  b^Omega_n := beta^B_n (q^{-2n} + x - 1 + x^2 q^{2n}) / (x(1+x))   (the definition, not the factored form)
Checks:
  P23   x K^Omega_k = K^flat_k - K_k                      (Proposition 2.3)            [DP; brute force for k <= kbrute]
  T32a  (xq)_inf (x^2Q;Q)_inf K^flat_k  = sum_{N>=M} lam^k t eps/(1+x)          (Theorem 3.2, first formula)
  T32b  (xq)_inf (x^2Q;Q)_inf K^Omega_k = sum_{N>=M} lam^k t (eps-1-x)/(x(1+x)) (Theorem 3.2, second formula)
  T32c  (xq)_inf (x^2Q;Q)_inf K_k       = sum_{N>=M} lam^k t                    (F-all-k, positive control)
  L41s  sum_{N>=M} lam^k t (eps-1-x)/(x(1+x)) = (sum_n lam_n^{k+1} b_n)(sum_m lam_m^{k+1} beta^J_m)  (Lemma 4.1 summed)
  B2    (x^2Q;Q)_inf J_k = sum_m lam_m^{k+1} beta^J_m
  TH4   (xq)_inf K^Omega_k = J_k * sum_n lam_n^{k+1} b^Omega_n               (Theorem 4)
  TH6a  (xq)_inf B^{(k)}_{k+1} = sum_n lam_n^{k+1} b^Omega_n                 (Theorem 6)
  TH6b  K^Omega_k C_1 = B^{(k)}_{k+1} C_k                                     (Theorem 6, division-free)
  Fk    K_k C_1 = B_{k+1} C_k                                                 (positive control, Theorem 1 of F-all-k)
Negative controls: perturbed b^Omega (x-1 -> x+1), Omega without its (1-q^{2(s_{i+1}-s_i)}) factors, eps_{N,M} with
q^{-N-M} -> q^{-N-M+1} for one mode, K^flat with q^{-r_k} instead of q^{-r_1}.
Usage: python3 r2_bivariate.py [L=40] [kmax=6] [kbrute=5]
"""
import sys, time, itertools
from flint import fmpz_poly
from r2_lib import mono, inv, poch, InvCache, BS, bpoch_x

L = int(sys.argv[1]) if len(sys.argv) > 1 else 40
kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
kbrute = int(sys.argv[3]) if len(sys.argv) > 3 else 5
D = L - 1
cache = InvCache(L, 4 * L + 10)
one = fmpz_poly([1])


def X(d):
    return BS.from_terms(D, L, [(d, 0, 1)])


inv1mx = BS.xseries(D, L, lambda d: 1)
inv1px = BS.xseries(D, L, lambda d: (-1) ** d)
inv1mx2 = BS.xseries(D, L, lambda d: 1 if d % 2 == 0 else 0)
xq_inf = bpoch_x(D, L, 1, 1, 1, 'inf')            # (xq;q)_inf
x2Q_inf = bpoch_x(D, L, 2, 2, 2, 'inf')           # (x^2 q^2; q^2)_inf


# ------------------------------------------------------------------ multisums
def states():
    # (R+S)^2+S^2-R < L: level 1 of K^flat carries the extra q^{-r_1}; at other levels the extra states have
    # w = 0 mod q^L and contribute nothing, so one superset is used everywhere.
    st = []
    for R in range(D + 1):
        for S in range(D + 1):
            if R + 2 * S <= D and (R + S) ** 2 + S * S - R < L:
                st.append((R, S))
    return st


ST = states()


def wexp(r, s):
    return (r + 2 * s, (r + s) ** 2 + s * s)


def K_dp(k, omega_variant='iv'):
    """returns (K_k, K^Omega_k, K^flat_k) by the level DP; Omega literal (variant 'nofactor' drops the
    (1-q^{2(s_{i+1}-s_i)}) factors: negative control)."""
    G0 = {p: BS.const(D, L, cache.q(p[0])) for p in ST}
    G1 = {p: BS.const(D, L, cache.q(p[0])).shift(0, p[0] + 2 * p[1]) for p in ST}
    kern, kins = {}, {}
    for i in range(k - 1, 0, -1):
        N0, N1 = {}, {}
        for (R, S) in ST:
            a0 = BS(D, L)
            a1 = BS(D, L)
            ai = BS(D, L)
            for (r, s) in ST:
                if r > R or s < S:
                    continue
                key = (R - r, s - S)
                if key not in kern:
                    kern[key] = cache.q(R - r).mul_low(cache.Q(s - S), L)
                    f = one if omega_variant == 'nofactor' else (one - mono(2 * (s - S), L))
                    kins[key] = kern[key].mul_low(f, L)
                xa, qb = wexp(r, s)
                g0 = G0[(r, s)].shift(xa, qb)
                g1 = G1[(r, s)].shift(xa, qb)
                a0 = a0 + g0.mulq(kern[key])
                a1 = a1 + g1.mulq(kern[key])
                ai = ai + g0.mulq(kins[key])
            N0[(R, S)] = a0
            N1[(R, S)] = a1 + ai.shift(0, R + 2 * S)
        G0, G1 = N0, N1
    K = BS(D, L)
    KO = BS(D, L)
    Kf = BS(D, L)
    for (R, S) in ST:
        xa, qb = wexp(R, S)
        K = K + G0[(R, S)].mulq(cache.Q(S)).shift(xa, qb)
        KO = KO + G1[(R, S)].mulq(cache.Q(S)).shift(xa, qb)
        Kf = Kf + G0[(R, S)].mulq(cache.Q(S)).shift(xa, qb - R)
    return K, KO, Kf


def K_brute(k, flat_level=1):
    """naive enumeration of all tuples; returns (K, K^Omega, K^flat) with K^flat using q^{-r_{flat_level}}."""
    K = [fmpz_poly() for _ in range(D + 1)]
    KO = [fmpz_poly() for _ in range(D + 1)]
    Kf = [fmpz_poly() for _ in range(D + 1)]
    Rm = int(L ** 0.5) + 2

    def rchains(n, top):
        if n == 0:
            yield ()
            return
        for r in range(top, -1, -1):
            for rest in rchains(n - 1, r):
                yield (r,) + rest

    def schains(n, bot):
        if n == 0:
            yield ()
            return
        for s in range(bot, Rm + 1):
            for rest in schains(n - 1, s):
                yield (s,) + rest

    cnt = 0
    for rs in rchains(k, Rm):
        if sum(r * r for r in rs) - rs[0] >= L:
            continue
        for ss in schains(k, 0):
            e = sum((rs[i] + ss[i]) ** 2 + ss[i] ** 2 for i in range(k))
            if e - rs[flat_level - 1] >= L:
                continue
            d = sum(rs[i] + 2 * ss[i] for i in range(k))
            if d > D:
                continue
            cnt += 1
            t = one
            for i in range(k):
                rn = rs[i + 1] if i + 1 < k else 0
                sp = ss[i - 1] if i > 0 else 0
                t = t.mul_low(cache.q(rs[i] - rn), L).mul_low(cache.Q(ss[i] - sp), L)
            om = fmpz_poly()
            for i in range(k - 1):
                om += mono(rs[i] + 2 * ss[i], L) * (one - mono(2 * ss[i + 1] - 2 * ss[i], L))
            om += mono(rs[k - 1] + 2 * ss[k - 1], L)
            base = t.mul_low(mono(e, L), L)
            K[d] += base
            KO[d] += base.mul_low(om, L)
            ef = e - rs[flat_level - 1]
            Kf[d] += t.mul_low(mono(ef, L), L)
    return BS(D, L, [p.truncate(L) for p in K]), BS(D, L, [p.truncate(L) for p in KO]), BS(D, L, [p.truncate(L) for p in Kf]), cnt


def K_operator_form(k, wrong=False):
    """Lemma 1.1 in operator form: sum_{i=1}^k omega(T^{i-1} O T^{k-i} u), O F(R,S) = q^{R+2S} F(R,S+1),
    T and omega implemented as operators on functions of (R,S).  wrong=True: negative control with the
    insertion placed one level off (omega(T^{i} O T^{k-i-1} u) for i < k, i.e. O applied before the bond)."""
    STset = set(ST)
    EXT = sorted(STset | {(R, S + 1) for (R, S) in ST})

    def Tapply(G, outs):
        out = {}
        for (R, S) in outs:
            acc = BS(D, L)
            for (r, s) in ST:
                if r > R or s < S or (r, s) not in G:
                    continue
                xa, qb = wexp(r, s)
                acc = acc + G[(r, s)].shift(xa, qb).mulq(cache.q(R - r).mul_low(cache.Q(s - S), L))
            out[(R, S)] = acc
        return out

    def Oapply(F):
        return {(R, S): (F[(R, S + 1)].shift(0, R + 2 * S) if (R, S + 1) in F else BS(D, L)) for (R, S) in ST}

    def omega(G):
        tot = BS(D, L)
        for (R, S) in ST:
            xa, qb = wexp(R, S)
            tot = tot + G[(R, S)].mulq(cache.Q(S)).shift(xa, qb)
        return tot

    u = {p: BS.const(D, L, cache.q(p[0])) for p in EXT}
    total = BS(D, L)
    for i in range(1, k + 1):
        a, b = (i - 1, k - i)
        if wrong and i < k:
            a, b = i, k - i - 1
        F = u
        for _ in range(b):
            F = Tapply(F, EXT)
        G = Oapply(F)
        for _ in range(a):
            G = Tapply(G, ST)
        total = total + omega(G)
    return total


def C_dp(k):
    f = {s: BS.from_terms(D, L, [(s, s * s, 1)]).mulq(cache.Q(s)) for s in range(D + 1) if s * s < L}
    for i in range(2, k + 1):
        g = {}
        for s2 in f:
            acc = BS(D, L)
            for s in f:
                if s <= s2:
                    acc = acc + f[s].mulq(cache.Q(s2 - s))
            g[s2] = acc.shift(s2, s2 * s2)
        f = g
    tot = BS(D, L)
    for s in f:
        tot = tot + f[s]
    return tot


def B_dp(k, lin_levels=()):
    """B^{(.)}_{k+1}: sum_{N_1>=...>=N_{k+1}} x^{sum N} q^{sum N^2 + sum_{j in lin} N_j} / ((q)_{N_1-N_2}...(Q;Q)_{N_{k+1}})."""
    Ns = [n for n in range(D + 1) if n * n < L]
    lin = lambda j, n: n if j in lin_levels else 0
    f = {n: BS.from_terms(D, L, [(n, n * n + lin(k + 1, n), 1)]).mulq(cache.Q(n)) for n in Ns}
    for j in range(k, 0, -1):
        g = {}
        for n2 in Ns:
            acc = BS(D, L)
            for n in Ns:
                if n <= n2:
                    acc = acc + f[n].mulq(cache.q(n2 - n))
            g[n2] = acc.shift(n2, n2 * n2 + lin(j, n2))
        f = g
    tot = BS(D, L)
    for n in Ns:
        tot = tot + f[n]
    return tot


# ------------------------------------------------------------------ spectral side
def PB(n):
    """(x^2;Q)_n (1 - x q^{2n}) / ((Q;Q)_n (1-x))  (beta^B_n without (-1)^n q^{n^2})."""
    return (bpoch_x(D, L, 2, 0, 2, n) * BS.from_terms(D, L, [(0, 0, 1), (1, 2 * n, -1)]) * inv1mx).mulq(cache.Q(n))


def PJ(m):
    """(x^2;Q)_m (1 - x^2 q^{4m}) / ((Q;Q)_m (1-x^2))  (beta^J_m without (-1)^m q^{m^2-m})."""
    return (bpoch_x(D, L, 2, 0, 2, m) * BS.from_terms(D, L, [(0, 0, 1), (2, 4 * m, -1)]) * inv1mx2).mulq(cache.Q(m))


PBc, PJc = {}, {}


def pb(n):
    if n not in PBc:
        PBc[n] = PB(n)
    return PBc[n]


def pj(m):
    if m not in PJc:
        PJc[m] = PJ(m)
    return PJc[m]


def place(series, a, b):
    """x^a q^b * series, allowing a = -1 (then series must be divisible by x)."""
    if a < 0:
        series = series.shift_down_x(-a)
        a = 0
    assert b >= 0, (a, b)
    return series.shift(a, b)


def spectral(k, kind, perturb=None):
    """sum_{N>=M} lam_{N,M}^k t_{N,M} * F, F = 1 ('t'), eps/(1+x) ('flat'), (eps-1-x)/(x(1+x)) ('omega')."""
    tot = BS(D, L)
    N = 0
    while (k + 1) * N - 1 <= D and (k + 2) * N * N - 2 * N < L:
        for M in range(0, N + 1):
            if (k + 1) * (N + M) - 1 > D:
                break
            prods = [(N, M)] if N == M else [(N, M), (M, N)]   # (n, m): beta^B_n beta^J_m
            for (n, m) in prods:
                sign = (-1) ** (N + M)
                a = (k + 1) * (N + M)
                b = (k + 1) * (N * N + M * M) + n * n + m * m - m
                ser = pb(n) * pj(m)
                if kind == 't':
                    pass
                else:
                    shiftq = N + M
                    if perturb == (N, M):
                        shiftq -= 1
                    epsnum = BS.from_terms(D, L, [(0, 0, 1), (1, 2 * N, 1)]) * BS.from_terms(D, L, [(0, 0, 1), (1, 2 * M, 1)])
                    if kind == 'flat':
                        ser = ser * epsnum * inv1px
                        b -= shiftq
                    else:  # omega: [(1+xq^{2N})(1+xq^{2M}) - q^{N+M}(1+x)] / x / (1+x)
                        br = epsnum - BS.from_terms(D, L, [(0, shiftq, 1), (1, shiftq, 1)])
                        ser = ser * br * inv1px
                        b -= shiftq
                        a -= 1
                tot = tot + place(ser, a, b).scale(sign)
        N += 1
    return tot


def nsum(k, perturb=False):
    """sum_n lam_n^{k+1} b^Omega_n, b^Omega_n = beta^B_n (q^{-2n} + x - 1 + x^2 q^{2n})/(x(1+x))."""
    tot = BS(D, L)
    n = 0
    while (k + 1) * n - 1 <= D and (k + 2) * n * n - 2 * n < L:
        c = 1 if (perturb and n >= 1) else -1
        # q^{-2n}(1 + (x + c) q^{2n} + x^2 q^{4n})
        fac = BS.from_terms(D, L, [(0, 0, 1), (1, 2 * n, 1), (0, 2 * n, c), (2, 4 * n, 1)])
        ser = pb(n) * fac * inv1px
        a = (k + 1) * n - 1
        b = (k + 1) * n * n + n * n - 2 * n
        tot = tot + place(ser, a, b).scale((-1) ** n)
        n += 1
    return tot


def jsum(k):
    tot = BS(D, L)
    m = 0
    while (k + 1) * m <= D and (k + 2) * m * m - m < L:
        tot = tot + place(pj(m), (k + 1) * m, (k + 2) * m * m - m).scale((-1) ** m)
        m += 1
    return tot


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    lines = []
    pr = lambda s: (print(s, flush=True), lines.append(s))
    pr("Referee 2, (2)+(5): bivariate checks mod (x^%d, q^%d) [= all powers of x mod q^%d], k = 1..%d, brute force k <= %d"
       % (D + 1, L, L, kmax, kbrute))
    allok = True

    def chk(name, A, B, expect=True):
        nonlocal allok
        fd = A.first_diff(B)
        ok = (fd is None) == expect
        allok &= ok
        if expect:
            pr("    %-6s %s" % (name, "OK" if fd is None else "FAIL first diff (x^%d q^%d: %d)" % fd))
        else:
            pr("    %-6s %s" % (name, "detected (good): first diff at x^%d q^%d" % fd[:2] if fd is not None else "NOT DETECTED (bad)"))

    C1 = C_dp(1)
    euler = bpoch_x(D, L, 1, 1, 2, 'inf', sign=-1)  # (-xq;q^2)_inf
    chk("C1=Eul", C1, euler)
    C1inv = C1.inverse()
    for k in range(1, kmax + 1):
        t1 = time.time()
        pr("  k = %d" % k)
        K, KO, Kf = K_dp(k)
        if k <= kbrute:
            Kb, KOb, Kfb, cnt = K_brute(k)
            pr("    brute force: %d tuples" % cnt)
            chk("K=bf", K, Kb)
            chk("KO=bf", KO, KOb)
            chk("Kf=bf", Kf, Kfb)
            chk("P23bf", KOb.shift(1, 0), Kfb - Kb)
        chk("P23", KO.shift(1, 0), Kf - K)
        chk("L11op", K_operator_form(k), KO)
        if k >= 2:
            chk("NC-L11", K_operator_form(k, wrong=True), KO, expect=False)
        Ck = C_dp(k)
        Jk = Ck * C1inv
        Bk1 = B_dp(k)
        Bkk1 = B_dp(k, lin_levels=(k, k + 1))
        ns = nsum(k)
        js = jsum(k)
        pref = xq_inf * x2Q_inf
        chk("T32a", pref * Kf, spectral(k, 'flat'))
        chk("T32b", pref * KO, spectral(k, 'omega'))
        chk("T32c", pref * K, spectral(k, 't'))
        chk("L41s", spectral(k, 'omega'), ns * js)
        chk("B2", x2Q_inf * Jk, js)
        chk("TH4", xq_inf * KO, Jk * ns)
        chk("TH6a", xq_inf * Bkk1, ns)
        chk("TH6b", KO * C1, Bkk1 * Ck)
        chk("Fk", K * C1, Bk1 * Ck)
        # negative controls
        chk("NC-b", xq_inf * KO, Jk * nsum(k, perturb=True), expect=False)
        if k >= 2:  # for k = 1 Omega = q^{r_1+2s_1} has no (1-q^..) factor: control vacuous
            _, KOnf, _ = K_dp(k, omega_variant='nofactor')
            chk("NC-Om", xq_inf * KOnf, Jk * ns, expect=False)
        chk("NC-eps", pref * Kf, spectral(k, 'flat', perturb=(1, 1)), expect=False)
        if k >= 2 and k <= kbrute:
            _, _, Kfb2, _ = K_brute(k, flat_level=k)
            chk("NC-fl", KO.shift(1, 0), Kfb2 - K, expect=False)
        pr("    (%.1fs)" % (time.time() - t1))
    pr("RESULT: %s   (%.1fs)" % ("ALL OK" if allok else "SOME CHECK FAILED", time.time() - t0))


if __name__ == '__main__':
    main()
