"""sw_certify_T1.py -- certificate / numerical checks for proofs/T1-BC.md and proofs/T1-A-even.md.

 (a) exact (polynomial) check of the three Bailey pairs (P1) = Slater G(1), (P2) = Slater C(1), (P3) = G(3) with
     q^{1/2} -> -q^{1/2}, for n <= NB;
 (b) exact check of the limiting Bailey step (S) [rho_1 -> infinity, rho_2 = -(aq)^{1/2}] on these pairs, n <= 6;
 (c) Theorem C  ((T1,C_r), Warnaar) and the companions (Wang-Wang Thm 4.1) as power series, r <= RC; the character
     form (Corollary C: Sun-Wang (3.6)-(3.13) for all r and all i);
 (d) Theorem B  ((T1,B_r)): even part, odd part, total, r <= RB; comparison with Du-Ji-Shen-Xu (1.7), (1.9);
     exact check that the Du-Ji-Shen-Xu matrix (1.4) is the defusion of A(T1,B_{2s})D;
 (e) Theorem A  ((T1,A_{2k})): exact fusion chain G_k -> C(A_{2k})^{-1} for k <= 10; identity (A.1) as power series;
     Lepowsky-Primc sector formula against the Kac-Peterson formula; final formula with string functions.
usage: python3 sw_certify_T1.py [N] [NB]
"""
import sys, math, time
from fractions import Fraction as Fr
from flint import fmpz_poly
from sw_lib import Ctx, Ser, compare, report
from fusion import inv, cartanA
from fusion2 import pair

N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 12
F = Fr
LOG = []


def rep(name, a, b):
    return report(name, a, b, LOG)


# ---------------------------------------------------------------- (a) Bailey pairs, exact, in p = q^{1/2}
X = fmpz_poly([0, 1])
ONE = fmpz_poly([1])


def mono(e):
    assert e >= 0
    return fmpz_poly([0] * e + [1])


def poch(e0, step, n, sign=-1):
    """prod_{j=0}^{n-1} (1 + sign * p^{e0 + j*step})   (sign=-1: (p^{e0}; p^{step})_n ; sign=+1: (-p^{e0};p^step)_n)"""
    r = ONE
    for j in range(n):
        r = r * (ONE + sign * mono(e0 + j * step))
    return r


def check_pair(name, a_exp, alpha, beta_den, beta_num=ONE, alpha_den=ONE, nmax=NB):
    """Bailey pair relative to a = q^{a_exp} = p^{2 a_exp} in base q = p^2:
       beta_n = sum_{r<=n} alpha_r / ((q;q)_{n-r} (aq;q)_{n+r});  alpha(r) -> polynomial numerator (in p),
       common denominator alpha_den;  beta_n = beta_num(n)/beta_den(n)."""
    ok = True
    for n in range(nmax + 1):
        # D = (q)_n (aq)_{2n}
        D = poch(2, 2, n) * poch(2 * a_exp + 2, 2, 2 * n)
        lhs = fmpz_poly([0])
        for r in range(n + 1):
            num = alpha(r)
            # D / ((q)_{n-r} (aq)_{n+r}) = (q^{n-r+1};q)_r * (a q^{n+r+1};q)_{n-r}
            co = poch(2 * (n - r + 1), 2, r) * poch(2 * a_exp + 2 * (n + r + 1), 2, n - r)
            lhs = lhs + num * co
        # claim: lhs/(alpha_den * D) = beta_num/beta_den
        bn = beta_num(n) if callable(beta_num) else beta_num
        good = (lhs * beta_den(n) == bn * alpha_den * D)
        ok = ok and good
    print("(a) Bailey pair %-44s exact for n <= %d : %s" % (name, nmax, ok), flush=True)
    assert ok
    return ok


def sgnp(e):
    return 1 if e % 2 == 0 else -1


# (P1) G(1): a = 1; alpha_r = (-1)^r q^{r(3r-1)/4}(1+q^{r/2}) = (-1)^r p^{r(3r-1)/2}(1+p^r);  beta_n = 1/((-p;p^2)_n (p^4;p^4)_n)
check_pair("(P1) [Slater G(1)], a=1",
           0, lambda r: ONE if r == 0 else sgnp(r) * mono(r * (3 * r - 1) // 2) * (ONE + mono(r)),
           lambda n: poch(1, 2, n, +1) * poch(4, 4, n))
# (P2) C(1): a = 1; alpha_{2m} = (-1)^m q^{3m^2-m}(1+q^{2m}), alpha_odd = 0;  beta_n = 1/((q;q^2)_n (q;q)_n)
check_pair("(P2) [Slater C(1)], a=1",
           0, lambda r: ONE if r == 0 else (fmpz_poly([0]) if r % 2 else sgnp(r // 2) * mono(2 * (3 * (r // 2) ** 2 - r // 2)) * (ONE + mono(2 * r))),
           lambda n: poch(2, 4, n) * poch(2, 2, n))
# (P3): a = q; alpha_r = (-1)^{floor((r+1)/2)} q^{(3r^2+r)/4}(1+q^{r+1/2})/(1-q);  beta_n = 1/((q^{1/2};q)_{n+1} (q^2;q^2)_n)
check_pair("(P3) [G(3), q^{1/2} -> -q^{1/2}], a=q",
           1, lambda r: sgnp((r + 1) // 2) * mono((3 * r * r + r) // 2) * (ONE + mono(2 * r + 1)),
           lambda n: poch(1, 2, n + 1) * poch(4, 4, n), alpha_den=ONE - mono(2))

# ---------------------------------------------------------------- (b) the half step (S) on Bailey pairs (exact, small n)
# (S): alpha'_n = a^{n/2} q^{n^2/2} alpha_n, beta'_n = (-(aq)^{1/2};q)_n^{-1} sum_k (-(aq)^{1/2};q)_k a^{k/2} q^{k^2/2} beta_k/(q)_{n-k}
# check that (alpha', beta') is again a Bailey pair, for (P1) (a=1) and (P3) (a=q), as power series mod p^{PS}
PS = 400


def ser(poly):
    return fmpz_poly(poly.coeffs()[:PS]) if poly.length() > PS else poly


def inv_series(poly):
    """power series inverse mod p^PS (constant term 1)"""
    c = [int(v) for v in poly.coeffs()] + [0] * PS
    assert c[0] == 1
    out = [0] * PS
    out[0] = 1
    for i in range(1, PS):
        s = 0
        for j in range(1, min(i, poly.length() - 1) + 1):
            if c[j]:
                s += c[j] * out[i - j]
        out[i] = -s
    return fmpz_poly(out)


def mul_s(a, b):
    return a.mul_low(b, PS)


for (nm, a_exp, alpha_num, alpha_den, beta_den) in [
        ("(P1)", 0, lambda r: ONE if r == 0 else sgnp(r) * mono(r * (3 * r - 1) // 2) * (ONE + mono(r)), ONE,
         lambda n: poch(1, 2, n, +1) * poch(4, 4, n)),
        ("(P3)", 1, lambda r: sgnp((r + 1) // 2) * mono((3 * r * r + r) // 2) * (ONE + mono(2 * r + 1)), ONE - mono(2),
         lambda n: poch(1, 2, n + 1) * poch(4, 4, n))]:
    okk = True
    ad_inv = inv_series(alpha_den)
    for n in range(0, 7):
        # (aq)^{1/2} = p^{a_exp+1}
        e = a_exp + 1
        bp = fmpz_poly([0])
        for k in range(n + 1):
            t = mul_s(poch(e, 2, k, +1), mono(a_exp * k + k * k))
            t = mul_s(t, inv_series(beta_den(k)))
            t = mul_s(t, inv_series(poch(2, 2, n - k)))
            bp = bp + t
        bp = mul_s(bp, inv_series(poch(e, 2, n, +1)))
        rhs = fmpz_poly([0])
        for r in range(n + 1):
            t = mul_s(mono(a_exp * r + r * r), alpha_num(r))
            t = mul_s(t, ad_inv)
            t = mul_s(t, inv_series(poch(2, 2, n - r) * poch(2 * a_exp + 2, 2, n + r)))
            rhs = rhs + t
        okk = okk and (ser(bp) == ser(rhs))
    print("(b) half step (S) applied to %s is again a Bailey pair (power series mod p^%d, n <= 6): %s" % (nm, PS, okk), flush=True)
    assert okk


def minmat(r):
    return [[min(i, j) + 1 for j in range(r)] for i in range(r)]


# ---------------------------------------------------------------- (c) Theorem C and companions
RC = 9
print("\n(c) (T1,C_r):  S_r^{(i)} = sum q^{(N_1^2+..+N_r^2)/2 + N_i+N_{i+2}+...}/((q)_{n_1}..(q)_{n_{r-1}}(q^2;q^2)_{n_r})", flush=True)
for r in range(1, RC + 1):
    NN = N if r <= 6 else max(30, N * 3 // 5)
    den = 16 * (2 * r + 3)
    c = Ctx(den, NN)
    pref = c.ns_pref()                                  # (-q^{1/2};q)_inf/(q)_inf
    d = [1] * (r - 1) + [2]
    assert pair(('T', 1), ('C', r)) == ([[F(x) for x in row] for row in minmat(r)], d) if r >= 2 else True
    for i in range(1, r + 2):
        B = [sum(1 for jj in range(i, r + 1, 2) if jj <= m) for m in range(1, r + 1)]
        j = 2 * r + 3 - 2 * i
        f = c.nahm(minmat(r), B, d)
        th = c.theta(lambda n: F((2 * r + 3) * n * n - j * n, 4), lambda n: 1 if n % 2 == 0 else -1)
        ok = rep("  r=%d i=%d B=%s: = (-q^1/2)_inf/(q)_inf * sum (-1)^n q^{((2r+3)n^2-%dn)/4}" % (r, i, B, j), f, pref * th)
        # character form
        p, pp = 4, 4 * r + 6
        Cc = F(j * j - (2 * r + 3), 16 * (2 * r + 3))
        h1 = F(j * j - 1, 16 * (2 * r + 3))
        h2 = F((4 * r + 6 + (1 if i % 2 else -1) * j) ** 2 - 1, 16 * (2 * r + 3))
        if i % 2:
            l1, l2 = (1, i), (1, 4 * r + 6 - i)
        else:
            l1, l2 = (1, 2 * r + 3 - i), (1, 2 * r + 3 + i)
        chars = c.sm(p, pp, *l1) - c.sm(p, pp, *l2)
        W = c.sm_weights(p, pp)

        def canon(l):
            return min(l, (p - l[0], pp - l[1]))
        assert W[('NS', h1)] == [canon(l1)] and W[('NS', h2)] == [canon(l2)], (W.get(('NS', h1)), l1, W.get(('NS', h2)), l2)
        rep("       = q^{-(%s)} (chi^NS_{1,%d} - chi^NS_{1,%d}) of SM(%d,4), eff. weights %s, %s" % (Cc, l1[1], l2[1], pp, h1, h2), f.shift(Cc), chars)
    # product form of Sun-Wang (2.8) for i = r+1 checked in sw_audit.py

# ---------------------------------------------------------------- (d) Theorem B
RB = 9
print("\n(d) (T1,B_r):  sum q^{(M_1^2+..+M_r^2)/4}/((q^2;q^2)_{n_1}..(q^2;q^2)_{n_{r-1}}(q)_{n_r}),  M_i = 2n_i+..+2n_{r-1}+n_r", flush=True)
for r in range(1, RB + 1):
    NN = N if r <= 6 else max(30, N * 3 // 5)
    c = Ctx(4, NN)
    if r >= 2:
        M, d = pair(('T', 1), ('B', r))
    else:
        M, d = [[F(1, 2)]], [1]
    # printed quadratic form (1/4) sum M_i^2
    Q = [[F(0)] * r for _ in range(r)]
    for i in range(r):
        v = [0] * r
        for jj in range(i, r - 1):
            v[jj] = 2
        v[r - 1] = 1
        for a in range(r):
            for b in range(r):
                Q[a][b] += F(v[a] * v[b], 2)
    assert Q == [[F(x) for x in row] for row in M], "quadratic form"
    f = c.nahm(M, [0] * r, d)
    cs = f.coeffs()
    fe = c.from_dict({F(k, 4): v for k, v in cs.items() if (k - r * 0) % 4 == 0}) if False else None
    # even / odd part in n_r: computed directly with the weights option replaced by two sums over n_r parity:
    # n_r = 2 mu: variables (n_1..n_{r-1}, mu) with 1/(q;q)_{2mu}: not a Nahm sum; separate by the exponent class instead:
    # exponent = sum_{j<r}(N'_j + n_r/2)^2 + n_r^2/4 = r n_r^2/4 mod 1  -> classes n_r even: integer ; n_r odd: r/4 mod 1 (for r = 0 mod 4 both integer)
    pe = c.prod(lambda n: (1 if n % 4 == 2 else 0) - (1 if n % 2 == 1 else 0) - (1 if n % 2 == 0 else 0))     # (-q;q^2)/(q^2;q^2)
    po = c.prod(lambda n: (1 if n % 4 == 0 else 0) - 2 * (1 if n % 2 == 0 else 0))                             # (-q^2;q^2)/(q^2;q^2)
    te = c.theta(lambda m: (4 * r + 6) * m * m - 2 * m, lambda m: 1 if m % 2 == 0 else -1)
    to = c.theta(lambda n: F((2 * r + 3) * n * n + (2 * r + 1) * n, 2), lambda n: 1 if ((n + 1) // 2) % 2 == 0 else -1).shift(F(r, 4))
    even, odd = pe * te, po * to
    rep("  r=%d: total = even + odd" % r, f, even + odd)
    if r % 4:
        fe = c.from_dict({F(k, 4): v for k, v in cs.items() if k % 4 == 0})
        fo = c.from_dict({F(k, 4): v for k, v in cs.items() if k % 4 != 0})
        rep("       even part (n_r even) = (-q;q^2)/(q^2;q^2) * sum (-1)^m q^{(4r+6)m^2-2m}", fe, even)
        rep("       odd part  (n_r odd)", fo, odd)
    # theta-product form of the odd part
    m8 = 8 * r + 12
    o2 = (c.prod(lambda n: 1 if n % m8 in (0, 2 * r + 5, 6 * r + 7) else 0) + c.prod(lambda n: 1 if n % m8 in (0, 2 * r + 1, 6 * r + 11) else 0).shift(1)).shift(F(r, 4))
    rep("       odd = q^{r/4}(-q^2;q^2)/(q^2;q^2) [(q^{2r+5},q^{6r+7},q^{8r+12};q^{8r+12}) + q (q^{2r+1},q^{6r+11},q^{8r+12};q^{8r+12})]", odd, po * o2)

# Du-Ji-Shen-Xu arXiv:2602.14483, Theorem 1.5 (1) [(1.7)] and (3) [(1.9)]
print("\n    Du-Ji-Shen-Xu (1.7), (1.9) versus (T1,B_{2R-2}), (T1,B_{2R-3}):", flush=True)


def djsx_AD(R, last):
    M = [[F(0)] * R for _ in range(R)]
    M[0][0] = F(2); M[0][R - 1] = M[R - 1][0] = F(1)
    for i in range(1, R - 1):
        for jj in range(1, R - 1):
            M[i][jj] = F(4 * min(i, jj))
        M[i][R - 1] = M[R - 1][i] = F(2 * i)
    M[R - 1][R - 1] = F(last)
    return M


for R in range(2, 6):
    c = Ctx(4, min(N, 60))
    d = [2] * (R - 1) + [1]
    # (1.7)
    fA = c.nahm(djsx_AD(R, R - 1), [0] * R, d)
    m16 = 16 * R - 4
    rhs = c.prod(lambda n: (1 if n % m16 in (0, 8 * R - 4, 8 * R) else 0) - (1 if n % 4 in (0, 1, 3) else 0))
    # second term: q^{(R-1)/2} (-q, q^{4R-2}, -q^{4R-1}; -q^{4R-1})_inf / (q^2,q^2,q^4;q^4)_inf ; triple product with base Q = -q^{4R-1}:
    #   (z, Q/z, Q; Q)_inf = sum_n (-1)^n z^n Q^{n(n-1)/2},  z = -q
    th2 = c.theta(lambda n: n + (4 * R - 1) * n * (n - 1) // 2, lambda n: 1 if (n * (n - 1) // 2) % 2 == 0 else -1)
    rhs2 = (th2 * c.prod(lambda n: -2 * (1 if n % 4 == 2 else 0) - (1 if n % 4 == 0 else 0))).shift(F(R - 1, 2))
    rep("    DJSX (1.7) R=%d (as printed)" % R, fA, rhs + rhs2)
    s = 2 * R - 2
    g = c.nahm(*((pair(('T', 1), ('B', s))[0], [0] * s, pair(('T', 1), ('B', s))[1]))) if s >= 2 else None
    if g is not None:
        rep("    (T1,B_%d) = DJSX A-family, R=%d, b_0" % (s, R), g, fA)
    # (1.9)
    fB = c.nahm(djsx_AD(R, F(2 * R - 1, 2)), [0] * (R - 1) + [F(-1, 2)], d)
    m16 = 16 * R - 12
    rhs = c.prod(lambda n: (1 if n % m16 in (0, 8 * R - 8, 8 * R - 4) else 0) - (1 if n % 4 in (0, 1, 3) else 0))
    th2 = c.theta(lambda n: n + (4 * R - 3) * n * (n - 1) // 2, lambda n: 1 if (n * (n - 1) // 2) % 2 == 0 else -1)
    rhs2 = (th2 * c.prod(lambda n: -2 * (1 if n % 4 == 2 else 0) - (1 if n % 4 == 0 else 0))).shift(F(2 * R - 3, 4))
    rep("    DJSX (1.9) R=%d (as printed)" % R, fB, rhs + rhs2)
    s = 2 * R - 3
    if s >= 2:
        Ms, ds = pair(('T', 1), ('B', s))
        rep("    (T1,B_%d) = DJSX B-family, R=%d, b_0" % (s, R), c.nahm(Ms, [0] * s, ds), fB)
    else:
        rep("    (T1,B_1) = sum q^{n^2/4}/(q)_n = DJSX B-family, R=2, b_0", c.nahm([[F(1, 2)]], [0]), fB)

# exact: DJSX matrix (1.4) of rank s+1 fuses to A(T1,B_{2s})D:  node 2j+1 := node 1 + node 2j (j = s-1,...,1), base q^2
okf = True
for s in range(2, 9):
    r = 2 * s
    M, d = pair(('T', 1), ('B', r))
    names = [1] + [2 * j for j in range(1, s)] + [r]
    G = djsx_AD(s + 1, s)
    A = {x: {y: G[i][jj] for jj, y in enumerate(names)} for i, x in enumerate(names)}
    for j in range(s - 1, 0, -1):
        a, b, new = 1, 2 * j, 2 * j + 1
        row = {kk: A[a][kk] + A[b][kk] for kk in names}
        diag = A[a][a] + A[b][b] + 2 * A[a][b]
        A[a][b] += 2; A[b][a] += 2           # base q^2: raise by d = 2
        for kk in names:
            A[kk][new] = row[kk]
        A[new] = dict(row); A[new][new] = diag
        names.append(new)
    okf = okf and all(A[i + 1][jj + 1] == M[i][jj] for i in range(r) for jj in range(r))
print("    exact: DJSX matrix (1.4) of rank s+1 --(s-1 fusions in base q^2)--> A(T1,B_{2s})D, s = 2..8 :", okf, flush=True)
assert okf

# ---------------------------------------------------------------- (e) Theorem A: (T1,A_{2k})
print("\n(e) (T1,A_{2k})", flush=True)


def Gk(k):
    """rank k+1 matrix on (b, m_1..m_{k-1}, a) <-> nodes (1, 2, 4, ..., 2k-2, 2k) of A_{2k}"""
    n = 2 * k + 1
    names = [1] + [2 * j for j in range(1, k)] + [2 * k]
    Ci = inv(cartanA(2 * k))
    A = {x: {y: Ci[x - 1][y - 1] for y in names} for x in names}
    for j in range(1, k):
        A[1][2 * j] -= 1; A[2 * j][1] -= 1
    return names, A


okf = True
for k in range(1, 11):
    names, A = Gk(k)
    Ci = inv(cartanA(2 * k))
    # explicit form: 2*Q = 2ab + 2k(a-b)^2/(2k+1) + 4L(a-b)/(2k+1) + 2 sum M_j^2 - 4L^2/(2k+1)
    for x in names:
        for y in names:
            def vec(z):
                # (a, b, m_1..m_{k-1}) coordinates of node z
                v = {'a': 0, 'b': 0}
                if z == 1: v['b'] = 1
                elif z == 2 * k: v['a'] = 1
                else: v[z // 2] = 1
                return v
            vx, vy = vec(x), vec(y)

            def Qf(v):
                a, b = v.get('a', 0), v.get('b', 0)
                m = [v.get(j, 0) for j in range(1, k)]
                L = sum((j + 1) * m[j] for j in range(k - 1))
                Ms = [sum(m[j:]) for j in range(k - 1)]
                return F(a * b) + F(k * (a - b) ** 2 + 2 * L * (a - b) - 2 * L * L, 2 * k + 1) + sum(x_ * x_ for x_ in Ms)
            vs = {kk: vx.get(kk, 0) + vy.get(kk, 0) for kk in set(vx) | set(vy)}
            bil = Qf(vs) - Qf(vx) - Qf(vy)
            okf = okf and (bil == A[x][y])
    for j in range(k - 1, 0, -1):
        a, b, new = 1, 2 * j, 2 * j + 1
        row = {kk: A[a][kk] + A[b][kk] for kk in names}
        diag = A[a][a] + A[b][b] + 2 * A[a][b]
        A[a][b] += 1; A[b][a] += 1
        for kk in names:
            A[kk][new] = row[kk]
        A[new] = dict(row); A[new][new] = diag
        names.append(new)
    okf = okf and all(A[i + 1][jj + 1] == Ci[i][jj] for i in range(2 * k) for jj in range(2 * k))
print("    exact: G_k = matrix of ab + (k(a-b)^2 + 2L(a-b) - 2L^2)/(2k+1) + sum M_j^2, and G_k --(k-1 fusions)--> C(A_{2k})^{-1}, k = 1..10 :", okf, flush=True)
assert okf

for k in range(1, 6):
    r = 2 * k
    K = k * (2 * k + 1)
    cc = F(3 * r, r + 4)
    den = math.lcm(2 * (r + 1), (cc / 24).denominator, 24 * (k + 2), 4 * k, K, 24)
    NN = min(N, 60) if k <= 3 else (30 if k == 4 else 20)
    c = Ctx(den, NN)
    f = c.nahm(inv(cartanA(r)), [0] * r)
    # (A.1): (1/(q)_inf) sum_{t in Z} sum_{m in N^{k-1}} q^{sum M_j^2 + (k t^2 + 2 L t - 2 L^2)/(2k+1)} / prod (q)_{m_j}
    # computed as sum over residue classes of L: for each m the t-sum is a theta series; use Nahm evaluator with weights replaced
    # by direct evaluation:  sum_t q^{(k t^2 + 2Lt)/(2k+1)} depends on L mod (2k+1)... (2Lt -> t -> t - s shifts L by s(…)); evaluate by brute force:
    Ck = inv(cartanA(k - 1)) if k >= 2 else []
    # PF sector sums:  P_Q = sum_{m: L = Q mod k} q^{m C^{-1} m}/(q)_m   (L = sum j m_j)
    tot = 0
    lp_ok = True
    for Qc in range(k):
        if k >= 2:
            # sector via a roots-of-unity-free method: restrict by enumerating with weights
            from nahm_flint import nahm as nahm_w
            # use Ctx evaluator on the lattice m, grouping by L mod k through exponent classes: exponent = sum M^2 - L^2/k = L(k-L)/k mod 1
            pass
        # theta: sum_{N = Q mod k} q^{N^2/K}
        th = c.theta(lambda n: F((k * n + Qc) ** 2, K))
        tot = tot + c.pf(k, 0, (2 * Qc) % (2 * k)) * th
    rhs = (tot * c.euler_inv()).shift(F(-1, 24))
    rep("    k=%d: q^{-c/24} Nahm(C(A_%d)^{-1}) = sum_{Q mod k} c^0_{2Q}(tau) sum_{N = Q (k)} q^{N^2/(k(2k+1))},  c=%s" % (k, r, cc), f.shift(-cc / 24), rhs)
    # (A.1) directly, independent of string functions: rank k+1 Nahm sum with matrix G_k
    names, A = Gk(k)
    Gm = [[A[x][y] for y in names] for x in names]
    g = c.nahm(Gm, [0] * (k + 1))
    rep("          Nahm(C(A_%d)^{-1}) = Nahm(G_%d) [rank %d]" % (r, k, k + 1), f, g)

# odd rank (not solved): exact check of the rank-k terminal form of C(A_{2k-1})^{-1}
#   Q'(b,m) = b^2/2 + sum M_j^2 - (L + b/2)^2/k   on nodes (1, 2, 4, ..., 2k-2), fusions node 2j+1 := node 1 + node 2j (j = k-1..1)
okf = True
for k in range(2, 9):
    r = 2 * k - 1
    Ci = inv(cartanA(r))
    names = [1] + [2 * j for j in range(1, k)]
    A = {x: {y: Ci[x - 1][y - 1] for y in names} for x in names}
    for j in range(1, k):
        A[1][2 * j] -= 1; A[2 * j][1] -= 1

    def Qp(v):
        b = v.get('b', 0)
        m = [v.get(j, 0) for j in range(1, k)]
        L = sum((j + 1) * m[j] for j in range(k - 1))
        Ms = [sum(m[j:]) for j in range(k - 1)]
        return F(b * b, 2) + sum(x_ * x_ for x_ in Ms) - F((2 * L + b) ** 2, 4 * k)

    def vec(z):
        return {'b': 1} if z == 1 else {z // 2: 1}
    for x in names:
        for y in names:
            vx, vy = vec(x), vec(y)
            vs = {kk: vx.get(kk, 0) + vy.get(kk, 0) for kk in set(vx) | set(vy)}
            okf = okf and (Qp(vs) - Qp(vx) - Qp(vy) == A[x][y])
    for j in range(k - 1, 0, -1):
        a, b, new = 1, 2 * j, 2 * j + 1
        row = {kk: A[a][kk] + A[b][kk] for kk in names}
        diag = A[a][a] + A[b][b] + 2 * A[a][b]
        A[a][b] += 1; A[b][a] += 1
        for kk in names:
            A[kk][new] = row[kk]
        A[new] = dict(row); A[new][new] = diag
        names.append(new)
    okf = okf and all(A[i + 1][jj + 1] == Ci[i][jj] for i in range(r) for jj in range(r))
print("    exact (odd rank, for the record): matrix of b^2/2 + sum M_j^2 - (L+b/2)^2/k --(k-1 fusions)--> C(A_{2k-1})^{-1}, k = 2..8 :", okf, flush=True)
assert okf

# Lepowsky-Primc, sector by sector, against Kac-Peterson:
#   sum_{m in N^{k-1}, L = Q mod k} q^{m C(A_{k-1})^{-1} m - c_PF/24}/(q)_m = eta(tau) c^0_{2Q}(tau),   L = sum j m_j
print("    Lepowsky-Primc sector formula vs Kac-Peterson string functions:", flush=True)
from nahm_flint import nahm as nahm_w
for k in range(2, 8):
    cpf = F(2 * (k - 1), k + 2)
    den = math.lcm(24 * (k + 2) * k, 2 * k)
    NN = min(N, 60) if k <= 5 else 30
    c = Ctx(den, NN)
    xs = nahm_w([[2 * x for x in row] for row in inv(cartanA(k - 1))], [0] * (k - 1), [1] * (k - 1), NN, k, weights=list(range(1, k)))
    okk = True
    for Qc in range(k):
        d_ = {}
        for L, arr in xs.items():
            if L % k == Qc:
                for e, v in enumerate(arr):
                    if v:
                        d_[F(e, k) - cpf / 24] = d_.get(F(e, k) - cpf / 24, 0) + v
        part = c.from_dict(d_)
        ok_, order, bad = compare(part, c.pf(k, 0, 2 * Qc))
        okk = okk and ok_
    print("      k=%d: for every Q mod k:  sum_{L = Q (k)} q^{mC^{-1}m - c/24}/(q)_m = eta c^0_{2Q}   mod q^%d : %s" % (k, NN - 1, okk), flush=True)
    LOG.append(("LP k=%d" % k, okk, NN, None))

nbad = sum(1 for x in LOG if not x[1])
print("\nSUMMARY: %d power-series comparisons, %d failed" % (len(LOG), nbad))
