"""Certificate for proofs/KRTW-general.md  (first KRTW Andrews-Gordon-type identity, all k).

Checks, in the order of the write-up (every check is exact: polynomial identities in Z[x,q] with python-flint
fmpz_mpoly, or power series in q with coefficients in Z[x] modulo q^N):

 [L1]  Lemma 1  (q-binomial theorem, used form):  sum_i [t;i] (-1)^i q^{i(i-1)/2} z^i = (z;q)_t              exact
 [L2]  Lemma 2:  sum_{t+2b=L} q^{t(t-1)/2} [R;t]/(q^2;q^2)_b = sum_{j=0}^L (-1)^j q^{Rj}/((q)_j (q)_{L-j})        exact
 [L3]  Lemma 3:  sum_r [R;r] y^r q^{r^2} (y q^{d+r+1};q)_{R-r} = sum_t [R;t][d;t](q)_t y^t q^{t^2}              exact
 [L4]  Lemma 4 (the key finite identity (*)), denominators cleared, 0 <= R, M <= RM                         exact
 [L4s] Lemma 4 with the infinite products (the form (**) used in the induction), and its M = infinity form   mod q^N
 [T1]  Theorem 1 (R-refined identity): left side summed directly from the definition (brute force over
       r_1>=...>=r_j, s_1<=...<=s_j) against the right side, all powers of x                                 mod q^N
 [T2]  Theorem 2 (E_k): KRTW_k(x) = Mel_k(x), left side by brute force (k<=3) and by the level-by-level
       transfer recursion with x (k<=KT), all powers of x                                                    mod q^N
 [C5]  Corollary 5 (extra variable y marking r_1) at y = 0, -1, q, q^2, q^3, -q, k <= 3                       mod q^N
 [B]   Bailey step (S2) with a = x^2 q: the pairs (alpha^(i), beta^(i)) of section 4 are Bailey pairs
       (checked for n <= 4, i <= 3), and Proposition 6: Mel_k(x) = very-well-poised single sum               mod q^N
 [T3]  Theorem 3: Mel_k(1) = (q,q^{k+1},q^{k+2};q^{k+2})_inf/((q;q^2)_inf (q)_inf) (k <= 20, mod q^NP) and the
       2k-fold KRTW sum (i) against the same product (transfer recursion, k <= 8, mod q^NP)

 [X1]  numerical remark (section 8): x = 1/q, sum of the k+1 theta quotients (k <= 6, mod q^80)
 [X2]  the second product form of (iii) quoted in section 8 (k <= 8, mod q^200)

usage: python3 krtw_general_certify.py            (about 15 seconds)
"""
import sys, itertools, time
from flint import fmpz_poly, fmpz_mpoly_ctx, Ordering
import krtw_ser as K
from krtw_ser import S, tr, mono, iqpoch, qbin, xpoch_inf, xpoch, xpoch_inv
from krtw_four_test import multisum, prod_side, first_diff

t0 = time.time()
OK = True
def report(name, ok, extra=""):
    global OK
    OK &= bool(ok)
    print("[%s] %s  %s" % (name, "OK " if ok else "FAIL", extra), flush=True)

# ---------------------------------------------------------------- exact polynomial arithmetic in Z[x,q]
ctx = fmpz_mpoly_ctx.get(("x", "q"), Ordering.lex)
X, Q = ctx.gens()
ONE = ctx.constant(1) if hasattr(ctx, "constant") else X ** 0
ZERO = ONE - ONE

def P(n, base=1):
    """(q^base;q^base)_n"""
    r = ONE
    for i in range(1, n + 1):
        r = r * (1 - Q ** (base * i))
    return r
_gb = {}
def gbin(n, k):
    """Gaussian binomial [n;k]_q as polynomial (q-Pascal recurrence)"""
    if k < 0 or k > n:
        return ZERO
    if k == 0 or k == n:
        return ONE
    if (n, k) not in _gb:
        _gb[(n, k)] = gbin(n - 1, k - 1) + Q ** k * gbin(n - 1, k)
    return _gb[(n, k)]
def poch(a, n):
    """(a;q)_n for a polynomial a"""
    r = ONE
    for i in range(n):
        r = r * (1 - a * Q ** i)
    return r
def Prange(lo, hi, base=1):
    """prod_{i=lo}^{hi} (1-q^{base*i})"""
    r = ONE
    for i in range(lo, hi + 1):
        r = r * (1 - Q ** (base * i))
    return r

# [L1]
ok = True
for t in range(0, 9):
    lhs = ZERO
    for i in range(t + 1):
        lhs = lhs + gbin(t, i) * (-1) ** i * Q ** (i * (i - 1) // 2) * X ** i
    ok &= (lhs == poch(X, t))
report("L1", ok, "q-binomial theorem, t <= 8 (z = x)")

# [L2]  multiply by D = (q)_L (q^2;q^2)_{floor(L/2)}
ok = True
for R in range(0, 9):
    for L in range(0, 11):
        h = L // 2
        lhs = ZERO
        for b in range(0, h + 1):
            t = L - 2 * b
            if t > R:
                continue
            lhs = lhs + Q ** (t * (t - 1) // 2) * gbin(R, t) * P(L) * Prange(b + 1, h, 2)
        rhs = ZERO
        for j in range(0, L + 1):
            rhs = rhs + (-1) ** j * Q ** (R * j) * gbin(L, j) * P(h, 2)
        ok &= (lhs == rhs)
report("L2", ok, "0 <= R <= 8, 0 <= L <= 10, exact")

# [L3]
ok = True
for R in range(0, 8):
    for d in range(0, 8):
        lhs = ZERO
        for r in range(R + 1):
            lhs = lhs + gbin(R, r) * X ** r * Q ** (r * r) * poch(X * Q ** (d + r + 1), R - r)
        rhs = ZERO
        for t in range(0, min(R, d) + 1):
            rhs = rhs + gbin(R, t) * gbin(d, t) * P(t) * X ** t * Q ** (t * t)
        ok &= (lhs == rhs)
report("L3", ok, "0 <= R, d <= 7, exact (y = x)")

# [L4]  (*) multiplied by D = (q)_M (q^2;q^2)_{floor(M/2)}
RM = 7
ok = True
for R in range(0, RM + 1):
    for M in range(0, RM + 1):
        h = M // 2
        lhs = ZERO
        for r in range(R + 1):
            for b in range(h + 1):
                n = r + 2 * b
                lhs = lhs + (gbin(R, r) * X ** n * Q ** (n * (n + 1) // 2 + r * (r - 1) // 2)
                             * poch(X * Q ** (r + M + 1), R - r) * Prange(M - 2 * b + 1, M) * Prange(b + 1, h, 2))
        rhs = ZERO
        for L in range(M + 1):
            rhs = rhs + X ** L * Q ** (L * (L + 1) // 2) * poch(X * Q ** (R + L + 1), M - L) * gbin(M, L) * P(h, 2)
        ok &= (lhs == rhs)
report("L4", ok, "key identity (*), 0 <= R, M <= %d, exact in Z[x,q]" % RM)

# ---------------------------------------------------------------- series checks
N = 36
K.setN(N)

def lhs_star(R, M):
    """left side of (**); M = None means M = infinity (factor 1/(q)_{M-2b} and the product dropped as in Lemma 4')"""
    s = S()
    for r in range(R + 1):
        b = 0
        while True:
            if M is not None and 2 * b > M:
                break
            n = r + 2 * b
            e = n * (n + 1) // 2 + r * (r - 1) // 2
            if e >= N:
                break
            c = tr(qbin(R, r) * iqpoch(b, 2))
            if M is None:
                s = s + S.xq(n, e) * c
            else:
                s = s + (xpoch_inf(-1, 1, r + M + 1) * tr(c * iqpoch(M - 2 * b))).shift(n, e)
            b += 1
    return s
def rhs_star(R, M):
    s = S()
    L = 0
    while True:
        if M is not None and L > M:
            break
        e = L * (L + 1) // 2
        if e >= N:
            break
        c = iqpoch(L) if M is None else tr(iqpoch(M - L) * iqpoch(L))
        s = s + (xpoch_inf(-1, 1, R + L + 1) * c).shift(L, e)
        L += 1
    return s
ok = all(lhs_star(R, M) == rhs_star(R, M) for R in range(0, 9) for M in list(range(0, 9)) + [None])
report("L4s", ok, "(**) for 0 <= R, M <= 8 and M = infinity, all powers of x, mod q^%d" % N)

def chain(j, A):
    """sum_{M_1>=...>=M_j>=0} x^{sum M} q^{sum M_i(M_i+1)/2} A(M_j) / prod_{i=1}^j (q)_{M_i-M_{i+1}}  (M_{j+1}=0)"""
    s = S()
    def rec(Ms, e, c):
        nonlocal s
        lev = len(Ms)
        if lev == j:
            s = s + (A(Ms[-1]) * tr(c * iqpoch(Ms[-1]))).shift(sum(Ms), e)
            return
        M = 0
        while lev == 0 or M <= Ms[-1]:
            ee = e + M * (M + 1) // 2
            if ee >= N:
                break
            cc = c if lev == 0 else tr(c * iqpoch(Ms[-1] - M))
            rec(Ms + [M], ee, cc)
            M += 1
    if j == 0:
        return A(None)
    rec([], 0, fmpz_poly([1]))
    return s

def brute_T1(j, R):
    """(q)_R * left side of Theorem 1, summed directly from the definition; R = None: no bound on r_1 and no
    factor (q)_R/(q)_{R-r_1} (= the KRTW sum of Theorem 2)."""
    s = S()
    rmax = 0
    while rmax * rmax < N:
        rmax += 1
    top = rmax if R is None else min(R, rmax)
    def rec_s(rs, ss, e, c):
        nonlocal s
        i = len(ss)
        if i == j:
            s = s + S.xq(sum(rs) + 2 * sum(ss), e) * c
            return
        sv = ss[-1] if ss else 0
        while True:
            ee = e + (rs[i] + sv) ** 2 + sv * sv + sv
            if ee >= N:
                break
            rec_s(rs, ss + [sv], ee, tr(c * iqpoch(sv - (ss[-1] if ss else 0), 2)))
            sv += 1
    for rs in itertools.product(range(top + 1), repeat=j):
        if any(rs[i] < rs[i + 1] for i in range(j - 1)):
            continue
        if sum(r * r for r in rs) >= N:
            continue
        c = fmpz_poly([1])
        for i in range(j):
            c = tr(c * iqpoch(rs[i] - (rs[i + 1] if i < j - 1 else 0)))
        if R is not None:
            c = tr(c * tr(K.qpoch(R) * iqpoch(R - rs[0])))
        rec_s(list(rs), [], 0, c)
    return s

ok = True
cases = [(1, R) for R in range(0, 7)] + [(2, R) for R in range(0, 6)] + [(3, R) for R in range(0, 5)] + [(4, R) for R in range(0, 4)]
for (j, R) in cases:
    rhs = chain(j, lambda Mj: xpoch_inf(-1, 1, R + Mj + 1))
    ok &= (brute_T1(j, R) == rhs)
report("T1", ok, "Theorem 1 for (j,R) in j<=4, R<=6,5,4,3 (brute force), all powers of x, mod q^%d" % N)

def transfer_x(k):
    """KRTW_k(x) by the level-by-level recursion (definition), with x^{sum (r_i+2 s_i)}"""
    rmax = 0
    while rmax * rmax < N:
        rmax += 1
    Phi = {(r, 0): S.one() for r in range(rmax + 1)}
    for i in range(1, k + 1):
        new = {}
        for s in range(rmax + 1):
            g = {}
            for r in range(rmax + 1):
                e = (r + s) ** 2 + s * s + s
                if e >= N:
                    break
                acc = S()
                for sp in range(s + 1):
                    p = Phi.get((r, sp))
                    if p is not None:
                        acc = acc + p * iqpoch(s - sp, 2)
                if not acc.iszero():
                    g[r] = acc.shift(r + 2 * s, e)
            if not g:
                continue
            for R in range(rmax + 1):
                acc = S()
                for r, v in g.items():
                    if r >= R:
                        acc = acc + v * iqpoch(r - R)
                if not acc.iszero():
                    new[(R, s)] = acc
        Phi = new
    tot = S()
    for (R, s), v in Phi.items():
        if R == 0:
            tot = tot + v
    return tot

def Mel(k):
    return chain(k, lambda Mj: S.one())

ok = True
for k in (1, 2, 3):
    ok &= (brute_T1(k, None) == Mel(k))
report("T2a", ok, "Theorem 2, k <= 3, left side by brute force, all powers of x, mod q^%d" % N)
KT = 7
ok = True
mel = {}
for k in range(1, KT + 1):
    mel[k] = Mel(k)
    ok &= (transfer_x(k) == mel[k])
report("T2b", ok, "Theorem 2, k <= %d, left side by transfer recursion, all powers of x, mod q^%d" % (KT, N))

# [C5] Corollary 5 (y marks r_1), checked at y = 0, -1, q, q^2, q^3
def brute_y(k, ysign, ypow):
    """sum y^{r_1} x^{sum N_i} W^{(i)},  y = ysign * q^ypow  (ysign in {0,1,-1})"""
    s = S()
    rmax = 0
    while rmax * rmax < N:
        rmax += 1
    def rec_s(rs, ss, e, c):
        nonlocal s
        i = len(ss)
        if i == k:
            s = s + S.xq(sum(rs) + 2 * sum(ss), e) * c
            return
        sv = ss[-1] if ss else 0
        while True:
            ee = e + (rs[i] + sv) ** 2 + sv * sv + sv
            if ee >= N:
                break
            rec_s(rs, ss + [sv], ee, tr(c * iqpoch(sv - (ss[-1] if ss else 0), 2)))
            sv += 1
    for rs in itertools.product(range(rmax + 1), repeat=k):
        if any(rs[i] < rs[i + 1] for i in range(k - 1)) or sum(r * r for r in rs) >= N:
            continue
        if ysign == 0 and rs[0] > 0:
            continue
        c = fmpz_poly([1])
        for i in range(k):
            c = tr(c * iqpoch(rs[i] - (rs[i + 1] if i < k - 1 else 0)))
        if ysign != 0:
            c = tr(c * mono(ypow * rs[0], ysign ** rs[0]))
        rec_s(list(rs), [], 0, c)
    return s
def rhs_y(k, ysign, ypow):
    def A(Mk):
        t = S()
        l = 0
        yl = fmpz_poly([1])          # (y;q)_l
        while True:
            e = l * (l + 1) // 2 + Mk * l
            if e >= N:
                break
            t = t + S.xq(l, e, (-1) ** l) * tr(yl * iqpoch(l))
            # (y;q)_{l+1} = (y;q)_l (1 - y q^l)
            yl = tr(yl * (fmpz_poly([1]) - (mono(ypow + l, ysign) if ysign != 0 else fmpz_poly([0]))))
            l += 1
        return t
    return chain(k, A)
ok = True
for k in (1, 2, 3):
    for (ys, yp) in ((0, 0), (-1, 0), (1, 1), (1, 2), (1, 3), (-1, 1)):
        ok &= (brute_y(k, ys, yp) == rhs_y(k, ys, yp))
report("C5", ok, "Corollary 5 (y^{r_1}-refinement) for k <= 3 at y = 0, -1, q, q^2, q^3, -q, all powers of x, mod q^%d" % N)

# [B] Bailey pairs relative to a = x^2 q
def alpha(i, r):
    """alpha^(i)_r = x^{ir} q^{i r(r+1)/2} * (-1)^r q^{r(r-1)/2} (x^2 q;q)_r (1-x^2 q^{2r+1}) / ((q)_r (1-x^2 q))"""
    e = i * r * (r + 1) // 2 + r * (r - 1) // 2
    if e >= N:
        return S()
    if r == 0:
        return S.one()
    # (x^2 q)_r/(1-x^2 q) = (x^2 q^2;q)_{r-1}
    a = xpoch(-1, 2, 2, r - 1) * (S.one() - S.xq(2, 2 * r + 1))
    return (a * iqpoch(r)).shift(i * r, e) * ((-1) ** r)
def beta(i, n):
    """beta^(i)_n = (1/(-xq)_n) sum_{n>=M_2>=...>=M_i>=0} prod_{l=2}^i x^{M_l} q^{M_l(M_l+1)/2} / ((q)_{n-M_2}...(q)_{M_i}),
    beta^(0)_n = delta_{n,0}, beta^(1)_n = 1/((q)_n (-xq)_n)"""
    if i == 0:
        return S.one() if n == 0 else S()
    s = S()
    def rec(Ms, e, c):
        nonlocal s
        if len(Ms) == i:
            s = s + S.xq(sum(Ms[1:]), e) * tr(c * iqpoch(Ms[-1]))
            return
        M = 0
        while M <= Ms[-1]:
            ee = e + M * (M + 1) // 2
            if ee >= N:
                break
            rec(Ms + [M], ee, tr(c * iqpoch(Ms[-1] - M)))
            M += 1
    rec([n], 0, fmpz_poly([1]))
    return s * xpoch_inv(1, 1, 1, n)
ok = True
for i in range(0, 4):
    for n in range(0, 5):
        rhs = S()
        for r in range(0, n + 1):
            # 1/((q)_{n-r} (x^2 q^2;q)_{n+r})
            rhs = rhs + alpha(i, r) * xpoch_inv(-1, 2, 2, n + r) * iqpoch(n - r)
        ok &= (beta(i, n) == rhs)
report("B1", ok, "(alpha^(i),beta^(i)) is a Bailey pair relative to a = x^2 q for i <= 3, n <= 4, mod q^%d" % N)

def vwp(k):
    s = S()
    r = 0
    while True:
        a = alpha(k + 1, r)       # x^{(k+1)r} q^{(k+1)r(r+1)/2} alpha_r
        if (k + 1) * r * (r + 1) // 2 + r * (r - 1) // 2 >= N:
            break
        s = s + a
        r += 1
    return s * xpoch_inf(1, 1, 1) * xpoch_inv(-1, 2, 2, N)
ok = True
for k in range(0, KT + 1):
    m = mel[k] if k >= 1 else S.one()
    ok &= (m == vwp(k))
report("B2", ok, "Proposition 6: Mel_k(x) = very-well-poised sum, 0 <= k <= %d, all powers of x, mod q^%d" % (KT, N))

# [T3]
NP = 300
K.setN(NP)
N = NP
def mel1(k):
    """Mel_k(1) mod q^NP by the chain recursion in one variable: f_i(M) = sum_{M'<=M} q^{M'(M'+1)/2} f_{i-1}(M')/(q)_{M-M'}"""
    mmax = 0
    while mmax * (mmax + 1) // 2 < NP:
        mmax += 1
    f = [iqpoch(M) for M in range(mmax + 1)]          # innermost: 1/(q)_{M_k}
    for i in range(k - 1):
        g = []
        for M in range(mmax + 1):
            acc = fmpz_poly([0])
            for Mp in range(M + 1):
                e = Mp * (Mp + 1) // 2
                if e >= NP:
                    break
                acc = acc + tr(tr(f[Mp] * mono(e)) * iqpoch(M - Mp))
            g.append(acc)
        f = g
    tot = fmpz_poly([0])
    for M in range(mmax + 1):
        e = M * (M + 1) // 2
        if e < NP:
            tot = tot + tr(f[M] * mono(e))
    return tot
ok = True
for k in range(1, 21):
    pr = prod_side([1, k + 1, k + 2], k + 2, NP, [(1, 2), (1, 1)])
    ok &= (first_diff(mel1(k), pr, NP) is None)
report("T3a", ok, "Mel_k(1) = (q,q^{k+1},q^{k+2};q^{k+2})_inf/((q;q^2)_inf (q)_inf), 1 <= k <= 20, mod q^%d" % NP)
ok = True
for k in range(1, 9):
    pr = prod_side([1, k + 1, k + 2], k + 2, NP, [(1, 2), (1, 1)])
    ok &= (first_diff(multisum(k, (0, 1), NP), pr, NP) is None)
report("T3b", ok, "KRTW identity (i): 2k-fold sum = product, 1 <= k <= 8, mod q^%d" % NP)

# [X1] remark of section 8: x = 1/q  (numerical observation; not used in any proof)
NX = 80
K.setN(NX)
N = NX
def mel_lin(k, lin):
    """sum_{M_1>=...>=M_k>=0} q^{sum M_i(M_i+1)/2 + lin*sum M_i}/prod (q)_{M_i-M_{i+1}}   (lin = -1: x = 1/q)"""
    mmax = 2
    while mmax * (mmax - 1) // 2 < NX:
        mmax += 1
    f = [iqpoch(M) for M in range(mmax + 1)]
    for i in range(k - 1):
        g = []
        for M in range(mmax + 1):
            acc = fmpz_poly([0])
            for Mp in range(M + 1):
                e = Mp * (Mp + 1) // 2 + lin * Mp
                if e < NX:
                    acc = acc + tr(tr(f[Mp] * mono(e)) * iqpoch(M - Mp))
            g.append(acc)
        f = g
    tot = fmpz_poly([0])
    for M in range(mmax + 1):
        e = M * (M + 1) // 2 + lin * M
        if e < NX:
            tot = tot + tr(f[M] * mono(e))
    return tot
ok = True
for k in range(1, 7):
    th = fmpz_poly([0])
    for i in range(1, k + 2):
        th = th + prod_side([i, k + 2 - i, k + 2], k + 2, NX, [])
    rhs = tr(th * prod_side([], 1, NX, [(1, 2), (1, 1)]))       # (-q)_inf/(q)_inf = 1/((q;q^2)_inf (q)_inf)
    m = mel_lin(k, -1)
    ok &= (first_diff(m, rhs, NX) is None)
    ok &= (first_diff(multisum(k, (-1, -1), NX), m, NX) is None)   # Theorem 2 at x = 1/q
report("X1", ok, "x = 1/q: KRTW_k(1/q) = Mel_k(1/q) = ((-q)_inf/(q)_inf) sum_{i=1}^{k+1} (q^i,q^{k+2-i},q^{k+2};q^{k+2})_inf, k <= 6, mod q^%d" % NX)

# [X2] the two product forms of (iii) quoted in section 8
NX = 200
K.setN(NX)
N = NX
def schain(k):
    """sum_{0<=s_1<=...<=s_k} q^{sum s_i^2}/prod (q^2;q^2)_{s_i-s_{i-1}}"""
    smax = 1
    while smax * smax < NX:
        smax += 1
    f = [iqpoch(s, 2) for s in range(smax + 1)]
    for i in range(k - 1):
        g = []
        for s in range(smax + 1):
            acc = fmpz_poly([0])
            for sp in range(s + 1):
                if sp * sp < NX:
                    acc = acc + tr(tr(f[sp] * mono(sp * sp)) * iqpoch(s - sp, 2))
            g.append(acc)
        f = g
    tot = fmpz_poly([0])
    for s in range(smax + 1):
        if s * s < NX:
            tot = tot + tr(f[s] * mono(s * s))
    return tot
ok = True
for k in range(1, 9):
    m = 2 * k + 4
    p3 = prod_side([k + 1, k + 2, k + 2, k + 3, m, m], m, NX, [(2, 2), (1, 1)])
    # (-q^2;q^2)_inf = 1/(q^2;q^4)_inf
    second = tr(schain(k) * prod_side([k + 2, k + 2, m], m, NX, [(2, 4), (2, 2)]))
    ok &= (first_diff(p3, second, NX) is None)
report("X2", ok, "(iii) product = [s-chain of (iii)] * (-q^2;q^2)_inf (q^{k+2},q^{k+2},q^{2k+4};q^{2k+4})_inf/(q^2;q^2)_inf, k <= 8, mod q^%d" % NX)

print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED", " (%.0f s)" % (time.time() - t0))
