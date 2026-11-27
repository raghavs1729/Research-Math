"""sw_lib.py -- exact q-series toolkit for the audit of Sun-Wang arXiv:2604.00847 (SUN-WANG-AUDIT.md).

A series is  sum_j c_j q^{(v+j)/den}  with integer coefficients, known for exponents (v+j) < P (in units 1/den).
Everything is exact integer arithmetic (python-flint fmpz_poly).

Contents
  Ctx(den, N)         context: all series are truncated at q^N, exponents in (1/den)Z
  ctx.nahm(M,B,d,C)   q^C * sum_n q^{n^T M n/2 + B.n}/prod (q^{d_i};q^{d_i})_{n_i}      (M symmetric pos. def.)
  ctx.vir(p,pp,r,s)   Virasoro minimal model character (Rocha-Caridi), incl. q^{h-c/24}
  ctx.sm(p,pp,r,s)    N=1 super-Virasoro character of SM(p,pp): NS if r+s even, R (sector Tr q^{L0-c/24},
                      prefactor prod (1+q^n)/(1-q^n), NO extra factor 2) if r+s odd
  ctx.u1(k,j)         U(1)_k character  eta^{-1} sum_n q^{k (n + j/(2k))^2},  j mod 2k,  h = j^2/(4k)
  ctx.pf(k,l,m)       Z_k parafermion character eta * c^l_m (Kac-Peterson string function)
"""
from fractions import Fraction as Fr
from flint import fmpz_poly
import math, sys
sys.setrecursionlimit(10000)


def _inv(M):
    n = len(M)
    A = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


class Ser:
    __slots__ = ("den", "v", "p", "P")

    def __init__(s, den, v, p, P):
        # normalise: strip to length P-v
        if not isinstance(p, fmpz_poly):
            p = fmpz_poly(p)
        L = P - v
        if L <= 0:
            p = fmpz_poly([0]);
        elif p.length() > L:
            p = p.truncate(L)
        s.den, s.v, s.p, s.P = den, v, p, P

    def _al(s, o):
        assert s.den == o.den
        v = min(s.v, o.v)
        return v, s.p.left_shift(s.v - v), o.p.left_shift(o.v - v), min(s.P, o.P)

    def __add__(s, o):
        if isinstance(o, int) and o == 0:
            return s
        v, a, b, P = s._al(o)
        return Ser(s.den, v, a + b, P)
    __radd__ = __add__

    def __sub__(s, o):
        v, a, b, P = s._al(o)
        return Ser(s.den, v, a - b, P)

    def __neg__(s):
        return Ser(s.den, s.v, -s.p, s.P)

    def __mul__(s, o):
        if isinstance(o, int):
            return Ser(s.den, s.v, s.p * o, s.P)
        assert s.den == o.den
        v = s.v + o.v
        P = min(s.P + o.v, o.P + s.v)
        L = max(P - v, 1)
        return Ser(s.den, v, s.p.mul_low(o.p, L), P)
    __rmul__ = __mul__

    def shift(s, a):
        """multiply by q^a"""
        k = Fr(a) * s.den
        assert k.denominator == 1, "shift not in (1/den)Z"
        return Ser(s.den, s.v + int(k), s.p, s.P + int(k))

    def cut(s, P):
        return Ser(s.den, s.v, s.p, min(P, s.P))

    def coeffs(s):
        """dict exponent(int, units 1/den) -> coeff for nonzero coefficients below P"""
        return {s.v + j: int(c) for j, c in enumerate(s.p.coeffs()) if c != 0 and s.v + j < s.P}

    def lead(s, n=6):
        d = s.coeffs()
        ks = sorted(d)[:n]
        return " + ".join("%d q^%s" % (d[k], Fr(k, s.den)) for k in ks)


def compare(a, b):
    """returns (ok, order, first_bad): a == b for all exponents < order (a Fraction, in q);
    first_bad = smallest exponent where they differ (Fraction) or None"""
    assert a.den == b.den
    P = min(a.P, b.P)
    da, db = a.coeffs(), b.coeffs()
    bad = [k for k in set(da) | set(db) if k < P and da.get(k, 0) != db.get(k, 0)]
    if bad:
        k = min(bad)
        return False, Fr(P, a.den), (Fr(k, a.den), da.get(k, 0), db.get(k, 0))
    return True, Fr(P, a.den), None


class Ctx:
    def __init__(c, den, N):
        c.den, c.N = den, N          # truncation: exponents < N (in q)
        c.P = N * den
        c._cache = {}

    # ---------- constructors
    def from_dict(c, d):
        """d: {Fraction exponent: int coeff}"""
        ks = [Fr(k) * c.den for k in d]
        for k in ks:
            assert k.denominator == 1, "exponent not in (1/den)Z: den=%d" % c.den
        if not d:
            return Ser(c.den, 0, fmpz_poly([0]), c.P)
        v = int(min(ks))
        L = c.P - v
        arr = [0] * max(L, 1)
        for k, val in d.items():
            j = int(Fr(k) * c.den) - v
            if j < L:
                arr[j] += val
        return Ser(c.den, v, fmpz_poly(arr), c.P)

    def ps(c, coeffs, step=1):
        """power series sum coeffs[i] q^{i*step}, step Fraction"""
        st = Fr(step) * c.den
        assert st.denominator == 1
        st = int(st)
        arr = [0] * c.P
        for i, a in enumerate(coeffs):
            if i * st < c.P:
                arr[i * st] = a
        return Ser(c.den, 0, fmpz_poly(arr), c.P)

    def one(c):
        return c.ps([1])

    def prod(c, f, step=1):
        """prod_{n>=1} (1 - x^n)^{f(n)}, x = q^step, f(n) integer (can be negative); as a series"""
        st = Fr(step)
        L = int(math.ceil(c.N / st)) + 1
        a = [0] * L
        a[0] = 1
        for n in range(1, L):
            e = f(n)
            if e > 0:
                for _ in range(e):
                    for j in range(L - 1, n - 1, -1):
                        a[j] -= a[j - n]
            elif e < 0:
                for _ in range(-e):
                    for j in range(n, L):
                        a[j] += a[j - n]
        return c.ps(a, st)

    def euler_inv(c):
        if "ei" not in c._cache:
            c._cache["ei"] = c.prod(lambda n: -1)
        return c._cache["ei"]

    def euler(c):
        if "e" not in c._cache:
            c._cache["e"] = c.prod(lambda n: 1)
        return c._cache["e"]

    def eta_inv(c):
        return c.euler_inv().shift(Fr(-1, 24))

    def theta(c, quad, sign=lambda k: 1):
        """sum_{k in Z} sign(k) q^{quad(k)}  (quad(k) -> Fraction, must grow quadratically)"""
        d = {}
        k = 0
        miss = 0
        while True:
            hit = False
            for kk in ((k, -k) if k else (0,)):
                e = Fr(quad(kk))
                if e < c.N:
                    d[e] = d.get(e, 0) + sign(kk)
                    hit = True
            if not hit:
                miss += 1
                if miss > 50:
                    break
            else:
                miss = 0
            k += 1
        return c.from_dict(d)

    # ---------- characters
    def vir(c, p, pp, r, s):
        """Virasoro minimal model M(p,pp) character chi_{r,s}, 1<=r<=p-1, 1<=s<=pp-1; exponent h-c/24 included"""
        assert math.gcd(p, pp) == 1 and 1 <= r < p and 1 <= s < pp
        num = c.theta(lambda k: k * (k * p * pp + r * pp - s * p)) - c.theta(lambda k: (k * p + r) * (k * pp + s))
        ex = Fr((r * pp - s * p) ** 2, 4 * p * pp) - Fr(1, 24)
        return (num * c.euler_inv()).shift(ex)

    def vir_h(c, p, pp, r, s):
        return Fr((r * pp - s * p) ** 2 - (p - pp) ** 2, 4 * p * pp)

    def ns_pref(c):
        if "ns" not in c._cache:
            # prod (1+q^{n-1/2})/(1-q^n) = prod_{m>=1} (1 - x^m)^{e(m)}, x = q^{1/2}:
            #   (1+x^{odd}) = (1-x^{2 odd})/(1-x^{odd});  1/(1-x^{even})
            def e(m):
                return (-1 if m % 2 == 1 else 0) + (1 if m % 4 == 2 else 0) + (-1 if m % 2 == 0 else 0)
            c._cache["ns"] = c.prod(e, Fr(1, 2))
        return c._cache["ns"]

    def r_pref(c):
        if "r" not in c._cache:
            # prod (1+q^n)/(1-q^n) = prod (1-q^{2n})/(1-q^n)^2
            c._cache["r"] = c.prod(lambda n: -2 + (1 if n % 2 == 0 else 0))
        return c._cache["r"]

    def sm(c, p, pp, r, s):
        """N=1 minimal model SM(p,pp) (pp-p even, gcd(p,(pp-p)/2)=1): NS character if r+s even, R if odd.
        R: q^{h-c/24} prod (1+q^n)/(1-q^n) * sum_k (...)   [no factor 2]."""
        assert (pp - p) % 2 == 0 and 1 <= r < p and 1 <= s < pp
        num = c.theta(lambda k: Fr(k * (k * p * pp + r * pp - s * p), 2)) - c.theta(lambda k: Fr((k * p + r) * (k * pp + s), 2))
        if (r + s) % 2 == 0:
            ex = Fr((r * pp - s * p) ** 2, 8 * p * pp) - Fr(1, 16)
            return (num * c.ns_pref()).shift(ex)
        ex = Fr((r * pp - s * p) ** 2, 8 * p * pp)
        return (num * c.r_pref()).shift(ex)

    def sm_weights(c, p, pp):
        """dict: ('NS'|'R', effective weight) -> list of labels (r,s) (one representative per orbit (r,s)~(p-r,pp-s)).
        effective weight = h - h_min, h_min = minimum over both sectors."""
        hs = {}
        for r in range(1, p):
            for s in range(1, pp):
                if (r, s) > (p - r, pp - s):
                    continue
                sec = 'NS' if (r + s) % 2 == 0 else 'R'
                h = Fr((r * pp - s * p) ** 2 - (pp - p) ** 2, 8 * p * pp) + (Fr(1, 16) if sec == 'R' else 0)
                hs[(sec, r, s)] = h
        hmin = min(hs.values())
        out = {}
        for (sec, r, s), h in hs.items():
            out.setdefault((sec, h - hmin), []).append((r, s))
        return out

    def sm_eff(c, p, pp, sec, w):
        """character of SM(p,pp) in sector sec with effective weight w (must be unique)"""
        W = c.sm_weights(p, pp)
        labs = W.get((sec, Fr(w)))
        assert labs is not None, "no %s weight %s in SM(%d,%d); have %s" % (sec, w, p, pp, sorted(k[1] for k in W if k[0] == sec))
        assert len(labs) == 1, "weight not unique: %s" % labs
        return c.sm(p, pp, *labs[0])

    def vir_eff(c, p, pp, w):
        hs = {}
        for r in range(1, p):
            for s in range(1, pp):
                if (r, s) > (p - r, pp - s):
                    continue
                hs[(r, s)] = c.vir_h(p, pp, r, s)
        hmin = min(hs.values())
        labs = [k for k, h in hs.items() if h - hmin == Fr(w)]
        assert len(labs) == 1, "Virasoro eff weight %s: labels %s" % (w, labs)
        return c.vir(p, pp, *labs[0])

    def u1(c, k, j):
        return (c.theta(lambda n: Fr((2 * k * n + j) ** 2, 4 * k)) * c.euler_inv()).shift(Fr(-1, 24))

    def pf(c, k, l, m):
        """Z_k parafermion character eta(tau) c^l_m(tau), l = 0..k, m = l mod 2 (m mod 2k).
        Kac-Peterson: c^l_m = eta^{-3} sum sign(x) q^{(k+2)x^2 - k y^2} over (x,y) with -|x| < y <= |x| and
        (x,y) or (1/2-x, 1/2+y) in ((l+1)/(2(k+2)), m/(2k)) + Z^2."""
        assert (l - m) % 2 == 0
        a, b = Fr(l + 1, 2 * (k + 2)), Fr(m, 2 * k)
        d = {}
        # enumerate x in a+Z or (1/2-a)+Z ; y accordingly
        X = int(math.sqrt(c.N + 2) * 3) + 5
        for (x0, y0) in ((a, b), (Fr(1, 2) - a, b - Fr(1, 2))):
            # (x,y) in (x0,y0)+Z^2 ; for the second family: 1/2 - x in a+Z <=> x in 1/2-a+Z ; 1/2+y in b+Z <=> y in b-1/2+Z
            x0f = x0 - math.floor(x0)
            y0f = y0 - math.floor(y0)
            for i in range(-X, X + 1):
                x = x0f + i
                ax = abs(x)
                # y in y0f+Z with -ax < y <= ax
                jlo = math.floor(-ax - y0f) + 1
                jhi = math.floor(ax - y0f)
                for j in range(jlo, jhi + 1):
                    y = y0f + j
                    if not (-ax < y <= ax):
                        continue
                    e = (k + 2) * x * x - k * y * y
                    if e < c.N + 1:
                        d[e] = d.get(e, 0) + (1 if x > 0 else -1)
        th = c.from_dict(d)
        e2 = c.euler_inv()
        return (th * e2 * e2).shift(Fr(-2, 24))

    # ---------- Nahm sums
    def nahm(c, M, B, d=None, C=0, count=False):
        """q^C sum_n q^{n^T M n/2 + B.n} / prod_i (q^{d_i};q^{d_i})_{n_i}, truncated at q^N.  M symmetric pos. def."""
        r = len(B)
        den, N = c.den, c.N
        if d is None:
            d = [1] * r
        M = [[Fr(x) for x in row] for row in M]
        B = [Fr(x) for x in B]
        C = Fr(C)
        assert all(M[i][j] == M[j][i] for i in range(r) for j in range(r)), "M not symmetric"
        M2 = [[x * den for x in row] for row in M]
        B2 = [2 * den * x for x in B]
        C1 = C * den
        for row in M2:
            for x in row:
                assert x.denominator == 1, "den too small for M"
        for x in B2:
            assert x.denominator == 1, "den too small for B"
        assert C1.denominator == 1, "den too small for C"
        M2 = [[int(x) for x in row] for row in M2]
        B2 = [int(x) for x in B2]
        C1 = int(C1)
        Mi = _inv(M)
        qmin = -sum(B[i] * Mi[i][j] * B[j] for i in range(r) for j in range(r)) / 2 + C     # unconstrained minimum
        amin = math.floor(qmin) - 1
        Lq = N - amin + 1                     # number of q-coefficients kept in the inner polynomials
        Nloc = float(N - C)                    # prune when (exponent without C) >= Nloc
        # Schur complement lower bounds
        S = [None] * (r + 1); L = [None] * (r + 1); c0 = [0.0] * (r + 1)
        for k in range(1, r + 1):
            if k < r:
                Myy = [[M[i][j] for j in range(k, r)] for i in range(k, r)]
                Miy = _inv(Myy)
                Myx = [[M[i][j] for j in range(k)] for i in range(k, r)]
                By = [B[i] for i in range(k, r)]
                m = r - k
                Sk = [[M[i][j] - sum(Myx[a][i] * Miy[a][b] * Myx[b][j] for a in range(m) for b in range(m)) for j in range(k)] for i in range(k)]
                Lk = [B[i] - sum(Myx[a][i] * Miy[a][b] * By[b] for a in range(m) for b in range(m)) for i in range(k)]
                ck = -sum(By[a] * Miy[a][b] * By[b] for a in range(m) for b in range(m)) / 2
            else:
                Sk = [row[:] for row in M]; Lk = B[:]; ck = Fr(0)
            S[k] = [[float(x) for x in row] for row in Sk]; L[k] = [float(x) for x in Lk]; c0[k] = float(ck)
        E = {}
        for dd in set(d):
            tab = [fmpz_poly([1])]
            n = 0
            while dd * n < Lq:
                n += 1
                geo = [0] * Lq
                for j in range(0, Lq, dd * n):
                    geo[j] = 1
                tab.append(tab[-1].mul_low(fmpz_poly(geo), Lq))
            E[dd] = tab            # tab[n] = 1/(q^dd;q^dd)_n mod q^Lq ; for n >= len-1 it equals tab[-1]
        acc = [fmpz_poly([0]) for _ in range(den)]
        x = [0] * r
        cnt = [0]
        eps = 1e-6
        nonneg = all(v >= 0 for row in M for v in row) and all(v >= 0 for v in B)
        Mf = [[float(v) for v in row] for row in M]; Bf = [float(v) for v in B]

        def lbu(k):
            Sk = S[k]; Lk = L[k]; v = c0[k]
            for i in range(k):
                xi = x[i]
                if xi:
                    row = Sk[i]; t = 0.0
                    for j in range(k):
                        t += row[j] * x[j]
                    v += xi * (Lk[i] + 0.5 * t)
            return v

        def lb0(k):
            v = 0.0
            for i in range(k):
                xi = x[i]
                if xi:
                    row = Mf[i]; t = 0.0
                    for j in range(k):
                        t += row[j] * x[j]
                    v += xi * (Bf[i] + 0.5 * t)
            return v

        def rec(k, cur):
            if k == r:
                e2 = sum(x[i] * M2[i][j] * x[j] for i in range(r) for j in range(r)) + sum(B2[i] * x[i] for i in range(r))
                assert e2 % 2 == 0, "exponent not in (1/den)Z"
                e = e2 // 2 + C1
                if e < N * den:
                    cnt[0] += 1
                    a, b = divmod(e, den)
                    assert a >= amin
                    acc[b] = acc[b] + cur.left_shift(a - amin)
                return
            tab = E[d[k]]; i = 0; prev = -1e300; lt = len(tab) - 1
            while True:
                x[k] = i
                if nonneg:
                    b = lb0(k + 1)
                    if b >= Nloc + eps:
                        break
                    rec(k + 1, cur.mul_low(tab[min(i, lt)], Lq))
                else:
                    b = lbu(k + 1)
                    if b >= Nloc + eps:
                        if i > 0 and b >= prev:
                            break
                    else:
                        rec(k + 1, cur.mul_low(tab[min(i, lt)], Lq))
                    prev = b
                i += 1
            x[k] = 0

        rec(0, fmpz_poly([1]))
        if count:
            print("   [nahm terms: %d]" % cnt[0], file=sys.stderr)
        # assemble
        arr = [0] * ((N - amin) * den)
        for b in range(den):
            cs = acc[b].coeffs()
            for j, cf in enumerate(cs):
                idx = j * den + b
                if idx < len(arr):
                    arr[idx] = int(cf)
        return Ser(den, amin * den, fmpz_poly(arr), N * den)


def report(name, lhs, rhs, log=None):
    ok, order, bad = compare(lhs, rhs)
    if ok:
        msg = "%-34s HOLDS  mod q^%s   [%s]" % (name, order, lhs.lead(4))
    else:
        msg = "%-34s FAILS  first difference at q^%s: lhs %d, rhs %d   (compared mod q^%s)\n      lhs: %s\n      rhs: %s" % (
            name, bad[0], bad[1], bad[2], order, lhs.lead(6), rhs.lead(6))
    print(msg, flush=True)
    if log is not None:
        log.append((name, ok, order, bad))
    return ok
