# Referee test of the level-m C_2 Rogers-Selberg identity with GENERIC (y1,y2):
#   sum_{lambda_1<=m} q^{|lambda|} P'_{2 lambda}(y1,y2;q)  vs  L^{(0)}_m(y;q)  with / without
#   1/( (q y1^2)_inf (q y2^2)_inf (q y1 y2)_inf ).
# LHS computed from Kostka-Foulkes polynomials (charge), independent of Bartlett-Warnaar's multisum.
import sys, itertools
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq

T = int(sys.argv[1]); N = int(sys.argv[2])   # t-degree <= T, q-degree < N ; y = (a t, b t)
a = Fr(2); b = Fr(-3, 5)

class TS:
    def __init__(self, c=None):
        self.c = c if c is not None else [fmpq_poly(0) for _ in range(T + 1)]
    def __add__(s, o): return TS([x + y for x, y in zip(s.c, o.c)])
    def __sub__(s, o): return TS([x - y for x, y in zip(s.c, o.c)])
    def __mul__(s, o):
        r = [fmpq_poly(0) for _ in range(T + 1)]
        for i, x in enumerate(s.c):
            if x.is_zero(): continue
            for j in range(T + 1 - i):
                if o.c[j].is_zero(): continue
                r[i + j] += x.mul_low(o.c[j], N)
        return TS(r)
    def scale(s, c): return TS([x * fmpq(c.numerator, c.denominator) for x in s.c])
    def inv(s):
        a0 = s.c[0]
        i0 = qinv(a0)
        r = [fmpq_poly(0) for _ in range(T + 1)]; r[0] = i0
        for d in range(1, T + 1):
            acc = fmpq_poly(0)
            for i in range(1, d + 1): acc += s.c[i].mul_low(r[d - i], N)
            r[d] = -acc.mul_low(i0, N)
        return TS(r)
    def __eq__(s, o): return all(x == y for x, y in zip(s.c, o.c))
def qinv(p):
    c = [Fr(int(v.p), int(v.q)) for v in p.coeffs()] + [Fr(0)] * N
    r = [Fr(0)] * N; r[0] = 1 / c[0]
    for n in range(1, N):
        r[n] = -sum(c[i] * r[n - i] for i in range(1, n + 1)) / c[0]
    return fmpq_poly([fmpq(v.numerator, v.denominator) for v in r])
def mono(td, qe, coef=Fr(1)):
    s = TS()
    if td <= T and qe < N:
        assert qe >= 0
        s.c[td] = fmpq_poly([0] * qe + [fmpq(coef.numerator, coef.denominator)])
    return s
ONE = mono(0, 0)
def fac(coef, td, qe):  # (1 - coef t^td q^qe)
    return ONE - mono(td, qe, coef)
def poch(coef, td, qe, n):  # (coef t^td q^qe; q)_n
    r = ONE
    for j in range(n): r = r * fac(coef, td, qe + j)
    return r
def poch_inf(coef, td, qe):
    r = ONE; j = 0
    while qe + j < N: r = r * fac(coef, td, qe + j); j += 1
    return r

# ---------------- Kostka-Foulkes via charge
def ssyt_2row(nu, mu):
    """SSYT of shape nu (<=2 rows) and content mu; returned as rows (lists)."""
    out = []
    def rec(i, rows):
        if i == len(mu):
            if [len(r) for r in rows] == list(nu) + [0] * (2 - len(nu)): out.append([list(r) for r in rows])
            return
        letter = i + 1; k = mu[i]
        l1, l2 = len(rows[0]), len(rows[1])
        # add a horizontal strip of size k: j1 boxes to row1, j2 to row2
        for j2 in range(0, k + 1):
            j1 = k - j2
            n1, n2 = l1 + j1, l2 + j2
            if n2 > l1: continue          # horizontal strip: row2 new boxes lie under old row1 boxes
            if n1 > nu[0] or n2 > (nu[1] if len(nu) > 1 else 0): continue
            rec(i + 1, [rows[0] + [letter] * j1, rows[1] + [letter] * j2])
    rec(0, [[], []])
    return out
def charge_word(w):
    w = list(w); total = 0
    while w:
        n = len(w); used = []
        pos = None; idx = 0; letter = 1; maxl = max(w)
        # find 1 scanning from the right
        p = n - 1
        while w[p] != 1: p -= 1
        used.append(p); cur = 0; pos = p
        while True:
            letter += 1
            if letter not in [w[i] for i in range(n) if i not in used]: break
            # search leftwards cyclically from pos for letter
            p = pos - 1; wrapped = False
            while True:
                if p < 0: p = n - 1; wrapped = True
                if p not in used and w[p] == letter: break
                p -= 1
            if wrapped: cur += 1
            total += cur; used.append(p); pos = p
        w = [w[i] for i in range(n) if i not in used]
    return total
def KF(nu, mu):
    """Kostka-Foulkes polynomial as dict power->coef"""
    res = {}
    for Tb in ssyt_2row(nu, mu):
        word = Tb[1] + Tb[0]  # bottom row first
        c = charge_word(word)
        res[c] = res.get(c, 0) + 1
    return res
# sanity table
assert KF((2,), (1, 1)) == {1: 1}
assert KF((1, 1), (1, 1)) == {0: 1}
assert KF((2, 1), (1, 1, 1)) == {1: 1, 2: 1}
assert KF((3,), (1, 1, 1)) == {3: 1}
assert KF((3, 1), (2, 1, 1)) == {1: 1, 2: 1}
assert KF((2, 2), (2, 1, 1)) == {1: 1}
assert KF((4,), (2, 1, 1)) == {3: 1}
assert KF((2, 2), (2, 2)) == {0: 1} and KF((3, 1), (2, 2)) == {1: 1} and KF((4,), (2, 2)) == {2: 1}
assert KF((4, 2), (2, 2, 2)) == {2: 1, 3: 1, 4: 1}; assert KF((2, 2), (1, 1, 1, 1)) == {2: 1, 4: 1}
assert KF((3,3),(2,2,2)) == {3: 1}

def schur2(nu):  # s_nu(a,b) for <=2 rows
    n1 = nu[0]; n2 = nu[1] if len(nu) > 1 else 0
    return (a * b) ** n2 * sum(a ** j * b ** (n1 - n2 - j) for j in range(n1 - n2 + 1))

def qpoch_q(n):  # (q;q)_n as fmpq_poly in q
    p = fmpq_poly([1])
    for j in range(1, n + 1): p = p * fmpq_poly([1] + [0] * (j - 1) + [-1])
    return p
def partitions_maxpart(n, m):
    def rec(n, mx):
        if n == 0: yield (); return
        for p in range(min(n, mx), 0, -1):
            for rest in rec(n - p, p): yield (p,) + rest
    return list(rec(n, m))

def LHS(m):
    s = TS()
    for size in range(0, T // 2 + 1):
        for lam in partitions_maxpart(size, m):
            mu = tuple(2 * p for p in lam)
            # b_mu(q) = prod_i (q)_{m_i(mu)}
            bmu = fmpq_poly([1])
            for part in set(mu): bmu = bmu * qpoch_q(mu.count(part))
            ibmu = qinv(bmu) if size else fmpq_poly([1])
            tot = fmpq_poly(0)
            for n1 in range(2 * size, -1, -1):
                n2 = 2 * size - n1
                if n2 > n1: break
                nu = (n1, n2) if n2 else (n1,)
                if size == 0: nu = ()
                kf = KF(nu, mu) if size else {0: 1}
                sv = schur2(nu) if size else Fr(1)
                for pw, cf in kf.items():
                    tot += fmpq_poly([0] * pw + [fmpq((cf * sv).numerator, (cf * sv).denominator)])
            tot = tot.mul_low(ibmu, N)
            s.c[2 * size] += fmpq_poly([0] * size + [1]).mul_low(tot, N) if size < N else 0
    return s

def L0(m):
    s = TS()
    for r1 in range(0, T + 1):
        for r2 in range(0, T + 1):
            td = 2 * (m + 1) * (r1 + r2)
            qe = (m + 1) * (r1 * r1 + r2 * r2) + r1 * (r1 - 1) + r2 * (r2 - 1)
            if td > T or qe >= N: continue
            coef = a ** (2 * (m + 1) * r1) * b ** (2 * (m + 1) * r2) * (a / b) ** r1 * (b / a) ** r2
            term = mono(td, qe, coef)
            # Delta ratio
            term = term * fac(a * a, 2, 2 * r1) * fac(b * b, 2, 2 * r2) * (fac(a * a, 2, 0) * fac(b * b, 2, 0)).inv()
            term = term * (mono(0, r1, a / (a - b)) - mono(0, r2, b / (a - b)))
            term = term * fac(a * b, 2, r1 + r2) * fac(a * b, 2, 0).inv()
            term = term * poch(a * a, 2, 0, r1) * poch(a * b, 2, 0, r1) * (poch(Fr(1), 0, 1, r1) * poch(a / b, 0, 1, r1)).inv()
            term = term * poch(b * b, 2, 0, r2) * poch(a * b, 2, 0, r2) * (poch(Fr(1), 0, 1, r2) * poch(b / a, 0, 1, r2)).inv()
            s = s + term
    return s

D = poch_inf(a * a, 2, 1) * poch_inf(b * b, 2, 1) * poch_inf(a * b, 2, 1)
for m in (1, 2, 3):
    lhs = LHS(m); l0 = L0(m)
    print("m=%d  T=%d N=%d:  LHS == L0/prod : %s   LHS == L0 : %s" % (m, T, N, lhs == l0 * D.inv(), lhs == l0))
    if m == 1:
        dif = lhs - l0
        print("   first nonzero coefficient of LHS - L0 at t-degree", next(i for i in range(T + 1) if not dif.c[i].is_zero()),
              ":", dif.c[next(i for i in range(T + 1) if not dif.c[i].is_zero())])
