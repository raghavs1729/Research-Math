# Referee: does T commute with H (section 6, item 2)?  Exact arithmetic at rational (q,x), F = finite-support
# functions on Z^2 (point masses and random combinations), evaluated at many points (R,S).
from fractions import Fraction as Fr
import random, itertools
q = Fr(3, 7); x = Fr(5, 11); Q = q * q
def poch(a, p, n):
    if n < 0: return None  # 1/(p;p)_n = 0
    r = Fr(1)
    for i in range(n): r *= (1 - a * p ** i)
    return r
def inv_qq(n):  # 1/(q;q)_n, 0 for n<0
    return Fr(0) if n < 0 else 1 / poch(q, q, n)
def inv_QQ(n):
    return Fr(0) if n < 0 else 1 / poch(Q, Q, n)
def T(F):  # F: dict (r,s)->value, finite support; returns function
    def TF(R, S):
        tot = Fr(0)
        for (r, s), v in F.items():
            if r <= R and s >= S:
                tot += x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s) * v * inv_qq(R - r) * inv_QQ(s - S)
        return tot
    return TF
def H_apply(G):  # G: function (R,S)->value ; returns function
    def HG(R, S):
        a = G(R, S) - G(R - 1, S) / q
        return (q ** (-R) * Q ** (-S) * a + q ** 3 * x * q ** (-R) * a + x / q * q ** R * G(R, S)
                + q ** 2 * x ** 2 * q ** R * Q ** S * (G(R, S) - G(R, S + 1) / Q))
    return HG
def H_dict(F, pts):  # H applied to finite-support F, as dict on its (finite) support
    G = lambda R, S: F.get((R, S), Fr(0))
    HG = H_apply(G)
    sup = set()
    for (r, s) in F:
        for d in ((0, 0), (1, 0), (0, -1)): sup.add((r + d[0], s + d[1]))
    return {p: HG(*p) for p in sup}
random.seed(1)
bad = 0; tested = 0
tests = [{(r0, s0): Fr(1)} for r0 in range(-2, 4) for s0 in range(-2, 4)]
for _ in range(5):
    tests.append({(random.randint(-2, 3), random.randint(-2, 3)): Fr(random.randint(-5, 5), random.randint(1, 5)) for _ in range(4)})
for F in tests:
    HT = H_apply(T(F)); TH = T(H_dict(F, None))
    for R in range(-3, 6):
        for S in range(-3, 6):
            tested += 1
            if HT(R, S) != TH(R, S): bad += 1
print("H T = T H : points tested %d, mismatches %d" % (tested, bad))
# restricted to functions supported on N^2 and points in N^2
bad2 = 0; t2 = 0
for F in tests:
    F = {k: v for k, v in F.items() if k[0] >= 0 and k[1] >= 0}
    if not F: continue
    HT = H_apply(T(F)); TH = T(H_dict(F, None))
    for R in range(0, 6):
        for S in range(0, 6):
            t2 += 1
            if HT(R, S) != TH(R, S): bad2 += 1
print("restricted to N^2: points %d, mismatches %d" % (t2, bad2))
# example of a mismatch
F = {(0, 0): Fr(1)}
HT = H_apply(T(F)); TH = T(H_dict(F, None))
print("F=delta_(0,0): (HT F)(0,0) =", HT(0, 0), "  (TH F)(0,0) =", TH(0, 0))
