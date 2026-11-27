# Durfee-rectangle identity of Lemma M2, Euler split of Lemma M1, closed form (3.1) of K_2, and the
# "direct reading" K_2 = (-xq;q)_inf J^(1)_1 J^(2)_1.
from flint import fmpz_poly
from ser import S, qpoch_inv_table, monomial, poch_inf, poch
import direct
N = 80
for d in (1, 2):
    t = qpoch_inv_table(d, 12, N); t2 = qpoch_inv_table(2 * d, 12, N)
    okD = all(t[A].mul_low(t[B], N) == sum((t[c].mul_low(t[A - c], N).mul_low(t[B - c], N).left_shift(d * (A - c) * (B - c)).truncate(N)
              for c in range(min(A, B) + 1)), fmpz_poly(0)) for A in range(9) for B in range(9))
    okE = all(t[n] == sum((t2[u].mul_low(t2[n - u], N).left_shift(d * (n - u)).truncate(N) for u in range(n + 1)), fmpz_poly(0)) for n in range(12))
    print("d=%d Durfee identity A,B<=8: %s ; Euler split n<=11: %s" % (d, okD, okE))
X = 12
K2 = direct.K(2, X, N)
tQ = qpoch_inv_table(2, X + 2, N)
s = S(X, N)
for p in range(0, X + 1):
    for u in range(0, X + 1):
        for v in range(0, X + 1):
            for al in range(0, X + 1):
                for ga in range(0, X + 1):
                    w = 4 * p + 2 * u + 2 * v + al + ga
                    if w > X: break
                    e = p * p + (p + u + v) ** 2 + (p + u + al) ** 2 + (p + v + ga) ** 2 + ga
                    if e >= N: continue
                    c = tQ[p].mul_low(tQ[u], N).mul_low(tQ[v], N).mul_low(tQ[al], N).mul_low(tQ[ga], N)
                    s.c[w] += c.left_shift(e).truncate(N)
print("(3.1) == K_2 mod x^%d q^%d:" % (X + 1, N), s == K2)
def fJ(e_):
    out = S(X, N)
    for sg in range(0, X // 2 + 1):
        t = monomial(X, N, 2 * sg, 2 * sg * sg) * poch(X, N, [(-1, 1, e_ + 2 * j) for j in range(sg)]).inv()
        t.c = [c.mul_low(tQ[sg], N) for c in t.c]
        out = out + t
    return out
print("K_2 == (-xq;q)_inf J1 J2:", K2 == poch_inf(X, N, -1, 1, 1, 1) * fJ(1) * fJ(2))
