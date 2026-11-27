# Second commuting operator found in commutant.py (null vector 2, coefficients identified by hand):
#   G = P^{-2}Y^{-1}(1-S)(1-qS) + Y^{-1} - x(1+1/q) S + x^2 Y (1-U) + x^2 P^2 Y (1-U)
from fractions import Fraction as Fr
for (q, x) in ((Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))):
    Q = q * q
    def inv_poch(p, n):
        if n < 0: return Fr(0)
        r = Fr(1)
        for i in range(1, n + 1): r /= (1 - p ** i)
        return r
    def T(F):
        return lambda R, S: sum((x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s) * v * inv_poch(q, R - r) * inv_poch(Q, s - S)
                                 for (r, s), v in F.items() if r <= R and s >= S), Fr(0))
    def Gop(F):
        return lambda R, S: (q ** (-2 * R) * Q ** (-S) * (F(R, S) - (1 + q) * F(R - 1, S) + q * F(R - 2, S))
                             + Q ** (-S) * F(R, S) - x * (1 + 1 / q) * F(R - 1, S)
                             + x * x * Q ** S * (F(R, S) - F(R, S + 1)) + x * x * q ** (2 * R) * Q ** S * (F(R, S) - F(R, S + 1)))
    bad = tot = 0
    for r0 in range(-3, 5):
        for s0 in range(-3, 5):
            F = {(r0, s0): Fr(1)}
            Gf = Gop(lambda R, S: F.get((R, S), Fr(0)))
            GF = {(r, s): Gf(r, s) for r in range(r0, r0 + 3) for s in range(s0 - 1, s0 + 1)}
            A = Gop(T(F)); B = T(GF)
            for R in range(-4, 7):
                for S in range(-4, 7):
                    tot += 1; bad += (A(R, S) != B(R, S))
    print("q=%s x=%s  G : [G,T]F != 0 at %d of %d" % (q, x, bad, tot))
