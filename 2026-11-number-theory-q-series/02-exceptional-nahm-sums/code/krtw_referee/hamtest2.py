# Referee: H (document) vs H' (found by commutant.py) at two rational points; exact arithmetic.
#   H' = P^{-1}Y^{-1}(1 - S) + x P^{-1}(1 - S) + x P + x^2 P Y (1 - U)
from fractions import Fraction as Fr
import sys
for (q, x) in ((Fr(3, 7), Fr(5, 11)), (Fr(-2, 9), Fr(13, 4))):
    Q = q * q
    def inv_poch(p, n):
        if n < 0: return Fr(0)
        r = Fr(1)
        for i in range(1, n + 1): r /= (1 - p ** i)
        return r
    def T(F):
        def TF(R, S):
            return sum((x ** (r + 2 * s) * q ** ((r + s) ** 2 + s * s) * v * inv_poch(q, R - r) * inv_poch(Q, s - S)
                        for (r, s), v in F.items() if r <= R and s >= S), Fr(0))
        return TF
    def mk(cS, cP1, cP, cPY, cU, cS2):
        # c0*P^-1Y^-1(1 - cS*Ssh) + cP1*P^-1(1 - cS2*Ssh) + cP*P + cPY*PY(1 - cU*Ush)
        def Hop(G):
            return lambda R, S: (q ** (-R) * Q ** (-S) * (G(R, S) - cS * G(R - 1, S))
                                 + cP1 * q ** (-R) * (G(R, S) - cS2 * G(R - 1, S))
                                 + cP * q ** R * G(R, S) + cPY * q ** R * Q ** S * (G(R, S) - cU * G(R, S + 1)))
        return Hop
    Hdoc = mk(1 / q, q ** 3 * x, x / q, q ** 2 * x * x, q ** -2, 1 / q)
    Href = mk(Fr(1), x, x, x * x, Fr(1), Fr(1))
    for name, Hop in (("document H", Hdoc), ("H'", Href)):
        bad = tot = 0
        for r0 in range(-3, 5):
            for s0 in range(-3, 5):
                F = {(r0, s0): Fr(1)}
                HF = {}
                G = lambda R, S, F=F: F.get((R, S), Fr(0))
                HG = Hop(G)
                for (r, s) in [(r0, s0), (r0 + 1, s0), (r0, s0 - 1)]: HF[(r, s)] = HG(r, s)
                A = Hop(T(F)); B = T(HF)
                for R in range(-4, 7):
                    for S in range(-4, 7):
                        tot += 1; bad += (A(R, S) != B(R, S))
        print("q=%s x=%s  %-10s : [H,T]F != 0 at %d of %d (point mass F, sample points)" % (q, x, name, bad, tot))
