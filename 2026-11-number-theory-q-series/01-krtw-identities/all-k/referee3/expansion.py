#!/usr/bin/env python3
"""Referee 3: exact x-adic check of the spectral expansion of u (F-all-k Theorem 6.4), q fixed rational.

Claim: u(R,S) = 1/(q)_R = sum_{N>=M>=0} u~_{N,M} phi_{N,M}(R,S),  u~ = V#(x^2q^{2L+2};Q)_{2M}(xq^{2M+1};q)_{2L}/(x^2Q;Q)_inf.
Multiply by (x^2Q;Q)_inf = sum_n (-1)^n Q^{n(n+1)/2} x^{2n}/(Q;Q)_n (Euler) and compare mod x^D, exactly, in Q[[x]];
phi_{N,M}(R,S) is the position-space eigenvector over Q(x) (valuations.py), V# the closed form of Cor. 5.2.
Only modes with L <= R, M >= S and L + 2(M - S) < D can contribute mod x^D (Lemma 6.2(a), ord V# = L);
we include all modes with L <= R and M <= S + D//2 + 1 and check that the extra ones are O(x^D).
Negative control: drop the mode (N,M) = (S+1, S) (or (S,S)) from the sum.
"""
import sys
from fractions import Fraction as Fr
import valuations as VL
from valuations import RF, Ctx, ONE, ZERO, xpow, lin

D = int(sys.argv[1]) if len(sys.argv) > 1 else 10


def series(rf, D):
    """power series of a rational function (den(0) != 0) mod x^D, as a list of Fractions"""
    num = [Fr(int(c.p), int(c.q)) for c in rf.n.coeffs()] + [Fr(0)] * D
    den = [Fr(int(c.p), int(c.q)) for c in rf.d.coeffs()] + [Fr(0)] * D
    assert den[0] != 0, "pole at x = 0"
    out = []
    for n in range(D):
        v = num[n] - sum(out[j] * den[n - j] for j in range(n))
        out.append(v / den[0])
    return out


def Vs_rf(C, N, M):
    q, Q = C.q, C.Q
    L = N - M
    r = xpow(L, (-1) ** L * q ** (L * (2 * L - 1) + 2 * M * L))
    num = Fr(1)
    for i in range(L):
        num *= 1 - q ** (2 * i + 1)
    den = C.qq(L)
    for i in range(L):
        den *= 1 - q ** (L + i)
    r = r * RF.c(num / den)
    r = r * C.xpoch(2, q * q, Q, L)
    r = r / (C.xpoch(1, q ** (2 * M + 1), q, L) * C.xpoch(1, q ** (2 * M + L), q, L))
    return r


def run(q):
    C = Ctx(q)
    Q = C.Q
    ok = True
    lines = []
    # (x^2Q;Q)_inf mod x^D
    e2 = [Fr(0)] * D
    for n in range(D):
        if 2 * n < D:
            e2[2 * n] = (-1) ** n * Q ** (n * (n + 1) // 2) / C.qq(n, Q)
    cache = {}
    for R in range(3):
        for S in range(3):
            target = [Fr(0)] * D
            inv = 1 / C.qq(R)
            target = [c * inv for c in e2]
            total = [Fr(0)] * D
            per_mode = {}
            for L in range(R + 1):
                for M in range(S, S + D // 2 + 2):
                    N = L + M
                    key = (N, M)
                    if key not in cache:
                        cache[key] = C.eig(N, M, 2)
                    ph = cache[key][(R, S)]
                    ut = Vs_rf(C, N, M) * C.xpoch(2, q ** (2 * L + 2), Q, 2 * M) * C.xpoch(1, q ** (2 * M + 1), q, 2 * L)
                    term = series(ut * ph, D) if not (ut * ph).iszero() else [Fr(0)] * D
                    per_mode[key] = term
                    total = [a + b for a, b in zip(total, term)]
                    # modes beyond the bound must be O(x^D)
                    if L + 2 * (M - S) >= D and any(term):
                        ok = False
                        lines.append(f"  FAIL: mode {key} at {(R, S)} not O(x^D)")
            good = total == target
            ok &= good
            # negative control
            drop = (S + 1, S) if R >= 1 else (S, S)
            nc = [a - b for a, b in zip(total, per_mode[drop])] == target
            lines.append(f"(R,S)=({R},{S}): sum_modes u~ phi == u mod x^{D}: {good};  without mode {drop}: equal = {nc}")
    return ok, lines


if __name__ == "__main__":
    allok = True
    for q in (Fr(2, 7), Fr(-3, 5)):
        ok, lines = run(q)
        print(f"=== q = {q}, exact in Q[[x]] mod x^{D} ===")
        for l in lines:
            print(l)
        allok &= ok
    print("THEOREM 6.4 EXPANSION CHECK " + ("PASSED" if allok else "FAILED"))
