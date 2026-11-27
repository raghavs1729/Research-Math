#!/usr/bin/env python3
"""Referee 3: (a) iv-all-k Lemma 1.1 on its own: K^Om_k (Omega literal) = sum_{i=1}^k omega(T^{i-1} O T^{k-i} u),
O F(R,S) = q^{R+2S} F(R,S+1), bivariate in Z[x][[q]] mod q^N (all powers of x), k <= KMAX;
(b) the x = 1 evaluations used in iv-all-k Theorem 5, mod q^NT, k <= 10:
    sum_n lam_n^{k+1} bOm_n |_{x=1} = (q^k, q^{k+4}, q^{2k+4}; q^{2k+4})_inf   (Jacobi triple product)
    (q^2;q^2)_inf J_k(1) = (q^{k+1}, q^{k+3}, q^{2k+4}; q^{2k+4})_inf,  J_k(1) = C_k(1)/C_1(1) from the C_k multisum.
"""
import sys
from flint import fmpz_poly

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 40
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 300

sys.argv = [sys.argv[0], str(NB), str(KMAX)]
import bivar as B

# enlarged state space so that O (which reads S+1) is exact at every state that T reads
B.STATES = [(r, s) for s in range(NB) for r in range(2 * NB) if (r + s) ** 2 + s * s - r < 2 * NB]


def O_apply(F):
    out = {}
    for (R, S) in B.STATES:
        src = F.get((R, S + 1))
        out[(R, S)] = B.shift(src, 0, R + 2 * S) if src is not None else B.zero()
    return out


def lemma11(k):
    u = {(r, s): B.mulq(B.ONE, B.Ainv[r]) for (r, s) in B.STATES}
    total = B.zero()
    for i in range(1, k + 1):
        F = u
        for _ in range(k - i):
            F = B.T_apply(F)
        F = O_apply(F)
        for _ in range(i - 1):
            F = B.T_apply(F)
        total = B.add(total, B.omega(F))
    return total


print(f"# (a) iv-all-k Lemma 1.1, bivariate mod q^{NB}, all powers of x")
ok = True
for k in range(1, KMAX + 1):
    lhs = B.KOm(k)
    rhs = lemma11(k)
    good = B.is_zero(B.sub(lhs, rhs))
    nc = B.first_bad(B.sub(B.KOm(k, drop_last=True), rhs))
    ok &= good
    print(f"k={k}: K^Om (Omega literal) == sum_i omega(T^(i-1) O T^(k-i) u): {good}; negative control (last Omega "
          f"term dropped) first difference at (x^d, q^n) = {nc}", flush=True)

# (b) theta evaluations at x = 1
N = NT


def trunc(c):
    return fmpz_poly(c[:N])


def prod(exps, mod):
    c = [1] + [0] * (N - 1)
    for a in exps:
        e = a
        while e < N:
            for n in range(N - 1, e - 1, -1):
                c[n] -= c[n - e]
            e += mod
    return c


def inv_poly_series(c):
    """1/c mod q^N, c[0] = 1"""
    out = [0] * N
    out[0] = 1
    for n in range(1, N):
        out[n] = -sum(c[j] * out[n - j] for j in range(1, n + 1))
    return out


Qinv = [[1] + [0] * (N - 1)]
for j in range(1, N):
    c = list(Qinv[-1])
    for n in range(2 * j, N):
        c[n] += c[n - 2 * j]
    Qinv.append(c)


def Ck1(k):
    """C_k(1) = sum_{0<=s_1<=...<=s_k} q^{sum s^2}/prod (Q;Q)_{s_i-s_{i-1}}"""
    smax = int(N ** 0.5) + 1
    E = {s: fmpz_poly([0] * (s * s) + [1]) if s * s < N else fmpz_poly([]) for s in range(smax + 1)}
    for _ in range(k - 1):
        newE = {}
        for s in range(smax + 1):
            acc = fmpz_poly([])
            for t in range(s, smax + 1):
                acc += E[t].mul_low(fmpz_poly(Qinv[t - s]), N)
            newE[s] = acc.left_shift(s * s).truncate(N)
        E = newE
    acc = fmpz_poly([])
    for s in range(smax + 1):
        acc += E[s].mul_low(fmpz_poly(Qinv[s]), N)
    return acc


print(f"# (b) x = 1 evaluations mod q^{N}")
C1 = Ck1(1)
C1inv = fmpz_poly(inv_poly_series(C1.coeffs() + [0] * N))
QQinf = fmpz_poly(prod([2], 2))
for k in range(1, 11):
    s = fmpz_poly([1])
    n = 1
    while (k + 2) * n * n - 2 * n < N:
        for e in ((k + 2) * n * n - 2 * n, (k + 2) * n * n + 2 * n):
            if e < N:
                s += fmpz_poly([0] * e + [(-1) ** n])
        n += 1
    jtp = (s - trunc(prod([k, k + 4, 2 * k + 4], 2 * k + 4))).is_zero()
    Jk1 = Ck1(k).mul_low(C1inv, N)
    thJ = (QQinf.mul_low(Jk1, N) - trunc(prod([k + 1, k + 3, 2 * k + 4], 2 * k + 4))).is_zero()
    ok &= jtp and thJ
    print(f"k={k}: sum lam^(k+1) bOm at x=1 = (q^k,q^(k+4),q^(2k+4);q^(2k+4)): {jtp};  (Q;Q)_inf J_k(1) = "
          f"(q^(k+1),q^(k+3),q^(2k+4);q^(2k+4)): {thJ}")
print("ALL CHECKS " + ("PASSED" if ok else "FAILED"))
