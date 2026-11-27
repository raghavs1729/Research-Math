"""Referee check (iv): bivariate power series in (x, q), exact integer coefficients, truncated at x^DX, q^DQ.

  (A) K_k directly from its definition (transfer operator T = E_R E_S W, own code)
  (B) the spectral formula of Thm 6.4:  sum_{N>=M} lam^k V V# (x^2 q^{2L+2};Q)_{2M} (xq)_{2N} / ((xq)_inf (x^2Q;Q)_inf)
  (C) sum_{N>=M} lam^k t_{N,M} / ((xq)_inf (x^2Q;Q)_inf)                      (Lemma 7.1 + Sec. 7)
  (D) B_{k+1} C_k / C_1 computed from the multisum definitions of B, C       (Theorem 1)
All four are compared coefficientwise.  A series is a list over x-degree of flint.fmpz_poly in q (or fmpq_poly
where needed)."""
import sys
from flint import fmpz_poly, fmpq_poly, fmpq

DX = 12
DQ = 40


def zero():
    return [fmpq_poly([]) for _ in range(DX + 1)]


def one():
    s = zero(); s[0] = fmpq_poly([1]); return s


def trunc(p):
    return fmpq_poly(p.coeffs()[:DQ]) if p.degree() >= DQ else p


def add(a, b):
    return [trunc(a[i] + b[i]) for i in range(DX + 1)]


def scal(a, c):
    return [a[i] * c for i in range(DX + 1)]


def mul(a, b):
    out = zero()
    for i in range(DX + 1):
        if a[i] == 0:
            continue
        for j in range(DX + 1 - i):
            if b[j] == 0:
                continue
            out[i + j] += a[i] * b[j]
    return [trunc(p) for p in out]


def mono(c, a, b):
    """c x^a q^b as a series (b >= 0)."""
    s = zero()
    if a <= DX and b < DQ:
        s[a] = fmpq_poly([0] * b + [c])
    return s


def mul_mono(F, c, a, b):
    out = zero()
    if b >= DQ:
        return out
    sh = fmpq_poly([0] * b + [c])
    for i in range(DX + 1 - a):
        if F[i] != 0:
            out[i + a] = trunc(F[i] * sh)
    return out


def qinv_series(p):
    """inverse of a q-polynomial with constant term +-1, mod q^DQ"""
    c = p.coeffs() + [0] * DQ
    c0 = fmpq(c[0])
    inv = [fmpq(0)] * DQ
    inv[0] = 1 / c0
    for n in range(1, DQ):
        s = fmpq(0)
        for k in range(1, n + 1):
            if c[k] != 0:
                s += fmpq(c[k]) * inv[n - k]
        inv[n] = -s / c0
    return fmpq_poly(inv)


def mul_factor(F, c, a, b):
    """F * (1 - c x^a q^b)"""
    return add(F, mul_mono(F, -c, a, b))


def div_factor(F, c, a, b):
    """F / (1 - c x^a q^b), (a,b) != (0,0)"""
    assert (a, b) != (0, 0)
    if a == 0:
        inv = qinv_series(fmpq_poly([1] + [0] * (b - 1) + [-c]) if b > 0 else None)
        return [trunc(F[i] * inv) for i in range(DX + 1)]
    G = zero()
    sh = fmpq_poly([0] * b + [c]) if b < DQ else fmpq_poly([])
    for i in range(DX + 1):
        g = F[i]
        if i - a >= 0 and b < DQ:
            g = g + G[i - a] * sh
        G[i] = trunc(g)
    return G


def poch_mul(F, c, a, b, step, n, inverse=False):
    """F * prod_{i<n} (1 - c x^a q^{b+i*step})  (or divided by it)"""
    for i in range(n):
        F = div_factor(F, c, a, b + i * step) if inverse else mul_factor(F, c, a, b + i * step)
    return F


# pure-q helpers
def qpoch_inv(step, start, n):
    """1/prod_{i<n}(1-q^{start+i*step}) as fmpq_poly mod q^DQ (start >= 1)"""
    p = fmpq_poly([1])
    for i in range(n):
        e = start + i * step
        if e < DQ:
            p = trunc(p * qinv_series(fmpq_poly([1] + [0] * (e - 1) + [-1])))
    return p


IQ = [qpoch_inv(1, 1, n) for n in range(DX + 2)]   # 1/(q)_n
IQQ = [qpoch_inv(2, 2, n) for n in range(DX + 2)]  # 1/(Q;Q)_n


def times_q(F, p):
    return [trunc(F[i] * p) for i in range(DX + 1)]


# ---------------- (A) K_k directly ----------------
def K_direct(k):
    Rm, Sm = DX, DX // 2
    grid = {(R, S): times_q(one(), IQ[R]) for R in range(Rm + 1) for S in range(Sm + 1)}   # u = 1/(q)_R

    def T(F):
        WF = {(r, s): mul_mono(F[(r, s)], 1, r + 2 * s, (r + s)**2 + s**2) for (r, s) in F}
        # E_S: sum_b u^b/(Q;Q)_b : (E_S G)(R,S) = sum_{s>=S} G(R,s)/(Q;Q)_{s-S}
        ES = {}
        for (R, S) in F:
            acc = zero()
            for s in range(S, Sm + 1):
                acc = add(acc, times_q(WF[(R, s)], IQQ[s - S]))
            ES[(R, S)] = acc
        out = {}
        for (R, S) in F:
            acc = zero()
            for r in range(R + 1):
                acc = add(acc, times_q(ES[(r, S)], IQ[R - r]))
            out[(R, S)] = acc
        return out
    F = grid
    for _ in range(k - 1):
        F = T(F)
    tot = zero()
    for (r, s) in F:
        tot = add(tot, times_q(mul_mono(F[(r, s)], 1, r + 2 * s, (r + s)**2 + s**2), IQQ[s]))
    return tot


# ---------------- spectral pieces ----------------
def xq_inf_inv():
    F = one()
    for i in range(1, DQ):
        F = div_factor(F, 1, 1, i)
    return F


def x2Q_inf_inv():
    F = one()
    for i in range(1, DQ):
        if 2 * i < DQ:
            F = div_factor(F, 1, 2, 2 * i)
    return F


def VVs_term(N, M):
    """V_{N,M} V#_{N,M} (x^2 q^{2L+2};Q)_{2M} (xq;q)_{2N}  as a series"""
    L = N - M
    e = 2 * M * N + M * (2 * M - 1) + L * (2 * L - 1) + 2 * M * L
    F = mono((-1)**L, 2 * M + L, e)
    if 2 * M + L > DX or e >= DQ:
        return zero()
    # V
    F = poch_mul(F, -1, 1, 0, 1, M)                       # (-x;q)_M
    F = times_q(F, qpoch_inv(2, 2, M))                    # 1/(Q;Q)_M
    F = times_q(F, qpoch_inv(2, 2 * (L + 1), M))          # 1/(Q^{L+1};Q)_M
    F = poch_mul(F, 1, 1, M, 1, M, inverse=True)          # 1/(xq^M;q)_M
    F = poch_mul(F, 1, 2, 2 * N, 2, M, inverse=True)      # 1/(x^2 q^{2N};Q)_M
    # V#: (q;q^2)_L/((q)_L (q^L;q)_L)
    p = fmpq_poly([1])
    for i in range(L):
        p = trunc(p * fmpq_poly([1] + [0] * (2 * i) + [-1]))   # (1 - q^{2i+1})
    p = trunc(p * qpoch_inv(1, 1, L))
    if L >= 1:
        p = trunc(p * qpoch_inv(1, L, L))
    F = times_q(F, p)
    F = poch_mul(F, 1, 2, 2, 2, L)                        # (x^2 Q;Q)_L
    F = poch_mul(F, 1, 1, 2 * M + 1, 1, L, inverse=True)  # 1/(xq^{2M+1};q)_L
    F = poch_mul(F, 1, 1, 2 * M + L, 1, L, inverse=True)  # 1/(xq^{2M+L};q)_L
    # extra factors
    F = poch_mul(F, 1, 2, 2 * L + 2, 2, 2 * M)            # (x^2 q^{2L+2};Q)_{2M}
    F = poch_mul(F, 1, 1, 1, 1, 2 * N)                    # (xq;q)_{2N}
    return F


def t_term(N, M):
    """t_{N,M} = lam (-1)^{N+M} q^{N^2+M^2-N-M}(q^N+q^M)(1+xq^{N+M})(1-xq^{2N})(1-xq^{2M})
                  (x^2;Q)_N (x^2;Q)_M / ((Q;Q)_N (Q;Q)_M (1-x)(1-x^2)(1+delta))"""
    e = 2 * (N * N + M * M) - N - M
    if N + M > DX or e >= DQ:
        return zero()
    F = mono((-1)**(N + M), N + M, e)
    F = add(mul_mono(F, 1, 0, N), mul_mono(F, 1, 0, M))   # (q^N + q^M)
    F = mul_factor(F, -1, 1, N + M)
    # (1-xq^{2N})(1-xq^{2M})(x^2;Q)_N(x^2;Q)_M/((1-x)(1-x^2)): cancel explicitly
    # (x^2;Q)_n = (1-x^2)(x^2 Q;Q)_{n-1} for n>=1
    num = []
    den_x = 2  # we must divide by (1-x)(1-x^2)
    if N >= 1:
        F = poch_mul(F, 1, 2, 2, 2, N - 1)   # (x^2Q;Q)_{N-1}, the (1-x^2) cancels 1/(1-x^2)
        cancelled_x2 = True
    else:
        cancelled_x2 = False
    if M >= 1:
        F = poch_mul(F, 1, 2, 2, 2, M - 1)
        Mfac = True  # an extra (1-x^2) = (1-x)(1+x) remains
    else:
        Mfac = False
    # remaining: (1-xq^{2N})(1-xq^{2M}) * [(1-x^2) if M>=1] / ((1-x) [*(1-x^2) if not cancelled])
    if not cancelled_x2:
        # N = 0 -> M = 0: t = lam * 2 * (1+x)(1-x)^2/((1-x)(1-x^2)) / 2 = 1  (handled directly)
        assert N == 0 and M == 0
        return one()
    if Mfac:
        F = mul_factor(F, -1, 1, 0)          # (1+x): (1-x^2)/(1-x)
        F = mul_factor(F, 1, 1, 2 * N)       # (1-xq^{2N})
        F = mul_factor(F, 1, 1, 2 * M)       # (1-xq^{2M})
    else:
        # M = 0: (1-xq^{2M}) = (1-x) cancels 1/(1-x)
        F = mul_factor(F, 1, 1, 2 * N)
    F = times_q(F, qpoch_inv(2, 2, N))
    F = times_q(F, qpoch_inv(2, 2, M))
    if N == M:
        F = scal(F, fmpq(1, 2))
    return F


def spectral(k, which):
    tot = zero()
    for N in range(0, DX + 1):
        for M in range(0, N + 1):
            if (N + M) * k > DX or (N * N + M * M) * k >= DQ:
                continue
            term = VVs_term(N, M) if which == "VV" else t_term(N, M)
            tot = add(tot, mul_mono(term, 1, (N + M) * k, (N * N + M * M) * k))
    return mul(mul(tot, xq_inf_inv()), x2Q_inf_inv())


# ---------------- (D) B_{k+1} C_k / C_1 ----------------
def B_series(m):
    # f_m(N) = x^N q^{N^2}/(Q;Q)_N ; f_j(N) = x^N q^{N^2} sum_{N'<=N} f_{j+1}(N')/(q)_{N-N'}
    f = [times_q(mono(1, n, n * n), IQQ[n]) for n in range(DX + 1)]
    for _ in range(m - 1):
        g = []
        for n in range(DX + 1):
            acc = zero()
            for n2 in range(n + 1):
                acc = add(acc, times_q(f[n2], IQ[n - n2]))
            g.append(mul_mono(acc, 1, n, n * n))
        f = g
    tot = zero()
    for n in range(DX + 1):
        tot = add(tot, f[n])
    return tot


def C_series(k):
    h = [times_q(mono(1, s, s * s), IQQ[s]) for s in range(DX + 1)]
    for _ in range(k - 1):
        g = []
        for s in range(DX + 1):
            acc = zero()
            for s2 in range(s + 1):
                acc = add(acc, times_q(h[s2], IQQ[s - s2]))
            g.append(mul_mono(acc, 1, s, s * s))
        h = g
    tot = zero()
    for s in range(DX + 1):
        tot = add(tot, h[s])
    return tot


def C1_inv():
    # 1/(-xq;Q)_inf
    F = one()
    for i in range(DQ):
        e = 2 * i + 1
        if e < DQ:
            F = div_factor(F, -1, 1, e)
    return F


def equal(a, b):
    return all(a[i] == b[i] for i in range(DX + 1))


if __name__ == "__main__":
    if len(sys.argv) > 2:
        DX, DQ = int(sys.argv[1]), int(sys.argv[2])
        IQ[:] = [qpoch_inv(1, 1, n) for n in range(DX + 2)]
        IQQ[:] = [qpoch_inv(2, 2, n) for n in range(DX + 2)]
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    print(f"truncation x^{DX}, q^{DQ}")
    c1i = C1_inv()
    for k in range(1, kmax + 1):
        A = K_direct(k)
        Bv = spectral(k, "VV")
        Cv = spectral(k, "t")
        D = mul(mul(B_series(k + 1), C_series(k)), c1i)
        print(f"k={k}: K_direct==spectral(V V#): {equal(A, Bv)}   spectral==sum lam^k t: {equal(Bv, Cv)}"
              f"   K_direct==B_(k+1) C_k/C_1: {equal(A, D)}")
    # negative control: perturb one spectral term
    A = K_direct(2)
    tot = spectral(2, "VV")
    pert = add(tot, mul(mul(mono(1, 2 * 2 + 1, 2 * 2 + 1 + 3), xq_inf_inv()), x2Q_inf_inv()))
    print("negative control (perturbed (1,1)... term) detected:", not equal(A, pert))
