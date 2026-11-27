#!/usr/bin/env python3
"""Referee 3, check (2): bivariate identities in Z[x][[q]] mod q^N, all powers of x.

Every series below has only monomials x^a q^b with a <= b, so mod q^N only x^0..x^{N-1} occur and the
truncation is exact.  Series are lists (index = power of x) of flint fmpz_poly in q truncated at q^N.

Checked (k = 1..KMAX):
  (F_k)      K_k C_1 = B_{k+1} C_k                                   [F-all-k Theorem 1]
  (F^Om_k)   K^Om_k C_1 = B^(k)_{k+1} C_k                            [iv-all-k Theorem 6, second part]
  (B1)       (xq)_inf B_{k+1} = sum_n lam_n^{k+1} beta^B_n
  (B2)       (x^2Q;Q)_inf C_k = C_1 sum_m lam_m^{k+1} beta^J_m
  Thm 4      (xq)_inf K^Om_k C_1 = C_k sum_n lam_n^{k+1} b^Om_n
  Thm 6      (xq)_inf B^(k)_{k+1} = sum_n lam_n^{k+1} b^Om_n
  Prop 2.3   x K^Om_k = K^flat_k - K_k   (K^flat: extra factor q^{-r_1})
  Thm 6.4    (xq)_inf (x^2Q;Q)_inf K_k = sum_{N>=M} lam_{N,M}^k V_{N,M} Vs_{N,M} (x^2 q^{2L+2};Q)_{2M} (xq;q)_{2N}
             with V, Vs the closed forms of Cor. 4.2 / Cor. 5.2 (NOT using Lemma 7.1)
  Thm 3.2    (xq)(x^2Q;Q) K^flat_k (1+x) = sum lam^k t eps ,  (xq)(x^2Q;Q) K^Om_k x(1+x) = sum lam^k t (eps-1-x)
K_k, K^Om_k (Omega literal), B, B^(k), C_k are computed from their multisum definitions by nested summation.
Negative controls: B_{k+2} in place of B_{k+1}; one spectral weight perturbed; Omega's last term dropped.
Usage: python3 bivar.py N KMAX
"""
import sys
import time
from flint import fmpz_poly

NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 40
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
ND = NQ  # x-degrees 0..NQ-1


def zero():
    return [fmpz_poly([]) for _ in range(ND)]


def qmono(e, c=1):
    if e >= NQ:
        return fmpz_poly([])
    assert e >= 0, e
    return fmpz_poly([0] * e + [c])


def bmono(a, b, c=1):
    """c x^a q^b"""
    z = zero()
    if a < ND and b < NQ:
        z[a] = qmono(b, c)
    return z


def add(A, B):
    return [a + b for a, b in zip(A, B)]


def sub(A, B):
    return [a - b for a, b in zip(A, B)]


def scal(A, c):
    return [a * c for a in A]


def mul(A, B):
    C = zero()
    for i, a in enumerate(A):
        if a.is_zero():
            continue
        for j in range(ND - i):
            b = B[j]
            if not b.is_zero():
                C[i + j] += a.mul_low(b, NQ)
    return C


def mulq(A, p):
    """multiply by a univariate q-series p"""
    return [a.mul_low(p, NQ) for a in A]


def shift(A, a, b, c=1):
    """multiply by c x^a q^b (b may be negative if the result has no negative powers)"""
    C = zero()
    for i in range(ND - a):
        p = A[i]
        if p.is_zero():
            continue
        cs = p.coeffs()
        if b >= 0:
            C[i + a] = (p.left_shift(b) * c).truncate(NQ)
        else:
            assert all(v == 0 for v in cs[:(-b)]), "negative q-power"
            C[i + a] = (p.right_shift(-b) * c).truncate(NQ)
    return C


def inv1(A, a, b, s=1):
    """multiply by 1/(1 - s x^a q^b), (a, b) != (0, 0), b >= 1"""
    assert b >= 1
    C = [p for p in A]
    for d in range(ND):
        if d - a >= 0:
            C[d] = C[d] + (C[d - a].left_shift(b) * s).truncate(NQ)
    # the recursion above is correct only if we process degrees in increasing order AND, for a = 0,
    # the q-direction recursion; handle a = 0 separately
    if a == 0:
        C = [p for p in A]
        for d in range(ND):
            cs = C[d].coeffs() + [0] * NQ
            cs = cs[:NQ]
            for n in range(b, NQ):
                cs[n] += s * cs[n - b]
            C[d] = fmpz_poly(cs).truncate(NQ)
    return C


def bin1(A, a, b, s=1):
    """multiply by (1 - s x^a q^b)"""
    return sub(A, shift(A, a, b, s))


def is_zero(A):
    return all(p.is_zero() for p in A)


def first_bad(A):
    best = None
    for d, p in enumerate(A):
        cs = p.coeffs()
        for n, v in enumerate(cs):
            if v != 0:
                if best is None or n < best[1]:
                    best = (d, n)
                break
    return best


# univariate 1/(q;q)_j, 1/(q^2;q^2)_j
def _div_geom(c, j):
    c = list(c)
    for n in range(j, NQ):
        c[n] += c[n - j]
    return c


_a = [[1] + [0] * (NQ - 1)]
_b = [[1] + [0] * (NQ - 1)]
for j in range(1, NQ + 2):
    _a.append(_div_geom(_a[-1], j))
    _b.append(_div_geom(_b[-1], 2 * j))
Ainv = [fmpz_poly(c) for c in _a]     # 1/(q;q)_j
Binv = [fmpz_poly(c) for c in _b]     # 1/(q^2;q^2)_j

ONE = bmono(0, 0)

# ---------------- multisums by nested summation ----------------
# superset of the states needed: K^flat has the weight q^{(r+s)^2+s^2-r} at level 1
STATES = [(r, s) for s in range(NQ) for r in range(NQ) if (r + s) ** 2 + s * s - r < NQ]


def wmul(F, r, s, extra_q=0):
    return shift(F, r + 2 * s, (r + s) ** 2 + s * s + extra_q)


def T_apply(F, extra=None):
    """(TF)(R,S) = sum_{r<=R, s>=S} w(r,s) q^{extra(r,s)} F(r,s) / ((q)_{R-r} (Q;Q)_{s-S}) on STATES"""
    G = {}
    for (r, s) in STATES:
        e = (r + s) ** 2 + s * s + (extra(r, s) if extra else 0)
        if e < NQ:
            G[(r, s)] = wmul(F[(r, s)], r, s, (extra(r, s) if extra else 0))
    out = {}
    for (R, S) in STATES:
        acc = zero()
        for (r, s), g in G.items():
            if r <= R and s >= S:
                acc = add(acc, mulq(mulq(g, Ainv[R - r]), Binv[s - S]))
        out[(R, S)] = acc
    return out


def omega(F, flat=False):
    acc = zero()
    for (r, s) in STATES:
        g = wmul(F[(r, s)], r, s, -r if flat else 0)
        acc = add(acc, mulq(g, Binv[s]))
    return acc


def K(k, flat=False):
    F = {(r, s): mulq(ONE, Ainv[r]) for (r, s) in STATES}
    for _ in range(k - 1):
        F = T_apply(F)
    return omega(F, flat)


def KOm(k, drop_last=False):
    F = {(r, s): mulq(ONE, Ainv[r]) for (r, s) in STATES}
    if drop_last:
        FO = {(r, s): zero() for (r, s) in STATES}
    else:
        FO = {(r, s): shift(mulq(ONE, Ainv[r]), 0, r + 2 * s) for (r, s) in STATES}
    for _ in range(k - 1):
        TF = T_apply(F)
        TYF = T_apply(F, extra=lambda r, s: 2 * s)
        TFO = T_apply(FO)
        FO = {(R, S): sub(add(TFO[(R, S)], shift(TF[(R, S)], 0, R + 2 * S)), shift(TYF[(R, S)], 0, R))
              for (R, S) in STATES}
        F = TF
    return omega(FO)


def Bsum(k, variant=False):
    """B_{k+1} (variant=False) or B^(k)_{k+1} (variant=True): indices N_1>=...>=N_{k+1}>=0"""
    nmax = int(NQ ** 0.5) + 1
    # D(m) at level k+1: 1/(Q;Q)_m
    D = {m: mulq(ONE, Binv[m]) for m in range(nmax + 1)}
    # levels k+1 -> 1: D_i(n) = sum_{m<=n} x^m q^{m^2 (+m if i+1 in {k,k+1})} D_{i+1}(m)/(q)_{n-m}
    for i in range(k, 0, -1):
        lin = 1 if (variant and (i + 1) in (k, k + 1)) else 0
        newD = {}
        for n in range(nmax + 1):
            acc = zero()
            for m in range(n + 1):
                if m * m + lin * m < NQ:
                    acc = add(acc, mulq(shift(D[m], m, m * m + lin * m), Ainv[n - m]))
            newD[n] = acc
        D = newD
    lin1 = 1 if (variant and 1 in (k, k + 1)) else 0
    acc = zero()
    for n in range(nmax + 1):
        if n * n + lin1 * n < NQ:
            acc = add(acc, shift(D[n], n, n * n + lin1 * n))
    return acc


def Csum(k):
    smax = int(NQ ** 0.5) + 1
    E = {s: bmono(s, s * s) for s in range(smax + 1)}          # level k
    for _ in range(k - 1):
        newE = {}
        for s in range(smax + 1):
            acc = zero()
            for t in range(s, smax + 1):
                acc = add(acc, mulq(E[t], Binv[t - s]))
            newE[s] = shift(acc, s, s * s)
        E = newE
    acc = zero()
    for s in range(smax + 1):
        acc = add(acc, mulq(E[s], Binv[s]))
    return acc


def poch(a, b, step, n, s=1):
    """prod_{i<n} (1 - s x^a q^{b+i*step}) as a bivariate series"""
    A = ONE
    for i in range(n):
        A = bin1(A, a, b + i * step, s)
    return A


def poch_inf(a, b, step, s=1):
    A = ONE
    i = 0
    while b + i * step < NQ:
        A = bin1(A, a, b + i * step, s)
        i += 1
    return A


def inv_poch(A, a, b, step, n, s=1):
    for i in range(n):
        A = inv1(A, a, b + i * step, s)
    return A


XQinf = poch_inf(1, 1, 1)          # (xq;q)_inf
X2Qinf = poch_inf(2, 2, 2)         # (x^2 q^2; q^2)_inf


def betaB(n):
    """beta^B_n = (-1)^n q^{n^2} (x^2;Q)_n (1-xq^{2n}) / ((Q;Q)_n (1-x)) as a series (polynomial in x)"""
    if n == 0:
        return ONE
    A = add(ONE, bmono(1, 0))                  # (1+x)
    A = mul(A, poch(2, 2, 2, n - 1))           # (x^2 Q;Q)_{n-1}
    A = bin1(A, 1, 2 * n)                      # (1 - x q^{2n})
    A = mulq(A, Binv[n])
    return shift(A, 0, n * n, (-1) ** n)


def betaJ(m):
    if m == 0:
        return ONE
    A = poch(2, 2, 2, m - 1)
    A = bin1(A, 2, 4 * m)
    A = mulq(A, Binv[m])
    return shift(A, 0, m * m - m, (-1) ** m)


def bOm_times_lam(n, k):
    """lam_n^{k+1} b^Om_n, b^Om_n = (-1)^n q^{n^2-2n} (x^2Q;Q)_{n-1}(1-xq^{2n})(1-(1-x)q^{2n}+x^2q^{4n})/(x (Q;Q)_n)"""
    if n == 0:
        return ONE
    A = poch(2, 2, 2, n - 1)
    A = bin1(A, 1, 2 * n)
    poly = sub(ONE, sub(bmono(0, 2 * n), bmono(1, 2 * n)))
    poly = add(poly, bmono(2, 4 * n))
    A = mul(A, poly)
    A = mulq(A, Binv[n])
    # times x^{(k+1)n - 1} q^{(k+1)n^2 + n^2 - 2n}
    e = (k + 2) * n * n - 2 * n
    return shift(A, (k + 1) * n - 1, e, (-1) ** n)


def lam_sum(fn, k):
    acc = zero()
    n = 0
    while True:
        if n > 0 and (k + 2) * n * n - 2 * n >= NQ:
            break
        if fn is bOm_times_lam:
            acc = add(acc, bOm_times_lam(n, k))
        else:
            term = fn(n)
            acc = add(acc, shift(term, (k + 1) * n, (k + 1) * n * n))
        n += 1
    return acc


def lamNM(N_, M_, k):
    return ((N_ + M_) * k, (N_ * N_ + M_ * M_) * k)


def t_series(N_, M_):
    """t_{N,M} without the lam factor:  (bB_N bJ_M + bB_M bJ_N)/(1+delta)"""
    A = add(mul(betaB(N_), betaJ(M_)), mul(betaB(M_), betaJ(N_)))
    if N_ == M_:
        A = [p * 1 for p in A]
        A = [fmpz_poly([c // 2 for c in p.coeffs()]) for p in A]
    return A


def eps_mul(A, N_, M_):
    """multiply by eps_{N,M} = q^{-N-M} + x q^{N-M} + x q^{M-N} + x^2 q^{N+M}  (negative powers allowed
    transiently: we multiply first by q^{N+M} and shift back at the end via the lam factor)"""
    # return q^{N+M} eps A
    return add(add(A, shift(A, 1, 2 * N_)), add(shift(A, 1, 2 * M_), shift(A, 2, 2 * N_ + 2 * M_)))


def V_series(N_, M_):
    """V_{N,M} = x^{2M} q^{2MN+M(2M-1)} (-x;q)_M / ((Q;Q)_M (Q^{N-M+1};Q)_M (xq^M;q)_M (x^2q^{2N};Q)_M)"""
    A = bmono(2 * M_, 2 * M_ * N_ + M_ * (2 * M_ - 1))
    if 2 * M_ >= ND or 2 * M_ * N_ + M_ * (2 * M_ - 1) >= NQ:
        return zero()
    A = mul(A, poch(1, 0, 1, M_, s=-1))                 # (-x;q)_M
    A = mulq(A, Binv[M_])                                 # 1/(Q;Q)_M
    for i in range(M_):                                   # 1/(Q^{L+1};Q)_M
        A = inv1(A, 0, 2 * (N_ - M_ + 1 + i))
    A = inv_poch(A, 1, M_, 1, M_)                         # 1/(xq^M;q)_M
    A = inv_poch(A, 2, 2 * N_, 2, M_)                     # 1/(x^2q^{2N};Q)_M
    return A


def Vs_series(N_, M_):
    """V#_{N,M} = (-x)^L q^{L(2L-1)+2ML} (q;q^2)_L/((q)_L (q^L;q)_L) (x^2Q;Q)_L/((xq^{2M+1};q)_L (xq^{2M+L};q)_L)"""
    L = N_ - M_
    e = L * (2 * L - 1) + 2 * M_ * L
    if L >= ND or e >= NQ:
        return zero()
    A = bmono(L, e, (-1) ** L)
    for i in range(L):
        A = bin1(A, 0, 2 * i + 1)        # (q;q^2)_L
    A = mulq(A, Ainv[L])                 # 1/(q)_L
    for i in range(L):
        A = inv1(A, 0, L + i)            # 1/(q^L;q)_L   (L >= 1 so exponents >= 1)
    A = mul(A, poch(2, 2, 2, L))         # (x^2Q;Q)_L
    A = inv_poch(A, 1, 2 * M_ + 1, 1, L)
    A = inv_poch(A, 1, 2 * M_ + L, 1, L)
    return A


def spectral_K(k, perturb=False):
    """sum_{N>=M} lam_{N,M}^k V Vs (x^2 q^{2L+2};Q)_{2M} (xq;q)_{2N}  (Theorem 6.4 x (xq)(x^2Q;Q))"""
    acc = zero()
    for N_ in range(0, NQ):
        for M_ in range(0, N_ + 1):
            a, b = lamNM(N_, M_, k)
            if b >= NQ or a >= ND:
                continue
            L = N_ - M_
            A = mul(V_series(N_, M_), Vs_series(N_, M_))
            A = mul(A, poch(2, 2 * L + 2, 2, 2 * M_))
            A = mul(A, poch(1, 1, 1, 2 * N_))
            if perturb and (N_, M_) == (1, 0):
                A = add(A, bmono(1, 1))
            acc = add(acc, shift(A, a, b))
    return acc


def spectral_t(k, which):
    """which = 'flat': sum lam^k t eps ;  'om': sum lam^k t (eps - 1 - x)"""
    acc = zero()
    for N_ in range(0, NQ):
        for M_ in range(0, N_ + 1):
            a, b = lamNM(N_, M_, k)
            a += N_ + M_
            b += N_ * N_ + M_ * M_          # t contains one more lam factor
            if b - N_ - M_ >= NQ or a >= ND:
                continue
            tt = t_series(N_, M_)
            E = eps_mul(tt, N_, M_)               # q^{N+M} eps t
            if which == 'flat':
                A = E
            else:
                A = sub(E, shift(add(tt, shift(tt, 1, 0)), 0, N_ + M_))   # q^{N+M}(eps - 1 - x) t
            acc = add(acc, shift(A, a, b - N_ - M_))
    return acc


if __name__ == "__main__":
    t0 = time.time()
    print(f"# referee3 bivar.py: all identities in Z[x][[q]] mod q^{NQ}, all powers of x (x^0..x^{ND-1})")
    C1 = Csum(1)
    allok = True
    sumBB = {}
    for k in range(1, KMAX + 1):
        Kk = K(k)
        Kf = K(k, flat=True)
        KO = KOm(k)
        Bk = Bsum(k)
        Bv = Bsum(k, variant=True)
        Ck = Csum(k)
        res = {}
        res['(F_k)  K C1 - B_{k+1} C_k'] = sub(mul(Kk, C1), mul(Bk, Ck))
        res['(F^Om) KOm C1 - B^(k) C_k'] = sub(mul(KO, C1), mul(Bv, Ck))
        sB = lam_sum(betaB, k)
        sJ = lam_sum(betaJ, k)
        sO = lam_sum(bOm_times_lam, k)
        res['(B1)'] = sub(mul(XQinf, Bk), sB)
        res['(B2)'] = sub(mul(X2Qinf, Ck), mul(C1, sJ))
        res['Thm 4'] = sub(mul(mul(XQinf, KO), C1), mul(Ck, sO))
        res['Thm 6'] = sub(mul(XQinf, Bv), sO)
        res['Prop 2.3'] = sub(shift(KO, 1, 0), sub(Kf, Kk))
        lhsK = mul(mul(XQinf, X2Qinf), Kk)
        res['Thm 6.4 (V,V# closed forms)'] = sub(lhsK, spectral_K(k))
        res['Thm 3.2 flat'] = sub(shift(add(mul(mul(XQinf, X2Qinf), Kf), shift(mul(mul(XQinf, X2Qinf), Kf), 1, 0)), 0, 0),
                                  spectral_t(k, 'flat'))
        lhsO = mul(mul(XQinf, X2Qinf), KO)
        res['Thm 3.2 Om'] = sub(add(shift(lhsO, 1, 0), shift(lhsO, 2, 0)), spectral_t(k, 'om'))
        line = []
        for name, R in res.items():
            ok = is_zero(R)
            allok &= ok
            line.append(f"{name}: {'OK' if ok else 'FAIL at ' + str(first_bad(R))}")
        # negative controls
        ncs = []
        ncs.append(f"K C1 - B_(k+2) C_k first nonzero (x^d,q^n) = {first_bad(sub(mul(Kk, C1), mul(Bsum(k + 1), Ck)))}")
        ncs.append(f"spectral with perturbed (1,0) weight: {first_bad(sub(lhsK, spectral_K(k, perturb=True)))}")
        ncs.append(f"Omega last term dropped: {first_bad(sub(mul(KOm(k, drop_last=True), C1), mul(Bv, Ck)))}")
        print(f"k={k}: " + "; ".join(line), flush=True)
        print(f"      negative controls: " + "; ".join(ncs), flush=True)
    print(f"# ALL BIVARIATE CHECKS {'PASSED' if allok else 'FAILED'} (time {time.time() - t0:.1f}s)")
