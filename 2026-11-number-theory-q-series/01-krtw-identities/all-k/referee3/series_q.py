#!/usr/bin/env python3
"""Referee 3, check (1): KRTW identities (i)-(iv) as q-series, directly from the printed multisums.

Written from scratch (no code from krtw4/code, referee, referee2).

The multisum  sum over r_1..r_k, s_1..s_k >= 0 of  prod_i q^{e(r_i,s_i)} / ((q;q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}})
(r_{k+1} = s_0 = 0, 1/(..)_m = 0 for m < 0) is evaluated by nested summation from the innermost level i = k
outwards: F_k(r,s) = 1/(q;q)_r (the bond r_k - r_{k+1} = r_k), and
    F_i(R,S) = sum_{r <= R, s >= S} q^{e(r,s)} F_{i+1}(r,s) / ((q;q)_{R-r} (q^2;q^2)_{s-S}),
with (R,S) = (r_i,s_i), (r,s) = (r_{i+1},s_{i+1}); finally  LHS = sum_{r,s} q^{e(r,s)} F_1(r,s)/(q^2;q^2)_s.
This is just the definition of the sum, reorganised; everything is exact integer arithmetic mod q^N.

Omega for (iv) is expanded literally: the i-th term (i < k) q^{r_i+2s_i}(1 - q^{2s_{i+1}-2s_i})
= q^{r_i+2s_i} - q^{r_i+2s_{i+1}} is inserted between levels i and i+1, the last term q^{r_k+2s_k} at level k.

Usage: python3 series_q.py N KMAX
"""
import sys
import time
from flint import fmpz_poly

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10


def P(c):
    return fmpz_poly(list(c)).truncate(N) if c else fmpz_poly([])


def mul(a, b):
    return a.mul_low(b, N)


def mono(e):
    if e >= N:
        return fmpz_poly([])
    return fmpz_poly([0] * e + [1])


def div_geom(c, j):
    """multiply the coefficient list c by 1/(1-q^j) mod q^N (j >= 1)"""
    c = list(c)
    for n in range(j, N):
        c[n] += c[n - j]
    return c


# 1/(q;q)_j and 1/(q^2;q^2)_b
A_list = [[1] + [0] * (N - 1)]
B_list = [[1] + [0] * (N - 1)]
for j in range(1, N + 1):
    A_list.append(div_geom(A_list[-1], j))
    B_list.append(div_geom(B_list[-1], 2 * j))
A = [P(c) for c in A_list]
B = [P(c) for c in B_list]


def states_for(ex):
    """all (r,s) with ex(r,s) < N (ex increasing in r and s)"""
    st = []
    s = 0
    while ex(0, s) < N:
        r = 0
        while ex(r, s) < N:
            st.append((r, s))
            r += 1
        s += 1
    return st


def T_apply(F, ex, states, extra=None):
    """(TF)(R,S) = sum_{r<=R, s>=S} q^{ex(r,s)+extra(r,s)} F(r,s) A_{R-r} B_{s-S}, on the states."""
    G = {}
    for (r, s) in states:
        e = ex(r, s) + (extra(r, s) if extra else 0)
        if e < N and (r, s) in F and not F[(r, s)].is_zero():
            G[(r, s)] = F[(r, s)].left_shift(e).truncate(N)
    smax = max(s for (_, s) in states)
    rmax = max(r for (r, _) in states)
    # U(R,s) = sum_{r<=R} A_{R-r} G(r,s)
    U = {}
    for s in range(smax + 1):
        rs = [r for (r, ss) in G if ss == s]
        for R in range(rmax + 1):
            acc = fmpz_poly([])
            for r in rs:
                if r <= R:
                    acc += mul(A[R - r], G[(r, s)])
            U[(R, s)] = acc
    out = {}
    for (R, S) in states:
        acc = fmpz_poly([])
        for s in range(S, smax + 1):
            u = U.get((R, s))
            if u is not None and not u.is_zero():
                acc += mul(B[s - S], u)
        out[(R, S)] = acc
    return out


def omega(F, ex, states):
    acc = fmpz_poly([])
    for (r, s) in states:
        e = ex(r, s)
        if e < N and not F[(r, s)].is_zero():
            acc += mul(B[s], F[(r, s)].left_shift(e).truncate(N))
    return acc


def multisum(k, ex):
    st = states_for(ex)
    F = {(r, s): A[r] for (r, s) in st}
    for _ in range(k - 1):
        F = T_apply(F, ex, st)
    return omega(F, ex, st)


def multisum_omega(k, ex, drop_minus=False):
    """sum Omega * prod(...), Omega expanded literally (two monomials per i < k, one for i = k)."""
    st = states_for(ex)
    F = {(r, s): A[r] for (r, s) in st}                       # unmarked chain, level k
    FO = {(r, s): A[r].left_shift(r + 2 * s).truncate(N) for (r, s) in st}   # level k with q^{r_k+2s_k}
    for _ in range(k - 1):
        TF = T_apply(F, ex, st)                              # sum w F bond
        TYF = T_apply(F, ex, st, extra=lambda r, s: 2 * s)   # sum w q^{2 s_{i+1}} F bond
        TFO = T_apply(FO, ex, st)
        newFO = {}
        for (R, S) in st:
            t = TFO[(R, S)] + TF[(R, S)].left_shift(R + 2 * S).truncate(N)
            if not drop_minus:
                t -= TYF[(R, S)].left_shift(R).truncate(N)
            newFO[(R, S)] = t
        FO = newFO
        F = TF
    return omega(FO, ex, st)


def qprod(exps, mod):
    """prod over a in exps of (q^a; q^mod)_inf mod q^N"""
    c = [1] + [0] * (N - 1)
    for a in exps:
        e = a
        while e < N:
            # multiply by (1 - q^e)
            for n in range(N - 1, e - 1, -1):
                c[n] -= c[n - e]
            e += mod
    return c


def div_by(c, exps, mod):
    """divide by prod (q^a;q^mod)_inf"""
    for a in exps:
        e = a
        while e < N:
            c = div_geom(c, e)
            e += mod
    return c


def rhs_i(k):
    c = qprod([1, k + 1, k + 2], k + 2)
    c = div_by(c, [1], 2)       # (q;q^2)_inf
    c = div_by(c, [1], 1)       # (q;q)_inf
    return P(c)


def rhs_ii(k):
    m = 2 * k + 4
    c = qprod([1, 2, 2 * k + 2, 2 * k + 3, 2 * k + 4, 2 * k + 4], m)
    c = div_by(c, [2], 2)
    c = div_by(c, [1], 1)
    return P(c)


def rhs_iii(k):
    m = 2 * k + 4
    c = qprod([k + 1, k + 2, k + 2, k + 3, 2 * k + 4, 2 * k + 4], m)
    c = div_by(c, [2], 2)
    c = div_by(c, [1], 1)
    return P(c)


def rhs_iv(k):
    m = 2 * k + 4
    c = qprod([k, k + 1, k + 3, k + 4, 2 * k + 4, 2 * k + 4], m)
    c = div_by(c, [2], 2)
    c = div_by(c, [1], 1)
    return P(c)


def first_diff(a, b):
    d = a - b
    cs = d.coeffs()
    for i, v in enumerate(cs):
        if v != 0:
            return i
    return None


ex_iii = lambda r, s: (r + s) ** 2 + s * s
ex_ii = lambda r, s: (r + s) ** 2 + s * s + r + 2 * s
ex_i = lambda r, s: (r + s) ** 2 + s * s + s

if __name__ == "__main__":
    t0 = time.time()
    print(f"# referee3 series_q.py: N = {N} (all checks mod q^{N}), k = 1..{KMAX}")
    allok = True
    for k in range(1, KMAX + 1):
        L2 = multisum(k, ex_ii)
        L3 = multisum(k, ex_iii)
        L4 = multisum_omega(k, ex_iii)
        L1 = multisum(k, ex_i)
        ok1 = (L1 - rhs_i(k)).is_zero()
        ok2 = (L2 - rhs_ii(k)).is_zero()
        ok3 = (L3 - rhs_iii(k)).is_zero()
        ok4 = (L4 - rhs_iv(k)).is_zero()
        allok &= ok2 and ok3 and ok4
        # negative controls
        nc1 = first_diff(L4, rhs_iii(k))                       # (iv) sum vs (iii) product
        nc2 = first_diff(multisum_omega(k, ex_iii, drop_minus=True), rhs_iv(k)) if k >= 2 else "n/a(k=1)"
        nc3 = first_diff(L3, rhs_iii(k + 1))                   # wrong modulus
        nc4 = first_diff(L1, rhs_iii(k))                        # (i) sum vs (iii) product
        print(f"k={k:2d}: (ii) {'OK' if ok2 else 'FAIL'}  (iii) {'OK' if ok3 else 'FAIL'}  (iv) {'OK' if ok4 else 'FAIL'}"
              f"  [(i), not claimed: {'OK' if ok1 else 'FAIL'}]"
              f"  | neg. controls, first differing power: sum(iv) vs prod(iii): {nc1};"
              f" Omega without the '-q^(2s_(i+1)-2s_i)' part vs prod(iv): {nc2};"
              f" sum(iii) vs prod(iii) at k+1: {nc3}; sum(i) vs prod(iii): {nc4}", flush=True)
        if k <= 3:
            print(f"      first coefficients of sum(iv): {L4.coeffs()[:16]}")
    print(f"# ALL (ii),(iii),(iv) CHECKS {'PASSED' if allok else 'FAILED'}; time {time.time() - t0:.1f}s")
