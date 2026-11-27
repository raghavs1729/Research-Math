# Direct evaluation of K_k, B_m, C_k, AG_j from their definitions (no data / moves).
from ser import S, qpoch_inv_table, monomial
from flint import fmpz_poly
import itertools

def _tab(N, X):
    return qpoch_inv_table(1, X + 2, N), qpoch_inv_table(2, X + 2, N)

def K(k, X, N, xpow=None):
    """K_k(x) = sum x^{sum(r_i+2s_i)} prod q^{(r_i+s_i)^2+s_i^2}/((q)_{r_i-r_{i+1}}(Q;Q)_{s_i-s_{i-1}})"""
    t1, t2 = _tab(N, X)
    out = S(X, N)
    # enumerate r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k with weight <= X
    def rec_r(i, rs, wleft):
        if i == k:
            rec_s(0, rs, [], wleft); return
        top = rs[-1] if rs else X
        for r in range(0, top + 1):
            if r > wleft: break
            rec_r(i + 1, rs + [r], wleft - r)
    def rec_s(i, rs, ss, wleft):
        if i == k:
            e = sum((rs[j] + ss[j]) ** 2 + ss[j] ** 2 for j in range(k))
            if e >= N: return
            p = fmpz_poly([1])
            for j in range(k):
                rn = rs[j + 1] if j + 1 < k else 0
                sp = ss[j - 1] if j > 0 else 0
                p = p.mul_low(t1[rs[j] - rn], N).mul_low(t2[ss[j] - sp], N)
            out.c[X - wleft] += p.left_shift(e).truncate(N)
            return
        lo = ss[-1] if ss else 0
        s = lo
        while 2 * s * (k - i) <= wleft:
            rec_s(i + 1, rs, ss + [s], wleft - 2 * s)
            s += 1
    rec_r(0, [], X)
    return out

def B(m, X, N):
    t1, t2 = _tab(N, X)
    out = S(X, N)
    def rec(i, Ns, wleft):
        if i == m:
            e = sum(v * v for v in Ns)
            if e >= N: return
            p = fmpz_poly([1])
            for j in range(m - 1):
                p = p.mul_low(t1[Ns[j] - Ns[j + 1]], N)
            p = p.mul_low(t2[Ns[-1]], N)
            out.c[X - wleft] += p.left_shift(e).truncate(N)
            return
        top = Ns[-1] if Ns else X
        for v in range(0, top + 1):
            if v > wleft: break
            rec(i + 1, Ns + [v], wleft - v)
    rec(0, [], X)
    return out

def C(k, X, N):
    t1, t2 = _tab(N, X)
    out = S(X, N)
    def rec(i, ss, wleft):
        if i == k:
            e = sum(v * v for v in ss)
            if e >= N: return
            p = fmpz_poly([1])
            for j in range(k):
                p = p.mul_low(t2[ss[j] - (ss[j - 1] if j else 0)], N)
            out.c[X - wleft] += p.left_shift(e).truncate(N)
            return
        lo = ss[-1] if ss else 0
        v = lo
        while v * (k - i) <= wleft:
            rec(i + 1, ss + [v], wleft - v); v += 1
    rec(0, [], X)
    return out

def AG1sq(X, N):
    """sum x^{2s} q^{2s^2}/(Q;Q)_s"""
    t1, t2 = _tab(N, X)
    out = S(X, N)
    s = 0
    while 2 * s <= X:
        e = 2 * s * s
        if e < N:
            out.c[2 * s] += t2[s].left_shift(e).truncate(N)
        s += 1
    return out
