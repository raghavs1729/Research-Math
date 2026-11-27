#!/usr/bin/env python3
"""Independent computation (plain Python integers, no Lean) of the first coefficients of
  LHS(k)  = sum_{r_1>=...>=r_k>=0, 0<=s_1<=...<=s_k} prod_i q^{(r_i+s_i)^2+s_i^2+s_i}
                / ((q;q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}})          (the printed KRTW sum (i))
  MEL(k)  = sum_{M_1>=...>=M_k>=0} prod_i q^{M_i(M_i+1)/2}/(q;q)_{M_i-M_{i+1}}
  PROD(k) = (q,q^{k+1},q^{k+2};q^{k+2})_inf / ((q;q^2)_inf (q;q)_inf)
for k = 1, 2, 3, modulo q^31.  The sums are enumerated by a depth-first search over the chains,
straight from the printed formulas.  Output: the lists used in scripts/Sanity.lean.
Usage: python3 scripts/sanity_independent.py"""
N = 30

def mul(a, b):
    c = [0] * (N + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(N + 1 - i):
                c[i + j] += ai * b[j]
    return c

def inv_poch(n, d):          # 1/(q^d;q^d)_n mod q^(N+1)
    c = [1] + [0] * N
    for m in range(1, n + 1):
        for i in range(d * m, N + 1):
            c[i] += c[i - d * m]
    return c

def shift(a, e):
    return [0] * e + a[:N + 1 - e] if e <= N else [0] * (N + 1)

def lhs(k):
    total = [0] * (N + 1)
    def rec(i, r_prev, s_prev, e, rs, ss):
        # choose (r_i, s_i) with r_i <= r_prev, s_i >= s_prev
        if i == k:
            term = [1] + [0] * N
            for t in range(k):
                rn = rs[t + 1] if t + 1 < k else 0
                sp = ss[t - 1] if t > 0 else 0
                term = mul(term, inv_poch(rs[t] - rn, 1))
                term = mul(term, inv_poch(ss[t] - sp, 2))
            term = shift(term, e)
            for n in range(N + 1):
                total[n] += term[n]
            return
        rmax = r_prev if r_prev is not None else N
        for r in range(rmax + 1):
            s = s_prev
            while True:
                w = (r + s) ** 2 + s * s + s
                if e + w > N:
                    break
                rec(i + 1, r, s, e + w, rs + [r], ss + [s])
                s += 1
    rec(0, None, 0, 0, [], [])
    return total

def mel(k):
    total = [0] * (N + 1)
    def rec(i, m_prev, e, ms):
        if i == k:
            term = [1] + [0] * N
            for t in range(k):
                mn = ms[t + 1] if t + 1 < k else 0
                term = mul(term, inv_poch(ms[t] - mn, 1))
            term = shift(term, e)
            for n in range(N + 1):
                total[n] += term[n]
            return
        mmax = m_prev if m_prev is not None else N
        for m in range(mmax + 1):
            w = m * (m + 1) // 2
            if e + w > N:
                break
            rec(i + 1, m, e + w, ms + [m])
    rec(0, None, 0, [])
    return total

def qinf(a, d):              # (q^a;q^d)_inf mod q^(N+1)
    c = [1] + [0] * N
    e = a
    while e <= N:
        c = [c[i] - (c[i - e] if i >= e else 0) for i in range(N + 1)]
        e += d
    return c

def inv(a):                  # 1/a mod q^(N+1), a[0] = 1
    b = [0] * (N + 1)
    b[0] = 1
    for n in range(1, N + 1):
        b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1))
    return b

def prod(k):
    num = mul(mul(qinf(1, k + 2), qinf(k + 1, k + 2)), qinf(k + 2, k + 2))
    den = mul(qinf(1, 2), qinf(1, 1))
    return mul(num, inv(den))

if __name__ == '__main__':
    for k in (1, 2, 3):
        a, b, c = lhs(k), mel(k), prod(k)
        assert a == b == c, (k, a, b, c)
        print(f"k = {k}: LHS = MEL = PROD mod q^{N+1}:")
        print("  ", a)
