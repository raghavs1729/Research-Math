#!/usr/bin/env python3
"""Independent computation (plain Python, no Lean) of the first coefficients of both sides of the
identities of lean-nahm, for scripts/Sanity.lean.

* Matrices: the Cartan matrices are built from the PAPER's Dynkin diagrams (E8: chain 1-...-7, node 8
  on node 5; E7: chain 1-...-6, node 7 on node 4), inverted with exact fractions; the linear terms are
  the paper's B^{(i)} (paper labelling).  The Lean side uses Mathlib's Bourbaki labelling, so this also
  checks the relabelling.
* Sums: depth-first enumeration of n with exponent <= N (all entries of C^{-1} are positive, so the
  exponent is monotone in each n_i), summand q^e / prod (q;q)_{n_i} expanded as a truncated series.
* Products: prod over n <= N, n mod 2k+1 not in {0, i, -i} of (1-q^n)^{-1}; theta series directly.
Series are lists of integers: coefficients of X^0..X^N (X = q for AG/E8, X = q^{1/4} for E7)."""
from fractions import Fraction as F

N = 40

def mul(a, b):
    c = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                c[i + j] += x * b[j]
    return c

def inv(a):
    assert a[0] == 1
    b = [0] * (N + 1); b[0] = 1
    for n in range(1, N + 1):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1))
    return b

def mono(e):
    c = [0] * (N + 1)
    if e <= N: c[e] = 1
    return c

def poch(n, d=1):  # (q;q)_n with q = X^d
    c = mono(0)
    for i in range(1, n + 1):
        f = mono(0); 
        if d * i <= N: f[d * i] -= 1
        c = mul(c, f)
    return c

def cartan(nodes, edges):
    C = [[0] * nodes for _ in range(nodes)]
    for i in range(nodes): C[i][i] = 2
    for a, b in edges: C[a - 1][b - 1] = C[b - 1][a - 1] = -1
    return C

def matinv(M):
    n = len(M); A = [[F(x) for x in r] + [F(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]

def nahm(A, B, scale, d):
    """sum over n in N^r of X^{scale*(n^T A n/2 + B.n)} / prod (X^d;X^d)_{n_i}"""
    r = len(A); tot = [0] * (N + 1)
    def expo(n):
        e = F(0)
        for i in range(r):
            for j in range(r):
                e += A[i][j] * n[i] * n[j]
        e = e / 2 + sum(B[i] * n[i] for i in range(r))
        e = e * scale
        assert e.denominator == 1
        return int(e)
    def rec(n):
        nonlocal tot
        if len(n) == r:
            t = mono(expo(n))
            for x in n: t = mul(t, inv(poch(x, d)))
            tot = [a + b for a, b in zip(tot, t)]
            return
        x = 0
        while True:
            m = n + [x] + [0] * (r - len(n) - 1)
            if expo(m) > N: break
            rec(n + [x]); x += 1
    rec([])
    return tot

def ag_matrix(r):  # (M_1^2+...+M_r^2) = n^T A n / 2 with A = (2 min(i,j))
    return [[2 * min(i, j) for j in range(1, r + 1)] for i in range(1, r + 1)]

def ag_lin(k, i):  # N_i + ... + N_{k-1} = sum_j max(0, j - i + 1) n_j
    return [max(0, j - i + 1) for j in range(1, k)]

def ag_product(k, i):
    D = 2 * k + 1; c = mono(0)
    for n in range(1, N + 1):
        if n % D not in (0, i, D - i):
            f = mono(0); f[n] -= 1; c = mul(c, f)
    return inv(c)

out = {}
for k, i in [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 2), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5)]:
    s = nahm(ag_matrix(k - 1), ag_lin(k, i), 1, 1); p = ag_product(k, i)
    assert s == p, (k, i)
    out[f'ag_{k}_{i}'] = s
# E8, paper labelling
C8 = cartan(8, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8)])
A8 = matinv(C8)
B8 = {5: [0] * 8, 4: [0, 0, 0, 1, 1, 1, 0, 1], 3: [0, 1, 1, 2, 2, 2, 0, 2], 2: [0, 2, 2, 4, 4, 3, 1, 3],
      1: [1, 3, 4, 6, 7, 5, 2, 4]}
for i in range(1, 6):
    s = nahm(A8, B8[i], 1, 1)
    assert s == ag_product(5, i), i
    out[f'e8_{i}'] = s
# E7, paper labelling, X = q^{1/4}
C7 = cartan(7, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (4, 7)])
A7 = matinv(C7)
s7 = nahm(A7, [0] * 7, 4, 4)
pinf4 = poch(N, 4)          # (X^4;X^4)_N = (X^4;X^4)_oo mod X^{N+1}
th = [0] * (N + 1)
for k in range(0, N + 1):
    if 3 * k * k <= N: th[3 * k * k] += 1
    if k >= 1 and 7 * k * k <= N: th[7 * k * k] += 1
rhs7 = mul(inv(pinf4), th)
assert s7 == rhs7
bil = [0] * (N + 1)
for k in range(-N, N + 1):
    for c in (3, 7):
        if c * k * k <= N: bil[c * k * k] += 1
assert [2 * x for x in s7] == mul(inv(pinf4), bil)
out['e7'] = s7; out['e7_theta_bilateral'] = bil
for key, v in out.items():
    print(key, '=', v)
print('all identities agree mod X^%d (independent computation)' % (N + 1))
