"""Referee re-computation for R05 (arXiv:2502.21050): Hankel determinants H_n(M(x)^r), H_n = det(a_{i+j})_{0<=i,j<n}.
Own code: Motzkin numbers from M_n = M_{n-1} + sum_{k} M_k M_{n-2-k} (functional equation M = 1 + xM + x^2 M^2),
checked against M_n = sum_k binom(n,2k) C_k (the source's definition); determinants by (1) fraction-free Bareiss on
Python integers and (2) exact rational Gaussian elimination (fractions) -- two different eliminations."""
from fractions import Fraction
from math import comb
def motzkin(D):
    M = [1]
    for n in range(1, D):
        M.append(M[n - 1] + sum(M[k] * M[n - 2 - k] for k in range(0, n - 1)))
    return M
D = 80
M = motzkin(D)
assert all(M[n] == sum(comb(n, 2 * k) * comb(2 * k, k) // (k + 1) for k in range(n // 2 + 1)) for n in range(D))
assert M[:8] == [1, 1, 2, 4, 9, 21, 51, 127]
def mul(a, b): return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(D)]
POW = {0: [1] + [0] * (D - 1)}
def mpow(r):
    if r not in POW: POW[r] = mul(mpow(r - 1), M)
    return POW[r]
def det_bareiss(A):
    n = len(A)
    if n == 0: return 1
    A = [row[:] for row in A]; sign = 1; prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            p = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if p is None: return 0
            A[k], A[p] = A[p], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]
def det_frac(A):
    n = len(A); A = [[Fraction(x) for x in row] for row in A]; d = Fraction(1)
    for k in range(n):
        p = next((i for i in range(k, n) if A[i][k] != 0), None)
        if p is None: return 0
        if p != k: A[k], A[p] = A[p], A[k]; d = -d
        d *= A[k][k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            if f:
                for j in range(k, n): A[i][j] -= f * A[k][j]
    assert d.denominator == 1
    return int(d)
def H(r, n, check=False):
    a = mpow(r); assert 2 * n - 2 < D
    A = [[a[i + j] for j in range(n)] for i in range(n)]
    d = det_bareiss(A)
    if check: assert d == det_frac(A)
    return d
# sanity: r = 1, 2, 3 (source's theorems for r <= 3)
assert [H(1, n) for n in range(12)] == [1] * 12
assert [H(2, n) for n in range(24)] == [1, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 0] * 2
assert all(H(3, 3 * n) == H(3, 3 * n + 1) == (-1) ** n * (n + 1) ** 2 and H(3, 3 * n + 2) == 0 for n in range(8))
assert all(H(6, 6 * n) == H(6, 6 * n + 1) == (-1) ** n * (n + 1) ** 5 for n in range(5))
print('sanity: source theorems for r = 1, 2, 3 and H_{6n}=H_{6n+1} for r = 6 reproduced')
print('H_9(M^9), H_10(M^9) =', H(9, 9, True), H(9, 10, True), '   printed conjecture: (-1)^1 * 2^8 =', -2 ** 8)
print('H_7(M^4), H_10(M^4) =', H(4, 7, True), H(4, 10, True), '   sum =', H(4, 7) + H(4, 10), '   printed conjecture: +-(2*4*1)^2 = +-64;  source Thm (r=4): both -64')
print('H_0..H_13(M^4) =', [H(4, n, True) for n in range(14)])
thm4 = lambda n: [1, 1, 2 * (2 * n + 1) * (8 * n - 1), -2 * (8 * n + 9) * (2 * n + 1), -1, -1, 4 * (n + 1) * (8 * n + 5), -64 * (n + 1) ** 2, 0, 0, -64 * (n + 1) ** 2, -4 * (n + 1) * (8 * n + 11)]
for n in range(3):
    got = [H(4, 12 * n + j) for j in range(12)]
    print(f'   source Thm r=4, n={n}: residues j where printed formula != true value:', [(j, thm4(n)[j], got[j]) for j in range(12) if thm4(n)[j] != got[j]])
print('H_9(M^5), H_12(M^5) =', H(5, 9, True), H(5, 12, True), '  (source Thm r=5 prints +1000 for both)')
print('H_12(M^12), H_13(M^12) =', H(12, 12), H(12, 13), '  printed: ', -2 ** 11)
print('r = 0 mod 3, sign of H_{rn} relative to (-1)^n (n+1)^(r-1), n = 1, 2, 3:')
for r in range(3, 22, 3):
    row = []
    for n in (1, 2, 3):
        if 2 * (r * n + 1) - 2 < D:
            h0, h1 = H(r, r * n), H(r, r * n + 1); pred = (-1) ** n * (n + 1) ** (r - 1)
            row.append((n, 'ok' if h0 == h1 == pred else ('SIGN' if h0 == h1 == -pred else f'OTHER {h0},{h1}')))
    print('   r =', r, row)
print('second conjecture, first line (r=4, n=0): H_0, H_1, H_4, H_5 =', H(4, 0), H(4, 1), H(4, 4), H(4, 5))
print('second conjecture, line 3, n=0: H_{2r-1} + H_{2r+2} for r = 4,5,7,8,10:', [(r, H(r, 2 * r - 1), H(r, 2 * r + 2), (2 * r) ** (r - 2)) for r in (4, 5, 7, 8, 10)])
print('r = 21, 24, 27 (n = 1): H_r, H_{r+1} and printed prediction -(2^(r-1)):', [(r, H(r, r), H(r, r + 1), -2 ** (r - 1)) for r in (21, 24, 27)])
print('second conjecture, line 3, n = 0 and 1, all r <= 13 with 3 !| r: sum H_{3rn+2r-1} + H_{3rn+2r+2} =', [(r, n, H(r, 3 * r * n + 2 * r - 1) + H(r, 3 * r * n + 2 * r + 2)) for r in (4, 5, 7, 8, 10, 11, 13) for n in (0, 1) if 2 * (3 * r * n + 2 * r + 2) - 2 < D])
