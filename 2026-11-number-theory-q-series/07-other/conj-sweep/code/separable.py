# arXiv:2404.18517 (Chen, Kitaev, Zhang), unimodality conjectures 1-3 for rmax on separable permutations.
# A = all separable, P = (+)-indecomposable ("irreducible"), M = (-)-indecomposable; x marks size, t marks rmax.
#   rmax(a (+) b) = rmax(b),  rmax(a (-) b) = rmax(a) + rmax(b).
#   P = x t + M^2/(1-M),   M = x t + P * P1/(1-P1)  (P1 = P at t=1),   A = P + P*P1/(1-P1).
import sys
from flint import fmpz_poly as Z
N = int(sys.argv[1])
zero = Z([0])
P = [zero]*(N+1); M = [zero]*(N+1); R = [zero]*(N+1)   # R = M/(1-M)
P1 = [0]*(N+1); S1 = [0]*(N+1)                           # S1 = P1/(1-P1) (integers)
PS = [zero]*(N+1)                                        # P * S1
for n in range(1, N+1):
    # terms of order n of M^2/(1-M) = M*R and of P*S1 only involve lower orders
    Q = zero
    for i in range(1, n): Q = Q + M[i]*R[n-i]
    PSn = zero
    for i in range(1, n): PSn = PSn + P[i]*S1[n-i]
    base = Z([0, 1]) if n == 1 else zero
    P[n] = base + Q; M[n] = base + PSn; PS[n] = PSn
    R[n] = M[n] + Q                                      # R = M + M*R
    P1[n] = int(P[n](1))
    S1[n] = P1[n] + sum(P1[i]*S1[n-i] for i in range(1, n))
A = [P[n] + PS[n] for n in range(N+1)]
def coeffs(p, n): 
    c = [int(x) for x in p.coeffs()]; c += [0]*(n+1-len(c)); return c[1:n+1]      # k = 1..n
assert coeffs(A[8], 8) == [1806, 2240, 2046, 1475, 745, 217, 28, 1], coeffs(A[8], 8)
assert coeffs(P[8], 8) == [0, 1092, 1288, 1069, 607, 195, 27, 1], coeffs(P[8], 8)
def unimodal_peak(c):
    pk = max(range(len(c)), key=lambda i: c[i])
    ok = all(c[i] <= c[i+1] for i in range(pk)) and all(c[i] >= c[i+1] for i in range(pk, len(c)-1))
    strict_unique = sum(1 for x in c if x == c[pk]) == 1
    return ok, pk+1, strict_unique
bad1 = []; bad2 = []; bad3 = []
for n in range(1, N+1):
    ok, pk, u = unimodal_peak(coeffs(A[n], n))
    if n >= 3 and not (ok and pk == 2): bad1.append((n, ok, pk))
    ok, pk, u = unimodal_peak(coeffs(P[n], n))
    if n >= 5 and not (ok and pk == 3): bad2.append((n, ok, pk))
    red = [a-b for a, b in zip(coeffs(A[n], n), coeffs(P[n], n))]
    ok, pk, u = unimodal_peak(red) if n >= 2 else (True, 1, True)
    if n >= 2 and not (ok and pk == 1): bad3.append((n, ok, pk))
print(f'n <= {N}: Conj 1 (all separable, peak k=2, n>=3) failures (n, unimodal, peak): {bad1[:8]}')
print(f'n <= {N}: Conj 2 (rmax on irreducible, peak k=3, n>=5) failures: {bad2[:8]}')
print(f'n <= {N}: Conj 3 (rmax on reducible, peak k=1) failures: {bad3[:8]}')
