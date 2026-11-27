"""V3 (exact rational arithmetic).  With A_k = C(2k,k)^4/256^k, lam_k = 4 sum_{j<k} (1/(j+1/2) - 1/(j+1)) (so that the
Frobenius solutions of L4 at 0 are y0 = sum A_k x^k, y1 = y0 log x + sum A_k lam_k x^k):
(a) psi0 := y0*theta(y1) - y1*theta(y0) is a power series sum w_n x^n with
        w_n = sum_k A_k A_{n-k} (1 + (2k-n) lam_k) = sum_k A_k A_{n-k} (1 + (2k-n)(lam_k - lam_{n-k})/2)  >=  sum_k A_k A_{n-k} > 0;
(b) (n+1)^5 w_{n+1} + n^5 w_{n-1} = atil(n) w_n (n >= 0), atil(n) = a_n/8, and p_{n-1} = 8^n n!^5 w_n;
(c) v_n (v_0 = 0, v_1 = 1/8, same recurrence for n >= 1) gives q_{n-1} = 8^n n!^5 v_n, and
        v_{n+1}/w_{n+1} - v_n/w_n = 1/(8 (n+1)^5 w_n w_{n+1})."""
from fractions import Fraction as Fr
from math import comb, factorial
N = 60
a = lambda n: 8 * (n**5 + (n + 1)**5) - 12 * (n**3 + (n + 1)**3) + 7 * (2 * n + 1)
b = lambda n: -64 * n**10
A = [Fr(comb(2 * k, k)**4, 256**k) for k in range(N + 2)]
lam = [Fr(0)]
for k in range(1, N + 2):
    j = k - 1
    lam.append(lam[-1] + 4 * (Fr(2, 2 * j + 1) - Fr(1, j + 1)))
assert all(lam[k + 1] > lam[k] for k in range(N))
w = []
for n in range(N + 1):
    s1 = sum(A[k] * A[n - k] * (1 + (2 * k - n) * lam[k]) for k in range(n + 1))
    s2 = sum(A[k] * A[n - k] * (1 + Fr(2 * k - n, 2) * (lam[k] - lam[n - k])) for k in range(n + 1))
    assert s1 == s2 and s1 >= sum(A[k] * A[n - k] for k in range(n + 1)) > 0
    w.append(s1)
for n in range(0, N):
    assert (n + 1)**5 * w[n + 1] + (n**5 * w[n - 1] if n else 0) == Fr(a(n), 8) * w[n], n
v = [Fr(0), Fr(1, 8)]
for n in range(1, N):
    v.append((Fr(a(n), 8) * v[n] - n**5 * v[n - 1]) / (n + 1)**5)
p = [1, a(0)]; q = [0, 1]
for n in range(1, N):
    p.append(a(n) * p[-1] + b(n) * p[-2]); q.append(a(n) * q[-1] + b(n) * q[-2])
for n in range(0, N + 1):
    assert p[n] == 8**n * factorial(n)**5 * w[n] and q[n] == 8**n * factorial(n)**5 * v[n], n
for n in range(0, N):
    assert v[n + 1] / w[n + 1] - v[n] / w[n] == 1 / (8 * (n + 1)**5 * w[n] * w[n + 1])
print("V3 PASS: w_n formula and positivity, recurrence, p_{n-1} = 8^n n!^5 w_n, q_{n-1} = 8^n n!^5 v_n, Casoratian (n <= %d)" % N)
print("   v_n/w_n at n = 10, 30, 60:", float(v[10] / w[10]), float(v[30] / w[30]), float(v[60] / w[60]))
