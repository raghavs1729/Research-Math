"""Audit 2: spot-check of the Sec. 7 search claim for (A,k,B,C,E) = (8,10,-7,1,4):
s(n) = max_{n-3<=m<=n} (ln GCD_m - ln m!)/m at n = 40, 80, 160, 320, 640 (claimed -2.04, -2.58, -2.84, -3.08 at the last four?)."""
from math import gcd, lgamma, log
A, k, B, C, E = 8, 10, -7, 1, 4
a = lambda n: A*n + k; b = lambda n: B*n*n + C*n + E
P, Q = [1, a(0)], [0, 1]
for i in range(1, 641):
    P.append(a(i)*P[-1] + b(i)*P[-2]); Q.append(a(i)*Q[-1] + b(i)*Q[-2])
p, q = P[1:], Q[1:]
s = lambda n: max((log(gcd(p[m], q[m])) - lgamma(m+1))/m for m in range(n-3, n+1))
print([round(s(n), 3) for n in (40, 80, 160, 320, 640)])
