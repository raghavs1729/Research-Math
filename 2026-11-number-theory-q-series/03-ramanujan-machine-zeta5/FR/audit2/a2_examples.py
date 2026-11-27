"""Audit 2: exact checks of the examples of THEOREM.md Sec. 6, and of two further fractions:

 (N1) PCF[n, 2(n+1)(n+2)]: a_n = A n + k with k = 0, D = A^2+4B = 9 a nonzero square, omega != 0.
      Printed Conjecture 1.4 requires k in Q\\{0} for (Eq. 6), so a_n = n is in neither printed family
      (Eq. 7 would need a_n = k'(2n+4)).  Non-terminating (b_n > 0).  Is GCD_n (and P_n) factorial?
      Also the paper's own example PCF[7n + k, 8n^2 - 2] at k = 0 ("any k in Z", paper Sec. 2.3).
 (N2) PCF[1, (n+1)^2]: A = 0 (deg a = 0), D = 4 a nonzero square: P_n factorial although outside
      Theorem D's hypotheses -- relevant to the Sec. 7 claim "outside the families ... P_n <= C^n n!^(1/2)".
"""
from math import gcd, lgamma, log


def conv(a, b, n):
    P, Q = [1, a(0)], [0, 1]           # index shift: P[i+1] = p_i
    for i in range(1, n + 1):
        P.append(a(i) * P[-1] + b(i) * P[-2])
        Q.append(a(i) * Q[-1] + b(i) * Q[-2])
    return P[1:], Q[1:]


def report(name, a, b, n=400, pts=(50, 100, 200, 400)):
    p, q = conv(a, b, n)
    G = [gcd(x, y) for x, y in zip(p, q)]
    print(name)
    for m in pts:
        Pm = gcd(G[m], G[m + 1]) if m + 1 <= n else None
        s = (log(G[m]) - lgamma(m + 1)) / m if G[m] else float('-inf')
        sp = (log(Pm) - lgamma(m + 1)) / m if Pm else float('-inf')
        print(f"   n={m:4d}  (ln GCD_n - ln n!)/n = {s:8.4f}   (ln P_n - ln n!)/n = {sp:8.4f}")
    return p, q, G


# Sec. 6 examples
p, q, G = report("Sec6.1 PCF[n+1, n(n-1)]", lambda n: n + 1, lambda n: n * (n - 1), n=200, pts=(50, 100, 199))
from math import factorial
print("   p_n == q_n for n>=0:", all(x == y for x, y in zip(p, q)),
      " GCD_n >= (n+1)! for n<=200:", all(G[m] >= factorial(m + 1) for m in range(201)))
p, q, G = report("Sec6.2 PCF[n^2+1, (n-1)^2(n^2-2)]", lambda n: n * n + 1, lambda n: (n - 1) ** 2 * (n * n - 2), n=120, pts=(50, 100, 119))
print("   p_n == q_n:", all(x == y for x, y in zip(p, q)), " GCD_n >= n!^2:", all(G[m] >= factorial(m) ** 2 for m in range(121)))
for name, a, b in [("Sec6.3a PCF[n-3, 2(n-1)(n-3)]", lambda n: n - 3, lambda n: 2 * (n - 1) * (n - 3)),
                   ("Sec6.3b PCF[n-1, -(n-1)(n-2)]", lambda n: n - 1, lambda n: -(n - 1) * (n - 2))]:
    p, q = conv(a, b, 60)
    print(name, "first n with p_n=q_n=0:", next(i for i in range(61) if p[i] == 0 == q[i]),
          " all zero after:", all(p[i] == 0 == q[i] for i in range(next(i for i in range(61) if p[i] == 0 == q[i]), 61)))

# Family-membership arithmetic for Sec. 6.3a: A = 1 = B/m - m with B = 2 -> m in {1, -2}
print("Sec6.3a: m with 2/m - m = 1:", [m for m in (1, -2) if 2 / m - m == 1], " sqrt|B| =", 2 ** 0.5)

# New: k = 0 members of (Eq. 6)
report("N1 PCF[n, 2(n+1)(n+2)]   (k=0, D=9, omega=-8)", lambda n: n, lambda n: 2 * (n + 1) * (n + 2), n=800, pts=(100, 200, 400, 799))
report("N1' PCF[7n, 8n^2-2]       (paper's Eq.6 example, m=1, k=0)", lambda n: 7 * n, lambda n: 8 * n * n - 2, n=800, pts=(100, 200, 400, 799))
report("control PCF[n+1, 2(n+1)(n+2)] (D=9, k=1: in Eq.6)", lambda n: n + 1, lambda n: 2 * (n + 1) * (n + 2), n=800, pts=(100, 200, 400, 799))
report("control outside: PCF[n, 3(n+1)(n+2)] (D=13)", lambda n: n, lambda n: 3 * (n + 1) * (n + 2), n=800, pts=(100, 200, 400, 799))
# New: A = 0, B square
report("N2 PCF[1, (n+1)^2]  (A=0, D=4)", lambda n: 1, lambda n: (n + 1) ** 2, n=800, pts=(100, 200, 400, 799))
report("N2' PCF[3, 4n^2+4n+1]=PCF[3,(2n+1)^2] (A=0, D=16)", lambda n: 3, lambda n: (2 * n + 1) ** 2, n=800, pts=(100, 200, 400, 799))
