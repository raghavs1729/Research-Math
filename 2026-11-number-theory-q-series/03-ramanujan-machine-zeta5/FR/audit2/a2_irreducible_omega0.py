"""Audit 2: Theorem D (1) when omega = 0 (b irreducible, a_n = A n + A(B+C)/(2B)).
The written proof excludes primes p | 2B*omega, i.e. ALL primes when omega = 0; check numerically that
P_n <= C^n n!^(1/2) nevertheless (the root-mod-p argument does not need p not dividing omega)."""
from math import gcd, lgamma, log
def conv(a, b, n):
    P, Q = [1, a(0)], [0, 1]
    for i in range(1, n + 1):
        P.append(a(i) * P[-1] + b(i) * P[-2]); Q.append(a(i) * Q[-1] + b(i) * Q[-2])
    return P[1:], Q[1:]
cases = [("PCF[2n+1, n^2+1]      D=8", lambda n: 2*n+1, lambda n: n*n+1),
         ("PCF[2n+1, -(n^2+1)]   D=0 (Euler shape, roots +-i)", lambda n: 2*n+1, lambda n: -(n*n+1)),
         ("PCF[2n+1, n^2-2]      D=8", lambda n: 2*n+1, lambda n: n*n-2),
         ("PCF[4n+2, -(4n^2+3)]  D=0", lambda n: 4*n+2, lambda n: -(4*n*n+3)),
         ("PCF[3n+2, n^2+n+1]... omega=2*2*1-3*2=-2 (control, omega!=0)", lambda n: 3*n+2, lambda n: n*n+n+1)]
for name, a, b in cases:
    p, q = conv(a, b, 1601)
    G = [gcd(x, y) for x, y in zip(p, q)]
    out = []
    for m in (100, 200, 400, 800, 1600):
        P = gcd(G[m], G[m+1])
        out.append(f"{(log(P) - 0.5*lgamma(m+1))/m:7.3f}/{(log(G[m]) - lgamma(m+1))/m:7.3f}")
    print(f"{name:60s} (lnP_n - .5 ln n!)/n / (lnGCD_n - ln n!)/n at n=100..1600:", " ".join(out))
