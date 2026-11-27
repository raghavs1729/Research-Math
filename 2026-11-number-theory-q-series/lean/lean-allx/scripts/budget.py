#!/usr/bin/env python3
"""Budget of the analytic inequality (floating point; planning aid only).
Theta1(t) >= (2/3) log2 t^{3/2} - Sbar (t+1)^{3/2} - losses(t);  prints the margin / t^{3/2}."""
import math
g = 0.5752; logpi = 1.1436; log2 = 0.693147
def chi(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def Lam(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0: n //= p
            return math.log(p) if n == 1 else 0.0
        p += 1
    return math.log(n) if n > 1 else 0.0
def Ahead(M, k):  # sum_{n<=M} Lam(n)(-chi(n))/n^k
    return sum(Lam(n) * (-chi(n)) / n ** k for n in range(2, M + 1))
def Sbar(M):
    # M = 3 mod 4: last head index; tail over n = 3 mod 4, n >= M+4, of log n / n^2 <= (log M + 1)/(4M)
    A = Ahead(M, 2) + (math.log(M) + 1) / (4 * M)
    # B lower bound: head to 7, minus tail over n = 1 mod 4, n >= 9 : log n / n^3 <= (2 log 5 + 1)/(16*25)
    B = Ahead(7, 3) - (2 * math.log(5) + 1) / (16 * 25)
    E = 2 * A + 2 - g - logpi
    Emaj = 3.5 - g - 3 * log2 - logpi + 5 * A - 3 * B
    return E / 3 + 5 / 12 * Emaj, A, B, Emaj
def margin(t, S, SM):
    L = math.log(t)
    main = 2 / 3 * log2 * t ** 1.5
    loss = (t * (log2 + math.log(math.sqrt(t) + 2)) + t / 2 * L + 0.75 * 2 * 0.6931472 * t ** (4 / 3)
            + 4 / 7 * t ** (7 / 6) * L + 5 / 6 * 2 * 0.6931472 * t ** 1.2 + 4 / 11 * t ** 1.1 * L)
    psi1 = S * (t + 1) ** 1.5 + (3 + 3 / 7 * SM) * t + SM / 5 + L + 3 + math.log(t + 1) / 3
    return (main - loss - psi1) / t ** 1.5
for M in (31, 47, 71, 103, 199):
    S, A, B, Emaj = Sbar(M)
    print("head", M, "A<=%.5f B>=%.5f Emaj<=%.4f Sbar<=%.5f" % (A, B, Emaj, S),
          " margins:", " ".join("%g:%.4f" % (t, margin(t, S, Emaj)) for t in (1e6, 2e6, 4e6, 1e7, 1e8)))
