#!/usr/bin/env python3
"""
explicit_failure.py q a k [k2 ...] -- Arb check of inequality (3) of Akbary-Francis at
explicit primorials: with N = N_k = product of the first k primes = a mod q (prime q),

   lhs(k) = N / ( phi(N) * (log(phi(q) log N))^(1/phi(q)) ) = prod p/(p-1) / (log((q-1) sum log p))^(1/(q-1)),
   (3):  lhs(k) > 1/C(q,a).

Prints lhs(k), 1/C(q,a) and whether (3) is certified to HOLD or to FAIL.
C(q,a) from mertens_const.logC (Arb; for a = 1 its value does not depend on the
branch check noted there: only real parts of log L(1,chi) enter).
"""
import sys
from flint import arb, ctx
from mertens_const import logC


def primes_in_class(q, a, k):
    out, n = [], a if a > 1 else q + 1
    while len(out) < k:
        if n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1)):
            out.append(n)
        n += q
    return out


if __name__ == "__main__":
    ctx.prec = 160
    q, a = int(sys.argv[1]), int(sys.argv[2])
    ks = [int(x) for x in sys.argv[3:]]
    lc = logC(q, a)
    rhs = (-lc).exp()
    ps = primes_in_class(q, a, max(ks))
    for k in ks:
        prod, th = arb(1), arb(0)
        for p in ps[:k]:
            prod *= arb(p) / (p - 1)
            th += arb(p).log()
        lhs = prod / ((q - 1) * th).log() ** (arb(1) / (q - 1))
        verdict = "HOLDS" if lhs > rhs else ("FAILS" if lhs < rhs else "undecided")
        print("q=%d a=%d k=%d p_k=%d : lhs = %s  1/C = %s  lhs - 1/C = %s  => (3) %s"
              % (q, a, k, ps[k - 1], lhs.str(18, radius=False), rhs.str(18, radius=False), (lhs - rhs).str(6, radius=False), verdict))
