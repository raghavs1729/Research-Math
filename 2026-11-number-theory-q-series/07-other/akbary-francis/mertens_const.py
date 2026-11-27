#!/usr/bin/env python3
"""
mertens_const.py -- the constant C(q, a) of Mertens' theorem in progressions
(Akbary-Francis Thm 1.8, after Languasco-Zaccagnini):
    C(q,a)^phi(q) = e^(-gamma) prod_p (1 - 1/p)^alpha(p;q,a),
    alpha = phi(q) - 1 if p = a mod q, -1 otherwise,
for a prime modulus q, in Arb.

    phi(q) log C(q,a) = -gamma - log(1 - 1/q) + sum_{chi != chi_0} conj(chi(a)) S(chi),
    S(chi) = sum_p chi(p) log(1 - 1/p)
           = -log L(1, chi) + sum_{k>=2} (1/k) [ P(k, chi^k) - P(k, chi) ],
    P(s, psi) = sum_p psi(p) p^(-s) = sum_{m>=1} mu(m)/m log L(ms, psi^m)   (s >= 2),
with log L(1, chi) = sum_p sum_k chi(p)^k/(k p^k) = the branch continuous on [1, oo)
vanishing at +oo; we check Re L(sigma, chi) > 0 at 257 sample points of [1, 2] (and
|L - 1| < 1 for sigma >= 2), so it is the principal logarithm (sampled check only).
Truncation: k <= KMAX, m*k <= KMAX; the neglected terms are bounded by
sum_{j > KMAX} d(j) * 2 * 2^(-j) * 2 < 2^(-KMAX + 12) (crude), added as a ball.
Used only for the numerical sanity check with primes (step 5), not in the proof.
"""
import sys
from flint import arb, acb, ctx, dirichlet_char

KMAX = 200


def mobius(n):
    r, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            r = -r
        p += 1
    return -r if n > 1 else r


def logL(q, n, s):
    """principal log of L(s, chi_q(n,.)) for integer s >= 2 (|L-1| < 1); n may be 1 (principal)."""
    if n % q == 1:
        return acb((arb(s).zeta() * (1 - arb(q) ** (-s))).log())
    return dirichlet_char(q, n).l(acb(s)).log()


def S(q, n):
    chi = dirichlet_char(q, n)
    # Re L > 0 on [1,2]
    # (sampled, not a covering: this script is a sanity-check tool, not part of the proof)
    N = 256
    for i in range(N + 1):
        assert chi.l(acb(arb(1) + arb(i) / N)).real > 0, "Re L(sigma, chi) > 0 fails"
    cache = {}

    def lg(s, nn):
        key = (s, nn % q)
        if key not in cache:
            cache[key] = logL(q, nn % q, s)
        return cache[key]

    def P(s, nn):
        tot = acb(0)
        m = 1
        while m * s <= KMAX:
            mu = mobius(m)
            if mu:
                tot += arb(mu) / m * lg(m * s, pow(nn, m, q))
            m += 1
        return tot

    tot = -chi.l(acb(1)).log()
    for k in range(2, KMAX + 1):
        tot += (P(k, pow(n, k, q)) - P(k, n)) / k
    err = arb(2) ** (-(KMAX - 12))
    return tot + acb(arb(0, err.upper()), arb(0, err.upper()))


def logC(q, a):
    phi = q - 1
    tot = acb(-arb.const_euler() - (1 - arb(1) / q).log())
    for n in range(2, q):
        chi = dirichlet_char(q, n)
        tot += chi(a).conjugate() * S(q, n)
    assert tot.imag.contains(0)
    return tot.real / phi


if __name__ == "__main__":
    ctx.prec = 160
    q, a = int(sys.argv[1]), int(sys.argv[2])
    lc = logC(q, a)
    print("log C(%d,%d) = %s" % (q, a, lc.str(35)))
    print("1/C(%d,%d)  = %s" % (q, a, (-lc).exp().str(35)))
