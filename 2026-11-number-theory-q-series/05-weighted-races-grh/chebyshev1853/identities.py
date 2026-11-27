#!/usr/bin/env python3
"""
identities.py -- numerical checks (mpmath, 30 digits; NOT certified) of the exact identities used in CHEBYSHEV1853.md.

 1. Explicit formula (Lemma 3):
      psi_e(y,chi) := sum_n chi(n) Lambda(n) e^{-n/y}
                    = - sum_rho Gamma(rho) y^rho - b (1 + 1/y) e^{-1/y} + [(1+b0) e^{-1/y} + E_1(1/y)]/y + R_1(y) - R_4(y),
      R_1(y) = y^{-2} sum_rho (rho(rho+1))^{-1} int_0^1 t^{rho+1} e^{-t/y} dt,   |R_1| <= (2/5) sigma_0^2 / y^2  (GRH),
      R_4(y) = y^{-2} int_1^oo sum_{n>=1} t^{-2n}/(2n(2n+1)) e^{-t/y} dt,          0 <= R_4 <= (log 2 - 1/2)/y^2.
    (zero sums over the listed zeros only: the neglected part of R_1 is O(1e-8/y^2), that of sum Gamma(rho) y^rho < 1e-400.)
 2. Laplace representation (Lemma 5):  1/log t = int_0^oo e^{-lambda t} mu(lambda) dlambda (t > 1),
      mu(lambda) = int_0^oo lambda^{u-1}/Gamma(u) du;   and the lower bound  mu(lambda) >= (1 - lambda + lambda log lambda)/(lambda log^2 lambda).
 3. Mixture identity (Lemma 5):  -C(x) = int_0^oo W(1/(lambda + 1/x)) mu(lambda) dlambda.
 4. Kernel identity (Lemma 1):   -C(x) = int_{p0}^oo Theta_1(t) F''(t) dt,  F(t) = e^{-t/x}/log t;
                                 W(y)  = y^{-2} int_{p0}^oo Theta_1(t) e^{-t/y} dt.
Usage: python3 identities.py [q]
"""
import sys, os
import mpmath as mp
mp.mp.dps = 30
q = int(sys.argv[1]) if len(sys.argv) > 1 else 4
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, '..', '..', 'prototype', 'explore')
gam = [mp.mpf(l.split()[0]) for l in open(os.path.join(ZDIR, 'chi%d_zeros.txt' % q)) if l.strip() and l[0] != '#']


def chi(n):
    if q == 4:
        return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
    return 0 if n % 3 == 0 else (1 if n % 3 == 1 else -1)


def primes_upto(n):
    bs = bytearray([1]) * (n + 1); bs[0] = bs[1] = 0
    i = 2
    while i * i <= n:
        if bs[i]:
            bs[i * i::i] = bytearray(len(bs[i * i::i]))
        i += 1
    return [i for i in range(2, n + 1) if bs[i]]


PR = primes_upto(200000)
LOGP = {p: mp.log(p) for p in PR if p < 20000}


def chivec():
    return [0, 1, 0, -1] if q == 4 else [0, 1, -1]


def LL(s):
    return mp.dirichlet(s, chivec(), derivative=1) / mp.dirichlet(s, chivec())


g = mp.euler
LL1 = LL(1)
sig0 = mp.log(q / mp.pi) - g + 2 * LL1
b = -mp.log(q / mp.pi) + g + mp.log(2) - LL1
b0 = -mp.log(q / mp.pi) + g / 2 - mp.digamma(mp.mpf(3) / 2) / 2 - LL(2)


def psi_e_direct(y):
    tot = mp.mpf(0)
    lim = 80 * y
    for p in PR:
        if p > lim:
            break
        lp = mp.log(p); pk = p; k = 1
        while pk <= lim:
            tot += chi(pk) * lp * mp.exp(-mp.mpf(pk) / y)
            pk *= p; k += 1
    return tot


def psi_e_formula(y, nz=None):
    y = mp.mpf(y)
    zs = gam if nz is None else gam[:nz]
    Z = mp.mpf(0); R1 = mp.mpf(0)
    for gm in zs:
        rho = mp.mpc(0.5, gm)
        Z += 2 * mp.re(mp.gamma(rho) * y ** rho)
        R1 += 2 * mp.re(y ** rho * mp.gammainc(rho + 2, 0, 1 / y) / (rho * (rho + 1)))    # y^{-2} * y^{rho+2} gamma(rho+2, 1/y)
    # sum_{n>=1} u^{2n}/(2n(2n+1)) = 1 - (1/2) log(1-u^2) - artanh(u)/u   (u = 1/t < 1)
    ser = lambda t: 1 - mp.log(1 - 1 / (t * t)) / 2 - t * mp.atanh(1 / t)
    R4 = mp.quad(lambda t: ser(t) * mp.exp(-t / y), [1, 1.01, 1.1, 2, 5, 20, mp.inf]) / y ** 2
    main = -Z - b * (1 + 1 / y) * mp.exp(-1 / y) + ((1 + b0) * mp.exp(-1 / y) + mp.e1(1 / y)) / y
    return main + R1 - R4, Z, R1, R4


def W(y):
    y = mp.mpf(y)
    tot = mp.mpf(0)
    for p in PR:
        if p > 80 * y:
            break
        c = chi(p)
        if c:
            tot -= c * LOGP[p] * mp.exp(-p / y)
    return tot


def C(x):
    x = mp.mpf(x)
    tot = mp.mpf(0)
    for p in PR:
        if p > 80 * x:
            break
        tot += chi(p) * mp.exp(-p / x)
    return tot


def mu(lam):
    lam = mp.mpf(lam)
    return mp.quad(lambda u: lam ** (u - 1) * mp.rgamma(u), [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128])


if __name__ == '__main__':
    ONLY1 = len(sys.argv) > 2 and sys.argv[2] == 'part1'
    print("q = %d: sigma_0^2 = %s, b = %s, b0 = %s, zeros used: %d" % (q, mp.nstr(sig0, 15), mp.nstr(b, 15), mp.nstr(b0, 15), len(gam)))
    print("\n1. explicit formula for psi_e(y,chi):")
    print("   %8s %22s %22s %10s %12s %12s %12s" % ("y", "direct", "formula", "diff", "zero sum", "R_1", "R_4"))
    for y in (1, 2, 5, 20, 100, 1000):
        d = psi_e_direct(y)
        f, Z, R1, R4 = psi_e_formula(y)
        print("   %8g %22s %22s %10s %12s %12s %12s  [bounds: |R_1| <= %s, R_4 <= %s]"
              % (y, mp.nstr(d, 18), mp.nstr(f, 18), mp.nstr(d - f, 3), mp.nstr(Z, 6), mp.nstr(R1, 6), mp.nstr(R4, 6),
                 mp.nstr(mp.mpf(2) / 5 * sig0 / mp.mpf(y) ** 2, 4), mp.nstr((mp.log(2) - mp.mpf(1) / 2) / mp.mpf(y) ** 2, 4)))
    if ONLY1:
        sys.exit(0)
    print("\n2. mu(lambda) and its Laplace transform:")
    for lam in ('1e-8', '1e-4', '0.01', '0.3', '0.9'):
        l = mp.mpf(lam)
        lowb = (1 - l + l * mp.log(l)) / (l * mp.log(l) ** 2)
        a = mp.log(1 / l)
        print("   lambda = %s: mu = %s, lower bound = %s, ratio = %s ; (1 + 2 gamma/a) = %s" % (lam, mp.nstr(mu(l), 12), mp.nstr(lowb, 12), mp.nstr(mu(l) / lowb, 8), mp.nstr(1 + 2 * g / a, 8)))
    mp.mp.dps = 15
    for t in (3, 10):
        v = mp.quad(lambda a: mp.exp(-mp.exp(-a) * t) * mu(mp.exp(-a)) * mp.exp(-a), [0, 2, 5, 10, 20, 40, 80, 200, 1000, 10000]) \
            + mp.quad(lambda l: mp.exp(-l * t) * mu(l), [1, 2, 4, 8, 16, 40])
        print("   t = %d: int_0^oo e^{-lambda t} mu dlambda = %s   vs 1/log t = %s   (lambda < e^-10000 neglected: < 1e-4)" % (t, mp.nstr(v, 10), mp.nstr(1 / mp.log(t), 10)))
    print("\n3. mixture identity  -C(x) = int_0^oo W(1/(lambda+1/x)) mu(lambda) dlambda:")
    for x in (2, 10):
        f = lambda l: W(1 / (l + 1 / mp.mpf(x))) * mu(l)
        # lambda in (0,1) via lambda = e^{-a}; the part a > A is bounded by W(x)/A
        v = mp.quad(lambda a: f(mp.exp(-a)) * mp.exp(-a), [0, 1, 3, 6, 12, 25, 50, 100, 200, 400]) + mp.quad(f, [1, 2, 4, 8, 16, 30])
        tail = W(x) / 400
        print("   x = %d: integral (lambda > e^-400) = %s, + tail (~ W(x)/400 = %s) = %s ;  -C(x) = %s" % (x, mp.nstr(v, 10), mp.nstr(tail, 4), mp.nstr(v + tail, 10), mp.nstr(-C(x), 10)))
    mp.mp.dps = 30
    print("\n4. kernel identities:")
    p0 = 3 if q == 4 else 2
    for x in (1, 7.5, 40):
        x = mp.mpf(x)
        lim = int(80 * x) + 10
        prs = [p for p in PR if p <= lim and chi(p)]
        # Theta_1 piecewise linear: on [p_i, p_{i+1}], Theta_1(t) = Theta(p_i) t - Wsum(p_i)
        Fpp = lambda t: mp.exp(-t / x) * (1 / (x * x * mp.log(t)) + 2 / (x * t * mp.log(t) ** 2) + 1 / (t * t * mp.log(t) ** 2) + 2 / (t * t * mp.log(t) ** 3))
        tot = mp.mpf(0); tot2 = mp.mpf(0); Th = mp.mpf(0); Ws = mp.mpf(0)
        for i, p in enumerate(prs):
            Th += -chi(p) * LOGP[p]; Ws += -chi(p) * p * LOGP[p]
            hi = prs[i + 1] if i + 1 < len(prs) else lim + 50
            tot += mp.quad(lambda t: (Th * t - Ws) * Fpp(t), [p, hi])
            tot2 += mp.quad(lambda t: (Th * t - Ws) * mp.exp(-t / x), [p, hi])
        print("   x = %s: int Theta_1 F'' = %s ;  -C(x) = %s  |  x^-2 int Theta_1 e^{-t/x} = %s ;  W(x) = %s"
              % (mp.nstr(x, 4), mp.nstr(tot, 15), mp.nstr(-C(x), 15), mp.nstr(tot2 / x ** 2, 15), mp.nstr(W(x), 15)))
