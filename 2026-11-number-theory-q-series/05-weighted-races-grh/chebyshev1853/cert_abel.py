#!/usr/bin/env python3
"""
cert_abel.py -- certified (Arb ball arithmetic, python-flint) inputs for CHEBYSHEV1853.md.
No floating-point decision in this file: every comparison is a comparison of balls.

Usage: python3 cert_abel.py [q ...]      (default 4 3)

 A. Finite prime data, all race primes p <= T = 2^22 (256-bit balls):
      A0. the endpoint u = T itself (T is not prime; the last prime is 4 194 301): Theta_1(T) > 0 and
          (T - vo(T))/sqrt(T) <= c1, so that A1 and A2 hold on the whole interval (0, T]     (added after the referee report)
      A1. Theta_1(p) > 0 for every race prime p0 < p <= T            (recomputed; = cert_small.py part B)
      A2. y - c1 sqrt(y) <= vo(y) < y on (0, T]                      (recomputed; = cert_small.py part C)
      A3. c3 := sup_{0<u<=T} |Theta(u)|/sqrt(u), certified enclosure (new; used for odd prime powers k >= 3)
      A4. unconditional positivity of W(y) for 0 < y <= T/60:
              (c_theta/2) y (T^2 + 2 T y + 2 y^2) e^{-(T-a)/y} < A   at y = T/60
          where A e^{-a/y} is a lower bound for int_{p0}^{T} Theta_1(t) e^{-t/y} dt  (see the note, Lemma 7).
 B. Zeros of L(s,chi) (GRH is assumed in the statements that use this part):
      B1. each listed ordinate gamma_j: Hardy Z changes sign on (gamma_j - 1e-20, gamma_j + 1e-20)  => a zero there
      B2. d := sigma_0^2 - sum_j 2/(1/4 + gamma_j^2) > 0 (enclosure), sigma_0^2 = log(q/pi) - gamma + 2 L'/L(1,chi)
          Under GRH sum_rho 1/|rho|^2 = sigma_0^2, so the zeros NOT accounted for by B1 (unlisted zeros and extra
          multiplicities) satisfy sum 1/|rho|^2 = d, come in conjugate pairs, hence each has |gamma| >= g := sqrt(2/d - 1/4).
      B3. S_Gamma := sum_rho |Gamma(rho)| = sum_rho sqrt(pi / cosh(pi gamma))  <=
              2 sum_j sqrt(pi/cosh(pi gamma_j^-)) + d * (1/4 + g^2) sqrt(pi/cosh(pi g))
          (the function (1/4+t^2) sqrt(pi/cosh(pi t)) is decreasing for t >= 2).
          Also the same with Gamma replaced by Gamma(rho/alpha)/alpha-type weights is NOT done here (see fujii.py).
 C. Closed-form constants as balls: sqrt(pi)/2, Gamma(1/4)/4, Gamma(1+1/k), b = L'/L(0,chi), b0.
"""
import sys, os, time
from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, '..', '..', 'prototype', 'explore')
T = 1 << 22
C1 = {4: '2.071193', 3: '2.081942'}


def sieve(n):
    bs = bytearray([1]) * (n + 1)
    bs[0] = bs[1] = 0
    i = 2
    while i * i <= n:
        if bs[i]:
            bs[i * i::i] = bytearray(len(bs[i * i::i]))
        i += 1
    return [i for i in range(2, n + 1) if bs[i]]


def s_of(p, q):
    if q == 4:
        return 0 if p == 2 else (1 if p % 4 == 3 else -1)
    return 0 if p == 3 else (1 if p % 3 == 2 else -1)


def LL1(q):
    """L'/L(1,chi) closed forms (Lerch).  chi_4: gamma + 2 log 2 + 3 log pi - 4 log Gamma(1/4);
       chi_3: gamma + log(2 pi) - 3 log(Gamma(1/3)/Gamma(2/3))."""
    g = arb.const_euler()
    pi = arb.pi()
    if q == 4:
        return g + 2 * arb(2).log() + 3 * pi.log() - 4 * (arb(1) / 4).lgamma()
    return g + (2 * pi).log() - 3 * ((arb(1) / 3).lgamma() - (arb(2) / 3).lgamma())


def partA(q, primes):
    ctx.prec = 256
    t0 = time.time()
    ell = 2 if q == 4 else 3
    p0 = 3 if q == 4 else 2
    c1 = arb(C1[q])
    Th = arb(0)      # Theta(p) = sum s log
    Wt = arb(0)      # sum s p log p
    vo = arb(0)      # theta of the principal character
    ok1 = ok2 = True
    c3 = arb(0); c3_at = 0
    minT1 = None; minT1_at = 0
    supc1 = arb(-10); supc1_at = 0
    th1_small = {}
    for p in primes:
        lp = arb(p).log()
        sp = arb(p).sqrt()
        # A2, left limit at p (for every prime, also p = ell): (p - vo(p^-))/sqrt p <= c1
        r = (arb(p) - vo) / sp
        if not (r <= c1):
            ok2 = False
        if r.mid() > supc1.mid():
            supc1 = r; supc1_at = p
        s = s_of(p, q)
        if s != 0:
            vo += lp
            Th += s * lp
            Wt += s * p * lp
            th1 = p * Th - Wt                  # Theta_1(p)
            if p <= 13:
                th1_small[p] = th1
            if p > p0:
                if not (th1 > 0):
                    ok1 = False
                if minT1 is None or th1.mid() < minT1.mid():
                    minT1 = th1; minT1_at = p
            # A3: sup over [p, p') of |Theta(u)|/sqrt u is |Theta(p)|/sqrt p
            r3 = abs(Th) / sp
            if r3.mid() > c3.mid():
                c3 = r3; c3_at = p
        # A2, value at p: vo(p) < p
        if not (vo < p):
            ok2 = False
    # Endpoint u = T = 2^22 (referee item REP-1): the loop above tests only at primes, and the last prime is
    # 4 194 301 < T.  On (p_last, T]:  Theta_1 is linear with slope Theta(p_last), so positivity at p_last and at T
    # gives positivity on the segment;  (u - vo(u))/sqrt(u) is increasing (vo constant), so its sup is at u = T;
    # vo(u) < u holds since vo is constant and vo(p_last) < p_last;  |Theta(u)|/sqrt(u) is decreasing (A3 needs nothing).
    p_last = primes[-1]
    th1_T = T * Th - Wt                        # Theta_1(T) = sum_{p<=T} s(p) log p (T - p)
    rT = (arb(T) - vo) / arb(T).sqrt()
    okT1 = bool(th1_T > 0)
    okT2 = bool(rT <= c1) and bool(vo < T)
    ok1 = ok1 and okT1
    ok2 = ok2 and okT2
    print("A0. q=%d: last prime <= 2^22 is %d;  endpoint T = 2^22:  Theta_1(T) = %s > 0: %s ;  (T - vo(T))/sqrt(T) = %s <= c1: %s"
          % (q, p_last, th1_T.str(15), okT1, rT.str(15), okT2))
    print("A1. q=%d: Theta_1(p) > 0 certified at all race primes p0 < p <= 2^22: %s ;  min = %s at p = %d"
          % (q, ok1, minT1.str(15), minT1_at))
    print("A2. q=%d: y - c1 sqrt(y) <= vo(y) < y on (0, 2^22] with c1 = %s: %s ;  sup = %s as y -> %d^-"
          % (q, C1[q], ok2, supc1.str(15), supc1_at))
    print("A3. q=%d: sup_{u<=2^22} |Theta(u)|/sqrt(u) = %s attained at u = %d" % (q, c3.str(15), c3_at))
    # rigorous: verify |Theta(p)|/sqrt p <= c3bound at every prime in a second pass
    c3bound = arb(c3.upper()) + arb('1e-12')
    Th = arb(0); ok3 = True
    for p in primes:
        s = s_of(p, q)
        if s:
            Th += s * arb(p).log()
            if not (abs(Th) / arb(p).sqrt() < c3bound):
                ok3 = False
    print("      certified: |Theta(u)| < %s sqrt(u) for 0 < u <= 2^22: %s" % (c3bound.str(12), ok3))
    # A4
    eps_a = arb('0.2') / (arb(T).log() ** 2)
    cth = 1 + eps_a
    if q == 4:
        # Theta_1 is linear on [5,7], Theta_1(5) = 2 log 3, Theta_1(7) = 4 log 3 - 2 log 5
        lo = min(th1_small[5].lower(), th1_small[7].lower()); a = 7
    else:
        # Theta_1 is linear on [2,5]: Theta_1(t) = (t-2) log 2; on [3,5] it is >= log 2
        lo = arb(2).log().lower(); a = 5
    A = 2 * arb(lo)
    y = arb(T) / 60
    lhs = cth / 2 * y * (arb(T) ** 2 + 2 * T * y + 2 * y * y) * (-(T - a) / y).exp()
    ok4 = bool(lhs < A)
    print("A4. q=%d: eps_a = 0.2/log^2(2^22) = %s, c_theta = %s;  A = %s, a = %d;  tail/A at y = T/60 = %s : %s"
          % (q, eps_a.str(8), cth.str(10), A.str(10), a, (lhs / A).str(5), ok4))
    print("      => W(y) > 0 and C(x) < 0 unconditionally for 0 < x, y <= T/60 = %.2f   [%.1f s]" % (T / 60, time.time() - t0))
    return ok1 and ok2 and ok3 and ok4, c3bound


def partB(q):
    ctx.prec = 160
    t0 = time.time()
    fn = os.path.join(ZDIR, 'chi%d_zeros.txt' % q)
    gam = [ln.split()[0] for ln in open(fn) if ln.strip() and not ln.startswith('#')]
    ch = dirichlet_char(q, q - 1)
    eps = arb('1e-20')
    ok = True
    prev_hi = None
    pi = arb.pi()
    sum_inv = arb(0)       # lower bound uses the upper end of gamma
    SG = arb(0)            # upper bound uses the lower end of gamma
    SG_first = None
    for gs in gam:
        g = arb(gs)
        lo, hi = g - eps, g + eps
        zl, zh = ch.hardy_z(acb(lo)).real, ch.hardy_z(acb(hi)).real      # Z is real on the real axis
        if not ((zl > 0 and zh < 0) or (zl < 0 and zh > 0)):
            ok = False
        if prev_hi is not None and not (lo > prev_hi):
            ok = False
        prev_hi = hi
        gball = arb(g.mid(), '1.01e-20')            # ball containing the zero
        sum_inv += 2 / (arb(1) / 4 + gball * gball)
        term = (pi / (pi * gball).cosh()).sqrt()
        SG += 2 * term
        if SG_first is None:
            SG_first = 2 * term
    sig0 = (arb(q) / pi).log() - arb.const_euler() + 2 * LL1(q)
    d = sig0 - sum_inv
    okd = bool(d > 0)
    g0 = (2 / d - arb(1) / 4).sqrt()
    glow = arb(g0.lower())
    rest = d * (arb(1) / 4 + glow * glow) * (pi / (pi * glow).cosh()).sqrt()
    okg = bool(glow > 2)
    SGtot = SG + arb(rest.upper())
    print("B1. q=%d: Hardy Z sign change certified at all %d listed ordinates (+-1e-20), intervals disjoint: %s" % (q, len(gam), ok))
    print("B2. q=%d: sigma_0^2 = %s ;  sum_j 2/(1/4+gamma_j^2) = %s ;  d = %s > 0: %s ;  g = %s"
          % (q, sig0.str(22), sum_inv.str(22), d.str(8), okd, g0.str(8)))
    print("B3. q=%d: 2 sum_j |Gamma(1/2+i gamma_j)| = %s (first zero alone: %s);  remainder <= %s"
          % (q, SG.str(20), SG_first.str(20), arb(rest.upper()).str(5)))
    print("      S_Gamma = sum_rho |Gamma(rho)| <= %s  (under GRH)   [%.1f s]" % (SGtot.str(20), time.time() - t0))
    return ok and okd and okg, SGtot, sig0


def partC(q, SG):
    ctx.prec = 128
    pi = arb.pi()
    g = arb.const_euler()
    ll1 = LL1(q)
    b = -(arb(q) / pi).log() + g + arb(2).log() - ll1
    print("C.  q=%d: sqrt(pi)/2 = %s ;  Gamma(1/4)/4 = %s ;  sqrt(pi)/2 - S_Gamma >= %s"
          % (q, (pi.sqrt() / 2).str(20), ((arb(1) / 4).gamma() / 4).str(20), arb((pi.sqrt() / 2 - SG).lower()).str(12)))
    print("      L'/L(1,chi) = %s ;  b = L'/L(0,chi) = %s" % (ll1.str(20), b.str(20)))
    print("      Gamma(1+1/k), k=3..7: " + ", ".join((1 + arb(1) / k).gamma().str(10) for k in range(3, 8)))
    print("      Gamma(1+1/(2k)), k=2..7: " + ", ".join((1 + arb(1) / (2 * k)).gamma().str(10) for k in range(2, 8)))


if __name__ == '__main__':
    qs = [int(a) for a in sys.argv[1:]] or [4, 3]
    primes = sieve(T)
    allok = True
    for q in qs:
        okA, c3b = partA(q, primes)
        okB, SG, sig0 = partB(q)
        partC(q, SG)
        allok = allok and okA and okB
    print("ALL CERTIFIED:", allok)
