#!/usr/bin/env python3
"""
fujii.py -- Fujii's weights e^{-(p/x)^alpha}: the Platt-Trudgian criterion and the small-y obstruction
(mpmath, 30 digits; numerical, not certified).

 1. Criterion (Platt-Trudgian (3)):  Sigma(alpha) := sum_rho |Gamma(rho/alpha)|  <  (1/2) Gamma(1/(2 alpha)),
    rho = 1/2 + i gamma over the listed zeros (both signs).  In the normalisation of the note:
        bias constant   beta(alpha) = Gamma(1 + 1/(2 alpha))            (prime squares: sum_p log p e^{-(p^2/y)^alpha} ~ beta sqrt y)
        zero amplitude  S(alpha)    = (1/alpha) sum_rho |Gamma(rho/alpha)|
    and the criterion is S(alpha) < beta(alpha).  Root alpha_* by bisection.
 2. W_alpha(y) = sum_p s(p) log p e^{-(p/y)^alpha} for small y: the least alpha for which min_y W_alpha(y) <= 0
    (an unconditional obstruction to "W_alpha > 0 for all y"; it comes from Theta(5) = log 3 - log 5 < 0).
Usage: python3 fujii.py [q]
"""
import sys, os
import mpmath as mp
mp.mp.dps = 30
q = int(sys.argv[1]) if len(sys.argv) > 1 else 4
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, '..', '..', 'prototype', 'explore')
gam = [mp.mpf(l.split()[0]) for l in open(os.path.join(ZDIR, 'chi%d_zeros.txt' % q)) if l.strip() and l[0] != '#']


def Sigma(alpha, nz=None):
    alpha = mp.mpf(alpha)
    zs = gam if nz is None else gam[:nz]
    tot = mp.mpf(0)
    for g in zs:
        t = abs(mp.gamma(mp.mpc(0.5, g) / alpha))
        tot += 2 * t
        if t < mp.mpf(10) ** -40:
            break
    return tot


def crit(alpha):
    return Sigma(alpha) - mp.gamma(1 / (2 * mp.mpf(alpha))) / 2


def s_of(p):
    if q == 4:
        return 0 if p == 2 else (1 if p % 4 == 3 else -1)
    return 0 if p == 3 else (1 if p % 3 == 2 else -1)


PR = [p for p in range(2, 400) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def Walpha(y, alpha):
    return sum(s_of(p) * mp.log(p) * mp.exp(-(mp.mpf(p) / y) ** alpha) for p in PR if (mp.mpf(p) / y) ** alpha < 80)


def minW(alpha, lo=2.0, hi=40.0, n=1500):
    best = None
    for i in range(n + 1):
        y = mp.mpf(lo) * (mp.mpf(hi) / lo) ** (mp.mpf(i) / n)
        w = Walpha(y, alpha)
        # normalise by the first term so that tiny y does not dominate
        p0 = 3 if q == 4 else 2
        r = w / (mp.log(p0) * mp.exp(-(mp.mpf(p0) / y) ** alpha))
        if best is None or r < best[0]:
            best = (r, y, w)
    return best


if __name__ == '__main__':
    print("q = %d, %d zeros (largest ordinate %s)" % (q, len(gam), mp.nstr(gam[-1], 8)))
    print("\n1. criterion  S(alpha) = (1/alpha) sum_rho |Gamma(rho/alpha)|  <  beta(alpha) = Gamma(1 + 1/(2 alpha)):")
    print("   %8s %14s %14s %10s %12s" % ("alpha", "beta", "S", "S/beta", "beta - S"))
    for a in (0.5, 1, 2, 3, 4, 4.19, 6, 8, 12, 16, 20, 20.40442, 20.40443, 24, 32, 64, 128, 256):
        S = Sigma(a) / a
        be = mp.gamma(1 + 1 / (2 * mp.mpf(a)))
        print("   %8g %14s %14s %10s %12s" % (a, mp.nstr(be, 10), mp.nstr(S, 10), mp.nstr(S / be, 6), mp.nstr(be - S, 6)))
    if q == 4:
        print("   Platt-Trudgian (1000 zeros): sum < 20.1276643 and Gamma(1/(2 alpha))/2 > 20.1276649 at alpha = 20.40442;")
        print("   here (511 zeros): sum = %s, Gamma(1/(2 alpha))/2 = %s" % (mp.nstr(Sigma('20.40442'), 12), mp.nstr(mp.gamma(1 / (2 * mp.mpf('20.40442'))) / 2, 12)))
    lo, hi = mp.mpf(4), mp.mpf(200)
    for _ in range(60):
        mid = (lo + hi) / 2
        if crit(mid) < 0:
            lo = mid
        else:
            hi = mid
    print("   alpha_* (root of S = beta, listed zeros only) = %s" % mp.nstr(lo, 12))
    # contribution of the zeros beyond the list at alpha_*: |Gamma((1/2+ig)/a)| ~ sqrt(2 pi) (g/a)^{1/(2a)-1/2} e^{-pi g/(2a)}
    print("   size of one term at the top of the list at alpha_*: %s" % mp.nstr(abs(mp.gamma(mp.mpc(0.5, gam[-1]) / lo)), 3))
    print("\n2. small y: min over y in [2, 40] of W_alpha(y) / (log p0 e^{-(p0/y)^alpha}):")
    for a in (1, 2, 3, 4, 5, 6, 8, 12, 16):
        r, y, w = minW(a)
        print("   alpha = %4g: min ratio = %s at y = %s  (W_alpha = %s)" % (a, mp.nstr(r, 6), mp.nstr(y, 6), mp.nstr(w, 6)))
    lo, hi = mp.mpf(1), mp.mpf(16)
    for _ in range(30):
        mid = (lo + hi) / 2
        if minW(mid, n=600)[0] > 0:
            lo = mid
        else:
            hi = mid
    print("   least alpha with min_y W_alpha(y) <= 0 on this grid: alpha_W ~ %s" % mp.nstr(hi, 6))
