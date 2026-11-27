#!/usr/bin/env python3
"""
constants.py -- closed-form constants for EXPLICIT.md and independent checks.

For chi in {chi_4, chi_3} (both real, primitive, odd):
  LL1  = L'/L(1,chi)
  sig0 = sum_rho 1/|rho|^2 (GRH) = log(q/pi) - gamma + 2 L'/L(1,chi)      [upper bound for S]
  b    = L'/L(0,chi)       = -log(q/pi) + gamma + log 2 - L'/L(1,chi)
  b0   = constant term of L'/L at s=-1 = -log(q/pi) + gamma/2 - psi(3/2)/2 - L'/L(2,chi)
  C    = 1 + b0  (constant term of the explicit formula for psi_1)
Checks:
  (a) b, b0 against direct numerical differentiation of L(s,chi) (mpmath Hurwitz zeta);
  (b) sig0 against the certified zero lists + Riemann-von Mangoldt tail;
  (c) sum_rho 1/(rho(rho+1)) = b0 - b + log 2 against the zero lists;
  (d) the explicit formula for psi_1(x,chi) against a direct prime-power sum, at small x
      (x = 1, 3, 5.5, 10, 20, 50, 200), where the truncation of the zero sum is small.
"""
import mpmath as mp, os, sys
mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = os.path.join(HERE, '..', '..', 'prototype', 'explore')

def chi_vec(q):
    return [0, 1, 0, -1] if q == 4 else [0, 1, -1]

def L(s, q):
    return mp.dirichlet(s, chi_vec(q))

def LL(s, q):
    # mpmath evaluates the derivative of the Dirichlet L-function analytically (Hurwitz zeta)
    return mp.dirichlet(s, chi_vec(q), derivative=1) / L(s, q)

def LL1_closed(q):
    """Independent closed forms (Lerch/Kronecker):  L'/L(1,chi_4) = gamma + 2 log 2 + 3 log pi - 4 log Gamma(1/4),
       L'/L(1,chi_3) = gamma + log(2 pi) - 6 log Gamma(1/3) + ... written via log(Gamma(1/3)/Gamma(2/3))."""
    g = mp.euler
    if q == 4:
        return g + 2 * mp.log(2) + 3 * mp.log(mp.pi) - 4 * mp.loggamma(mp.mpf(1) / 4)
    # L'(1,chi)/L(1,chi) = gamma + log(2 pi) - 3 log( Gamma(1/3)/Gamma(2/3) )   for chi_3
    return g + mp.log(2 * mp.pi) - 3 * mp.log(mp.gamma(mp.mpf(1) / 3) / mp.gamma(mp.mpf(2) / 3))

def consts(q):
    g = mp.euler
    LL1 = LL(1, q)
    sig0 = mp.log(q / mp.pi) - g + 2 * LL1
    b = -mp.log(q / mp.pi) + g + mp.log(2) - LL1
    b0 = -mp.log(q / mp.pi) + g / 2 - mp.digamma(mp.mpf(3) / 2) / 2 - LL(2, q)
    return dict(LL1=LL1, sig0=sig0, b=b, b0=b0, C=1 + b0)

def zeros(q):
    fn = os.path.join(ZDIR, 'chi%d_zeros.txt' % q)
    zs = []
    for line in open(fn):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        zs.append(mp.mpf(line.split()[-1] if len(line.split()) > 1 and False else line.split()[0]))
    return zs

def chi(n, q):
    return chi_vec(q)[n % q]

def psi1_direct(x, q):
    """psi_1(x,chi) = sum_{n<=x} chi(n) Lambda(n) (x-n), exact prime-power sum"""
    N = int(x)
    sieve = bytearray([1]) * (N + 1)
    tot = mp.mpf(0)
    for p in range(2, N + 1):
        if sieve[p]:
            for m in range(p * p, N + 1, p):
                sieve[m] = 0
            lp = mp.log(p)
            pk = p
            while pk <= N:
                tot += chi(pk, q) * lp * (x - pk)
                pk *= p
    return tot

if __name__ == '__main__':
    for q in (4, 3):
        c = consts(q)
        print('=== q = %d ===' % q)
        for k in ('LL1', 'sig0', 'b', 'b0', 'C'):
            print('  %-5s = %s' % (k, mp.nstr(c[k], 25)))
        print('  check: L\'/L(1) closed form (Lerch) = %s   diff %s' % (mp.nstr(LL1_closed(q), 25), mp.nstr(LL1_closed(q) - c['LL1'], 5)))
        # (a) direct checks of b and b0
        b_dir = LL(0, q)
        eps = mp.mpf(10) ** (-12)
        b0_dir = (LL(-1 + eps, q) - 1 / eps + LL(-1 - eps, q) + 1 / eps) / 2
        print('  check (a): L\'/L(0) direct = %s   diff %s' % (mp.nstr(b_dir, 20), mp.nstr(b_dir - c['b'], 5)))
        print('  check (a): b0 direct (symmetric eps=1e-12) = %s   diff %s' % (mp.nstr(b0_dir, 20), mp.nstr(b0_dir - c['b0'], 5)))
        zs = zeros(q)
        T = zs[-1]
        # (b) sum over zeros of 2/(1/4+g^2), tail by dN = (1/2pi) log(q t/2pi) dt
        part = sum(2 / (mp.mpf(1) / 4 + g * g) for g in zs)
        tail = mp.quad(lambda t: 2 / (mp.mpf(1) / 4 + t * t) * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
        print('  check (b): %d zeros up to %s: partial %s + RvM tail %s = %s  vs closed form %s'
              % (len(zs), mp.nstr(T, 10), mp.nstr(part, 12), mp.nstr(tail, 6), mp.nstr(part + tail, 12), mp.nstr(c['sig0'], 12)))
        # true S = sum_rho 1/|rho(rho+1)|
        Spart = sum(2 / mp.sqrt((mp.mpf(1) / 4 + g * g) * (mp.mpf(9) / 4 + g * g)) for g in zs)
        Stail = mp.quad(lambda t: 2 / mp.sqrt((mp.mpf(1) / 4 + t * t) * (mp.mpf(9) / 4 + t * t)) * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
        print('  S = sum_rho 1/|rho(rho+1)| ~ %s (zeros %s + tail %s)   [proof uses S <= sig0]' % (mp.nstr(Spart + Stail, 10), mp.nstr(Spart, 10), mp.nstr(Stail, 6)))
        # (c)
        re_part = sum(2 * (1 / ((mp.mpf(1) / 2 + 1j * g) * (mp.mpf(3) / 2 + 1j * g))).real for g in zs)
        re_tail = mp.quad(lambda t: 2 * (1 / ((mp.mpf(1) / 2 + 1j * t) * (mp.mpf(3) / 2 + 1j * t))).real * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
        print('  check (c): sum_rho 1/(rho(rho+1)) from zeros = %s ; closed form b0 - b + log2 = %s'
              % (mp.nstr(re_part + re_tail, 10), mp.nstr(c['b0'] - c['b'] + mp.log(2), 10)))
        # (d) explicit formula for psi_1
        # small x only: there the bound for the omitted zeros is small enough for the test to mean something
        # (at x = 10^3, 10^4 the omitted-zero bound is 10^2..10^3 and the comparison is vacuous)
        for x in ('1', '3', '5.5', '10', '20', '50', '200'):
            x = mp.mpf(x)
            direct = psi1_direct(x, q)
            zsum = sum(2 * (x ** (mp.mpf(3) / 2 + 1j * g) / ((mp.mpf(1) / 2 + 1j * g) * (mp.mpf(3) / 2 + 1j * g))).real for g in zs)
            triv = mp.nsum(lambda n: x ** (-2 * n) / (2 * n * (2 * n + 1)), [1, mp.inf])
            formula = -zsum - c['b'] * x + mp.log(x) + c['C'] - triv
            # absolute tail bound for the omitted zeros
            tb = x ** 1.5 * mp.quad(lambda t: 2 / (t * t) * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, mp.inf])
            print('  check (d): x=%-6s psi_1 direct = %-22s formula(%d zeros) = %-22s diff = %-12s (omitted-zero bound ~ %s)'
                  % (mp.nstr(x, 6), mp.nstr(direct, 15), len(zs), mp.nstr(formula, 15), mp.nstr(direct - formula, 5), mp.nstr(tb, 4)))
