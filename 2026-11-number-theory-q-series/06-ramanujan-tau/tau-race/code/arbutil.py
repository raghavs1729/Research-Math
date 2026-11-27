"""Helpers: outward-rounded decimal printing of Arb balls, primes."""
from flint import arb, fmpz
from fractions import Fraction

def _exact_fraction(a):
    """exact rational value of an exact arb (zero radius)"""
    m, e = a.man_exp()
    m = int(m); e = int(e)
    return Fraction(m * (1 << e)) if e >= 0 else Fraction(m, 1 << (-e))

def dec_down(x, k=15):
    """decimal string D with k fractional digits and D <= every point of the ball x"""
    lo = _exact_fraction(x.lower())
    n = (lo.numerator * 10**k) // lo.denominator        # floor
    return _fmt(n, k)

def dec_up(x, k=15):
    hi = _exact_fraction(x.upper())
    n = -((-hi.numerator * 10**k) // hi.denominator)     # ceil
    return _fmt(n, k)

def _fmt(n, k):
    s = "-" if n < 0 else ""
    n = abs(n)
    return "%s%d.%0*d" % (s, n // 10**k, k, n % 10**k)

def encl(x, k=15):
    return [dec_down(x, k), dec_up(x, k)]

def primes_upto(N):
    import numpy as np
    s = np.ones(N+1, dtype=bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]
