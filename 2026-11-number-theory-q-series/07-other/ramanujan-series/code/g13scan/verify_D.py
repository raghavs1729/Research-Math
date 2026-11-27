# target: the diagonal series of the G13 triangle pair, d_n ~ (3/4)^n (1/6)_n(5/6)_n/((1/4)_n(3/4)_n) * rational
from mpmath import mp, mpf, binomial, pi, zeta, log, sqrt, catalan, nsum, inf, rf
mp.dps = 60
def dterm(n):
    # d_1 = -74/315 and ratio d_{n+1}/d_n from g13_diag.py; closed form up to a constant:
    n = mpf(n)
    return rf(mpf(5)/6, n)*rf(mpf(7)/6, n)*rf(mpf(1), n)*rf(mpf(1)/2, n)/(rf(mpf(5)/4, n)*rf(mpf(7)/4, n)*rf(mpf(2), n)*rf(mpf(3)/2, n))*(12*n*n+18*n+7)*(mpf(3)/4)**n*n/(n+1)*0+ \
           rf(mpf(5)/6, n)*rf(mpf(7)/6, n)*rf(mpf(1), n)*rf(mpf(1)/2, n)/(rf(mpf(5)/4, n)*rf(mpf(7)/4, n)*rf(mpf(2), n)*rf(mpf(3)/2, n))*(12*n*n+18*n+7)*(mpf(3)/4)**n
E = [dict(id='D13', start=1, summand=dterm, rhs=lambda: mpf(0), src='diagonal series of the G13 pair')]
ENTRIES = E
