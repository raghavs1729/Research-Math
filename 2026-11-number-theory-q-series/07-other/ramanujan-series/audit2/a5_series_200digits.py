#!/usr/bin/env python3
"""
audit2 / a5: independent high-precision summation of all eleven series against the closed forms.
Direct summation of the series exactly as printed in Sun's papers (binomial coefficients via math.comb,
harmonic numbers by running sums), mpmath at 260 digits; the number of terms is chosen so that the
(geometric) tail is below 10^-215.  Prints the number of agreeing digits.
"""
import sys
from math import comb
from mpmath import mp, mpf, pi, sqrt, log, catalan, zeta, nstr

mp.dps = 260
K = (zeta(2, mpf(1)/3) - zeta(2, mpf(2)/3))/9          # L_{-3}(2)
G = catalan

class H:                                                # harmonic numbers H_m as mpf, on demand, increasing m
    def __init__(self): self.m, self.v = 0, mpf(0)
    def __call__(self, m):
        assert m >= self.m
        while self.m < m:
            self.m += 1; self.v += mpf(1)/self.m
        return self.v

def S1(N):  return sum(mpf((3*k+2)**2*16**k)/((2*k+1)*(6*k+1)*(6*k+5)*comb(6*k, 3*k)) for k in range(N))
def S2(N):  return sum(mpf((8*k-3)*comb(4*k, 2*k))/(k*(4*k-1)*9**k*comb(2*k, k)**2) for k in range(1, N))
def S3(N):  return sum(mpf((4*k-1)*3**k)/((2*k-1)*k**2*comb(2*k, k)*comb(3*k, k)) for k in range(1, N))
P4 = lambda k: 360*k**3+612*k**2+230*k+15
def w4(k):  return mpf(comb(3*k, k)*comb(6*k, 3*k)**2)/((3*k+1)*(3*k+2)*2**(15*k))
def S4(N):  return sum(P4(k)*w4(k) for k in range(N))
def S5(N):  return sum(mpf((77*k**3-86*k**2+29*k-3)*(-32)**k)/(k**3*(3*k-1)*(3*k-2)*comb(3*k, k)*comb(4*k, k)**2) for k in range(1, N))
def S6(N):
    h = H(); return sum(mpf((3*k+2)**2*16**k)/((2*k+1)*(6*k+1)*(6*k+5)*comb(6*k, 3*k))*h(k) for k in range(N))
def S7(N):  return sum(mpf(comb(2*k, k)*16**k*(1656*k**4-2604*k**3+1158*k**2-229*k+15))/(k**2*(2*k-1)**2*(6*k-1)*(6*k-5)*comb(6*k, 3*k)*comb(3*k, k)) for k in range(1, N))
def S8(N):
    h1, h2 = H(), H(); s = mpf(0)
    for k in range(N):
        a = h1(k); b = h2(2*k)
        s += w4(k)*(P4(k)*(b - a) - 180*k**2 + 36*k + 23)
    return s
def S9(N):
    h1, h2 = H(), H(); s = mpf(0)
    for k in range(1, N):
        a = h1(k-1); b = h2(2*k-1)
        s += mpf(comb(4*k, 2*k))*((8*k-3)*(5*b - 4*a) - 6)/(k*(4*k-1)*9**k*comb(2*k, k)**2)
    return s
def S10(N):
    h1, h3, h6 = H(), H(), H(); s = mpf(0)
    for k in range(N):
        a = h1(k); b = h3(3*k); c = h6(6*k)
        f = mpf(1296*k**3+1980*k**2+912*k+139)/((3*k+1)*(3*k+2))
        s += w4(k)*(P4(k)*(4*c - 3*b - a) + f)
    return s
def S11(N): return sum(mpf((1344*k**3+944*k**2+156*k+9)*comb(2*k, k)**2*comb(4*k, 2*k)**3)/((2*k+1)*2**(22*k)) for k in range(N))

JOBS = [  # name, function, number of terms, closed form, convergence rate
 ('T1 ', S1, 380,  pi/(2*sqrt(3)),                         '1/4'),
 ('T2 ', S2, 245,  sqrt(3)*pi/18,                          '1/9'),
 ('T3 ', S3, 245,  2*K,                                    '1/9'),
 ('T4 ', S4, 3100, 32*sqrt(2)/pi,                          '27/32'),
 ('T5 ', S5, 180,  -6*G,                                   '-27/512'),
 ('T6 ', S6, 385,  mpf(5)/2*K - pi*log(432)/(6*sqrt(3)),   '1/4'),
 ('T7 ', S7, 275,  mpf(0),                                 '4/27'),
 ('T8 ', S8, 3100, 48*sqrt(2)*log(2)/pi,                   '27/32'),
 ('T9 ', S9, 245,  mpf(3)/2*K,                             '1/9'),
 ('T10', S10, 3100, 192*sqrt(2)*log(2)/pi,                 '27/32'),
 ('T11', S11, 130, 64*sqrt(2)/pi**2,                       '1/64'),
]
ok = True
for name, f, N, val, rate in JOBS:
    s = f(N)
    err = abs(s - val)
    digits = int(-mp.log10(err)) if err != 0 else mp.dps
    ok &= digits >= 200
    print(f'{name} rate {rate:8s} terms {N:5d}  sum = {nstr(s, 30)}  closed form = {nstr(val, 30)}  agreement: {digits} digits')
print('ALL ELEVEN AGREE TO >= 200 DIGITS' if ok else 'SOME SERIES DISAGREES / BELOW 200 DIGITS')
sys.exit(0 if ok else 1)
