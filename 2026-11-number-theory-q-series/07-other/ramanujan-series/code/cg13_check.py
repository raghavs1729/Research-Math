"""Numerical check of the three-term identity (Au, arXiv:2312.14051, Example I) for generic a and of its limit a -> 1/2.
L1(a) + L2(a) = Rhs(a), with
 L1(a)  = sum_k (6a+4k+1)(a+1/2)_k^5 (4a+1)_k / ((1+2k-2a)(4a+k)(1/2-a)_k (2a+1)_k^5)
 L2(a)  = -sum_n (-1)^n (820(a+n)^2+180(a+n)+13) Gamma(1/2-a) Gamma(2a+1)^5 (a+1/2)_n^10 / (4096^a Gamma(a+1/2)^5 Gamma(4a+1) (2a+1)_{2n}^5)
 Rhs(a) = -sum_n P(a+n) (a+1/2)_n^16 (4a+1)_{4n} / (65536 a^5 (2a+1)_{2n}^10)."""
import mpmath as mp
mp.mp.dps = 50
rf, G = mp.rf, mp.gamma
P = lambda x: 43680*x**4 + 20632*x**3 + 4340*x*x + 466*x + 21
def L1(a): return mp.nsum(lambda k: (6*a+4*k+1)*rf(a+0.5, k)**5*rf(4*a+1, k)/((1+2*k-2*a)*(4*a+k)*rf(0.5-a, k)*rf(2*a+1, k)**5), [0, mp.inf])
def L2(a): return -mp.nsum(lambda n: (-1)**int(n)*(820*(a+n)**2+180*(a+n)+13)*G(0.5-a)*G(2*a+1)**5*rf(a+0.5, n)**10/(mp.mpf(4096)**a*G(a+0.5)**5*G(4*a+1)*rf(2*a+1, 2*n)**5), [0, mp.inf])
def Rhs(a): return -mp.nsum(lambda n: P(a+n)*rf(a+0.5, n)**16*rf(4*a+1, 4*n)/(65536*a**5*rf(2*a+1, 2*n)**10), [0, mp.inf])
for a in (mp.mpf('0.3'), mp.mpf('0.13'), mp.mpf('0.45'), mp.mpf('0.62')):
    print('a =', a, ' L1+L2-Rhs =', mp.nstr(L1(a) + L2(a) - Rhs(a), 8), ' Rhs =', mp.nstr(Rhs(a), 12))
half = mp.mpf(1)/2
print('Rhs(1/2) =', Rhs(half), '  -pi^4/30 =', -mp.pi**4/30)
for d in (mp.mpf(10)**-6, mp.mpf(10)**-9):
    print('delta', mp.nstr(d, 2), ' L1+L2 at 1/2+-delta:', mp.nstr(L1(half+d) + L2(half+d), 20), mp.nstr(L1(half-d) + L2(half-d), 20))
