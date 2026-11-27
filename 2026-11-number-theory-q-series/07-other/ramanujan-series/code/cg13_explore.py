"""Cohen-Guillera #13: Cullen's 1/pi^4 series shifted by 1/2, conjectured value (2048/15)/pi^4.
Explore Au's WZ pair (Dougall7F6(1/2-n,...;k+1+2n)) near n = 1/2."""
import mpmath as mp
from wz import *
from numerics import *
mp.mp.dps = 60
F = seed_instance('Dougall7F6', dict(a=-1,b=-1,c=-1,d=-1,e=-1), 2, consts=dict(a='1/2',b='1/2',c='1/2',d='1/2',e='1/2'), k0=0)
print('F =', F)
R = wz_mate(F); assert check_wz(F, R)
W = WZNum(F, R)
H = mp.mpf(1)/2
ss = [H]*7 + [mp.mpf(1)/4, mp.mpf(3)/4]
P = lambda n: 43680*n**4 + 20632*n**3 + 4340*n*n + 466*n + 21
def tshift(n, x=H):   # Cohen-Guillera summand of S(H,a,P;x), a = 4096
    v = P(n + x)/mp.mpf(4096)**(n + x)
    for s in ss: v *= mp.gamma(s + n + x)/mp.gamma(s)
    return v/mp.gamma(1 + n + x)**9
Sx = mp.nsum(lambda n: tshift(n), [0, mp.inf])
print('S(1/2) =', Sx, ' (2048/15)/pi^4 =', mp.mpf(2048)/15/mp.pi**4)
for eps in (mp.mpf(10)**-8, mp.mpf(10)**-10):
    a = H + eps
    print('eps', mp.nstr(eps, 3))
    print('  G(n+a,0)/tshift(n):', [mp.nstr(W.Gv(n + a, 0)/tshift(n), 12) for n in range(4)])
    print('  F(a,k) k=0..4:', [mp.nstr(W.Fv(a, k), 10) for k in range(5)])
    print('  g(n+a) ~ G(n+a,K=3000):', [mp.nstr(W.Gv(n + a, 3000), 10) for n in range(3)])
