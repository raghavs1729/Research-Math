"""Cohen-Guillera #12: Zhao's 1/pi^4 series shifted by 1/2, observed value 9216 zeta(3)/pi^7."""
import mpmath as mp
from wz import *
from numerics import *
mp.mp.dps = 50
F = seed_instance('Dougall7F6', dict(a=1, b=0, c=0, d=0, e=0), 1, consts=dict(a='1/2', b='1/2', c='1/2', d='1/2', e='1/2'), k0=0)
print('F =', str(F)[:900])
R = wz_mate(F); assert check_wz(F, R)
W = WZNum(F, R)
print('R(n,0) =', pari.factor(pari.numerator(pari.subst(R, 'k', 0))), '/', pari.factor(pari.denominator(pari.subst(R, 'k', 0))))
H = mp.mpf(1)/2
ss = [H]*5 + [mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4]
P = lambda n: 4528*n**4 + 3180*n**3 + 972*n*n + 147*n + 9
a0 = -mp.mpf(256)/27
def tshift(n, x=H):
    v = P(n + x)*(-1)**int(n)/abs(a0)**(n + x)
    for s in ss: v *= mp.gamma(s + n + x)/mp.gamma(s)
    return v/mp.gamma(1 + n + x)**9
Sx = mp.nsum(lambda n: tshift(n), [0, mp.inf])
print('S(1/2) =', mp.nstr(Sx, 30), ' 9216 zeta(3)/pi^7 =', mp.nstr(9216*mp.zeta(3)/mp.pi**7, 30))
for eps in (mp.mpf(10)**-8, mp.mpf(10)**-10):
    a = H + eps
    print('eps', mp.nstr(eps, 3))
    print('  G(n+a,0)/tshift(n):', [mp.nstr(W.Gv(n + a, 0)/tshift(n), 12) for n in range(4)])
    print('  F(a,k) k=0..4:', [mp.nstr(W.Fv(a, k), 10) for k in range(5)])
    print('  G(n+a,K=3000):', [mp.nstr(W.Gv(n + a, 3000), 10) for n in range(3)])
