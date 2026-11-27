"""Parameter derivatives of the Theorem-1 WZ pair: F(n,k;a,b) = seed(1/2 + a + n, k + b + n) (real version)."""
import mpmath as mp
from wz import *
from seeds_extra import load_ebisu
from numerics import pari_to_func
mp.mp.dps = 40
EB = load_ebisu()
name = [n for n in EB if n.startswith('Eb024')][0]
params, T = EB[name]
F = T.subs({'a': '1/2+a+n', 'k': 'k+b+n'})
import time; t0 = time.time()
R = wz_mate(F); assert check_wz(F, R); print('R(n,k,a,b) found in', round(time.time()-t0,1), 's')
print('R =', R)
Rf = pari_to_func(R, ['n','k','a','b'])
def Fv(n, k, a, b):
    # real version: replace (-1/3)^(k+b+n) by (-1)^(k+n) 3^-(k+b+n)
    vals = {'n': mp.mpf(n), 'k': mp.mpf(k), 'a': mp.mpf(a), 'b': mp.mpf(b)}
    r = mp.mpf(1)
    for base, L in F.bases:
        e = sum((mp.mpf(x.numerator)/x.denominator)*(1 if v == '1' else vals[v]) for v, x in L.items())
        if base < 0:
            e_int = int(n) + int(k)
            r *= (-1)**e_int * mp.power(-mp.mpf(base.numerator)/base.denominator, e)
        else:
            r *= mp.power(mp.mpf(base.numerator)/base.denominator, e)
    for L, p in F.gammas:
        xx = sum((mp.mpf(c.numerator)/c.denominator)*(1 if v == '1' else vals[v]) for v, c in L.items())
        r *= mp.gamma(xx)**p if p > 0 else mp.rgamma(xx)**(-p)
    return r
def SG(a, b): return mp.fsum(Rf(mp.mpf(n), mp.mpf(0), a, b)*Fv(n, 0, a, b) for n in range(0, 60))
def SF(a, b): return mp.fsum(Fv(0, k, a, b) for k in range(0, 90))
for (a, b) in ((0, 0), (0.05, 0), (0, 0.07), (0.03, -0.04)):
    print((a, b), mp.nstr(SG(a, b), 25), mp.nstr(SF(a, b), 25))
dGa = mp.diff(lambda a: SG(a, 0), 0); dGb = mp.diff(lambda b: SG(0, b), 0)
dFa = mp.diff(lambda a: SF(a, 0), 0); dFb = mp.diff(lambda b: SF(0, b), 0)
print('d/da:', mp.nstr(dGa, 20), mp.nstr(dFa, 20)); print('d/db:', mp.nstr(dGb, 20), mp.nstr(dFb, 20))
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
c = mp.sqrt(6)/mp.pi
print('pslq d/da SF:', mp.pslq([dFa/c, K, mp.pi*mp.log(2)/mp.sqrt(3), mp.pi*mp.log(3)/mp.sqrt(3), mp.pi/mp.sqrt(3)], maxcoeff=2000, maxsteps=20000))
print('pslq d/db SF:', mp.pslq([dFb/c, K, mp.pi*mp.log(2)/mp.sqrt(3), mp.pi*mp.log(3)/mp.sqrt(3), mp.pi/mp.sqrt(3)], maxcoeff=2000, maxsteps=20000))
