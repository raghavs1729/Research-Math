"""Parameter derivatives of the Theorem-2 pair: F(n,k;a,b,c) = Seed4(a-n, -1/2+b-2n, k+1/2+c+n)."""
import mpmath as mp, time
from wz import *
from numerics import pari_to_func
mp.mp.dps = 40
params, T = SEEDS['Seed4']
F = T.subs({'a': 'a-n', 'b': '-1/2+b-2*n', 'k': 'k+1/2+c+n'})
print('F =', F)
t0 = time.time(); R = wz_mate(F); assert check_wz(F, R); print('R found', round(time.time()-t0, 1), 's')
print('R(n,0) =', pari.factor(pari.subst(R, 'k', 0)))
Rf = pari_to_func(R, ['n', 'k', 'a', 'b', 'c'])
def Fv(n, k, a, b, c):
    return F.evalf({'n': mp.mpf(n), 'k': mp.mpf(k), 'a': mp.mpf(a), 'b': mp.mpf(b), 'c': mp.mpf(c)}, mp)
def SG(a, b, c): return mp.fsum(Rf(mp.mpf(n), mp.mpf(0), mp.mpf(a), mp.mpf(b), mp.mpf(c))*Fv(n, 0, a, b, c) for n in range(0, 90))
def SF(a, b, c): return mp.nsum(lambda k: Fv(0, k, a, b, c), [0, mp.inf])
print('SG(0)=', mp.nstr(SG(0,0,0), 20), 'SF(0)=', mp.nstr(SF(0,0,0), 20))
for pt in ((0.03, 0, 0), (0, 0.04, 0), (0, 0, 0.05), (0.02, -0.03, 0.01)):
    print(pt, 'SG', mp.nstr(SG(*pt), 18), 'SF', mp.nstr(SF(*pt), 18))
H = mp.harmonic
K = (mp.zeta(2, mp.mpf(1)/3) - mp.zeta(2, mp.mpf(2)/3))/9
e = mp.mpf('1e-10')
d = {'a': (SG(e,0,0)-SG(-e,0,0))/(2*e), 'b': (SG(0,e,0)-SG(0,-e,0))/(2*e), 'c': (SG(0,0,e)-SG(0,0,-e))/(2*e)}
print({k_: mp.nstr(v, 15) for k_, v in d.items()})
c0 = -2*mp.sqrt(3)/mp.pi   # G(n,0) = c0 t(n+1)
B35 = mp.fsum(mp.binomial(4*k,2*k)*((8*k-3)*(5*H(2*k-1)-4*H(k-1)) - 6)/(k*(4*k-1)*mp.mpf(9)**k*mp.binomial(2*k,k)**2) for k in range(1, 120))
print('B35 check:', mp.nstr(B35 - 3*K/2, 5))
base = SG(0,0,0)
for vec in ([d['a']/c0, K, mp.pi/mp.sqrt(3), mp.pi*mp.log(3)/mp.sqrt(3), mp.pi*mp.log(2)/mp.sqrt(3)], [d['b']/c0, K, mp.pi/mp.sqrt(3), mp.pi*mp.log(3)/mp.sqrt(3), mp.pi*mp.log(2)/mp.sqrt(3)], [d['c']/c0, K, mp.pi/mp.sqrt(3), mp.pi*mp.log(3)/mp.sqrt(3), mp.pi*mp.log(2)/mp.sqrt(3)]):
    print('pslq:', mp.pslq(vec, maxcoeff=5000, maxsteps=50000, tol=mp.mpf(10)**-18))
