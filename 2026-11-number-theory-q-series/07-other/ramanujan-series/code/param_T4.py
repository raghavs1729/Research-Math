"""Parameter derivatives of the Theorem-4 pair: F(n,k;a,b,c) = Seed1(1/2+a+n, 1+c+2n, k+b+n)."""
import mpmath as mp, time
from wz import *
from numerics import pari_to_func
mp.mp.dps = 40
params, T = SEEDS['Seed1']
F = T.subs({'a': '1/2+a+n', 'c': '1+c+2*n', 'k': 'k+b+n'})
print('F =', F)
t0 = time.time(); R = wz_mate(F); assert check_wz(F, R); print('R found', round(time.time()-t0, 1), 's')
print('R(n,0) =', pari.factor(pari.subst(R, 'k', 0)))
Rf = pari_to_func(R, ['n', 'k', 'a', 'b', 'c'])
def Fv(n, k, a, b, c):
    return F.evalf({'n': mp.mpf(n), 'k': mp.mpf(k), 'a': mp.mpf(a), 'b': mp.mpf(b), 'c': mp.mpf(c)}, mp)
def SG(a, b, c): return mp.fsum(Rf(mp.mpf(n), mp.mpf(0), mp.mpf(a), mp.mpf(b), mp.mpf(c))*Fv(n, 0, a, b, c) for n in range(0, 700))
H = mp.harmonic
t = lambda k: (360*k**3+612*k**2+230*k+15)*mp.binomial(3*k,k)*mp.binomial(6*k,3*k)**2/((3*k+1)*(3*k+2)*mp.mpf(2)**(15*k))
h = mp.mpf('1e-12')
base = SG(h, 0, 0)   # at exactly 0 the term n=0 is 0*inf
print('SG(~0) =', mp.nstr(base, 20), ' 1/(4 sqrt pi) =', mp.nstr(1/(4*mp.sqrt(mp.pi)), 20))
d = {}
for nm, f in (('a', lambda x: SG(x, 0, 0)), ('b', lambda x: SG(h, x, 0)), ('c', lambda x: SG(h, 0, x))):
    e = mp.mpf('1e-8')
    d[nm] = (f(2*e) - f(e))/e if nm == 'a' else (f(e) - f(-e))/(2*e)
    print('d/d%s SG =' % nm, mp.nstr(d[nm], 15))
c0 = (1/(4*mp.sqrt(mp.pi)))/(32*mp.sqrt(2)/mp.pi)
B46 = mp.fsum(mp.binomial(3*k,k)*mp.binomial(6*k,3*k)**2*((360*k**3+612*k**2+230*k+15)*(H(2*k)-H(k)) - 180*k**2+36*k+23)/((3*k+1)*(3*k+2)*mp.mpf(2)**(15*k)) for k in range(700))
print('B46 lhs*c0 =', mp.nstr(B46*c0, 15), ' check 48 sqrt2 log2/pi:', mp.nstr(B46 - 48*mp.sqrt(2)*mp.log(2)/mp.pi, 5))
print('pslq [B46*c0, da, db, dc, base, base*log2]:', mp.pslq([B46*c0, d['a'], d['b'], d['c'], base, base*mp.log(2)], maxcoeff=200, maxsteps=20000))
