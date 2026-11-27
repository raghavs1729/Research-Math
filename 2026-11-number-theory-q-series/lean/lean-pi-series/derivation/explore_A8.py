"""Explore all scan hits for Sun's 64 sqrt2/pi^2 series (A8 = N2): is there a pair with an easy limit?"""
import sys
sys.path.insert(0, '../../ramanujan-series/code')
import mpmath as mp
from explore import make_F, ALL
from wz import *
from numerics import WZNum
mp.mp.dps = 40
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
def t(k): return (1344*k**3+944*k**2+156*k+9)*mp.binomial(2*k,k)**2*mp.binomial(4*k,2*k)**3/((2*k+1)*mp.mpf(2)**(22*k))
HITS = [('Seed9',[-2,-1,-1],1,['1','1/4','3/4'],'0'), ('Seed9',[-2,-1,-1],1,['3/4','1/4','1/4'],'0'), ('Seed9',[-2,-1,-1],1,['3/4','1/4','1/4'],'1/4'),
        ('Seed9',[-2,-1,-1],1,['5/4','3/4','3/4'],'0'),
        ('Dixon6F5',[-1,0,0],1,['1/2','1/2','1/2'],'1/4'), ('Dixon6F5',[-1,0,0],1,['1/2','1/4','3/4'],'0'),
        ('Dixon6F5',[-1,0,0],1,['1/4','1/4','1/4'],'0'), ('Dixon6F5',[-1,0,0],1,['3/4','3/4','3/4'],'0')]
sel = [int(x) for x in sys.argv[1:]] or range(len(HITS))
for i in sel:
    name, sh, K, base, k0 = HITS[i]
    print('=' * 100); print(i, HITS[i])
    F = make_F(name, sh, K, base, k0)
    R = wz_mate(F); assert check_wz(F, R)
    print('F =', F)
    print('sigma =', show(F.ratio('k'))); print('rho =', show(F.ratio('n')))
    print('R den =', pari.factor(pari.denominator(R)), ' deg num (k,n):', pari.poldegree(pari.numerator(R),'k'), pari.poldegree(pari.numerator(R),'n'))
    W = WZNum(F, R)
    for e in (mp.mpf('1e-8'), mp.mpf('1e-12')):
        def sf(f, *a):
            try: return mp.nstr(f(*a), 10)
            except Exception as ex: return 'ERR'
        print(' eps=%s' % mp.nstr(e, 3))
        print('  G(n+e,0)/t(n):', [sf(lambda n: W.Gv(n+e,0)/t(n), n) for n in range(0,4)])
        print('  G(n+e,0)/t(n-1):', [sf(lambda n: W.Gv(n+e,0)/t(n-1), n) for n in range(1,4)])
        print('  F(e,k):', [sf(W.Fv, e, k) for k in range(5)]); print('  F(1+e,k):', [sf(W.Fv, 1+e, k) for k in range(5)]); print('  F(3+e,k):', [sf(W.Fv, 3+e, k) for k in range(7)])
        for n in (0,1,2): print('  G(%d+e,K): K=1e3,1e5:'%n, [sf(W.Gv, n+e, K) for K in (10**3,10**5)])
