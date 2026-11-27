import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 40
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
def t(k): return (77*k**3-86*k**2+29*k-3)*mp.mpf(-32)**k/(k**3*(3*k-1)*(3*k-2)*mp.binomial(3*k,k)*mp.binomial(4*k,k)**2)
for (sh, K, base, k0, j) in (([-1,1],1,['-1','2'],'3/2',2), ([-1,1],1,['0','1'],'1/2',1), ([-1,0],2,['0','1'],'1/2',1)):
    F = make_F('Seed1', sh, K, base, k0)
    print('F =', F)
    R = wz_mate(F); assert check_wz(F, R)
    print('R =', show(R))
    W = WZNum(F, R)
    e = mp.mpf('1e-15')
    print(' G(n+e,0)/(e t(n+j)):', [mp.nstr(W.Gv(n+e,0)/e/t(n+j), 15) for n in range(3)])
    print(' F(e,k)/e:', [mp.nstr(W.Fv(e,k)/e,12) for k in range(5)])
    print(' sum_k F(e,k)/e =', mp.nstr(hsum_k(W, e)/e, 15), '  sum_n G(n+e,0)/e =', mp.nstr(hsum_n(W, e)/e, 15), ' catalan', mp.catalan)
    print(' g(n+e)/e at K=1e5,1e7,1e9:', [[mp.nstr(W.Gv(n+e,KK)/e, 8) for KK in (10**5,10**7,10**9)] for n in range(3)])
    print(' L: sum_k F(N+e,k)/e:', [mp.nstr(hsum_k(W, N+e)/e, 8) for N in (1,2,5,10)])
