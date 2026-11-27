import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 40
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
for (sh, K) in (([1,2],1), ([-3,0],1)):
    F = make_F('Seed1', sh, K, ['1/2','1'], 0)
    print('F =', F)
    R = wz_mate(F); assert check_wz(F, R)
    print('R =', show(R))
    print('sigma =', show(F.ratio('k'))); print('rho =', show(F.ratio('n')))
    W = WZNum(F, R)
    def t(k): return (360*k**3+612*k**2+230*k+15)*mp.binomial(3*k,k)*mp.binomial(6*k,3*k)**2/((3*k+1)*(3*k+2)*mp.mpf(2)**(15*k))
    print('G(n,0)/t(n):', [mp.nstr(W.Gv(n,0)/t(n), 20) for n in range(4)])
    print('F(0,k):', [mp.nstr(W.Fv(0,k),10) for k in range(4)], 'F(1,k):', [mp.nstr(W.Fv(1,k),10) for k in range(4)])
    for n in (0,1,2,3): print(' g(%d): K=1e4,1e6,1e8:'%n, [mp.nstr(W.Gv(n,K),15) for K in (10**4,10**6,10**8)])
    gl = [W.Gv(n,mp.mpf(10)**12) for n in range(60)]
    print(' sum g ~', mp.nstr(mp.fsum(gl),15), ' ratios g(n+1)/g(n):', [mp.nstr(gl[i+1]/gl[i],12) for i in range(4)])
    print(' sum G(n,0) =', mp.nstr(hsum_n(W,0),20), ' 1/(4 sqrt pi)=', mp.nstr(1/(4*mp.sqrt(mp.pi)),20))
