import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 40
F = make_F('Seed4', [-1,-2], 1, ['0','-1/2'], '1/2')
print('F =', F)
R = wz_mate(F); assert check_wz(F, R)
print('R =', R); print(pari.factor(pari.numerator(R)), pari.factor(pari.denominator(R)))
print('sigma =', pari.factor(F.ratio('k'))); print('rho =', pari.factor(F.ratio('n')))
W = WZNum(F, R)
def t(k): return (8*k-3)*mp.binomial(4*k,2*k)/(k*(4*k-1)*mp.mpf(9)**k*mp.binomial(2*k,k)**2)
for n in range(4): print(n, mp.nstr(W.Gv(n,0)/t(n+1), 25), mp.nstr(-2*mp.sqrt(3)/mp.pi, 25))
print('F(0,k):', [mp.nstr(W.Fv(0,k)*(3*k+1)*(3*k+2),12) for k in range(4)], mp.nstr(-mp.sqrt(3)/mp.pi,12))
for n in (0,1,2): print('G(n,K)*sqrt(K) K=1e3,1e4,1e5:', [mp.nstr(W.Gv(n,K)*mp.sqrt(K),10) for K in (1000,10000,100000)], ' F(n,K)*K^?:', [mp.nstr(W.Fv(n,K)*K**2,8) for K in (1000,10000)])
print('sum_k F(N,k):', [mp.nstr(hsum_k(W,N),8) for N in (1,2,3,5,10)])
