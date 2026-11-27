import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 50
name = [n for n in ALL if n.startswith('Eb024')][0]
params, T = ALL[name]
print(name); print('seed:', T)
F = make_F(name, [1], 1, ['1/2'], 0)
print('F =', F)
R = wz_mate(F); assert check_wz(F, R)
print('R =', R); print(pari.factor(pari.numerator(R)), pari.factor(pari.denominator(R)))
W = WZNum(F, R)
def t(k): return (3*k+2)**2*mp.mpf(16)**k/((2*k+1)*(6*k+1)*(6*k+5)*mp.binomial(6*k,3*k))
for n in range(4): print(n, mp.nstr(W.Gv(n,0)/t(n), 30))
print('sumF0 =', hsum_k(W, 0), ' 1/sqrt2 =', 1/mp.sqrt(2))
print('F(0,k) k=0..3:', [mp.nstr(W.Fv(0,k),15) for k in range(4)])
print('sum G =', hsum_n(W, 0), ' target*c: ', mp.pi/(2*mp.sqrt(3)))
print('sum_k F(n,k) n=1..:', [mp.nstr(hsum_k(W,n),8) for n in (1,2,3,5,10,20)])
print('|F| sums:', [mp.nstr(mp.fsum(abs(W.Fv(n,k)) for k in range(0,400)),8) for n in (1,2,3,5,10,20)])
