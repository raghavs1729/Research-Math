import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 40
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
def t(k): return mp.binomial(2*k,k)*mp.mpf(16)**k*(1656*k**4-2604*k**3+1158*k**2-229*k+15)/(mp.mpf(k)**2*(2*k-1)**2*(6*k-1)*(6*k-5)*mp.binomial(6*k,3*k)*mp.binomial(3*k,k))
F = make_F('Seed3', [-2,-2,-3], 2, ['1/2','0','1/2'], '1/2')
print('F =', F)
R = wz_mate(F); assert check_wz(F, R)
print('R =', show(R)); print('sigma =', show(F.ratio('k'))); print('rho=', show(F.ratio('n')))
W = WZNum(F, R)
e = mp.mpf('1e-18')
print(' G(n+e,0)/(e t(n+1)):', [mp.nstr(W.Gv(n+e,0)/e/t(n+1), 18) for n in range(4)])
print(' F(e,k)/e:', [mp.nstr(W.Fv(e,k)/e,15) for k in range(6)])
print(' F(e,k)/e * (2k+1)^2:', [mp.nstr(W.Fv(e,k)/e*(2*k+1)**2,15) for k in range(6)])
print(' sum t =', mp.nsum(t,[1,mp.inf]))
