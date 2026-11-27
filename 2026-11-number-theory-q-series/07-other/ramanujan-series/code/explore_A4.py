import mpmath as mp
from wz import *
from seeds_extra import SEEDS, load_ebisu
from numerics import *
mp.mp.dps = 50
EB = load_ebisu()
name = [n for n in EB if n.startswith('Eb014')][0]
params, T = EB[name]
print(name); print(T)
F = T.subs({'a': '1/2+n', 'k': 'k-2*n'})
print(F)
R = wz_mate(F); assert check_wz(F, R)
print('R =', pari.factor(pari.numerator(R)), '/', pari.factor(pari.denominator(R)))
W = WZNum(F, R)
eps = mp.mpf(10)**(-30)
for n in range(0, 6):
    print(n, 'G(n,0)~', mp.nstr(W.Gv(n, eps), 20), ' g(n)~', mp.nstr(W.Gv(n, 4000), 12), mp.nstr(W.Gv(n, 16000), 12))
# target summand
def t(k): return (4*k-1)*mp.mpf(3)**k/((2*k-1)*k**2*mp.binomial(2*k,k)*mp.binomial(3*k,k))
for k in range(1,6): print(k, mp.nstr(t(k),20), mp.nstr(W.Gv(k-1, eps)/t(k), 20))
print('sum F(0,k)=', mp.nsum(lambda k: W.Fv(0,k), [0, mp.inf]))
print('sum F(3,k)=', mp.nsum(lambda k: W.Fv(3,k), [6, mp.inf]))
