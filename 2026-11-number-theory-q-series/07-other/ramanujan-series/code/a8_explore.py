import mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
mp.mp.dps = 40
F = make_F('Seed9', [-2,-1,-1], 1, ['1','1/4','3/4'], 0)
R = wz_mate(F); assert check_wz(F, R)
print('F =', F)
print('den R:', pari.factor(pari.denominator(R)))
print('deg num R in k,n:', pari.poldegree(pari.numerator(R),'k'), pari.poldegree(pari.numerator(R),'n'))
W = WZNum(F, R)
for a in ('0.3', '0.1', '0.01', '0.001'):
    a = mp.mpf(a)
    try: SF = hsum_k(W, a)
    except Exception as ex: SF = repr(ex)[:50]
    SG = hsum_n(W, a)
    gs = lambda K: mp.fsum(W.Gv(n + a, K) for n in range(0, 60))
    g1, g2 = gs(4000), gs(16000)
    SFN = [hsum_k(W, N + a) for N in (10, 30)]
    print('a =', a); print('  SF', SF); print('  SG', SG, ' SG*a^2', SG*a**2); print('  sum g', g1, g2); print('  SFN', SFN, 'SFN*a^2', SFN[-1]*a**2)
    if not isinstance(SF, str): print('  check', SF + g2 - SG - SFN[-1])
print('F(a,k) for a=0.3, k=0..4:', [mp.nstr(W.Fv(0.3,k),8) for k in range(5)], ' k=200,400:', mp.nstr(W.Fv(0.3,200),8), mp.nstr(W.Fv(0.3,400),8))
