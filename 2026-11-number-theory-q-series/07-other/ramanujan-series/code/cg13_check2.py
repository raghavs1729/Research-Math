import mpmath as mp
from wz import *
from numerics import *
mp.mp.dps = 40
Fs = seed_instance('Dougall7F6', dict(a=-1,b=-1,c=-1,d=-1,e=-1), 2, consts=dict(a='1/2',b='1/2',c='1/2',d='1/2',e='1/2'), k0=0)
R = wz_mate(Fs); assert check_wz(Fs, R)
Rf = pari_to_func(R, ['n', 'k'])
# leading coefficient r(n) = lim R/k^3
Rn, Rd = pari.numerator(R), pari.denominator(R)
dn, dd = pari.poldegree(Rn, 'k'), pari.poldegree(Rd, 'k')
r = pari.pollead(Rn, 'k')/pari.pollead(Rd, 'k')
print('deg_k R num/den:', dn, dd, ' r(n) =', pari.factor(r))
rfun = pari_to_func(r, ['n'])
G = mp.gamma
def E(n, k):
    return (2*k+3*n+mp.mpf(1)/2)*G(n+0.5)**10*G(k+n+0.5)**5*G(k+4*n)*mp.rgamma(k-n+1.5)/(G(2*n)**5*G(k+2*n+1)**5)
def v(n):   # lim_k R E = r(n) * 2 * Gamma(n+1/2)^10/Gamma(2n)^5   (E ~ 2k * k^-4 * ...)
    return rfun(n)*2*G(n+0.5)**10/G(2*n)**5
def T1(a): return mp.nsum(lambda k: E(a, k), [0, mp.inf])
def T2(a): return mp.nsum(lambda m: (-1)**int(m)*v(a+m), [0, mp.inf])
def T3(a): return mp.nsum(lambda m: (-1)**int(m)*Rf(a+m, mp.mpf(0))*E(a+m, 0), [0, mp.inf])
# sanity of v: compare with R*E at large k
print('v check:', mp.nstr(v(mp.mpf('0.8')), 12), mp.nstr(Rf(mp.mpf('0.8'), mp.mpf(10**6))*E(mp.mpf('0.8'), 10**6), 12))
for a in (mp.mpf('0.3'), mp.mpf('0.45'), mp.mpf('0.7')):
    print('a', a, 'T1+T2-T3 =', mp.nstr(T1(a) + T2(a) - T3(a), 8), ' T1,T2,T3 =', mp.nstr(T1(a), 10), mp.nstr(T2(a), 10), mp.nstr(T3(a), 10))
h = mp.mpf(1)/2
print('at 1/2: T1 =', mp.nstr(T1(h), 20), ' 2 zeta(3) =', mp.nstr(2*mp.zeta(3), 20), ' T2 =', mp.nstr(T2(h), 20), ' T3 =', mp.nstr(T3(h), 8))
d1 = mp.diff(T1, h); d2 = mp.diff(T2, h); d3 = mp.diff(T3, h)
print("T1' =", mp.nstr(d1, 15), ' -9 zeta(4) =', mp.nstr(-9*mp.zeta(4), 15))
print("T2' =", mp.nstr(d2, 15), ' pi^4/6 =', mp.nstr(mp.pi**4/6, 15))
print("T3' =", mp.nstr(d3, 15), ' pi^4/15 =', mp.nstr(mp.pi**4/15, 15))
