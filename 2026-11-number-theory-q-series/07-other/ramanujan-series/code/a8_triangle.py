import mpmath as mp
from explore import *
mp.mp.dps = 40
def show(x):
    num, den = pari.numerator(x), pari.denominator(x)
    c = pari.content(num)/pari.content(den)
    return str(c) + ' * ' + str(pari.factor(num/pari.content(num))) + ' / ' + str(pari.factor(den/pari.content(den)))
F = make_F('Seed9', [-2,-1,-1], 1, ['1','1/4','3/4'], 0)
R = wz_mate(F); assert check_wz(F, R)
print('sigma =', show(F.ratio('k'))); print('rho =', show(F.ratio('n')))
print('R den =', pari.factor(pari.denominator(R)), ' deg num (k,n):', pari.poldegree(pari.numerator(R),'k'), pari.poldegree(pari.numerator(R),'n'))
print('R(n,0) =', show(pari.subst(R,'k',0)))
W = WZNum(F, R)
e = mp.mpf('1e-12')
ph = mp.exp(-1j*mp.pi*(mp.mpf(3)/4))   # (-1)^(3/4-n-eps) ~ e^{i pi (3/4 - n)}
def Phi(N, k): return (W.Fv(N+e, k)*e**2)
def Gam(n, k): return (W.Gv(n+e, k)*e**2)
def t(k): return (1344*k**3+944*k**2+156*k+9)*mp.binomial(2*k,k)**2*mp.binomial(4*k,2*k)**3/((2*k+1)*mp.mpf(2)**(22*k))
print('Gam(n,0)/t(n):', [mp.nstr(Gam(n,0)/t(n), 12) for n in range(0,5)])
print('Gam(n,0)/t(n-1):', [mp.nstr(Gam(n,0)/t(n-1), 12) for n in range(1,5)])
for N in (1,2,3,6):
    lhs = mp.fsum(Gam(n,0) for n in range(N)); rhs = -mp.fsum(Phi(N,k) for k in range(N))
    print(N, 'sum Gam(n,0) =', mp.nstr(lhs, 12), ' -T_N =', mp.nstr(rhs, 12), ' Gam(N-1,N-1)=', mp.nstr(Gam(N-1,N-1),8), 'Gam(N-1,N)=', mp.nstr(Gam(N-1,N),8))
print('Phi(8,k)*e^{i..}:', [mp.nstr(Phi(8,k), 6) for k in range(9)])
print('64 sqrt2/pi^2 * ? :', mp.nstr(mp.fsum(Gam(n,0) for n in range(60)), 15), mp.nstr((1-1j)/(4*mp.pi**2), 15))
