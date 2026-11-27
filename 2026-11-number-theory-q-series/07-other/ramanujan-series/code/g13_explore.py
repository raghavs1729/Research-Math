import sys, json, mpmath as mp
from explore import *
mp.mp.dps = 40
full = [n for n in ALL if n.startswith('Eb019')][0]
F = make_F(full, [-1], -1, ['0'], '0')
R = wz_mate(F); assert check_wz(F, R)
sg, rh = F.ratio('k'), F.ratio('n')
print('sigma =', sg); print('rho =', rh); print('R =', R)
sig = pari_to_func(sg, ['n', 'k']); rho = pari_to_func(rh, ['n', 'k']); Rf = pari_to_func(R, ['n', 'k'])
# triangle function: F(n,k) for 0 <= k < n, normalised by F(1,0) = 1, built from the ratios
def Ftri(n, k):
    # F(1,0)=1; go up in n at k=0 (rho(n,0) finite for n>=1?), then right in k
    v = mp.mpf(1)
    for j in range(1, n):
        v *= rho(mp.mpf(j), mp.mpf(0))
    for j in range(0, k):
        v *= sig(mp.mpf(n), mp.mpf(j))
    return v
print('rho(n,0):', [mp.nstr(rho(mp.mpf(j), mp.mpf(0)), 8) for j in range(1, 6)])
print('F(n,k) triangle:', [[mp.nstr(Ftri(n, k), 6) for k in range(n)] for n in range(1, 6)])
# target summand
t = lambda j: mp.rf(1, j)*mp.rf(mp.mpf(5)/6, j)*mp.rf(mp.mpf(7)/6, j)/mp.rf(mp.mpf(3)/2, j)**3*(7 + 10*j)*(mp.mpf(-9)/16)**j
print('G(n,0)/t(n-1):', [mp.nstr(Rf(mp.mpf(n), mp.mpf(0))*Ftri(n, 0)/t(n-1), 15) for n in range(1, 7)])
# G(n,n): limit k -> n of R(n,k) F(n,k): F(n,k) = F(n,n-1)*sigma(n,n-1+eps..)  -> compute R(n,k)*sigma(n,k-1) at k -> n
def Gdiag(n):
    e = mp.mpf(10)**-20
    return Rf(mp.mpf(n), mp.mpf(n) + e)*sig(mp.mpf(n), mp.mpf(n) - 1 + e)*Ftri(n, n-1)
d = lambda n: Ftri(n+1, n) + Gdiag(n)
print('d_n:', [mp.nstr(d(n), 12) for n in range(1, 8)])
print('d_{n+1}/d_n:', [mp.nstr(d(n+1)/d(n), 12) for n in range(1, 8)])
T = lambda n: mp.fsum(Ftri(n, k) for k in range(n))
for N in (2, 4, 8, 16, 32):
    lhs = mp.fsum(Rf(mp.mpf(n), mp.mpf(0))*Ftri(n, 0) for n in range(1, N))
    rhs = T(1) - T(N) + mp.fsum(d(n) for n in range(1, N))
    print(N, mp.nstr(lhs, 20), mp.nstr(rhs, 20), 'T_N =', mp.nstr(T(N), 12))
print('G(n,0)/t(n):', [mp.nstr(Rf(mp.mpf(n), mp.mpf(0))*Ftri(n, 0)/t(n), 15) for n in range(1, 7)])
c0 = Rf(mp.mpf(1), mp.mpf(0))*Ftri(1, 0)/t(1)
S = mp.nsum(t, [0, mp.inf])
print('sum t =', S, ' target', mp.sqrt(3)/4*(3*mp.log(3)-2*mp.log(2)))
SG = mp.nsum(lambda n: Rf(mp.mpf(n), mp.mpf(0))*Ftri(int(n), 0), [1, mp.inf])
print('sum_{n>=1} G(n,0) =', SG, ' = c0*(S - t(0)) ?', c0*(S - t(0)))
D = mp.nsum(lambda n: d(int(n)), [1, mp.inf])
print('T_1 =', T(1), ' sum d_n =', D, ' T_1 + sum d =', T(1) + D)
# closed form guess for d_n: compare with (-1/3)^n-type or 3^n/4^n C(2n,n)-type terms
for n in range(1, 6):
    print(n, mp.nstr(d(n), 15), mp.identify(d(n)*mp.binomial(2*n, n)*(2*n+1)*(2*n+3)), mp.nstr(d(n)*4**n/3**n, 12))
