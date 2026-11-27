#!/usr/bin/env python3
"""
audit2 / a3: numerical tests of the limiting steps (not of the final identities, which a5 checks).

T4 : lim_K Q(0,K)h(0,K) = 1/(sqrt2 pi);  sum_k h(n,k) <= 7 h(n,0) (n>=2);  the identity after K -> oo:
     sum_{n<N} t(n)/64 + N sum_k h(N,k) = 1/(sqrt2 pi)  for N = 1..8.
T10: constancy in b of Phi(b) = sum_n Q(n,b) h(n,b)  (real second variable), |b| <= 1/10.
T8 : constancy in s of sum_n Q_s(n,0) h_s(n,0) = 1/(4 sqrt pi), |s| <= 1/10 (N_s read from ../certificates/T8.gp).
T5 : sum_k psi(0,k) = 2G; after K -> oo: sum_{n<N} gamma(n,0) = 2G - sum_k psi(N,k), N = 1..6; |sum_k psi(N,k)| <= 6*8^-N.
T7 : sum_{k<N} Phi(N,k) -> 0 geometrically.
T3 : G'(n,0) = -(3/(4 pi)) t(n+1);  D_0 = -3K/(2 pi);  I = -4K;  sum_{n<N} G'(n,0) = D_0 - D_N and the bound on D_N.
T2 : sum_k f(0,k) = pi/(3 sqrt3), g(n,0) = 2 t(n+1); T9: Phi(c) = (psi(c+2/3)-psi(c+1/3))/3 for |c| <= 1/10.
T6 : Phi(a,b) = sum_n G(n,0;a,b) = sum_k F(0,k;a,b) at a few (a,b); sum_k F(0,k;a,0) = 2^(-1/2).
"""
import os, re, sys
from math import comb
from mpmath import (mp, mpf, pi, sqrt, log, gamma, rgamma, rf, hyp3f2, catalan, zeta, quad, diff, nsum, inf,
                    psi as polygamma, nstr, binomial, loggamma, exp)
mp.dps = 40
OK = True
def check(cond, msg):
    global OK
    OK &= bool(cond)
    print(('  ok   ' if cond else '  FAIL ') + msg)
K3 = (zeta(2, mpf(1)/3) - zeta(2, mpf(2)/3))/9
HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------- T4 / T10
print('T4')
def h(n, x):      # real second variable; h(0,0) = 1
    return mpf(2)**(-2*x-5*n)*rf(mpf(1)/2, 3*n)*gamma(2*x+3*n+mpf(1)/2)/(gamma(mpf(1)/2)*rf(mpf(1)/2, n)*gamma(x+n+1)*gamma(n+1)*gamma(x+3*n+1))
def Q(n, k): n = mpf(n); k = mpf(k); return (64*k**3+(160*n+144)*k**2+(240*n**2+432*n+92)*k+360*n**3+612*n**2+230*n+15)/(64*(k+3*n+1)*(k+3*n+2))
def t4(n):   return mpf((360*n**3+612*n**2+230*n+15)*comb(3*n, n)*comb(6*n, 3*n)**2)/((3*n+1)*(3*n+2)*2**(15*n))
lim = 1/(sqrt(mpf(2))*pi)
check(abs(1/(gamma(mpf(1)/4)*gamma(mpf(3)/4)) - lim) < mpf(10)**(-35), '1/(Gamma(1/4)Gamma(3/4)) = 1/(sqrt2 pi)')
vals = [(Kk, Kk*(Q(0, mpf(Kk))*h(0, mpf(Kk)) - lim)) for Kk in (10**2, 10**4, 10**6, 10**8)]
print('      K (Q(0,K)h(0,K) - 1/(sqrt2 pi)) for K = 1e2,1e4,1e6,1e8:', [nstr(v, 8) for _, v in vals])
check(all(abs(v) < 1 for _, v in vals), 'Q(0,K) h(0,K) -> 1/(sqrt2 pi) with error O(1/K)')
def rowsum(n):    # sum_k h(n,k)/h(n,0) = 3F2(3n/2+1/4, 3n/2+3/4, 1; n+1, 3n+1; 1).
    # mpmath's 3F2 at z = 1 is reliable for small n (cross-checked below) but loses everything for n >= 40,
    # where the terms decay fast and direct summation is used instead (accuracy ~1e-15, enough for "<= 7").
    if n <= 10:
        return hyp3f2(mpf(3*n)/2+mpf(1)/4, mpf(3*n)/2+mpf(3)/4, 1, n+1, 3*n+1, 1)
    s_, term, k = mpf(0), mpf(1), 0
    while term > mpf(10)**(-15):
        s_ += term
        term *= mpf((4*k+6*n+1)*(4*k+6*n+3))/(16*(k+n+1)*(k+3*n+1)); k += 1
    return s_ + term*(k+3*n+1)/n
rs = {n: rowsum(n) for n in (1, 2, 3, 5, 10, 40, 200, 1000)}
def rowsum_direct(n, tol=mpf(10)**(-45)):
    s_, term, k = mpf(0), mpf(1), 0
    while term > tol and k < 4*10**6:
        s_ += term
        term *= mpf((4*k+6*n+1)*(4*k+6*n+3))/(16*(k+n+1)*(k+3*n+1)); k += 1
    return s_, term*(k+3*n+1)/n
d3, tail3 = rowsum_direct(10, mpf(10)**(-30))
check(abs(d3 - rowsum(10)) <= 2*tail3 + mpf(10)**(-30), 'cross-check of the two ways of computing sum_k h(10,k)/h(10,0) (3F2 at 1 versus direct summation with tail bound)')
print('      sum_k h(n,k)/h(n,0):', {n: nstr(v, 6) for n, v in rs.items()})
check(all(v <= 7 for n, v in rs.items() if n >= 2) and rs[1] <= 9, 'sum_k h(n,k) <= 7 h(n,0) for the tested n >= 2 (and <= 9 h(1,0)); values tend to 4')
ok = True
for N in range(1, 9):
    lhs = sum(t4(n) for n in range(N))/64 + N*h(N, 0)*rowsum(N)
    ok &= abs(lhs - lim) < mpf(10)**(-33)
check(ok, 'after K -> oo:  sum_{n<N} t(n)/64 + N sum_k h(N,k) = 1/(sqrt2 pi)  for N = 1..8 (33 digits)')
print('T10 (real second variable)')
for b in (mpf(-1)/10, mpf(-1)/20, mpf(1)/20, mpf(1)/10):
    s = sum(Q(n, b)*h(n, b) for n in range(0, 700))
    check(abs(s - lim) < mpf(10)**(-30), f'Phi({nstr(b,3)}) = sum_n Q(n,b)h(n,b) = 1/(sqrt2 pi)   (difference {nstr(s-lim,3)})')
# derivative identity of T10 at the level of the summand
def dlogh(n): return diff(lambda x: log(h(n, x)), 0)
H = lambda m: polygamma(0, m+1) + mp.euler
check(all(abs(dlogh(n) - ((4*H(6*n)-3*H(3*n)-H(n)) - 6*log(2))) < mpf(10)**(-25) for n in range(0, 6)),
      'd/dx log h(n,x)|_0 = 4H_{6n} - 3H_{3n} - H_n - 6 log 2')
P = lambda n: 360*n**3+612*n**2+230*n+15
f10 = lambda n: mpf(1296*n**3+1980*n**2+912*n+139)/((3*n+1)*(3*n+2))
check(all(abs(diff(lambda x: log(Q(n, x)), 0) - f10(n)/P(n)) < mpf(10)**(-25) for n in range(0, 6)), 'd/dx log Q(n,x)|_0 = f(n)/P(n)')

# ------------------------------------------------------------------------------------------- T8
print('T8 (parameter s)')
src = open(os.path.join(HERE, '..', 'certificates', 'T8.gp')).read()
Nsrc = re.search(r'^N\s*=\s*(.*?);', src, flags=re.S | re.M).group(1).replace('^', '**')
def Qs(n, k, s): return eval(Nsrc, {'__builtins__': {}}, {'n': mpf(n), 'k': mpf(k), 's': s})/(64*(2*n+2*s+1)*(k+3*n+1)*(k+3*n+2))
def hs(n, k, s): return mpf(2)**(-2*k-5*n+s-mpf(3)/2)*gamma(3*n+s+mpf(1)/2)*gamma(2*k+3*n-s+mpf(1)/2)/(gamma(n+s+mpf(1)/2)*gamma(k+n-s+1)*gamma(n+1)*gamma(k+3*n+1))
for s in (mpf(-1)/10, mpf(-1)/20, mpf(0), mpf(1)/20, mpf(1)/10):
    v = sum(Qs(n, 0, s)*hs(n, 0, s) for n in range(0, 700))
    check(abs(v - 1/(4*sqrt(pi))) < mpf(10)**(-30), f's = {nstr(s,3)}: sum_n Q_s(n,0)h_s(n,0) = 1/(4 sqrt pi)   (difference {nstr(v-1/(4*sqrt(pi)),3)})')
check(all(abs(diff(lambda s: log(hs(n, 0, s)), 0) - (3*log(2) - 2*(H(2*n)-H(n)))) < mpf(10)**(-25) for n in range(6)),
      'd/ds log h_s(n,0)|_0 = 3 log 2 - 2(H_{2n}-H_n)')
check(all(abs(diff(lambda s: log(Qs(n, 0, s)), 0) - mpf(360*n*n-72*n-46)/P(n)) < mpf(10)**(-25) for n in range(6)),
      'd/ds log Q_s(n,0)|_0 = (360n^2-72n-46)/P(n)')

# ------------------------------------------------------------------------------------------- T5
print('T5')
def sg5(n, k):  return mpf((2*k+n+1)*(2*k+n+2))/((2*k+2*n+3)*(2*k+4*n+3))
def rho50(n):   return -mpf((n+1)*(3*n+1)*(3*n+2)*(3*n+3)*(n+1))/((4*n+1)*(4*n+3)*(2*n+3)*(4*n+3)*(4*n+5))
def R5(n, k):   return mpf((16*n+12)*k**2+(66*n**2+86*n+28)*k+77*n**3+145*n**2+88*n+17)/((4*n+1)*(4*n+3)*(2*k+4*n+3))
psi0 = [mpf(1)]
for n in range(12): psi0.append(psi0[-1]*rho50(n))
def rowsum5(n):   # sum_k psi(n,k)/psi(n,0) = 3F2((n+1)/2, (n+2)/2, 1; n+3/2, 2n+3/2; 1)
    return hyp3f2(mpf(n+1)/2, mpf(n+2)/2, 1, n+mpf(3)/2, 2*n+mpf(3)/2, 1)
check(abs(rowsum5(0) - 2*catalan) < mpf(10)**(-30), 'sum_k psi(0,k) = sum 4^k/((2k+1)^2 C(2k,k)) = 2G')
check(abs(quad(lambda th: th/mp.sin(th), [0, pi/2]) - 2*catalan) < mpf(10)**(-30), 'int_0^{pi/2} theta/sin(theta) = 2G')
ok = True; ok2 = True
for N in range(1, 7):
    lhs = sum(R5(n, 0)*psi0[n] for n in range(N)); tail = psi0[N]*rowsum5(N)
    ok &= abs(lhs - (2*catalan - tail)) < mpf(10)**(-30)
    ok2 &= abs(psi0[N])*hyp3f2(mpf(N+1)/2, mpf(N+2)/2, 1, N+mpf(3)/2, 2*N+mpf(3)/2, 1) <= 6*mpf(8)**(-N)
check(ok, 'after K -> oo:  sum_{n<N} gamma(n,0) = 2G - sum_k psi(N,k),  N = 1..6 (30 digits)')
check(ok2, 'sum_k |psi(N,k)| <= 6 * 8^-N for N = 1..6')
g0K = [R5(0, Kk)*mpf(4)**0*exp(Kk*log(4) + 2*loggamma(Kk+1) - loggamma(2*Kk+1))/(2*Kk+1)**2*sqrt(mpf(Kk)) for Kk in (10**2, 10**4, 10**6)]
print('      gamma(0,K) sqrt(K), K = 1e2,1e4,1e6:', [nstr(v, 8) for v in g0K], ' (bounded, so gamma(0,K) = O(K^-1/2) -> 0)')

# ------------------------------------------------------------------------------------------- T7
print('T7')
def Phi7(n, k):
    return ((-1)**(k+1)*(4*k+4*n+3)*mpf(2)**(1-6*n)*gamma(k+1)*gamma(4*n+1)*gamma(2*n+1)
            / (pi*gamma(n-k)*(2*k+2*n+1)*(2*k+2*n+2)*gamma(k+2*n+mpf(3)/2)*gamma(k+3*n+mpf(3)/2)))
vals = [6*pi**2*sum(Phi7(N, k) for k in range(N)) for N in (1, 2, 5, 10, 20, 40)]
print('      6 pi^2 sum_{k<N} Phi(N,k), N = 1,2,5,10,20,40:', [nstr(v, 6) for v in vals])
check(abs(vals[0] + mpf(32)/75) < mpf(10)**(-30) and all(abs(vals[i+1]) < abs(vals[i]) for i in range(5)) and abs(vals[-1]) < mpf(10)**(-30),
      'equals the partial sum t(1) = -32/75 at N = 1 and tends to 0 geometrically (ratio about 4/27 per step)')

# ------------------------------------------------------------------------------------------- T3
print('T3')
def c3(n): return mpf(9)**(-n)*gamma(3*n+mpf(3)/2)*gamma(n+1)**2/(4*gamma(2*n+1)*gamma(3*n+1)*gamma(n+mpf(5)/6)*gamma(n+mpf(7)/6))
def F3(n, x): return c3(n)*(mpf(3)/4)**x*rf(-x, 2*n)*gamma(x+n+1)/gamma(x+n+mpf(3)/2)
def R3(n, x): n = mpf(n); return (8*n*n+6*n-(2*n+3)*x-x*x)/((3*n+1)*(3*n+2))
def G3(n, x): return R3(n, x)*F3(n, x)
t3 = lambda k: mpf((4*k-1)*3**k)/((2*k-1)*k**2*comb(2*k, k)*comb(3*k, k))
xs = [mpf('0.3'), mpf('1.7'), mpf('4.25')]
check(all(abs(F3(n+1, x) - F3(n, x) - (G3(n, x+1) - G3(n, x))) < mpf(10)**(-30) for n in range(4) for x in xs), 'WZ equation for real x (n = 0..3, three non-integer x)')
dG = [diff(lambda x: G3(n, x), 0) for n in range(7)]
check(all(abs(dG[n] + 3/(4*pi)*t3(n+1)) < mpf(10)**(-25) for n in range(7)), "G'(n,0) = -(3/(4 pi)) t(n+1), n = 0..6")
def D(n): return nsum(lambda k: diff(lambda x: F3(n, x), k), [0, inf])
D0 = D(0)
check(abs(D0 + 3*K3/(2*pi)) < mpf(10)**(-22), f'D_0 = sum_k F\'(0,k) = -3K/(2 pi)   (difference {nstr(D0 + 3*K3/(2*pi), 3)})')
I = quad(lambda t: log(3*t/4)/((1-3*t/4)*sqrt(1-t)), [0, mpf(1)/2, mpf(9)/10, mpf(99)/100, 1])
check(abs(I + 4*K3) < mpf(10)**(-20), f'I = int_0^1 log(3t/4)/((1-3t/4) sqrt(1-t)) dt = -4K   (difference {nstr(I + 4*K3, 3)})')
ok = True
for N in range(1, 5):
    DN = D(N)
    bound = c3(N)*gamma(2*N+1)*4*mp.beta(N, mpf(1)/2)/(3*sqrt(pi)*(2*N+1))
    ok &= abs(sum(dG[:N]) - (D0 - DN)) < mpf(10)**(-20) and abs(DN) <= bound
check(ok, "sum_{n<N} G'(n,0) = D_0 - D_N and |D_N| <= c_N (2N)! 4 B(N,1/2)/(3 sqrt(pi)(2N+1)),  N = 1..4")

# ------------------------------------------------------------------------------------------- T2 / T9
print('T2 / T9')
def f2(n, x):
    Fv = (2*x+n+1)*mpf(3)**(mpf(1)/2-2*n)*gamma(n+1)*gamma(x+mpf(1)/2)*gamma(3*x+n+1)/(gamma(-2*n-mpf(1)/2)*gamma(x+n+mpf(3)/2)*gamma(n+mpf(1)/2)*gamma(3*x+2*n+3))
    return -(pi/sqrt(mpf(3)))*Fv
def R2(n, x): return (3*x+2*n+2)*(6*x+8*n+5)/(9*(2*n+1)*(2*x+n+1))
t2 = lambda k: mpf((8*k-3)*comb(4*k, 2*k))/(k*(4*k-1)*9**k*comb(2*k, k)**2)
check(all(abs(f2(0, mpf(k)) - mpf(1)/((3*k+1)*(3*k+2))) < mpf(10)**(-30) for k in range(8)), 'f(0,k) = 1/((3k+1)(3k+2)) from the Gamma closed form')
check(all(abs(R2(mpf(n), mpf(0))*f2(n, mpf(0)) - 2*t2(n+1)) < mpf(10)**(-30) for n in range(8)), 'g(n,0) = 2 t(n+1)')
for c in (mpf(-1)/10, mpf(1)/20, mpf(1)/10):
    lhs = sum(R2(mpf(n), c)*f2(n, c) for n in range(0, 60))
    rhs = (polygamma(0, c+mpf(2)/3) - polygamma(0, c+mpf(1)/3))/3
    check(abs(lhs - rhs) < mpf(10)**(-30), f'T9: Phi({nstr(c,3)}) = sum_n g(n,c) = (psi(c+2/3)-psi(c+1/3))/3   (difference {nstr(lhs-rhs,3)})')
check(abs((polygamma(1, mpf(2)/3) - polygamma(1, mpf(1)/3))/3 + 3*K3) < mpf(10)**(-30), "Phi'(0) = (psi'(2/3)-psi'(1/3))/3 = -3K")

# ------------------------------------------------------------------------------------------- T6
print('T6')
def F6(n, k, a, b):
    u = k+3*n+2*a+b
    return ((-1)**(n+k)*mpf(3)**(-(k+n+b))*(mpf(3)/4)**(a+n+mpf(1)/2)*gamma(u+1)*gamma(k+b-a+mpf(1)/2)*gamma(2*n+2*a+mpf(3)/2)*gamma(n+a+1)**2
            / (gamma(2*n+2*a+1)*gamma(mpf(1)/2-a-n)*gamma(u+mpf(3)/2)*gamma(k+n+b+1)*gamma(n+a+mpf(3)/4)*gamma(n+a+mpf(5)/4)))
def R6(n, k, a, b):
    u = k+3*n+2*a+b
    return 3*(u+2)**2/(4*(u+2)**2-1)
for (a, b) in ((mpf(0), mpf(0)), (mpf(1)/20, mpf(-1)/10), (mpf(-1)/10, mpf(1)/10), (mpf(1)/10, mpf(0))):
    lhs = sum(R6(mpf(n), 0, a, b)*F6(n, 0, a, b) for n in range(80)); rhs = sum(F6(0, k, a, b) for k in range(90))
    extra = ''
    if b == 0:
        extra = f';  sum_k F(0,k;a,0) - 2^(-1/2) = {nstr(rhs - 1/sqrt(mpf(2)), 3)}'
        check(abs(rhs - 1/sqrt(mpf(2))) < mpf(10)**(-30), f'T6: sum_k F(0,k;{nstr(a,3)},0) = 2^(-1/2)')
    check(abs(lhs - rhs) < mpf(10)**(-30), f'T6: sum_n G(n,0;a,b) = sum_k F(0,k;a,b) at (a,b) = ({nstr(a,3)},{nstr(b,3)})   (difference {nstr(lhs-rhs,3)}){extra}')
t1 = lambda n: mpf((3*n+2)**2*16**n)/((2*n+1)*(6*n+1)*(6*n+5)*comb(6*n, 3*n))
check(all(abs(R6(mpf(n), 0, mpf(0), mpf(0))*F6(n, 0, mpf(0), mpf(0)) - sqrt(mpf(6))/pi*t1(n)) < mpf(10)**(-30) for n in range(8)), 'G(n,0;0,0) = (sqrt6/pi) t(n)')

print('ALL LIMIT-STEP CHECKS OK' if OK else 'SOME LIMIT-STEP CHECK FAILED')
sys.exit(0 if OK else 1)
