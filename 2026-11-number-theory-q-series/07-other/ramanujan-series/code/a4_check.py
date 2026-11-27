import mpmath as mp
mp.mp.dps = 60
pi = mp.pi
def c(n): return (mp.mpf(1)/9)**n * mp.gamma(3*n+mp.mpf(3)/2)*mp.factorial(n)**2/(4*mp.factorial(2*n)*mp.factorial(3*n)*mp.gamma(n+mp.mpf(5)/6)*mp.gamma(n+mp.mpf(7)/6))
def F(n,x): return c(n)*(mp.mpf(3)/4)**x*mp.rf(-x,2*n)*mp.gamma(x+n+1)/mp.gamma(x+n+mp.mpf(3)/2)
def R(n,x): return (-x**2-(2*n+3)*x+8*n**2+6*n)/((3*n+1)*(3*n+2))
def G(n,x): return R(n,x)*F(n,x)
# WZ check at random real x
for (n,x) in [(0,mp.mpf('0.3')),(2,mp.mpf('1.7')),(5,mp.mpf('7.25'))]:
    print('WZ', mp.nstr(F(n+1,x)-F(n,x)-(G(n,x+1)-G(n,x)), 5))
def t(k): return (4*k-1)*mp.mpf(3)**k/((2*k-1)*k**2*mp.binomial(2*k,k)*mp.binomial(3*k,k))
dG = lambda n: mp.diff(lambda x: G(n,x), 0)
for n in range(4): print(n, mp.nstr(dG(n)/t(n+1), 30), mp.nstr(-3/(4*pi),30))
K = (mp.zeta(2,mp.mpf(1)/3)-mp.zeta(2,mp.mpf(2)/3))/9
# D_N by direct sum and by the integral formula
def D_direct(N, kmax=600):
    return mp.fsum(mp.diff(lambda x: F(N,x), k) for k in range(0,kmax))
def Psi(N,y): return -mp.factorial(2*N)*y**(2*N)/(1-y)**(2*N+1)*mp.quad(lambda s: s**(2*N)/(1+s), [0,(1-y)/y])
def D_int(N): return c(N)/mp.sqrt(pi)*mp.quad(lambda t_: t_**N*(1-t_)**(-mp.mpf(1)/2)*Psi(N,3*t_/4), mp.linspace(0,1,9))
mp.mp.dps = 30
print('D0 direct', D_direct(0), ' -3K/(2pi) =', -3*K/(2*pi))
for N in (1,2,3): print(N, 'D_N direct', mp.nstr(D_direct(N),15), 'integral', mp.nstr(D_int(N),15))
print('sum t', mp.nsum(t,[1,mp.inf]), '2K', 2*K)
