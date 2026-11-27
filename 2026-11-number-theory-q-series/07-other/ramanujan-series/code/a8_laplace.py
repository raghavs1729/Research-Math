import mpmath as mp
mp.mp.dps = 30
def loga(N, k):
    N = mp.mpf(N); k = mp.mpf(k)
    lg = mp.loggamma
    # |Gamma(5/4-N)| / (|Gamma(1/4-N)| |Gamma(3/4-N)|) via reflection: Gamma(z)Gamma(1-z) = pi/sin(pi z)
    # |Gamma(5/4-N)| = pi/(|sin(pi(5/4-N))| Gamma(N-1/4)) etc.
    s = mp.sin(mp.pi/4)
    pref = (mp.log(mp.pi/s) - lg(N-mp.mpf(1)/4)) - (mp.log(mp.pi/s) - lg(N+mp.mpf(3)/4)) - (mp.log(mp.pi/s) - lg(N+mp.mpf(1)/4))
    val = (mp.log(2*k+1) + (2-4*N)*mp.log(2) + pref + lg(N+mp.mpf(3)/4) - lg(2*N-1)
           + lg(2*k+N+mp.mpf(1)/4) + lg(2*k+N+mp.mpf(3)/4) - 2*lg(N-k) - 2*lg(k+N+1)
           - mp.re(lg(2*k-N+mp.mpf(11)/4)) - mp.re(lg(2*k-N+mp.mpf(9)/4)))
    return val
target = mp.sqrt(2)/(4*mp.pi**2)
for N in (5, 10, 50, 200, 1000, 5000):
    T = mp.fsum(mp.exp(loga(N, k)) for k in range(0, N))
    kmax = max(range(N), key=lambda k: loga(N, k))
    print(N, mp.nstr(T, 18), mp.nstr(T - target, 5), 'argmax k/N =', kmax/N, ' a(kmax)*sqrt(N) =', mp.nstr(mp.exp(loga(N,kmax))*mp.sqrt(N), 12))
print('target sqrt2/(4 pi^2) =', mp.nstr(target, 18), ' sqrt3/2 =', mp.nstr(mp.sqrt(3)/2, 8))
