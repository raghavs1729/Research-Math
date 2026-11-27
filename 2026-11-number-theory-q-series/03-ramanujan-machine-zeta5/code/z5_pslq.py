from mpmath import mp, mpf, zeta, log, pi, pslq, catalan, sqrt, gamma
mp.dps=300
exec(open('z5_values.py').read())
G=catalan
beta4=(zeta(4,mpf(1)/4)-zeta(4,mpf(3)/4))/4**4
sets={
 'zetas':  [mpf(1), zeta(3), zeta(5)],
 'zetas+log': [mpf(1), zeta(3), zeta(5), log(2)],
 'zetas+pi2': [mpf(1), zeta(3), zeta(5), pi**2, pi**4],
 'all1': [mpf(1), zeta(3), zeta(5), log(2), pi**2, pi**4],
 'pi': [mpf(1), pi, pi**2, pi**3, pi**4, pi**5],
 'invpi': [mpf(1), 1/pi, 1/pi**2, 1/pi**3, 1/pi**4],
 'mix': [mpf(1), zeta(3), zeta(5), pi**2, pi**4, pi**2*zeta(3), log(2), pi**2*log(2)],
 'cat': [mpf(1), G, beta4, pi, pi**3, pi**2, zeta(3)],
 'gam': [mpf(1), gamma(mpf(1)/4)**4/pi**3, gamma(mpf(1)/4)**8/pi**6, 1/pi, 1/pi**2, gamma(mpf(1)/4)**4/pi, pi**2/gamma(mpf(1)/4)**4*1],
}
for name,v in (('U1',U1),('U2',U2),('U3',U3)):
    for sn,K in sets.items():
        vec=K+[v*k for k in K]
        r=pslq(vec,maxcoeff=10**12,maxsteps=200000)
        print(name,sn,r, flush=True)
