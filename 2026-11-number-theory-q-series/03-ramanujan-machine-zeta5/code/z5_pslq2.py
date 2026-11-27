from mpmath import mp, mpf, zeta, log, pi, pslq, catalan
mp.dps=300
exec(open('z5_values.py').read())
print('U1 check', U1-96/(372*zeta(5)-pi**5))
G=catalan
bases={
 'A':[zeta(5), pi**5],
 'B':[zeta(5), pi**5, pi**2*zeta(3), pi**4*log(2), zeta(3)],
 'C':[zeta(5), pi**2*zeta(3), pi**4, pi**2, zeta(3)],
 'D':[zeta(5), pi**4*log(2), pi**2*zeta(3), pi**4, pi**2, zeta(3), log(2)],
 'E':[zeta(5), pi**5, pi**3, pi, pi**2*zeta(3), zeta(3), pi**4, pi**2],
 'F':[zeta(5), pi**2*zeta(3), pi**4*log(2)],
 'G':[zeta(5), pi**2*zeta(3)],
 'H':[zeta(5), pi*G*0+pi**3, pi**5, pi],
}
for name,v in (('U2',U2),('U3',U3)):
    for bn,K in bases.items():
        vec=[mpf(1)]+[v*k for k in K]
        r=pslq(vec,maxcoeff=10**15,maxsteps=400000)
        print(name,bn,r,flush=True)
        vec=[mpf(1)]+K+[v]+[v*k for k in K]
        r=pslq(vec,maxcoeff=10**12,maxsteps=400000)
        print(name,bn,'mobius',r,flush=True)
