import mpmath as mp, sys
mp.mp.dps=30
q=int(sys.argv[1])
chi=[0,1,0,-1] if q==4 else [0,1,-1]
L=lambda s: mp.dirichlet(s,chi)
LL=lambda s: mp.diff(L,s)/L(s)
b=LL(0)
c=mp.limit(lambda s: LL(s)-1/(s+1), -1) if False else None
# constant term at -1 via symmetric difference
eps=mp.mpf('1e-8')
c=(LL(-1+eps)-1/eps + LL(-1-eps)+1/eps)/2
print("b=",b," b0(const term of L'/L at -1)=",c)
zs=[mp.mpf(l) for l in open('/home/user/weighted-prime-races/prototype/explore/chi%d_zeros.txt'%q)]
rhos=[mp.mpc(0.5,g) for g in zs]+[mp.mpc(0.5,-g) for g in zs]
s2=sum(1/abs(r)**2 for r in rhos)
print("partial sum 1/|rho|^2",s2," closed form",mp.log(q/mp.pi)-mp.euler+2*LL(1))
SG=sum(abs(mp.gamma(r)) for r in rhos); print("S_Gamma (listed)",SG)
# direct psi_e
N=60000
import math
def lam(n):
    # von Mangoldt
    for p in range(2,int(n**0.5)+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return math.log(p) if m==1 else 0.0
    return math.log(n)
sieve=[0.0]*(N+1)
isp=[True]*(N+1)
for p in range(2,N+1):
    if isp[p]:
        for m in range(p*p,N+1,p): isp[m]=False
        pk=p
        while pk<=N: sieve[pk]=math.log(p); pk*=p
for y in [1,2,5,20,100,1000]:
    y=mp.mpf(y)
    lim=min(N,int(60*y)+10)
    direct=mp.fsum(chi[n%q]*sieve[n]*mp.e**(-n/y) for n in range(2,lim+1) if sieve[n] and chi[n%q])
    Z=-mp.fsum(mp.gamma(r)*y**r for r in rhos)
    R1=y**-2*mp.fsum(1/(r*(r+1))*y**(r+2)*mp.gammainc(r+2,0,1/y) for r in rhos)
    R4=y**-2*mp.quad(lambda t: (-mp.log(1-t**-2)/2 - mp.atanh(1/t)*t + 1)*mp.e**(-t/y),[1,1.001,1.1,2,mp.inf])
    form=Z - b*(1+1/y)*mp.e**(-1/y) + ((1+c)*mp.e**(-1/y)+mp.e1(1/y))/y + R1 - R4
    print("y=%s direct=%s  formula=%s  diff=%s"%(y,mp.nstr(direct,15),mp.nstr(form,15),mp.nstr(direct-form,3)))
