import flint, sys, time
ell=19; p=17
def kron(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
D=int(sys.argv[1])
P=lambda c: flint.nmod_poly(c,ell)
phi=[0]*D; n=0
while n*n<D: phi[n*n]=1 if n==0 else 2; n+=1
F=[0]*D
for d in range(1,D,2):
    for m in range(d,D,2*d): F[m]=(F[m]+d)%ell
phin=P(phi); Fn=P(F)
pb=phin.inverse_series_trunc(D)          # sum pbar(n)(-q)^n
W=[int(pb[n]) if kron(n,ell)==1 else 0 for n in range(D)]
# is W in Lambda_{719/2} (x) F_19 ?  basis B_j = F^j phi^(719-4j), j=0..179
K=719; dim=K//4+1
t0=time.time()
phi4=phin.pow_trunc(4,D)
# powers phi^(719-4j) = phi^3 * (phi^4)^(179-j): build descending list via repeated multiplication
pw=[None]*dim; cur=phin.pow_trunc(3,D)
for j in range(dim-1,-1,-1):
    pw[j]=cur; cur=cur.mul_low(phi4,D)
Fp=P([1]); B=[]
for j in range(dim):
    B.append(Fp.mul_low(pw[j],D)); Fp=Fp.mul_low(Fn,D)
print('basis built',time.time()-t0,flush=True)
def expand(vals):
    rem=list(vals[:dim]); co=[]
    for j in range(dim):
        a=rem[j]%ell; co.append(a)
        if a:
            for n in range(j,dim): rem[n]=(rem[n]-a*int(B[j][n]))%ell
    return co
def combo(co):
    s=P([0])
    for j,a in enumerate(co):
        if a: s=s+a*B[j]
    return s
w0=expand(W)
w0s=combo(w0)
mism=[n for n in range(D) if int(w0s[n])!=W[n]]
print('W in lattice mod 19? mismatches:',len(mism),mism[:5], ' (D=%d)'%D)
print('w0 coords:',w0)
lam=359
eps=kron(-1,p)**lam
c1=pow(p,lam-1,ell); c2=pow(p,2*lam-1,ell)
def b(n):
    v=W[p*p*n]+eps*kron(n,p)*c1*W[n]
    if n%(p*p)==0: v+=c2*W[n//(p*p)]
    return v%ell
beta=expand([b(n) for n in range(dim)])
bs=combo(beta)
N=(D-1)//(p*p)
print('T(289) consistency for n <',N,':',all(int(bs[n])==b(n) for n in range(N)))
for s in (1,-1):
    h=[(int(bs[n])-s*eps*c1*W[n])%ell for n in range(D)]
    for mod,res in ((8,None),(32,12)):
        for r in range(mod):
            nz=sum(1 for n in range(r,D,mod) if h[n])
            if nz==0: print(f'(n/17)={s:+d}: h(n)=0 mod 19 for ALL n = {r} mod {mod}, n < {D}')
    print(f'(n/17)={s:+d}: h nonzero count at n=12 mod 32:',sum(1 for n in range(12,D,32) if h[n]), ' c1=',c1)
