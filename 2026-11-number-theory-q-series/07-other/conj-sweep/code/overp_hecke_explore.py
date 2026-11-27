import flint, sys
from math import gcd
def kron(a,p):  # Legendre symbol, p odd prime
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def explore(ell, p, lam, coeffs, res8, D):
    # c = sum_j coeffs[j] F^j phi^(2*lam+1-4j)  (weight lam+1/2 on Gamma_0(4)), integer series to q^D
    Z=flint.fmpz_poly
    phi=[0]*D; n=0
    while n*n<D: phi[n*n]=1 if n==0 else 2; n+=1
    phi=Z(phi)
    F=[0]*D
    for n in range(1,D,2):
        F[n]=sum(d for d in range(1,n+1) if n%d==0)
    F=Z(F)
    def tr(x): return Z(x.coeffs()[:D])
    def pw(x,e):
        r=Z([1])
        for _ in range(e): r=tr(r*x)
        return r
    k2=2*lam+1
    basis=[tr(pw(F,j)*pw(phi,k2-4*j)) for j in range(len(coeffs))]
    dim=k2//4+1
    fullbasis=[tr(pw(F,j)*pw(phi,k2-4*j)) for j in range(dim)]
    c=Z([0])
    for j,a in enumerate(coeffs): c=c+a*basis[j]
    cc=[int(x) for x in c.coeffs()]+[0]*D
    # check c(n) = pbar(ell n)(-1)^n mod ell
    Dp=min(D, 4000)
    phim=flint.nmod_poly([int(x)%ell for x in phi.coeffs()[:Dp]],ell)
    pb=phim.inverse_series_trunc(Dp)   # sum pbar(n)(-q)^n
    okc=all((int(pb[ell*n])-cc[n])%ell==0 for n in range(Dp//ell))
    print(f'c(n) = pbar({ell}n)(-1)^n mod {ell} for n < {Dp//ell}:',okc)
    eps = (-1)**(lam*((p-1)//2))   # (-1/p)^lam
    def b(n):
        v=cc[p*p*n]+eps*kron(n,p)*p**(lam-1)*cc[n]
        if n%(p*p)==0: v+=p**(2*lam-1)*cc[n//(p*p)]
        return v
    # express b in the basis F^j phi^(k2-4j), j=0..dim-1 (triangular, leading coefficient 1 at q^j)
    al=[]; rem=[b(n) for n in range(dim)]
    fb=[[int(x) for x in f.coeffs()]+[0]*D for f in fullbasis]
    for j in range(dim):
        a=rem[j]; al.append(a)
        for n in range(dim): rem[n]-=a*fb[j][n]
    print('T(p^2)c in basis:',al, 'mod ell:',[a%ell for a in al])
    N=(D-1)//(p*p)
    bser=[sum(al[j]*fb[j][n] for j in range(dim)) for n in range(N)]
    print('consistency b(n) = c(p^2 n)+... for n <',N,':',all(bser[n]==b(n) for n in range(N)))
    for sgn in (1,-1):
        # for (n/p)=sgn: c(p^2 n) = b(n) - eps*sgn*p^(lam-1) c(n)
        h=[(bser[n]-eps*sgn*p**(lam-1)*cc[n])%ell for n in range(N)]
        for r in range(8):
            z=[n for n in range(N) if n%8==r]
            allz=all(h[n]==0 for n in z)
            sub=[n for n in z if kron(n,p)==sgn]
            print(f'  (n/p)={sgn:+d}: h(n)=0 mod {ell} for all n={r} mod 8: {allz};  for those with (n/{p})={sgn:+d}: {all(h[n]==0 for n in sub)}')
if sys.argv[2]=="17": explore(17,11,7,[1,13,13],3,int(sys.argv[1]))
if sys.argv[2]=="23": explore(23,13,10,[1,9,5,14,17,20],5,int(sys.argv[1]))
