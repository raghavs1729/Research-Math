"""Numerical modularity test (NOT a proof) for the Sun--Wang pairs (X, A_{k-1}) of proofs/X-A-parafermion.md.

If f(tau) = q^{-c/24} F(q) is a modular function (weight 0, finite-order multiplier) then, as q = e^{-eps} -> 1,
log f(e^{-eps}) = alpha/eps + beta + (exponentially small), so D(eps) := eps*log f(e^{-eps}) is linear in eps up to
exponentially small terms.  A wrong prefactor q^{-(c+x)/24} adds -x*eps^2/24 to D.  The script prints the second
divided differences of D on eps = 0.5, 0.4, 0.3, 0.25 (40-digit arithmetic), and for contrast the same with c + 1/10.
alpha*6/pi^2 is the effective central charge.  Output: certificates/xa_modularity_numeric.txt
"""
import numpy as np, itertools, math, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40

def gram(kind, l):
    if kind == 'A':  G=[[2 if i==j else (-1 if abs(i-j)==1 else 0) for j in range(l)] for i in range(l)]; d=[1]*l; h=l+1
    if kind == 'D':
        G=[[2 if i==j else 0 for j in range(l)] for i in range(l)]
        for i in range(l-2): G[i][i+1]=G[i+1][i]=-1
        G[l-3][l-1]=G[l-1][l-3]=-1; d=[1]*l; h=2*l-2
    if kind == 'C':  # short,...,short,long
        G=[[0]*l for _ in range(l)]
        for i in range(l): G[i][i]=2 if i<l-1 else 4
        for i in range(l-1): G[i][i+1]=G[i+1][i]=-1 if i<l-2 else -2
        d=[1]*(l-1)+[2]; h=2*l
    if kind == 'B':  # long,...,long,short
        G=[[0]*l for _ in range(l)]
        for i in range(l): G[i][i]=4 if i<l-1 else 2
        for i in range(l-1): G[i][i+1]=G[i+1][i]=-2
        d=[2]*(l-1)+[1]; h=2*l
    if kind == 'G': G=[[2,-3],[-3,6]]; d=[1,3]; h=6
    if kind == 'F': G=[[2,-1,0,0],[-1,2,-2,0],[0,-2,4,-2],[0,0,-2,4]]; d=[1,1,2,2]; h=12
    if kind == 'T':
        G=[[0]*l for _ in range(l)]
        for i in range(l): G[i][i]=2 if i<l-1 else 1
        for i in range(l-1): G[i][i+1]=G[i+1][i]=-1
        d=[1]*l; h=2*l+1
    return G,d,h

def F_at(G,d,k,q,tol=mp.mpf(10)**(-30)):
    l=len(G)
    var=[(i,s) for i in range(l) for s in range(1,k)]
    R=len(var)
    M=[[mp.mpf(G[i][j])*(min(s,t)-mp.mpf(s*t)/k) for (j,t) in var] for (i,s) in var]
    dd=[d[i] for (i,s) in var]
    # 1/(q^d;q^d)_n tables
    lam=min(np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in M])))
    L=-mp.log(q)
    nmax=int(math.sqrt(2*90/ (float(L)*lam)))+2
    ip={}
    for dv in set(dd):
        t=[mp.mpf(1)]
        for n in range(1,nmax+1): t.append(t[-1]/(1-q**(dv*n)))
        ip[dv]=t
    # recursive summation with pruning by Cholesky bound
    Mf=np.array([[float(x) for x in r] for r in M])
    Lc=np.linalg.cholesky(Mf)
    total=mp.mpf(0)
    p=[0]*R
    bound=2*90/float(L)
    def rec(i,acc):
        nonlocal total
        if i<0:
            E=sum(M[a][b]*p[a]*p[b] for a in range(R) for b in range(R))/2
            term=q**E
            for a in range(R): term*=ip[dd[a]][p[a]]
            total+=term; return
        c=sum(Lc[j][i]*p[j] for j in range(i+1,R))
        v=0
        while v<=nmax:
            y=Lc[i][i]*v+c
            if acc+y*y<=bound:
                p[i]=v; rec(i-1,acc+y*y)
            elif y>0: break
            v+=1
        p[i]=0
    rec(R-1,0.0)
    return total

def test(kind,l,k,epss):
    G,d,h=gram(kind,l)
    c=mp.mpf(sum(d)*(k-1)*h)/(h+k)
    Ds=[]
    for e in epss:
        q=mp.e**(-mp.mpf(e))
        F=F_at(G,d,k,q)
        f=F*q**(-c/24)
        Ds.append(mp.mpf(e)*mp.log(f))
    # divided differences
    xs=[mp.mpf(e) for e in epss]
    dd1=[(Ds[i+1]-Ds[i])/(xs[i+1]-xs[i]) for i in range(len(xs)-1)]
    dd2=[(dd1[i+1]-dd1[i])/(xs[i+2]-xs[i]) for i in range(len(xs)-2)]
    # same with a wrong prefactor (c+1/10) for contrast
    Dw=[Ds[i]-xs[i]*xs[i]*(mp.mpf(1)/240) for i in range(len(xs))]
    w1=[(Dw[i+1]-Dw[i])/(xs[i+1]-xs[i]) for i in range(len(xs)-1)]
    w2=[(w1[i+1]-w1[i])/(xs[i+2]-xs[i]) for i in range(len(xs)-2)]
    ceff=Ds[0]*6/mp.pi**2  # crude
    print("(%s_%d, A_%d) c=%s: slope %s ; alpha*6/pi^2 ~ %s ; 2nd div.diff %s   [with c+1/10: %s]"%(kind,l,k-1,mp.nstr(c,8),mp.nstr(dd1[0],10),mp.nstr((Ds[0]-dd1[0]*xs[0])*6/mp.pi**2,10),[mp.nstr(x,3) for x in dd2],[mp.nstr(x,3) for x in w2]))
    sys.stdout.flush()

if __name__=="__main__":
    eps=[0.5,0.4,0.3,0.25]
    for case in [('A',1,2),('A',2,2),('C',2,2),('G',2,2),('C',3,2),('B',3,2),('T',1,2),('T',1,4),('T',3,2),('D',4,2),('T',1,6)]:
        test(*case,eps)
