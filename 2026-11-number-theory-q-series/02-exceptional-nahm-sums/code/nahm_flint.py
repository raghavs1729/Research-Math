"""Fast exact evaluation of (generalized) Nahm sums with python-flint.

nahm(M, B, d, N, den) = sum_{n in N^r} q^{(1/2) n^T M n + B.n} / prod_i (q^{d_i};q^{d_i})_{n_i}
returned as a list of N*den integer coefficients in t = q^(1/den)  (truncation O(q^N)).
M: symmetric positive definite (Fractions), B: Fractions, d: positive ints.
Pruning uses a float lower bound (unconstrained minimum over the remaining variables, convex in each
coordinate); it is only used with a safety margin, and every exponent is computed exactly.
"""
from fractions import Fraction as Fr
from flint import fmpz_poly
import sys
sys.setrecursionlimit(10000)

def _inv(M):
    n=len(M); A=[[Fr(x) for x in row]+[Fr(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        pv=A[c][c]; A[c]=[x/pv for x in A[c]]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return [row[n:] for row in A]

def trunc(p,P):
    if p.degree()>=P:
        return fmpz_poly(p.coeffs()[:P])
    return p

def nahm(M,B,d,N,den=1,weights=None,count=False):
    r=len(B); P=N*den
    M=[[Fr(x) for x in row] for row in M]; B=[Fr(x) for x in B]
    # exact integer data for exponent: 2*den*e = n^T (den M) n + 2 den B.n
    M2=[[x*den for x in row] for row in M]; B2=[2*den*x for x in B]
    for row in M2:
        for x in row: assert x.denominator==1, "den too small for M"
    for x in B2: assert x.denominator==1, "den too small for B"
    M2=[[int(x) for x in row] for row in M2]; B2=[int(x) for x in B2]
    # Schur-complement lower bounds: for prefix of length k:
    #   lb = (1/2) x^T S x + L.x + c0
    S=[None]*(r+1); L=[None]*(r+1); c0=[0.0]*(r+1)
    for k in range(1,r+1):
        if k<r:
            Myy=[[M[i][j] for j in range(k,r)] for i in range(k,r)]
            Mi=_inv(Myy)
            Myx=[[M[i][j] for j in range(k)] for i in range(k,r)]
            By=[B[i] for i in range(k,r)]
            m=r-k
            # w = Myx x + By ; lb = 1/2 x Mxx x + Bx.x - 1/2 w^T Mi w
            Sk=[[M[i][j]-sum(Myx[a][i]*Mi[a][b]*Myx[b][j] for a in range(m) for b in range(m)) for j in range(k)] for i in range(k)]
            Lk=[B[i]-sum(Myx[a][i]*Mi[a][b]*By[b] for a in range(m) for b in range(m)) for i in range(k)]
            ck=-sum(By[a]*Mi[a][b]*By[b] for a in range(m) for b in range(m))/2
        else:
            Sk=[row[:] for row in M]; Lk=B[:]; ck=Fr(0)
        S[k]=[[float(x) for x in row] for row in Sk]; L[k]=[float(x) for x in Lk]; c0[k]=float(ck)
    # 1/(q^d;q^d)_n tables, in t
    need=set(d); E={}
    for dd in need:
        tab=[fmpz_poly([1])]; n=0
        while True:
            n+=1
            step=dd*n*den
            if step>=P and n>1: break
            geo=[0]*P
            for j in range(0,P,step): geo[j]=1
            tab.append(trunc(tab[-1]*fmpz_poly(geo),P))
            if step>=P: break
        E[dd]=tab
    acc=[0]*P; x=[0]*r; cnt=[0]
    xacc={} if weights is not None else None
    eps=1e-6
    nonneg = all(v>=0 for row in M for v in row) and all(v>=0 for v in B)
    Mf=[[float(v) for v in row] for row in M]; Bf=[float(v) for v in B]
    def lbu(k):   # unconstrained minimum over remaining variables (convex in x)
        Sk=S[k]; Lk=L[k]; v=c0[k]
        for i in range(k):
            xi=x[i]
            if xi:
                row=Sk[i]; t=0.0
                for j in range(k): t+=row[j]*x[j]
                v+=xi*(Lk[i]+0.5*t)
        return v
    def lb0(k):   # value with remaining variables = 0 (exact constrained min when M,B >= 0; monotone)
        v=0.0
        for i in range(k):
            xi=x[i]
            if xi:
                row=Mf[i]; t=0.0
                for j in range(k): t+=row[j]*x[j]
                v+=xi*(Bf[i]+0.5*t)
        return v
    def rec(k,cur):
        if k==r:
            e2=sum(x[i]*M2[i][j]*x[j] for i in range(r) for j in range(r))+sum(B2[i]*x[i] for i in range(r))
            assert e2%2==0, "exponent not in (1/den)Z"
            e=e2//2
            if e<P:
                cnt[0]+=1
                cs=cur.coeffs()
                if weights is None:
                    for j,c in enumerate(cs[:P-e]): acc[e+j]+=int(c)
                else:
                    w=sum(weights[i]*x[i] for i in range(r))
                    a=xacc.setdefault(w,[0]*P)
                    for j,c in enumerate(cs[:P-e]): a[e+j]+=int(c)
            return
        tab=E[d[k]]; i=0; prev=-1e300
        while True:
            x[k]=i
            if nonneg:
                b=lb0(k+1)
                if b>=N+eps: break
                if i<len(tab): rec(k+1,trunc(cur*tab[i],P))
            else:
                b=lbu(k+1)
                if b>=N+eps:
                    if i>0 and b>=prev: break
                else:
                    if i<len(tab): rec(k+1,trunc(cur*tab[i],P))
                prev=b
            i+=1
        x[k]=0
    rec(0,fmpz_poly([1]))
    if count: print("terms:",cnt[0],file=sys.stderr)
    return acc if weights is None else xacc

def qprod_inv(residues,m,P):
    """coefficients of prod_{n>=1, n mod m in residues} 1/(1-t^n) mod t^P"""
    c=[0]*P; c[0]=1
    for n in range(1,P):
        if n%m in residues:
            for j in range(n,P): c[j]+=c[j-n]
    return c
def qprod(residues,m,P):
    """coefficients of prod_{n>=1, n mod m in residues} (1-t^n) mod t^P"""
    c=[0]*P; c[0]=1
    for n in range(1,P):
        if n%m in residues:
            for j in range(P-1,n-1,-1): c[j]-=c[j-n]
    return c
def mul(a,b,P):
    return [int(v) for v in (trunc(fmpz_poly(a)*fmpz_poly(b),P).coeffs()+[0]*P)[:P]]
