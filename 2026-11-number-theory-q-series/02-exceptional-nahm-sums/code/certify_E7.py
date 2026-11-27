"""Certificate / numerical checks for proofs/E7.md.
 (a) exact: three fusion steps turn the rank-4 matrix of the proof into C(E7)^{-1};
 (b) exact polynomial check of  c_n = 1 + q^{n^2}  for 1 <= n <= NMAX  [proved in general in the text];
 (c) formula (T) for T(a), a <= 5, as power series mod q^NT;
 (d) the E7 identity itself, as power series in t = q^(1/4), mod q^NF.
Usage: python3 certify_E7.py [NF] [NT] [NMAX]
"""
import sys
from fractions import Fraction as Fr
from flint import fmpz_poly, fmpq_poly
from nahm_flint import nahm, qprod, qprod_inv, mul
NF=int(sys.argv[1]) if len(sys.argv)>1 else 60
NT=int(sys.argv[2]) if len(sys.argv)>2 else 60
NMAX=int(sys.argv[3]) if len(sys.argv)>3 else 14

def inv(M):
    n=len(M); A=[[Fr(x) for x in row]+[Fr(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        pv=A[c][c]; A[c]=[x/pv for x in A[c]]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return [row[n:] for row in A]
def fuse(A,names,a,b,new):
    row={k:A[a][k]+A[b][k] for k in names}
    diag=A[a][a]+A[b][b]+2*A[a][b]
    A[a][b]+=1; A[b][a]+=1
    for k in names: A[k][new]=row[k]
    A[new]=dict(row); A[new][new]=diag
    names.append(new)
def show(A,names):
    for i in names: print("   ",i.rjust(3),[str(A[i][j]) for j in names])

# (a) E7: chain 1-2-3-4-5-6, node 7 attached to node 4
C7=[[0]*7 for _ in range(7)]
for i in range(7): C7[i][i]=2
for (i,j) in [(1,2),(2,3),(3,4),(4,5),(5,6),(4,7)]: C7[i-1][j-1]=C7[j-1][i-1]=-1
A7=inv(C7)
print("C(E7)^{-1} ="); [print("   ",[str(x) for x in row]) for row in A7]
names=['1','6','2','5']     # variables (a, n1, n2, n3) of the proof
start=[[Fr(3,2),0,1,2],[0,2,2,2],[1,2,4,4],[2,2,4,6]]
A={x:{} for x in names}
for i,x in enumerate(names):
    for j,y in enumerate(names): A[x][y]=Fr(start[i][j])
print("start (nodes 1,6,2,5 <-> a, n1, n2, n3):"); show(A,names)
for (a,b,new) in [('1','6','7'),('5','6','4'),('1','2','3')]:
    fuse(A,names,a,b,new); print("after fusion  node %s := node %s + node %s"%(new,a,b)); show(A,names)
ok=all(A[str(i+1)][str(j+1)]==A7[i][j] for i in range(7) for j in range(7))
print("E7 CERTIFICATE: final matrix == C(E7)^{-1} :",ok); assert ok

# (b) c_n = (1-q^n) sum_{r,j>=0, 2r+j<=n} (-1)^r q^{(3r^2-r)/2+j^2+2rj} (q)_{n-r-1}/((q)_r (q)_j (q)_{n-2r-j})
x=fmpq_poly([0,1])
def qp(n):
    p=fmpq_poly([1])
    for k in range(1,n+1): p=p*(1-x**k)
    return p
for n in range(1,NMAX+1):
    # common denominator D = (q)_n^3 is enough? use exact rational arithmetic via cross-multiplication with L = prod of all three
    num=fmpq_poly([0]); L=qp(n)*qp(n)*qp(n)
    for r in range(0,n//2+1):
        for j in range(0,n-2*r+1):
            den=qp(r)*qp(j)*qp(n-2*r-j)
            quo=L//den; assert quo*den==L
            num+= (-1)**r * x**((3*r*r-r)//2+j*j+2*r*j) * qp(n-r-1) * quo
    lhs=(1-x**n)*num; rhs=(1+x**(n*n))*L
    assert lhs==rhs, n
print("(b) c_n == 1 + q^(n^2) exactly for 1 <= n <= %d : True"%NMAX)

# (c) T(a) = sum q^{N1^2+N2^2+N3^2+a(N2+N3)}/((q)_{n1}(q)_{n2}(q)_{n3})  vs formula (T)
P=NT
def ser(poly): return [int(c) for c in (poly.coeffs()+[0]*P)[:P]]
def qpz(n):
    p=fmpz_poly([1])
    for k in range(1,n+1): p=p*fmpz_poly([1]+[0]*(k-1)+[-1])
    return p
pinf=qprod_inv({0},1,P)
def series_div(num,den):   # num/den as truncated series, den(0)=1
    a=ser(num); d=ser(den); out=[0]*P
    for i in range(P):
        s=a[i]-sum(d[k]*out[i-k] for k in range(1,i+1))
        out[i]=s
    return out
for a in range(0,6):
    lhs=nahm([[2,2,2],[2,4,4],[2,4,6]],[0,a,2*a],[1,1,1],P,1)
    tot=[0]*P
    r=0
    while 3*a*r+(9*r*r-r)//2<P:
        e=3*a*r+(9*r*r-r)//2
        # factor (1-q^{a+2r})(q)_{a+r-1}/(q)_r  (for a=r=0: 1)   times  sum_j [a;j] q^{j^2+2rj}
        if a==0 and r==0: fac=fmpz_poly([1]); dn=fmpz_poly([1])
        else: fac=fmpz_poly([1]+[0]*(a+2*r-1)+[-1])*qpz(a+r-1); dn=qpz(r)
        sj=fmpz_poly([0])
        for j in range(0,a+1):
            g=qpz(a); dd=qpz(j)*qpz(a-j); gq=g//dd; assert gq*dd==g
            sj+=gq*fmpz_poly([0]*(j*j+2*r*j)+[1])
        term=series_div(fac*sj,dn)
        for i in range(P-e): tot[e+i]+=(-1)**r*term[i]
        r+=1
    rhs=mul(tot,pinf,P)
    print("(c) a=%d: T(a) == formula (T) mod q^%d : %s"%(a,P,lhs==rhs)); assert lhs==rhs

# (d) E7 identity in t=q^(1/4)
den=4; PP=NF*den
f=nahm(A7,[0]*7,[1]*7,NF,den,count=True)
th=[0]*PP
k=0
while 3*k*k<PP:
    th[3*k*k]+=1
    if k>0 and 7*k*k<PP: th[7*k*k]+=1
    k+=1
pin=[0]*PP
for i,c in enumerate(qprod_inv({0},1,NF)): pin[4*i]=c
print("(d) sum_{n in N^7} q^{n C(E7)^-1 n/2}/(q)_n == (1/(q)_inf)(sum_{k>=0} q^{3k^2/4} + sum_{k>=1} q^{7k^2/4})  mod q^%d : %s"%(NF,f==mul(th,pin,PP)))
assert f==mul(th,pin,PP)
print("ALL CHECKS PASSED")
