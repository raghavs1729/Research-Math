"""Exact certificate for the fusion chains used in proofs/E8-E6.md, plus numerical confirmation.

Part 1 (exact, rational arithmetic): starting from the Andrews-Gordon matrix 2*min(i,j) (rank 4) and
applying four fusion steps (Lemma 2 of the proof) one obtains C(E8)^{-1} up to the stated relabelling;
starting from blockdiag(2 C(A2)^{-1}, (2)) and applying three fusion steps one obtains C(E6)^{-1}.
Part 2 (numerical): all identities of Corollaries 1 and 2 as power series to the requested order.
Usage: python3 certify_E8_E6.py [N_E8] [N_E6]
"""
import sys, time
from fractions import Fraction as Fr
from nahm_flint import nahm, qprod_inv, qprod, mul

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
    """Lemma 2: add species new = a+b.  A is a dict-of-dict keyed by names."""
    row={k:A[a][k]+A[b][k] for k in names}
    diag=A[a][a]+A[b][b]+2*A[a][b]
    A[a][b]+=1; A[b][a]+=1
    for k in names: A[k][new]=row[k]
    A[new]=dict(row); A[new][new]=diag
    names.append(new)

def show(A,names,order=None):
    order=order or names
    for i in order: print("   ",i.rjust(6),[str(A[i][j]) for j in order])

# ---------------------------------------------------------------- E8
# E8 Dynkin diagram: chain 1-2-3-4-5-6-7, node 8 attached to node 5.
C8=[[0]*8 for _ in range(8)]
for i in range(8): C8[i][i]=2
for (i,j) in [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(5,8)]:
    C8[i-1][j-1]=C8[j-1][i-1]=-1
A8=inv(C8)
print("C(E8)^{-1} ="); [print("   ",[str(x) for x in row]) for row in A8]
# start: AG species m1..m4 with matrix 2 min(i,j); identify m1=node1, m2=node7, m3=node2, m4=node8
names=['1','7','2','8']
A={x:{} for x in names}
for i,x in enumerate(names):
    for j,y in enumerate(names): A[x][y]=Fr(2*min(i+1,j+1))
print("\nstart (species = E8 nodes 1,7,2,8 <-> AG variables m1,m2,m3,m4):"); show(A,names)
for (a,b,new) in [('7','8','4'),('1','8','6'),('1','4','5'),('1','2','3')]:
    fuse(A,names,a,b,new)
    print("\nafter fusion  node %s := node %s + node %s :"%(new,a,b)); show(A,names)
ok=all(A[str(i+1)][str(j+1)]==A8[i][j] for i in range(8) for j in range(8))
print("\nE8 CERTIFICATE: final matrix == C(E8)^{-1} :",ok); assert ok
# composition matrix P (rows: E8 nodes, columns: m1..m4)
P8={'1':(1,0,0,0),'7':(0,1,0,0),'2':(0,0,1,0),'8':(0,0,0,1),'4':(0,1,0,1),'6':(1,0,0,1),'5':(1,1,0,1),'3':(1,0,1,0)}
P8=[P8[str(i+1)] for i in range(8)]
print("P (rows = nodes 1..8) =",P8)

# ---------------------------------------------------------------- E6
# E6 Dynkin diagram: chain 1-2-3-4-5, node 6 attached to node 3.
C6=[[0]*6 for _ in range(6)]
for i in range(6): C6[i][i]=2
for (i,j) in [(1,2),(2,3),(3,4),(4,5),(3,6)]:
    C6[i-1][j-1]=C6[j-1][i-1]=-1
A6=inv(C6)
print("\nC(E6)^{-1} ="); [print("   ",[str(x) for x in row]) for row in A6]
names=['1','5','6']
A={x:{} for x in names}
start=[[Fr(4,3),Fr(2,3),0],[Fr(2,3),Fr(4,3),0],[0,0,2]]
for i,x in enumerate(names):
    for j,y in enumerate(names): A[x][y]=Fr(start[i][j])
print("start (species = E6 nodes 1,5 <-> (i,j) of the A2 sum; node 6 <-> Rogers-Ramanujan variable):"); show(A,names)
for (a,b,new) in [('5','6','2'),('1','6','4'),('1','2','3')]:
    fuse(A,names,a,b,new)
    print("\nafter fusion  node %s := node %s + node %s :"%(new,a,b)); show(A,names)
ok=all(A[str(i+1)][str(j+1)]==A6[i][j] for i in range(6) for j in range(6))
print("\nE6 CERTIFICATE: final matrix == C(E6)^{-1} :",ok); assert ok
P6={'1':(1,0,0),'5':(0,1,0),'6':(0,0,1),'2':(0,1,1),'4':(1,0,1),'3':(1,1,1)}
P6=[P6[str(i+1)] for i in range(6)]
print("P (rows = nodes 1..6; columns = (i, j, m)) =",P6)

# ---------------------------------------------------------------- numerics
N8=int(sys.argv[1]) if len(sys.argv)>1 else 150
N6=int(sys.argv[2]) if len(sys.argv)>2 else 80
t0=time.time()
print("\n--- Corollary 1 (E8, modulus 11), power series mod q^%d ---"%N8)
for i in range(5,0,-1):
    b=[max(0,j-i+1) for j in range(1,5)]
    B=[sum(P8[n][j]*b[j] for j in range(4)) for n in range(8)]
    f=nahm(A8,B,[1]*8,N8,1)
    res=set(range(1,11))-{i,11-i}
    g=qprod_inv(res,11,N8)
    print(" i=%d  B=%s : E8 sum == prod_{n != 0,+-%d mod 11} 1/(1-q^n) : %s   (%.0fs)"%(i,B,i,f==g,time.time()-t0)); sys.stdout.flush()
    assert f==g
print("\n--- Corollary 2 (E6), power series in q^(1/3) mod q^%d ---"%N6)
den=3; P=N6*den
def lift(c,k):   # series in q^k -> series in t=q^(1/3) of length P : exponent multiplied by 3k
    out=[0]*P
    for i,v in enumerate(c):
        if 3*k*i<P: out[3*k*i]=v
    return out
def shift(c,s):  # multiply by t^s
    return ([0]*s+c)[:P] if s>=0 else c[-s:]+[0]*(-s)
def theta_J(a,m,PP):
    """J_{a,m} = (q^a,q^{m-a},q^m;q^m)_inf as series in q of length PP"""
    c=[0]*PP
    k=0
    while True:
        done=True
        for kk in ([k,-k] if k else [0]):
            e=(m*kk*kk+(m-2*a)*kk)//2     # sum (-1)^k q^{m k(k-1)/2 + a k} with k -> -k variants
            e=m*kk*(kk-1)//2+a*kk
            if e<PP: c[e]+=(-1)**(kk%2); done=False
        if done and k>0: break
        k+=1
    return c
G=qprod_inv({1,4},5,N6); H=qprod_inv({2,3},5,N6)
# Example 10 (Cao-Rosengren-Wang (1.6a),(1.6b)), written in q (their q^3 -> q):
#   S4: sum q^{(2/3)(i^2+ij+j^2)}/((q)_i(q)_j) = (J_{7,15} - q J_{2,15} + 2 q^{2/3} J_{3,15})/J_1
#   S3: sum q^{(2/3)(i^2+ij+j^2)-(2i+j)/3}/((q)_i(q)_j) = (2 J_{6,15} + q^{1/3} J_{4,15} + q^{4/3} J_{1,15})/J_1
J1inv=lift(qprod_inv({0},1,N6),1)
def comb(terms):
    tot=[0]*P
    for coef,sh,a in terms:
        s=shift(lift(theta_J(a,15,N6),1),sh)
        tot=[x+coef*y for x,y in zip(tot,s)]
    return mul(tot,J1inv,P)
S4=comb([(1,0,7),(-1,3,2),(2,2,3)])
S3=comb([(2,0,6),(1,1,4),(1,4,1)])
A2=[[Fr(4,3),Fr(2,3)],[Fr(2,3),Fr(4,3)]]
z0=nahm(A2,[0,0],[1,1],N6,den); z1=nahm(A2,[Fr(-2,3),Fr(-1,3)],[1,1],N6,den)
print(" rank-2 input (Zagier Ex. 10, proved by Cherednik-Feigin / Cao-Rosengren-Wang): B=0:",z0==S4,"  B=(-2/3,-1/3):",z1==S3)
assert z0==S4 and z1==S3
for (bi,bj,nm,ser) in [(Fr(0),Fr(0),'S4',S4),(Fr(-2,3),Fr(-1,3),'S3',S3),(Fr(-1,3),Fr(-2,3),'S3',S3)]:
    for (bm,rn,rs) in [(0,'G',G),(1,'H',H)]:
        b=(bi,bj,Fr(bm))
        B=[sum(P6[n][j]*b[j] for j in range(3)) for n in range(6)]
        f=nahm(A6,B,[1]*6,N6,den)
        g=mul(lift(rs,1),ser,P)
        print(" B=%s : E6 sum == %s(q) * %s : %s   (%.0fs)"%([str(x) for x in B],rn,nm,f==g,time.time()-t0)); sys.stdout.flush()
        assert f==g
# Sun-Wang (3.19): RHS in Virasoro characters of M(6,5):  G(q) * (chi_{1,1} + q^3 chi_{4,1} + 2 q^{2/3} chi_{1,3}) (normalised)
def vir(r,s,PP):
    c=[0]*PP
    for k in range(-PP,PP+1):
        e1=k*(30*k+6*r-5*s); e2=(5*k+r)*(6*k+s)
        if 0<=e1<PP: c[e1]+=1
        if 0<=e2<PP: c[e2]-=1
    return c
pinv=qprod_inv({0},1,N6)
from flint import fmpz_poly
def m1(a,b): 
    return [int(v) for v in ((fmpz_poly(a)*fmpz_poly(b)).coeffs()+[0]*N6)[:N6]]
ch11=m1(vir(1,1,N6),pinv); ch41=m1(vir(4,1,N6),pinv); ch13=m1(vir(1,3,N6),pinv)
tot=[0]*P
for coef,sh,c in [(1,0,ch11),(1,9,ch41),(2,2,ch13)]:
    s=shift(lift(c,1),sh); tot=[x+coef*y for x,y in zip(tot,s)]
f=nahm(A6,[0]*6,[1]*6,N6,den)
print(" Sun-Wang (3.19): E6 sum == G(q) (chi_0 + chi_3 + 2 chi_{2/3})^{M(6,5)} (normalised) :",f==mul(lift(G,1),tot,P))
assert f==mul(lift(G,1),tot,P)
print("ALL CHECKS PASSED")
