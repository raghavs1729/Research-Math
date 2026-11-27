"""Fusion / defusion of Nahm data.

Fusion identity (Durfee rectangle / q-Chu-Vandermonde limit), valid for a,b >= 0:
    1/((q)_a (q)_b) = sum_{c=0}^{min(a,b)} q^{(a-c)(b-c)} / ((q)_{a-c} (q)_{b-c} (q)_c).
Hence a rank-r Nahm sum with matrix A (species a,b) equals the rank-(r+1) Nahm sum with a new
species c whose row is row_a+row_b, with A_ab increased by 1.  'Defusion' is the inverse step.
"""
from fractions import Fraction as Fr
from itertools import combinations

def inv(M):
    n=len(M); A=[[Fr(x) for x in row]+[Fr(int(i==j)) for j in range(n)] for i,row in enumerate(M)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        pv=A[c][c]; A[c]=[x/pv for x in A[c]]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]; A[r]=[x-f*y for x,y in zip(A[r],A[c])]
    return [row[n:] for row in A]

def cartanA(r): return [[2 if i==j else (-1 if abs(i-j)==1 else 0) for j in range(r)] for i in range(r)]
def cartanT(r):
    C=cartanA(r); C[r-1][r-1]=1; return C
def cartanD(r):
    C=cartanA(r); C[r-2][r-1]=C[r-1][r-2]=0; C[r-3][r-1]=C[r-1][r-3]=-1; return C
def cartanE(n):
    C=cartanA(n); C[n-2][n-1]=C[n-1][n-2]=0; C[n-4][n-1]=C[n-1][n-4]=-1; return C

def defusions(A):
    """all (a,b,c): species c is the fusion of a and b"""
    r=len(A); out=[]
    for c in range(r):
        for a,b in combinations([i for i in range(r) if i!=c],2):
            g=A[a][b]-1
            if A[c][c]!=A[a][a]+A[b][b]+2*g: continue
            if A[a][c]!=A[a][a]+g or A[b][c]!=A[b][b]+g: continue
            if all(A[c][k]==A[a][k]+A[b][k] for k in range(r) if k not in (a,b,c)):
                out.append((a,b,c))
    return out

def defuse(A,labels,a,b,c):
    r=len(A); keep=[i for i in range(r) if i!=c]
    B=[[A[i][j] for j in keep] for i in keep]
    ia,ib=keep.index(a),keep.index(b)
    B[ia][ib]-=1; B[ib][ia]-=1
    return B,[labels[i] for i in keep]

def reduce_all(A,labels=None,seen=None,path=()):
    """DFS over all defusion sequences; returns dict of terminal forms -> path"""
    if labels is None: labels=list(range(1,len(A)+1))
    if seen is None: seen={}
    ds=defusions(A)
    key=(tuple(tuple(r) for r in A),tuple(labels))
    if not ds:
        seen.setdefault(key,path); return seen
    for (a,b,c) in ds:
        B,l2=defuse(A,labels,a,b,c)
        reduce_all(B,l2,seen,path+((labels[a],labels[b],labels[c]),))
    return seen

def show(name,A):
    res=reduce_all(A)
    print("==",name,"rank",len(A))
    forms={}
    for (M,lab),path in res.items():
        forms.setdefault(tuple(sorted(map(tuple,M))),(M,lab,path))
    for k,(M,lab,path) in forms.items():
        print(" terminal rank",len(M),"labels",lab,"path",path)
        for row in M: print("   ",[str(x) for x in row])
        try:
            print("    inverse:",[[str(x) for x in row] for row in inv(M)])
        except StopIteration: print("    singular")

if __name__=="__main__":
    for n in (6,7,8): show("E%d inverse"%n, inv(cartanE(n)))
    for r in (3,4,5,6): show("D%d inverse"%r, inv(cartanD(r)))
    for r in (2,3,4,5): show("A%d inverse"%r, inv(cartanA(r)))
    for r in (2,3,4): show("T%d inverse"%r, inv(cartanT(r)))
