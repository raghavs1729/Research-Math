"""Fusion classes of generalized Nahm data (M, d): exponent (1/2) n^T M n, denominators (q^{d_i};q^{d_i})_{n_i}.
Fusion of species a,b (d_a=d_b=d): new species c with row_c=row_a+row_b, M_ab += d  (q -> q^d in the Durfee identity).
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations
from fusion import inv, cartanA, cartanT, cartanD, cartanE

def cartanB(r):
    C=cartanA(r); C[r-2][r-1]=-2; return C      # Bourbaki: a_{r-1,r}=-2 ; long roots 1..r-1, short root r
def cartanC(r):
    C=cartanA(r); C[r-1][r-2]=-2; return C      # a_{r,r-1}=-2 ; short roots 1..r-1, long root r
def cartanF4():
    return [[2,-1,0,0],[-1,2,-2,0],[0,-1,2,-1],[0,0,-1,2]]   # long,long,short,short
def cartanG2():
    return [[2,-1],[-3,2]]                      # short, long   (alpha1 short)
def dvec(name,r):
    if name=='B': return [2]*(r-1)+[1]
    if name=='C': return [1]*(r-1)+[2]
    if name=='F': return [2,2,1,1]
    if name=='G': return [1,3]
    return [1]*r
def cartan(name,r):
    return {'A':cartanA,'T':cartanT,'D':cartanD,'E':cartanE,'B':cartanB,'C':cartanC}[name](r) if name not in 'FG' else (cartanF4() if name=='F' else cartanG2())

def kron(X,Y):
    return [[X[i][j]*Y[k][l] for j in range(len(X)) for l in range(len(Y))] for i in range(len(X)) for k in range(len(Y))]
def pair(X,Y):
    """returns (M=A D symmetric exponent matrix, d) for A=C(X) (x) C(Y)^{-1}, D=D(X)(x)D(Y)"""
    (nx,rx),(ny,ry)=X,Y
    A=kron([[Fr(v) for v in row] for row in cartan(nx,rx)], inv(cartan(ny,ry)))
    d=[a*b for a in dvec(nx,rx) for b in dvec(ny,ry)]
    M=[[A[i][j]*d[j] for j in range(len(d))] for i in range(len(d))]
    assert all(M[i][j]==M[j][i] for i in range(len(d)) for j in range(len(d))), "AD not symmetric"
    return M,d

def defusions(M,d):
    r=len(M); out=[]
    for c in range(r):
        for a,b in combinations([i for i in range(r) if i!=c],2):
            if not (d[a]==d[b]==d[c]): continue
            g=M[a][b]-d[a]
            if M[c][c]!=M[a][a]+M[b][b]+2*g: continue
            if M[a][c]!=M[a][a]+g or M[b][c]!=M[b][b]+g: continue
            if all(M[c][k]==M[a][k]+M[b][k] for k in range(r) if k not in (a,b,c)):
                out.append((a,b,c))
    return out
def defuse(M,d,lab,a,b,c):
    r=len(M); keep=[i for i in range(r) if i!=c]
    B=[[M[i][j] for j in keep] for i in keep]
    ia,ib=keep.index(a),keep.index(b)
    B[ia][ib]-=d[a]; B[ib][ia]-=d[a]
    return B,[d[i] for i in keep],[lab[i] for i in keep]
def canon(M,d):
    """canonical form under simultaneous permutation (brute force for small rank, else sorted-invariant key)"""
    r=len(M)
    if r<=7:
        best=None
        idx=list(range(r))
        # restrict permutations: sort by (d, diagonal) blocks
        keyf=lambda i:(d[i],M[i][i])
        order=sorted(idx,key=keyf)
        groups=[]
        for i in order:
            if groups and keyf(groups[-1][0])==keyf(i): groups[-1].append(i)
            else: groups.append([i])
        def rec(gi,cur):
            nonlocal best
            if gi==len(groups):
                t=tuple(tuple(M[i][j] for j in cur) for i in cur)
                if best is None or t<best: best=t
                return
            for p in permutations(groups[gi]): rec(gi+1,cur+list(p))
        rec(0,[])
        return (tuple(sorted(d)),best)
    return (tuple(sorted(d)),tuple(sorted(tuple(sorted(row)) for row in M)))
def terminals(M,d,lab=None,seen=None,memo=None):
    if lab is None: lab=list(range(1,len(M)+1))
    if seen is None: seen={}
    if memo is None: memo=set()
    key=(tuple(tuple(r) for r in M),tuple(lab))
    if key in memo: return seen
    memo.add(key)
    ds=defusions(M,d)
    if not ds:
        seen.setdefault(canon(M,d),(M,d,lab)); return seen
    for (a,b,c) in ds:
        B,d2,l2=defuse(M,d,lab,a,b,c)
        terminals(B,d2,l2,seen,memo)
    return seen
def blocks(M):
    r=len(M); comp=list(range(r))
    def find(i):
        while comp[i]!=i: i=comp[i]
        return i
    for i in range(r):
        for j in range(r):
            if M[i][j]!=0: comp[find(i)]=find(j)
    g={}
    for i in range(r): g.setdefault(find(i),[]).append(i)
    return list(g.values())

if __name__=="__main__":
    import sys
    types=[('A',r) for r in range(1,9)]+[('T',r) for r in range(1,9)]+[('D',r) for r in range(4,9)]+[('E',r) for r in (6,7,8)]+[('B',r) for r in range(2,6)]+[('C',r) for r in range(2,6)]+[('F',4),('G',2)]
    maxrank=int(sys.argv[1]) if len(sys.argv)>1 else 8
    table={}
    for X in types:
        for Y in types:
            if X[1]*Y[1]>maxrank: continue
            M,d=pair(X,Y)
            T=terminals(M,d)
            for key,(M2,d2,lab) in T.items():
                table.setdefault(key,[]).append((X,Y))
            name="(%s%d,%s%d)"%(X[0],X[1],Y[0],Y[1])
            descr=[]
            for key,(M2,d2,lab) in T.items():
                bl=blocks(M2)
                descr.append("rank %d%s"%(len(M2)," blocks "+str([len(b) for b in bl]) if len(bl)>1 else ""))
            print(name,"rank",len(M),"->",descr)
    print("\n== coincidences of terminal forms ==")
    for key,ps in table.items():
        if len(ps)>1:
            print(["(%s%d,%s%d)"%(X[0],X[1],Y[0],Y[1]) for X,Y in ps], "rank",len(key[1]), "d",key[0])
