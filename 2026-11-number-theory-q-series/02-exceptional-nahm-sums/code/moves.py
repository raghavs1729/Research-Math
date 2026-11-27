"""Rank-changing exact transformations of generalized Nahm data (M, B, d):
    F(M,B,d) = sum_{n in N^r} q^{(1/2) n^T M n + B.n} / prod_i (q^{d_i};q^{d_i})_{n_i}.

FUSION (Durfee rectangle, base q^d):  1/((Q)_a (Q)_b) = sum_c Q^{(a-c)(b-c)}/((Q)_{a-c}(Q)_{b-c}(Q)_c), Q=q^d.
   species a,b (d_a=d_b=d)  ->  a,b,c  with row_c=row_a+row_b, B_c=B_a+B_b, M_ab += d.
SELF-FUSION (1/(x;Q)_inf = (-x;Q)_inf/(x^2;Q^2)_inf):
   1/(Q)_a = sum_{m+2c=a} Q^{m(m-1)/2}/((Q)_m (Q^2;Q^2)_c).
   species a (d_a=d, diag alpha, linear beta) -> m (d, diag alpha+d, linear beta-d/2), c (2d, diag 4 alpha, linear 2 beta),
   M_mc = 2 alpha, row_c = 2 row_a on the other species.
The functions below find the inverse (rank-lowering) moves.
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations

def defusions(M,B,d):
    r=len(M); out=[]
    for c in range(r):
        for a,b in combinations([i for i in range(r) if i!=c],2):
            if not (d[a]==d[b]==d[c]): continue
            g=M[a][b]-d[a]
            if B[c]!=B[a]+B[b]: continue
            if M[c][c]!=M[a][a]+M[b][b]+2*g: continue
            if M[a][c]!=M[a][a]+g or M[b][c]!=M[b][b]+g: continue
            if all(M[c][k]==M[a][k]+M[b][k] for k in range(r) if k not in (a,b,c)):
                out.append(('F',a,b,c))
    return out
def selfdefusions(M,B,d):
    r=len(M); out=[]
    for m in range(r):
        for c in range(r):
            if c==m or d[c]!=2*d[m]: continue
            al=M[m][m]-d[m]
            if M[c][c]!=4*al or M[m][c]!=2*al: continue
            if B[c]!=2*(B[m]+Fr(d[m],2)): continue
            if all(M[c][k]==2*M[m][k] for k in range(r) if k not in (m,c)):
                out.append(('S',m,c))
    return out
def apply(M,B,d,lab,mv):
    r=len(M)
    if mv[0]=='F':
        _,a,b,c=mv; keep=[i for i in range(r) if i!=c]
        M2=[[M[i][j] for j in keep] for i in keep]; ia,ib=keep.index(a),keep.index(b)
        M2[ia][ib]-=d[a]; M2[ib][ia]-=d[a]
        return M2,[B[i] for i in keep],[d[i] for i in keep],[lab[i] for i in keep]
    else:
        _,m,c=mv; keep=[i for i in range(r) if i!=c]
        M2=[[M[i][j] for j in keep] for i in keep]; im=keep.index(m)
        M2[im][im]-=d[m]
        B2=[B[i] for i in keep]; B2[im]+=Fr(d[m],2)
        return M2,B2,[d[i] for i in keep],[lab[i] for i in keep]
def canon(M,B,d):
    r=len(M); idx=list(range(r))
    keyf=lambda i:(d[i],M[i][i],B[i],tuple(sorted(M[i])))
    order=sorted(idx,key=keyf); groups=[]
    for i in order:
        if groups and keyf(groups[-1][0])==keyf(i): groups[-1].append(i)
        else: groups.append([i])
    best=[None]
    def rec(gi,cur):
        if gi==len(groups):
            t=tuple(tuple(M[i][j] for j in cur) for i in cur)
            if best[0] is None or t<best[0][0]: best[0]=(t,tuple(cur))
            return
        for p in permutations(groups[gi]): rec(gi+1,cur+list(p))
    rec(0,[])
    cur=best[0][1]
    return (tuple(d[i] for i in cur),tuple(B[i] for i in cur),best[0][0])
def terminals(M,B,d,lab=None,seen=None,memo=None,path=()):
    if lab is None: lab=[str(i+1) for i in range(len(M))]
    if seen is None: seen={}
    if memo is None: memo=set()
    key=(tuple(tuple(r) for r in M),tuple(B),tuple(d))
    if key in memo: return seen
    memo.add(key)
    mvs=defusions(M,B,d)+selfdefusions(M,B,d)
    if not mvs:
        seen.setdefault(canon(M,B,d),(M,B,d,lab,path)); return seen
    for mv in mvs:
        M2,B2,d2,l2=apply(M,B,d,lab,mv)
        desc=(mv[0],)+tuple(lab[i] for i in mv[1:])
        terminals(M2,B2,d2,l2,seen,memo,path+(desc,))
    return seen
def show(name,M,B,d,lab=None):
    T=terminals(M,B,d,lab)
    print("==",name,"rank",len(M))
    for key,(M2,B2,d2,l2,path) in T.items():
        print("  terminal rank",len(M2),"labels",l2,"d",d2,"B",[str(x) for x in B2])
        print("    path",path)
        for row in M2: print("     ",[str(x) for x in row])
    return T

def krtw(k,lin):
    """KRTW sums in variables a_j=r_j-r_{j+1} (base q), b_l=s_l-s_{l-1} (base q^2):
       exponent sum_i (r_i+s_i)^2+s_i^2 + lin_r*sum r_i + lin_s*sum s_i"""
    n=2*k; M=[[Fr(0)]*n for _ in range(n)]; B=[Fr(0)]*n
    # r_i = sum_{j>=i} a_j ; s_i = sum_{l<=i} b_l ; indices a_j -> j-1, b_l -> k+l-1
    import itertools
    def rvec(i): return [1 if (j>=i) else 0 for j in range(1,k+1)]+[0]*k
    def svec(i): return [0]*k+[1 if (l<=i) else 0 for l in range(1,k+1)]
    for i in range(1,k+1):
        u=[x+y for x,y in zip(rvec(i),svec(i))]; v=svec(i)
        for p in range(n):
            for q in range(n):
                M[p][q]+=2*(u[p]*u[q]+v[p]*v[q])
        for p in range(n):
            B[p]+=lin[0]*rvec(i)[p]+lin[1]*svec(i)[p]
    d=[1]*k+[2]*k
    lab=["a%d"%j for j in range(1,k+1)]+["b%d"%l for l in range(1,k+1)]
    return M,B,d,lab
if __name__=="__main__":
    for k in (1,2,3):
        for name,lin in (("(i) +s",(0,1)),("(ii) +r+2s",(1,2)),("(iii) 0",(0,0))):
            M,B,d,lab=krtw(k,lin)
            show("KRTW k=%d %s"%(k,name),M,B,d,lab)
