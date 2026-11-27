"""BFS in the equivalence class of generalized Nahm data under exact moves:
 FUSE / DEFUSE, SELF-FUSE / SELF-DEFUSE, SPLIT / MERGE  (see moves.py).  SPLIT: 1/(x;Q)_inf = 1/((x;Q^2)_inf (xQ;Q^2)_inf):
   species n (d) -> s,t (2d,2d) with identical rows (M_ss=M_tt=M_st=M_nn), B_s=B_n, B_t=B_n+d.
"""
import sys
from fractions import Fraction as Fr
from itertools import combinations
from moves import defusions, selfdefusions, apply, canon, krtw

def merges(M,B,d):
    r=len(M); out=[]
    for s in range(r):
        for t in range(r):
            if s==t or d[s]!=d[t] or d[s]%2: continue
            if B[t]!=B[s]+Fr(d[s],2): continue
            if not (M[s][s]==M[t][t]==M[s][t]): continue
            if all(M[s][k]==M[t][k] for k in range(r) if k not in (s,t)):
                out.append(('M',s,t))
    return out
def up_moves(M,B,d,maxd):
    r=len(M); out=[]
    for a,b in combinations(range(r),2):
        if d[a]==d[b]: out.append(('f',a,b))
    for a in range(r):
        if 2*d[a]<=maxd: out.append(('s',a)); out.append(('p',a))
    return out
def do(M,B,d,mv):
    r=len(M)
    M=[row[:] for row in M]; B=B[:]; d=d[:]
    if mv[0] in 'FS':
        M2,B2,d2,_=apply(M,B,d,list(range(r)),mv); return M2,B2,d2
    if mv[0]=='M':
        _,s,t=mv; keep=[i for i in range(r) if i!=t]
        M2=[[M[i][j] for j in keep] for i in keep]; d2=[d[i] for i in keep]; B2=[B[i] for i in keep]
        d2[keep.index(s)]=d[s]//2
        return M2,B2,d2
    if mv[0]=='f':
        _,a,b=mv
        row=[M[a][k]+M[b][k] for k in range(r)]
        g=M[a][b]
        cc=M[a][a]+M[b][b]+2*g
        M[a][b]+=d[a]; M[b][a]+=d[a]
        row[a]=M[a][a]+g; row[b]=M[b][b]+g
        for k in range(r): M[k].append(row[k])
        M.append(row+[cc]); B.append(B[a]+B[b]); d.append(d[a]); return M,B,d
    if mv[0]=='s':
        _,a=mv; al=M[a][a]
        row=[2*M[a][k] for k in range(r)]; row[a]=2*al
        M[a][a]=al+d[a]
        for k in range(r): M[k].append(row[k])
        M.append(row+[4*al]); B.append(2*B[a]); B[a]-=Fr(d[a],2); d.append(2*d[a]); return M,B,d
    if mv[0]=='p':
        _,a=mv
        row=[M[a][k] for k in range(r)]
        for k in range(r): M[k].append(row[k])
        M.append(row+[M[a][a]]); B.append(B[a]+d[a]); d[a]*=2; d.append(d[a]); return M,B,d
def explore(M,B,d,maxrank,maxd=4,limit=200000):
    start=canon(M,B,d); seen={start:(None,None)}; frontier=[(M,B,d,start)]
    while frontier and len(seen)<limit:
        nxt=[]
        for (M,B,d,key) in frontier:
            mvs=defusions(M,B,d)+selfdefusions(M,B,d)+merges(M,B,d)
            if len(M)<maxrank: mvs+=up_moves(M,B,d,maxd)
            for mv in mvs:
                M2,B2,d2=do(M,B,d,mv)
                k2=canon(M2,B2,d2)
                if k2 not in seen:
                    seen[k2]=(key,mv); nxt.append((M2,B2,d2,k2))
        frontier=nxt
    return seen
def pathto(seen,key):
    p=[]
    while seen[key][0] is not None:
        p.append(seen[key][1]); key=seen[key][0]
    return p[::-1]
if __name__=="__main__":
    k=int(sys.argv[1]); R=int(sys.argv[2])
    # target class: start from KRTW (i)
    for name,lin in (("(i)",(0,1)),("(ii)",(1,2)),("(iii)",(0,0))):
        M,B,d,lab=krtw(k,lin)
        seen=explore(M,B,d,R)
        print("KRTW k=%d %s: class size (rank<=%d): %d"%(k,name,R,len(seen)))
        byrank={}
        for key in seen: byrank.setdefault(len(key[0]),[]).append(key)
        for r in sorted(byrank)[:2]:
            for key in byrank[r][:6]:
                print("  rank",r,"d",key[0],"B",[str(x) for x in key[1]],"M",[[str(x) for x in row] for row in key[2]], "path",pathto(seen,key))
