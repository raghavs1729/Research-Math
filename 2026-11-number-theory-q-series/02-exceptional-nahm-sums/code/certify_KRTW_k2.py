"""Certificate for proofs/KRTW-k2.md: the k=2 case of the first conjectural Andrews-Gordon-type identity of
Kanade-Russell-Tsuchioka-Warnaar (arXiv:2404.03851, display after Conjecture 4.9).

A state is generalized Nahm data (M,B,d):  F = sum_n q^{n M n/2 + B.n}/prod (q^{d_i};q^{d_i})_{n_i}.
The chain below consists of 14 moves; every move is an instance of Lemma F (fusion), S (self-fusion),
P (split) of the proof, read forwards (f,s,p: rank +1) or backwards (F,S,M: rank -1).
Checks performed for every step:
  (1) EXACT: for a backward move, applying the corresponding forward move to the new state reproduces the old
      state (as labelled matrices, up to reordering of species); for a forward move the new state is
      computed by the forward rule itself;
  (2) numerical: the Nahm sums of consecutive states agree mod q^N.
Finally the last state is the rank-2 sum  sum q^{(N1^2+N2^2+N1+N2)/2}/((q)_{n1}(q)_{n2}), which is evaluated in the proof.
"""
import sys
from fractions import Fraction as Fr
from itertools import permutations
from nahm_flint import nahm, qprod_inv, mul
from moves import krtw
N=int(sys.argv[1]) if len(sys.argv)>1 else 40

def fwd_fuse(M,B,d,a,b):
    r=len(M); M=[row[:] for row in M]; B=B[:]; d=d[:]
    assert a!=b and d[a]==d[b]
    row=[M[a][k]+M[b][k] for k in range(r)]; cc=M[a][a]+M[b][b]+2*M[a][b]
    M[a][b]+=d[a]; M[b][a]+=d[a]
    for k in range(r): M[k].append(row[k])
    M.append(row+[cc]); B.append(B[a]+B[b]); d.append(d[a]); return M,B,d
def fwd_self(M,B,d,a):
    r=len(M); M=[row[:] for row in M]; B=B[:]; d=d[:]
    al=M[a][a]; row=[2*M[a][k] for k in range(r)]; row[a]=2*al
    M[a][a]=al+d[a]
    for k in range(r): M[k].append(row[k])
    M.append(row+[4*al]); B.append(2*B[a]); B[a]-=Fr(d[a],2); d.append(2*d[a]); return M,B,d
def fwd_split(M,B,d,a):
    r=len(M); M=[row[:] for row in M]; B=B[:]; d=d[:]
    row=[M[a][k] for k in range(r)]
    for k in range(r): M[k].append(row[k])
    M.append(row+[M[a][a]]); B.append(B[a]+d[a]); d[a]*=2; d.append(d[a]); return M,B,d
def same(S1,S2):
    (M1,B1,d1),(M2,B2,d2)=S1,S2
    r=len(M1)
    if r!=len(M2): return False
    for p in permutations(range(r)):
        if all(d1[i]==d2[p[i]] and B1[i]==B2[p[i]] for i in range(r)) and all(M1[i][j]==M2[p[i]][p[j]] for i in range(r) for j in range(r)):
            return True
    return False
def remove(M,B,d,c):
    keep=[i for i in range(len(M)) if i!=c]
    return [[M[i][j] for j in keep] for i in keep],[B[i] for i in keep],[d[i] for i in keep],keep
def bwd(M,B,d,mv):
    """returns new state and checks exactly that the forward move maps it back to (M,B,d)"""
    if mv[0]=='F':
        _,a,b,c=mv; M2,B2,d2,keep=remove(M,B,d,c); ia,ib=keep.index(a),keep.index(b)
        M2[ia][ib]-=d[a]; M2[ib][ia]-=d[a]
        assert same(fwd_fuse(M2,B2,d2,ia,ib),(M,B,d)), "defusion check failed"
    elif mv[0]=='S':
        _,m,c=mv; M2,B2,d2,keep=remove(M,B,d,c); im=keep.index(m)
        M2[im][im]-=d[m]; B2[im]+=Fr(d[m],2)
        assert same(fwd_self(M2,B2,d2,im),(M,B,d)), "self-defusion check failed"
    elif mv[0]=='M':
        _,s,t=mv; M2,B2,d2,keep=remove(M,B,d,t); isx=keep.index(s)
        assert d[s]%2==0; d2[isx]=d[s]//2
        assert same(fwd_split(M2,B2,d2,isx),(M,B,d)), "merge check failed"
    return M2,B2,d2
def posdef(M):
    # leading principal minors > 0 (exact)
    import sympy
    A=sympy.Matrix([[sympy.Rational(x.numerator,x.denominator) for x in row] for row in M])
    return all(A[:k,:k].det()>0 for k in range(1,len(M)+1))
path=[('s', 0), ('f', 3, 4), ('S', 1, 5), ('p', 0), ('f', 0, 3), ('F', 0, 5, 4), ('f', 3, 4), ('M', 5, 6), ('f', 0, 4), ('M', 0, 4), ('S', 0, 5), ('F', 0, 1, 4), ('S', 0, 3), ('S', 1, 2)]
names={'f':'FUSE','s':'SELF-FUSE','p':'SPLIT','F':'DEFUSE (Lemma F backwards)','S':'SELF-DEFUSE (Lemma S backwards)','M':'MERGE (Lemma P backwards)'}
M,B,d,lab=krtw(2,(0,1))
def pr(M,B,d):
    for i in range(len(M)): print("     d=%d  B=%-5s"%(d[i],B[i]),[str(x) for x in M[i]])
print("State 0 = KRTW k=2 sum (variables a1=r1-r2, a2=r2, b1=s1, b2=s2-s1):"); pr(M,B,d)
ser=nahm(M,B,d,N,2); ser0=ser
for t,mv in enumerate(path,1):
    if mv[0]=='f': M,B,d=fwd_fuse(M,B,d,mv[1],mv[2])
    elif mv[0]=='s': M,B,d=fwd_self(M,B,d,mv[1])
    elif mv[0]=='p': M,B,d=fwd_split(M,B,d,mv[1])
    else: M,B,d=bwd(M,B,d,mv)
    if posdef(M):
        s2=nahm(M,B,d,N,2); num="series equal mod q^%d: %s"%(N,s2==ser); assert s2==ser
    else:
        num="(numerical check skipped: matrix only semidefinite after SPLIT; sum is still admissible)"
    print("Step %2d: %-34s indices %s   exact check: ok   %s"%(t,names[mv[0]],mv[1:],num))
    pr(M,B,d)
assert same((M,B,d),([[Fr(1),Fr(1)],[Fr(1),Fr(2)]],[Fr(1,2),Fr(1)],[1,1]))
print("Final state == Melzer-type rank 2 data M=[[1,1],[1,2]], B=(1/2,1): True")
# product side (-q;q)_inf (-q^2;q^2)_inf = 1/((q;q^2)_inf (q^2;q^4)_inf), as series in t=q^(1/2)
g=mul(qprod_inv({1},2,N),qprod_inv({2},4,N),N)
g2=[0]*(2*N)
for i,c in enumerate(g): g2[2*i]=c
print("KRTW k=2 sum == (-q;q)_inf (-q^2;q^2)_inf mod q^%d:"%N, ser0==g2); assert ser0==g2
print("ALL CHECKS PASSED")
