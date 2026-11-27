"""x-refined conjecture (E_k) and the full x-refined KRTW(i)_k = Melzer_k."""
import sys
from fractions import Fraction as Fr
from nahm_flint import nahm
from moves import krtw
def run(k,N):
    # full: KRTW(i)_k with x^{sum (r_i+2 s_i)}  vs  Melzer_k with x^{sum N_i}
    M,B,d,lab=krtw(k,(0,1))
    # weights: a_j contributes to r_i for i<=j : j ; b_l contributes 2 to s_i for i>=l : 2(k+1-l)
    w=[j for j in range(1,k+1)]+[2*(k+1-l) for l in range(1,k+1)]
    L=nahm(M,B,d,N,1,weights=w)
    Mm=[[Fr(min(i,j)) for j in range(1,k+1)] for i in range(1,k+1)]
    R=nahm(Mm,[Fr(j,2) for j in range(1,k+1)],[1]*k,N,1,weights=list(range(1,k+1)))
    print("k=%d  x-refined KRTW(i) == Melzer (x^{sum N_i}):"%k, L==R)
    # (E_k): drop b_1 (s_1=0)
    idx=[i for i in range(2*k) if i!=k]
    M2=[[M[i][j] for j in idx] for i in idx]; B2=[B[i] for i in idx]; d2=[d[i] for i in idx]; w2=[w[i] for i in idx]
    L2=nahm(M2,B2,d2,N,1,weights=w2)
    # RHS: Melzer with last level full: exponent sum_i N_i(N_i+1)/2 + N_k(N_k-1)/2 ; N_k = n_k
    M3=[row[:] for row in Mm]; M3[k-1][k-1]+=1; B3=[Fr(j,2) for j in range(1,k+1)]; B3[k-1]-=Fr(1,2)
    R2=nahm(M3,B3,[1]*k,N,1,weights=list(range(1,k+1)))
    print("k=%d  (E_k):"%k, L2==R2)
for k,N in ((2,40),(3,36),(4,24)): run(k,N)
