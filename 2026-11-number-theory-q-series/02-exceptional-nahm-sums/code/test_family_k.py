"""Family generalising Theorem 4 (E7 is k=3):
   sum_{a,n_1..n_k} q^{(k/4)a^2 + a(N_2+...+N_k) + N_1^2+...+N_k^2}/((q)_a prod (q)_{n_i})
      = (1/(q)_inf) ( sum_{m>=0} q^{k m^2/4} + sum_{m>=1} q^{(k+4) m^2/4} )."""
import sys
from fractions import Fraction as Fr
from nahm_flint import nahm, qprod_inv, mul
N=int(sys.argv[1]) if len(sys.argv)>1 else 60
for k in range(1,7):
    r=k+1
    M=[[Fr(0)]*r for _ in range(r)]
    M[0][0]=Fr(k,2)
    for i in range(1,k+1):
        M[0][i]=M[i][0]=Fr(i-1)
        for j in range(1,k+1): M[i][j]=Fr(2*min(i,j))
    den=4; P=N*den
    f=nahm(M,[0]*r,[1]*r,N,den)
    th=[0]*P; m=0
    while k*m*m<P:
        th[k*m*m]+=1
        if m>0 and (k+4)*m*m<P: th[(k+4)*m*m]+=1
        m+=1
    pin=[0]*P
    for i,c in enumerate(qprod_inv({0},1,N)): pin[4*i]=c
    print("k=%d: identity holds mod q^%d: %s"%(k,N,f==mul(th,pin,P)))
