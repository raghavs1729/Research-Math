from fractions import Fraction as Fr
from nahm_flint import nahm, qprod, qprod_inv, mul
from moves import krtw
N=60
def prodser(num_res, m, den_list, N):
    # prod_{n mod m in num_res}(1-q^n) / prod over (res set, mod) in den_list
    c=qprod(num_res,m,N)
    for (rs,mm) in den_list: c=mul(c,qprod_inv(rs,mm,N),N)
    return c
for k in (1,2,3,4):
    # Melzer-type: M=min(i,j), B=(j/2)_j  (sum q^{sum N_j(N_j+1)/2})
    M=[[Fr(min(i,j)) for j in range(1,k+1)] for i in range(1,k+1)]
    B=[Fr(j,2) for j in range(1,k+1)]
    f=nahm(M,B,[1]*k,N,1)
    # (-q)_inf (q,q^{k+1},q^{k+2};q^{k+2})/(q)_inf
    g=prodser({1%(k+2),(k+1)%(k+2),0},k+2,[({1},2),({0},1)],N)
    Mk,Bk,dk,_=krtw(k,(0,1))
    h=nahm(Mk,Bk,dk,N,1) if k<=3 else None
    print(k,"Melzer+half-linear == (-q)(q,q^{k+1},q^{k+2};q^{k+2})/(q):",f==g," KRTW(i)==same:",(h==g) if h else None)
