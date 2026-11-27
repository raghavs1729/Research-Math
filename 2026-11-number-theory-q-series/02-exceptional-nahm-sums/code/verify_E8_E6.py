import sys, time
from fractions import Fraction as Fr
from nahm_flint import nahm, qprod_inv, mul
from fusion import inv, cartanE, cartanT, cartanA
N=int(sys.argv[1]) if len(sys.argv)>1 else 100
t0=time.time()
A8=inv(cartanE(8))
f=nahm(A8,[0]*8,[1]*8,N,1,count=True)
g=qprod_inv({1,2,3,4,7,8,9,10},11,N)
print("E8: sum_{n in N^8} q^{n C(E8)^-1 n/2}/(q)_n == 1/(q,q^2,q^3,q^4,q^7,q^8,q^9,q^10;q^11)_inf  mod q^%d :"%N, f==g, " time %.1fs"%(time.time()-t0)); sys.stdout.flush()
# AG multisum itself (rank 4) for cross-check
T4=[[2*x for x in row] for row in inv(cartanT(4))]
h=nahm(T4,[0]*4,[1]*4,N,1,count=True)
print("AG rank 4 multisum == product:", h==g, " time %.1fs"%(time.time()-t0)); sys.stdout.flush()
# E6: sum q^{n C(E6)^-1 n /2}/(q)_n = RR(q) * sum_{a,b} q^{(2/3)(a^2+ab+b^2)}/((q)_a (q)_b), in t=q^(1/3)
N6=int(sys.argv[2]) if len(sys.argv)>2 else N
A6=inv(cartanE(6))
f6=nahm(A6,[0]*6,[1]*6,N6,3,count=True)
A2=[[2*x for x in row] for row in inv(cartanA(2))]
z10=nahm(A2,[0,0],[1,1],N6,3)
rr=[0]*(3*N6)
r1=qprod_inv({1,4},5,N6)
for i,c in enumerate(r1): rr[3*i]=c
print("E6: sum == RR * Zagier-ex10(B=0)  mod q^%d :"%N6, f6==mul(rr,z10,3*N6), " time %.1fs"%(time.time()-t0))
