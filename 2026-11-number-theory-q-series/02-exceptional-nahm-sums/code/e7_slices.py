from fractions import Fraction as Fr
from nahm_flint import nahm, qprod, mul
N=40
# T(a) = sum q^{N1^2+N2^2+N3^2 + a(N2+N3)}/..., variables (d,b,c): N1=d+b+c,N2=b+c,N3=c ; M=[[2,2,2],[2,4,4],[2,4,6]], B=a*(0,1,2)
M=[[2,2,2],[2,4,4],[2,4,6]]
eta=qprod({0},1,N)
qa=[1]+[0]*(N-1)
for a in range(0,7):
    T=nahm(M,[0,a,2*a],[1,1,1],N,1)
    if a>0:
        # divide by (1-q^a)
        inv=[0]*N
        for j in range(0,N,a): inv[j]=1
        qa=mul(qa,inv,N)
    g=mul(mul(T,qa,N),eta,N)
    print("a=%d  (q)_inf*T(a)/(q)_a ="%a, [(i,c) for i,c in enumerate(g) if c][:22])
