from fractions import Fraction as Fr
from nahm_flint import nahm, qprod, qprod_inv, mul
from fusion import inv, cartanE
N=40; den=4
A7=inv(cartanE(7))
for row in A7: print([str(x) for x in row])
f=nahm(A7,[0]*7,[1]*7,N,den,count=True)
print([ (i,c) for i,c in enumerate(f) if c][:60])
# multiply by (q;q)_inf
eta=[0]*(N*den)
e1=qprod(set(range(1)),1,N)
for i,c in enumerate(e1): eta[den*i]=c
g=mul(f,eta,N*den)
print("f*(q)_inf:",[ (i,c) for i,c in enumerate(g) if c])
