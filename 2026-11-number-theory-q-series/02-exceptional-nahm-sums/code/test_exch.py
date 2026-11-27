from fractions import Fraction as Fr
from nahm_flint import nahm
N=40
# LHS_3: species (a1,a2,b): M=[[2,2,0],[2,4,2],[0,2,4]], B=(0,0,1), d=(1,1,2), x-weights (1,2,2)
L=nahm([[2,2,0],[2,4,2],[0,2,4]],[0,0,1],[1,1,2],N,1,weights=[1,2,2])
# RHS': species (n,Y): M=[[1,1],[1,3]], B=(1/2,1/2), x-weights (1,2)
R=nahm([[1,1],[1,3]],[Fr(1,2),Fr(1,2)],[1,1],N,1,weights=[1,2])
print("x-refined rank3 == rank2 :", L==R, sorted(L.keys())[:8])
if L!=R:
    for w in sorted(set(L)|set(R)):
        a=L.get(w,[0]*N); b=R.get(w,[0]*N)
        if a!=b: print(w,[x-y for x,y in zip(a,b)][:25]); break
L0=nahm([[2,2,0],[2,4,2],[0,2,4]],[0,0,1],[1,1,2],N,1); R0=nahm([[1,1],[1,3]],[Fr(1,2),Fr(1,2)],[1,1],N,1)
print("x=1:",L0==R0)
