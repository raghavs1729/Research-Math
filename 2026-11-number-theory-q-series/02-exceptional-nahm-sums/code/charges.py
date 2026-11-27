"""Which linear charges (x^{w.n}) are preserved along the k=2 move chain KRTW(i)_2 -> Melzer_2 ?
Each species carries a symbolic charge = vector in Q^4 (coefficients on the initial species a1,a2,b1,b2).
Up-moves define the charge of new species; down-moves impose linear constraints."""
from fractions import Fraction as Fr
from classsearch import do
from moves import krtw
path=[('s', 0), ('f', 3, 4), ('S', 1, 5), ('p', 0), ('f', 0, 3), ('F', 0, 5, 4), ('f', 3, 4), ('M', 5, 6), ('f', 0, 4), ('M', 0, 4), ('S', 0, 5), ('F', 0, 1, 4), ('S', 0, 3), ('S', 1, 2)]
M,B,d,lab=krtw(2,(0,1))
W=[[Fr(int(i==j)) for j in range(4)] for i in range(4)]
cons=[]
add=lambda u,v:[x+y for x,y in zip(u,v)]
sub=lambda u,v:[x-y for x,y in zip(u,v)]
for mv in path:
    if mv[0]=='f': W.append(add(W[mv[1]],W[mv[2]]))
    elif mv[0]=='s': W.append([2*x for x in W[mv[1]]])
    elif mv[0]=='p': W.append(W[mv[1]][:])
    elif mv[0]=='F': a,b,c=mv[1:]; cons.append(sub(W[c],add(W[a],W[b]))); del W[c]
    elif mv[0]=='S': m,c=mv[1:]; cons.append(sub(W[c],[2*x for x in W[m]])); del W[c]
    elif mv[0]=='M': s,t=mv[1:]; cons.append(sub(W[s],W[t])); del W[t]
    M,B,d=do(M,B,d,mv)
print("constraints:"); [print("  ",[str(x) for x in c]) for c in cons if any(c)]
print("final species charges (in terms of a1,a2,b1,b2 weights):"); [print("  ",[str(x) for x in w]) for w in W]
# nullspace of constraints
import sympy
A=sympy.Matrix([[sympy.Rational(x.numerator,x.denominator) for x in c] for c in cons])
print("consistent charges (nullspace basis):",A.nullspace())
