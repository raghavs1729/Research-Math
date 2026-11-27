from fractions import Fraction as Fr
from classsearch import explore, pathto
from moves import canon
from replay import replay
M=[[Fr(2),Fr(2),Fr(0)],[Fr(2),Fr(4),Fr(2)],[Fr(0),Fr(2),Fr(4)]]; B=[Fr(0),Fr(0),Fr(1)]; d=[1,1,2]
tgt=canon([[Fr(1),Fr(1)],[Fr(1),Fr(3)]],[Fr(1,2),Fr(1,2)],[1,1])
for R in (4,5,6):
    seen=explore([r[:] for r in M],B[:],d[:],R,maxd=2,limit=300000)
    print("maxrank",R,"class size",len(seen),"target reached:",tgt in seen)
    if tgt in seen:
        p=pathto(seen,tgt); print(len(p),p)
        replay([r[:] for r in M],B[:],d[:],['a1','a2','b'],p)
        break
