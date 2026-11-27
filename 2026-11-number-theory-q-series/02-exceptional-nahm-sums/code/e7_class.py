import sys
from fractions import Fraction as Fr
from classsearch import explore, pathto
F=Fr
M=[[F(3,2),F(1),F(2),F(0)],[F(1),F(4),F(4),F(2)],[F(2),F(4),F(6),F(2)],[F(0),F(2),F(2),F(2)]]
B=[F(0)]*4; d=[1]*4
R=int(sys.argv[1]); maxd=int(sys.argv[2])
seen=explore(M,B,d,R,maxd=maxd,limit=250000)
byrank={}
for key in seen: byrank.setdefault(len(key[0]),[]).append(key)
print("class size",len(seen),{r:len(v) for r,v in sorted(byrank.items())})
for r in sorted(byrank):
    if r<=4:
        for key in byrank[r][:12]:
            print(" rank",r,"d",key[0],"B",[str(x) for x in key[1]],"M",[[str(x) for x in row] for row in key[2]])
