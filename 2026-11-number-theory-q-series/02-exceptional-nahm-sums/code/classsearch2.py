import sys
from fractions import Fraction as Fr
from classsearch import explore, pathto
from moves import krtw
k=int(sys.argv[1]); R=int(sys.argv[2]); maxd=int(sys.argv[3]) if len(sys.argv)>3 else 2
for name,lin in (("(i)",(0,1)),("(ii)",(1,2)),("(iii)",(0,0))):
    M,B,d,lab=krtw(k,lin)
    seen=explore(M,B,d,R,maxd=maxd,limit=400000)
    byrank={}
    for key in seen: byrank.setdefault(len(key[0]),[]).append(key)
    print("KRTW k=%d %s: class size (rank<=%d): %d"%(k,name,R,len(seen)), {r:len(v) for r,v in sorted(byrank.items())})
    for r in sorted(byrank):
        if r<2*k:
            for key in byrank[r][:10]:
                print("  rank",r,"d",key[0],"B",[str(x) for x in key[1]],"M",[[str(x) for x in row] for row in key[2]], "path",pathto(seen,key))
