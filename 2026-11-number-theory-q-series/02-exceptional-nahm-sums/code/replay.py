import sys
from fractions import Fraction as Fr
from classsearch import do
from moves import krtw
def lab_do(lab,d,mv,cnt):
    r=len(lab); lab=lab[:]
    if mv[0]=='F': a,b,c=mv[1:]; print("   DEFUSE  %s = %s + %s"%(lab[c],lab[a],lab[b])); del lab[c]
    elif mv[0]=='S': m,c=mv[1:]; print("   SELF-DEFUSE  (%s,%s) -> %s"%(lab[m],lab[c],lab[m]+"'")); lab[m]=lab[m]+"'"; del lab[c]
    elif mv[0]=='M': s,t=mv[1:]; print("   MERGE (%s,%s)"%(lab[s],lab[t])); lab[s]="<%s|%s>"%(lab[s],lab[t]); del lab[t]
    elif mv[0]=='f': a,b=mv[1:]; print("   FUSE new = %s + %s"%(lab[a],lab[b])); lab.append("(%s+%s)"%(lab[a],lab[b]))
    elif mv[0]=='s': a=mv[1]; print("   SELF-FUSE %s"%lab[a]); lab.append("2[%s]"%lab[a]); lab[a]="m[%s]"%lab[a]
    elif mv[0]=='p': a=mv[1]; print("   SPLIT %s"%lab[a]); lab.append("%s~"%lab[a]); 
    return lab
def replay(M,B,d,lab,path):
    print("start:"); show(M,B,d,lab)
    for mv in path:
        lab=lab_do(lab,d,mv,0)
        M,B,d=do(M,B,d,mv)
        show(M,B,d,lab)
def show(M,B,d,lab):
    for i in range(len(M)):
        print("      %-28s d=%d B=%-5s"%(lab[i],d[i],str(B[i])),[str(x) for x in M[i]])
if __name__=="__main__":
    path=[('s', 0), ('f', 3, 4), ('S', 1, 5), ('p', 0), ('f', 0, 3), ('F', 0, 5, 4), ('f', 3, 4), ('M', 5, 6), ('f', 0, 4), ('M', 0, 4), ('S', 0, 5), ('F', 0, 1, 4), ('S', 0, 3), ('S', 1, 2)]
    M,B,d,lab=krtw(2,(0,1))
    replay(M,B,d,lab,path)
