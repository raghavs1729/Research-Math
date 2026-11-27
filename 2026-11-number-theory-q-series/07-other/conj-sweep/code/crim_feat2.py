from crim_struct import *
from crim_feat import runs
import collections
def red(w, f):
    return tuple(f(l) for ch,l in runs(w))
maps = {
 'parity': lambda l: l%2,
 '1/even/odd>=3': lambda l: 1 if l==1 else (2 if l%2==0 else 3),
 'cap2': lambda l: min(l,2), 'cap3': lambda l: min(l,3),
 '1,2,odd>=3,even>=4': lambda l: l if l<=2 else (3 if l%2 else 4),
}
for L in (8,10,12,14,16):
    W=list(words(L)); 
    for name,f in maps.items():
        d={};ok=True
        for w in W:
            k=red(w,f); v=isP(w)
            if d.setdefault(k,v)!=v: ok=False;break
        print(L,name,ok, end=' | ')
    print()
# statistics by number of runs
for L in (10,12):
    c=collections.Counter(); t=collections.Counter()
    for w in words(L):
        k=len(runs(w)); t[k]+=1; c[k]+=isP(w)
    print(L,[(k,c[k],t[k]) for k in sorted(t)])
