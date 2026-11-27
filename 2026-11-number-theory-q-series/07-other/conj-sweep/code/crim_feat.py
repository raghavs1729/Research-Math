from crim_struct import *
import itertools
def runs(w):
    r=[];
    for ch in w:
        if r and r[-1][0]==ch: r[-1][1]+=1
        else: r.append([ch,1])
    return r
feats = {
 'pre2': lambda w: w[:2], 'suf2': lambda w: w[-2:], 'pre3': lambda w: w[:3], 'suf3': lambda w: w[-3:],
 'pE': lambda w: w.count('E')%2, 'nruns2': lambda w: len(runs(w))%2,
 'run1par': lambda w: runs(w)[0][1]%2, 'runlastpar': lambda w: runs(w)[-1][1]%2,
 'run1is1': lambda w: runs(w)[0][1]==1, 'runlastis1': lambda w: runs(w)[-1][1]==1,
 'oddruns': lambda w: sum(1 for ch,l in runs(w) if l%2)%2,
 'allrunpar': lambda w: tuple(l%2 for ch,l in runs(w)) if len(runs(w))<=4 else None,
 'run2par': lambda w: runs(w)[1][1]%2, 'runm2par': lambda w: runs(w)[-2][1]%2,
 'nruns_small': lambda w: min(len(runs(w)),6),
 'Eoddruns': lambda w: sum(1 for ch,l in runs(w) if l%2 and ch=='E')%2,
 'evenpos': lambda w: sum(1 for i,ch in enumerate(w) if ch=='E' and i%2==0)%2,
}
names=list(feats)
for L in (10,12,14):
    W=list(words(L)); val={w:isP(w) for w in W}
    found=[]
    for k in range(1,5):
        for comb in itertools.combinations(names,k):
            d={};ok=True
            for w in W:
                key=tuple(feats[n](w) for n in comb)
                if d.setdefault(key,val[w])!=val[w]: ok=False;break
            if ok: found.append(comb)
        if found: break
    print(L,found[:6])
