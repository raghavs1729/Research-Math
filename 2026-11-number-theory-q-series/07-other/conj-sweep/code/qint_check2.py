# Independent (plain Python integers) check of the counterexamples in R01.
import itertools
def mul(p,q):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r
def P(As,b,r):
    p=[1]
    for a in As: p=mul(p,[1]*a)
    q=[0]*(r*(b-1)+1)
    for j in range(b): q[r*j]=1
    return mul(p,q)
def uni(p):
    i=0
    while i+1<len(p) and p[i]<=p[i+1]: i+=1
    while i+1<len(p) and p[i]>=p[i+1]: i+=1
    return i==len(p)-1
c=P((2,)*6,2,3); print('(1+q)^6(1+q^3) =',c,'unimodal:',uni(c))
cnt=0; lst=[]
for As in itertools.combinations_with_replacement(range(1,14),6):
    if any(a%3==0 for a in As): continue
    s=1+sum(a//3 for a in As)
    for b in range(s+1,21):
        if uni(P(As,b,3)): cnt+=1; lst.append((As,b))
print('r=3,k=6,a<=13,b<=20: unimodal although the condition fails:',cnt); print(lst[:12])
