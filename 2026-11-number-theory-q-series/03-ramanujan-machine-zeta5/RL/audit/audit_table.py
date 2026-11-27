# Exact (p + q R2) arithmetic via the contiguous relations quoted in THEOREM.md s.3; integer table k=0..8.
from fractions import Fraction as F
from math import gcd
from mpmath import mp, mpf, sqrt, pi, e, erf, hyp1f1, nstr
mp.dps=50; R2=sqrt(pi*e/2)*erf(1/sqrt(2))
mu={0:(F(1),F(0)),1:(F(0),F(1))}   # value = p + q*R2
for a in range(1,12):
    # a mu_{a+1} = (3/2-a) mu_{a-1} + (2a-1) mu_a
    p=(F(3,2)-a)*mu[a-1][0]+(2*a-1)*mu[a][0]; q=(F(3,2)-a)*mu[a-1][1]+(2*a-1)*mu[a][1]
    mu[a+1]=(p/a,q/a)
M5={a:(3*(mu[a][0]-mu[a-1][0]),3*(mu[a][1]-mu[a-1][1])) for a in range(1,12)}
M5[0]=(F(1),F(0))
M7={a:(5*(M5[a][0]-M5[a-1][0]),5*(M5[a][1]-M5[a-1][1])) for a in range(1,12)}
val=lambda t: t[0]+t[1]*R2
for a in range(0,6): assert abs(val(mu[a])-hyp1f1(a,1.5,0.5))<mpf(10)**-45 and abs(val(M5[a])-hyp1f1(a,2.5,0.5))<mpf(10)**-45
for a in range(1,6): assert abs(val(M7[a])-hyp1f1(a,3.5,0.5))<mpf(10)**-45
for k in range(0,9):
    # C = 1 + (k+1)/5 * M7[k+2]/M5[k+1] = (num)/(den), num = 5*M5[k+1] + (k+1)*M7[k+2], den = 5*M5[k+1]
    n=(5*M5[k+1][0]+(k+1)*M7[k+2][0], 5*M5[k+1][1]+(k+1)*M7[k+2][1]); d=(5*M5[k+1][0],5*M5[k+1][1])
    import math
    L=1
    for x in n+d: L=L*x.denominator//gcd(L,x.denominator)
    ints=[int(x*L) for x in (n[1],n[0],d[1],d[0])]
    g=0
    for x in ints: g=gcd(g,x)
    al,be,ga,de=[x//g for x in ints]
    if ga<0 or (ga==0 and de<0): al,be,ga,de=-al,-be,-ga,-de
    print(f"k={k}: ({al} R2 + {be})/({ga} R2 + {de})  check vs formula: {nstr((al*R2+be)/(ga*R2+de) - (1+(k+1)*hyp1f1(k+2,3.5,0.5)/(5*hyp1f1(k+1,2.5,0.5))),3)}")
