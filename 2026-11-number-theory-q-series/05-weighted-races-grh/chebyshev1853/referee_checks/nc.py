import numpy as np, math, sys
P=int(6.2e7)
bs=np.ones(P+1,dtype=bool); bs[:2]=False
for i in range(2,int(P**0.5)+1):
    if bs[i]: bs[i*i::i]=False
pr=np.nonzero(bs)[0].astype(np.float64)
def C(x,q=4):
    p=pr[pr<=60*x]
    if q==4: p=p[p!=2]; chi=np.where(p%4==1,1.0,-1.0)
    else: p=p[p!=3]; chi=np.where(p%3==1,1.0,-1.0)
    e=np.exp(-p/x)
    c=math.fsum(chi*e); w=-math.fsum(chi*np.log(p)*e)
    return c,w
for q in (4,3):
  xs=[13,13.43,14,68.5,68.86,69,69.79,70,71,100,1e3,6309.6,1e4,1e5,1e6] if q==4 else [13,13.43,13.53,14,100,2089.3,1e3,1e6]
  for x in xs:
    c,w=C(x,q)
    print(q,x,"C=%.10f"%c,"N_C=%.6f"%(-c*math.log(x)/math.sqrt(x)),"W/sqrt=%.6f"%(w/math.sqrt(x)))
# fine scan near crossing and peak
import numpy as np
def nc(x): c,_=C(x); return -c*math.log(x)/math.sqrt(x)
lo,hi=68.0,70.0
for _ in range(40):
    m=(lo+hi)/2
    if nc(m)<0.88: lo=m
    else: hi=m
print("crossing q4",lo)
best=max((nc(x),x) for x in np.linspace(5800,6900,111)); print("peak",best)
