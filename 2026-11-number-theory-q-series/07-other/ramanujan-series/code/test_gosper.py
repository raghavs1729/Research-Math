import time
from wz import *
t=time.time()
for name,(params,T) in SEEDS.items():
    if name in ('Dougall5F4','Dougall7F6','Seed7','Dixon3F2','Watson3F2','Seed3','Seed9'): continue
    ok=[]
    for p in params:
        F = T.subs({p: p+'+n'})
        try:
            R = wz_mate(F)
            ok.append(check_wz(F,R))
        except NoMate as e:
            ok.append('NoMate:'+str(e))
    print(name, ok, round(time.time()-t,2), flush=True)
t=time.time()
F = seed_instance('Dougall7F6', dict(a=-1,b=-1,c=-1,d=-1,e=-1), 2, consts=dict(a='1/2',b='1/2',c='1/2',d='1/2',e='1/2'), k0=1)
R = wz_mate(F, verbose=True); print(check_wz(F,R), time.time()-t)
print(pari.factor(pari.denominator(R)))
num=pari.numerator(R); print(pari.poldegree(num,'k'), pari.poldegree(num,'n'))
print(pari.factor(num)[:, 0][0:3] if False else '')
print(num)
