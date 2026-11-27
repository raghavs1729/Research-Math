import time, sys
from wz import *
import mpmath as mp
mp.mp.dps = 30
# 1. numerically test each seed sums to 1 at random params (where convergent)
import random
random.seed(1)
tests = {
 'Gauss2F1': dict(a=0.3,b=0.4,c=1.9),
 'Dixon3F2': dict(a=0.9,b=0.3,c=0.2),
 'Watson3F2': dict(a=0.3,b=0.2,c=1.4),
 'Dougall5F4': dict(a=1.9,b=0.3,c=0.2,d=0.4),
 'Seed1': dict(a=0.3,c=1.9),
 'Seed3': dict(a=1.9,b=0.3,d=0.4),
 'Seed4': dict(a=1.9,b=0.3),
}
for name, vals in tests.items():
    params, T = SEEDS[name]
    v = {p: mp.mpf(x) for p, x in vals.items()}
    def f(kk):
        vv = dict(v); vv['k'] = mp.mpf(kk); vv['n']=mp.mpf(0)
        return T.evalf(vv, mp)
    s = mp.nsum(f, [0, mp.inf])
    print(name, s)
# terminating ones: Dougall7F6 with e=-N... balanced: e = -N? check Balanced3F2 with c=-N
for name, vals in {'Balanced3F2': dict(a=0.3,b=0.4,c=-5,d=1.7+0.01), 'Dougall7F6': dict(a=1.9,b=0.3,c=0.2,d=0.4,e=-6),
                   'Seed2': dict(a=0.3,c=-5,d=1.7), 'Seed7': dict(a=1.9,b=0.3,d=0.4,e=-6), 'Seed9': None, 'Seed10': None}.items():
    if vals is None: continue
    params, T = SEEDS[name]
    eps = mp.mpf(10)**(-12)
    v = {p: mp.mpf(x) for p, x in vals.items()}
    # perturb the terminating parameter to avoid 1/Gamma(-N) = 0 * inf
    tp = 'c' if name in ('Balanced3F2','Seed2') else 'e'
    N = -int(vals[tp])
    v[tp] = v[tp] + eps
    s = 0
    for kk in range(0, N+1):
        vv = dict(v); vv['k'] = mp.mpf(kk); vv['n']=mp.mpf(0)
        s += T.evalf(vv, mp)
    print(name, s)
