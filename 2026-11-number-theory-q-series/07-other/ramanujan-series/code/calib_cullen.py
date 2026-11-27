import time, mpmath as mp
from wz import *
from numerics import *
mp.mp.dps = 40
F = seed_instance('Dougall7F6', dict(a=-1,b=-1,c=-1,d=-1,e=-1), 2, consts=dict(a='1/2',b='1/2',c='1/2',d='1/2',e='1/2'), k0=1)
R = wz_mate(F); assert check_wz(F,R)
W = WZNum(F,R)
a = mp.mpf('0.13')
SF = W.sumF(a)
SG = W.sumG(a)
# g(n) = lim_k G(n,k)
def g(n, K=4000):
    return W.Gv(n, K)
print('sumF', SF); print('sumG', SG)
for K in (2000,4000,8000):
    sg = mp.nsum(lambda n: W.Gv(n+a, K), [0, mp.inf]); print(K, 'sum g ~', sg, 'residual', SF + sg - SG)
