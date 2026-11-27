# arXiv:2503.06228 (Nath, Sarma), Conjecture 1: ped(225n+t) = 0 mod 192, t in {43, 88, 133, 223}; ped: f_4/f_1.
from radu import *
out='../cert/radu_ped_mod192.txt'; open(out,'w').write('arXiv:2503.06228 Conjecture 1\n\n')
r={1:-1,4:1}
done=set()
for t in (43,88,133,223):
    if t in done: continue
    ok,info=certificate(f'ped(225n+{t}) = 0 mod 192',4,r,225,t,192,out=out,claimed=[43,88,133,223])
    done|=set(info[3])
