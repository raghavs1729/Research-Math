# arXiv:2603.08510, Conjecture 1: pbar(16(56n+k)) = 0 mod 7, k in {11,43,51}  (pbar: overpartitions, f_2/f_1^2)
from radu import *
out='../cert/radu_overpartition_mod7.txt'; open(out,'w').write('arXiv:2603.08510 Conjecture 1\n\n')
r={1:-2,2:1}
m=896
print([ (16*k, P_set(m,r,16*k)) for k in (11,43,51)])
ok,info=certificate('pbar(896n+176) = 0 mod 7 (and the other residues of P)',2,r,m,176,7,out=out,claimed=[176,688,816])
