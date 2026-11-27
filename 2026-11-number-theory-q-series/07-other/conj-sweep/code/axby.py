# Chen-Zhu arXiv:2506.03620 Conjectures 1,2: 13/66*pi(s) <= pi(a,b) <= 1/2*pi(s), s=ab-a-b, 3<=a<b coprime.
import sys, numpy as np
from math import gcd
def run(amax, bmax_of):
    Smax = max((a*bmax_of(a)) for a in range(3, amax+1)) + 10
    sv = np.ones(Smax+1, dtype=bool); sv[:2] = False
    for i in range(2, int(Smax**.5)+1):
        if sv[i]: sv[i*i::i] = False
    pc = np.cumsum(sv)
    lo = (1, 1, None); hi = (0, 1, None); eq_lo = []; eq_hi = []; viol = []
    cnt = 0
    for a in range(3, amax+1):
        for b in range(a+1, bmax_of(a)+1):
            if gcd(a, b) != 1: continue
            s = a*b - a - b
            # n in T iff n >= b*y with n = b*y mod a for the unique y in [0,a-1]; n<=s not in T iff s-n in T
            inT = np.zeros(s+1, dtype=bool)
            for y in range(a):
                if b*y > s: break
                inT[b*y::a] = True
            p = int(np.count_nonzero(inT & sv[:s+1])); P = int(pc[s]); cnt += 1
            # compare p/P with 13/66 and 1/2 exactly
            if 66*p < 13*P or 2*p > P: viol.append((a, b, p, P))
            if 66*p == 13*P: eq_lo.append((a, b, p, P))
            if 2*p == P: eq_hi.append((a, b, p, P))
            if p*lo[1] < lo[0]*P: lo = (p, P, (a, b))
            if p*hi[1] > hi[0]*P: hi = (p, P, (a, b))
    print(f'pairs tested {cnt}; violations {viol[:20]}; equality lower {eq_lo}; equality upper {eq_hi}')
    print('min ratio', lo, lo[0]/lo[1], 'max ratio', hi, hi[0]/hi[1])
if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'wide': run(400, lambda a: max(a+1, 3000))
    if mode == 'small_a': run(12, lambda a: 300000)
