# Test of Conjecture [conj:unimodalitylemma] of arXiv:2605.12822 (Connelly, Ito, Martinez, Shevchenko, Yang)
# P = [a_1]_q ... [a_k]_q [b]_{q^r};  claim: (r | a_i for some i) or b <= 1 + sum floor(a_i/r)  ==> unimodal;
# and if k<=3 or r<=3 the condition is necessary.
import sys, itertools
import numpy as np
def poly(As, b, r):
    p = np.zeros(1 + r*(b-1), dtype=object); p[::r] = 1
    p = np.array([int(x) for x in p], dtype=np.int64)
    for a in As:
        if a == 1: continue
        c = np.concatenate(([0], np.cumsum(np.concatenate((p, np.zeros(a-1, dtype=np.int64))))))
        p = c[a:] - c[:-a] if False else (c[1:] - np.concatenate((np.zeros(a, dtype=np.int64), c[1:-a])))
    return p
def unimodal(p):
    d = np.diff(p)
    # symmetric polynomial: unimodal iff nondecreasing up to the middle
    h = (len(p)) // 2
    return bool((d[:h] >= 0).all()) if len(d) else True
def run(k, rmax, amax, bmax):
    bad_suff = []; bad_nec = []; n = 0
    for r in range(2, rmax+1):
        for As in itertools.combinations_with_replacement(range(1, amax+1), k):
            cond_div = any(a % r == 0 for a in As)
            s = 1 + sum(a // r for a in As)
            for b in range(1, bmax+1):
                cond = cond_div or b <= s
                u = unimodal(poly(As, b, r)); n += 1
                if cond and not u: bad_suff.append((r, As, b))
                if (not cond) and u and (k <= 3 or r <= 3): bad_nec.append((r, As, b))
    return n, bad_suff, bad_nec
if __name__ == '__main__':
    k, rmax, amax, bmax = map(int, sys.argv[1:5])
    n, bs, bn = run(k, rmax, amax, bmax)
    print(f'k={k} r<={rmax} a<={amax} b<={bmax}: {n} cases; sufficiency failures {len(bs)}; necessity failures {len(bn)}')
    print('first sufficiency failures:', bs[:15]); print('first necessity failures:', bn[:15])
