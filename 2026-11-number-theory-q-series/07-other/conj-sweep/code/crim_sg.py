# CRIM: are the normal-play and misere Grundy functions regular (finite Nerode index) on boundary words?
import itertools, functools, sys
sys.setrecursionlimit(100000)
def canon(w): return w.lstrip('N').rstrip('E')
def mex(s):
    m = 0
    while m in s: m += 1
    return m
@functools.lru_cache(maxsize=None)
def G(w, misere):
    if w == '': return 1 if misere else 0
    return mex({G(canon(w[:i] + w[i+1:]), misere) for i in range(len(w))})
if __name__ == '__main__':
    SUF, PRE, mis = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    f = lambda w: G(canon(w), mis)
    sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
    sig = lambda u: tuple(f(u+v) for v in sufs)
    classes = {}; mx = 0
    for l in range(0, PRE+1):
        new = 0
        for u in itertools.product('EN', repeat=l):
            u = ''.join(u)
            if l and u[0] == 'N': continue
            s = sig(u); mx = max(mx, max(s))
            if s not in classes: classes[s] = u; new += 1
        print(f'misere={mis} prefix length {l}: new classes {new}, total {len(classes)}, max value seen {mx}', flush=True)
