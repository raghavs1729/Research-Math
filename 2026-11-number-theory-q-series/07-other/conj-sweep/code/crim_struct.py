# CRIM: P-positions as boundary words; look for a local characterisation.
import sys, itertools, functools
sys.setrecursionlimit(100000)
def canon(w):
    w = w.lstrip('N').rstrip('E')
    return w
@functools.lru_cache(maxsize=None)
def isP(w):
    if w == '': return True
    seen = set()
    for i in range(len(w)):
        c = canon(w[:i] + w[i+1:])
        if c in seen: continue
        seen.add(c)
        if isP(c): return False
    return True
def words(L):
    for mid in itertools.product('EN', repeat=L-2): yield 'E' + ''.join(mid) + 'N'
if __name__ == '__main__':
    for L in (4, 6, 8, 10, 12, 14):
        Pw = [w for w in words(L) if isP(w)]
        print(L, len(Pw))
        if L <= 6: print('  ', Pw)
        # is P-status determined by (prefix of length a, suffix of length b)?
        for a, b in ((2,2),(3,2),(2,3),(3,3),(4,3),(3,4),(4,4)):
            if a + b > L: continue
            d = {}
            ok = True
            for w in words(L):
                k = (w[:a], w[-b:]); v = isP(w)
                if d.setdefault(k, v) != v: ok = False
            if ok:
                print(f'  L={L}: determined by prefix {a}, suffix {b}:', sorted(k for k, v in d.items() if v)); break
