# arXiv:2506.04991 (Nim on integer partitions): P-NIM on one partition; a move removes a nonempty set of rows
# or a nonempty set of columns.  "heavy": SG value = length of the longest play = (#rows) + (#columns) - 1.
# Conjectures conj:shallowstair and conj:heavychoppedrect.
import sys, functools, itertools
sys.setrecursionlimit(1000000)
def runs(w):
    r = []
    for ch in w:
        if r and r[-1][0] == ch: r[-1][1] += 1
        else: r.append([ch, 1])
    return r
@functools.lru_cache(maxsize=None)
def sg(w):
    if not w: return 0
    R = runs(w); opts = set()
    for letter in 'EN':
        idx = [k for k, (ch, l) in enumerate(R) if ch == letter]
        for xs in itertools.product(*[range(R[k][1]+1) for k in idx]):
            if not any(xs): continue
            parts = []
            it = iter(xs); 
            for k, (ch, l) in enumerate(R):
                parts.append(ch * (l - next(it)) if ch == letter else ch * l)
            v = ''.join(parts).lstrip('N').rstrip('E')
            opts.add(sg(v))
    m = 0
    while m in opts: m += 1
    return m
def word(lam):
    w = ''; prev = 0
    for p in reversed(lam): w += 'E'*(p-prev) + 'N'; prev = p
    return w
def heavy(lam): return sg(word(lam)) == len(lam) + lam[0] - 1
if __name__ == '__main__':
    LMAX = int(sys.argv[1])
    # sanity: rectangle [c^r] heavy iff C(c+r-2, r-1) odd (their theorem)
    from math import comb
    assert all(heavy((c,)*r) == (comb(c+r-2, r-1) % 2 == 1) for r in range(1, 8) for c in range(1, 8))
    bad = []; cnt = 0
    for k in range(1, 12):
        for s in range(1, 12):
            for i in range(1, 12):
                lam = tuple(i + (k-1-j)*s for j in range(k))
                if len(lam) + lam[0] > LMAX: continue
                cnt += 1
                if not heavy(lam): bad.append((i, s, k, lam, sg(word(lam))))
    print(f'conj:shallowstair: {cnt} partitions with rows+columns <= {LMAX}; not heavy: {bad[:8]}', flush=True)
    bad2 = []; cnt2 = 0
    for a in range(0, 9):
        for b in range(0, a+1):
            if a + b + 2 > LMAX or comb(a+b, b) % 2 == 0: continue
            lo = tuple(a+1-j for j in range(b+1)); hi = (a+1,)*(b+1)
            # all partitions between lo and hi (b+1 parts, lo_j <= lam_j <= a+1)
            def gen(j, prev):
                if j == b+1: yield (); return
                for v in range(lo[j], min(prev, a+1)+1):
                    for rest in gen(j+1, v): yield (v,) + rest
            for lam in gen(0, a+1):
                cnt2 += 1
                if not heavy(lam): bad2.append((a, b, lam, sg(word(lam))))
    print(f'conj:heavychoppedrect: {cnt2} partitions tested (a+b+2 <= {LMAX}); not heavy: {bad2[:8]}', flush=True)
