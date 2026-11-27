# Independent implementation of P-NIM Sprague-Grundy values: partitions as tuples, subsets of rows/columns removed literally.
import sys, functools, itertools
sys.setrecursionlimit(1000000)
def conj(lam): return tuple(sum(1 for p in lam if p >= j) for j in range(1, lam[0]+1)) if lam else ()
def row_removals(lam):
    # all partitions obtained by deleting a nonempty set of parts (as multisets)
    out = set()
    from collections import Counter
    c = sorted(Counter(lam).items(), reverse=True)
    for xs in itertools.product(*[range(m+1) for _, m in c]):
        if not any(xs): continue
        out.add(tuple(v for (v, m), x in zip(c, xs) for _ in range(m - x)))
    return out
@functools.lru_cache(maxsize=None)
def sg(lam):
    if not lam: return 0
    opts = {sg(mu) for mu in row_removals(lam)} | {sg(conj(mu)) for mu in row_removals(conj(lam))}
    m = 0
    while m in opts: m += 1
    return m
if __name__ == '__main__':
    for lam in [(9,)*8, (9,8,7,6,5,4,3,2), (9,8,8,8,8,6,5,2), (9,9,7,7,6,5,4,2), (9,9,9,9,9,9,9,8), (9,9,9,9,9,9,9,2)]:
        print(lam, 'SG =', sg(lam), ' longest play =', len(lam) + lam[0] - 1, flush=True)
