# arXiv:2610.06925 "Removing m from m-pile Divisor Nim": conjecture g(P) <= 2*h_min (h_min = smallest nonempty heap).
# Move: remove d >= 1 stones from one heap, where d divides the size of every other nonempty heap.
import sys, functools, itertools
from math import gcd
sys.setrecursionlimit(1000000)
@functools.lru_cache(maxsize=None)
def g(P):                      # P sorted tuple of positive heap sizes
    if not P: return 0
    opts = set()
    for i in range(len(P)):
        if i and P[i] == P[i-1]: continue
        G = 0
        for j in range(len(P)):
            if j != i: G = gcd(G, P[j])
        h = P[i]
        ds = range(1, h+1) if G == 0 else [d for d in range(1, min(h, G)+1) if G % d == 0]
        for d in ds:
            Q = list(P); Q[i] = h - d
            opts.add(g(tuple(sorted(x for x in Q if x > 0))))
    m = 0
    while m in opts: m += 1
    return m
worst = {}
for m, N in ((2, int(sys.argv[1])), (3, int(sys.argv[2])), (4, int(sys.argv[3])), (5, int(sys.argv[4]))):
    bad = []; mx = (0, None)
    for P in itertools.combinations_with_replacement(range(1, N+1), m):
        v = g(P); hm = P[0]
        if v > 2*hm: bad.append((P, v))
        if v - 2*hm > mx[0] - 0 or mx[1] is None: 
            if mx[1] is None or v - 2*hm > mx[0]: mx = (v - 2*hm, P)
    print(f'{m} heaps, sizes <= {N}: violations of g <= 2 h_min: {len(bad)} {bad[:6]}; max of g - 2 h_min = {mx[0]} at {mx[1]}', flush=True)
