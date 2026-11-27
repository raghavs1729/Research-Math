# Cross-checks of the CRIM automaton (cert/crim_dfa.json) against published theorems of arXiv:2606.16828
# and against brute-force game search (independent of the learning code).
import json, itertools, functools, sys
D = json.load(open('../cert/crim_dfa.json')); delta = D['delta']; acc = D['acc']
def word(lam):
    lam = [x for x in lam if x > 0]
    w = ''; prev = 0
    for part in reversed(lam): w += 'E'*(part-prev) + 'N'; prev = part
    return w
def dfaP(lam):
    q = 0
    for x in word(lam): q = delta[f'{q}{x}']
    return acc[q]
# brute force directly on partitions (tuples), rows/columns removed literally
@functools.lru_cache(maxsize=None)
def bruteP(lam):
    if not lam: return True
    opts = set()
    for i in range(len(lam)): opts.add(lam[:i] + lam[i+1:])
    for j in range(1, lam[0]+1): opts.add(tuple(x for x in (p - (1 if p >= j else 0) for p in lam) if x > 0))
    return not any(bruteP(o) for o in opts)
def partitions(n, mx=None):
    mx = mx or n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n-k, k): yield (k,) + rest
bad = 0; cnt = 0
for n in range(0, int(sys.argv[1])+1):
    for lam in partitions(n):
        cnt += 1
        if dfaP(lam) != bruteP(lam): bad += 1
print(f'all {cnt} partitions of n <= {sys.argv[1]}: automaton vs brute force (row/column removal on tuples): mismatches {bad}')
ok2 = all(dfaP((a, b)) == (a % 2 == 0) for a in range(1, 80) for b in range(1, a+1))
ok3 = all(dfaP((a, b, c)) == (a % 2 == 1 and b % 2 == 1 and ((b == 1 and c == 1) or c > 1)) for a in range(1, 50) for b in range(1, a+1) for c in range(1, b+1))
print('Proposition thm:2parts ([l1,l2] in P iff l1 even) reproduced for l1 < 80:', ok2)
print('Theorem thm:3parts reproduced for l1 < 50:', ok3)
app = [(4,2,2,2),(6,2,2,2),(6,4,1,1),(6,3,2,2),(6,4,2,1),(4,4,3,3),(6,4,2,2),(6,6,1,1),(8,2,2,2),(8,4,1,1),(6,4,3,2),(6,5,3,1),(5,4,2,2,1),(5,3,3,2,2),(5,3,2,2,2),(5,3,3,2,1),(8,3,2,2),(8,4,2,1),(8,6,1,1)]
print('appendix "unexplained losing partitions" (first 19) accepted by the automaton:', all(dfaP(l) for l in app))
# counts by semi-perimeter L
def count(L):
    from collections import Counter
    c = Counter({int(delta['0E']): 1})
    for _ in range(L-1):
        c2 = Counter()
        for q, m in c.items():
            for x in 'EN': c2[delta[f'{q}{x}']] += m
        c = c2
    # words of length L starting with E and ending with N: redo with last letter forced N
    c = Counter({int(delta['0E']): 1})
    for _ in range(L-2):
        c2 = Counter()
        for q, m in c.items():
            for x in 'EN': c2[delta[f'{q}{x}']] += m
        c = c2
    return sum(m for q, m in c.items() if acc[delta[f'{q}N']])
print('number of P-positions with rows+columns = L, L=2..24:', [count(L) for L in range(2, 25)])
