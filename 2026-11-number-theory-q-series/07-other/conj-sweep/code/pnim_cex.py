# Counterexamples to conj:heavychoppedrect of arXiv:2506.04991, with sanity checks against the paper's theorems.
import sys
from math import comb
from pnim_check2 import sg as sg2            # tuples implementation
sys.argv = ['x', '0']
from pnim import sg as sg1, word             # boundary-word implementation
def partitions(n, mx=None):
    mx = mx or n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n-k, k): yield (k,) + rest
# sanity 1: Theorem thm:partition_interval: SG = 1 iff lam = [1] or [r,r,r-1,...,2] <= lam <= [r^r] for some r >= 2
def in_interval(lam):
    r = len(lam)
    if lam == (1,): return True
    if r < 2: return False
    lo = (r, r) + tuple(range(r-1, 1, -1))
    return all(lo[j] <= lam[j] <= r for j in range(r))
ok = True; cnt = 0
for n in range(1, 27):
    for lam in partitions(n):
        cnt += 1
        v = sg1(word(lam))
        if (v == 1) != in_interval(lam): ok = False
        if n <= 18 and v != sg2(lam): ok = False
print(f'sanity: SG = 1 characterisation (their Theorem) and agreement of the two implementations on all {cnt} partitions of n <= 26 (n <= 18 for the comparison):', ok)
print('sanity: rectangles [c^r], r, c <= 9: heavy iff C(c+r-2, r-1) odd:', all((sg1(word((c,)*r)) == r+c-1) == (comb(c+r-2, r-1) % 2 == 1) for r in range(1, 10) for c in range(1, 10)))
a, b = 8, 7
lo = tuple(a+1-j for j in range(b+1))
def gen(j, prev):
    if j == b+1: yield (); return
    for v in range(lo[j], min(prev, a+1)+1):
        for rest in gen(j+1, v): yield (v,) + rest
allp = list(gen(0, a+1)); bad = [(lam, sg1(word(lam)), sg2(lam)) for lam in allp if sg1(word(lam)) != a+b+1]
print(f'a = 8, b = 7: rectangle [9^8] has SG {sg2((9,)*8)} (heavy, C(15,7) = {comb(15,7)} odd); interval [9,8,...,2] <= lam <= [9^8]: {len(allp)} partitions, {len(bad)} NOT heavy')
for lam, v1, v2 in sorted(bad, key=lambda t: sum(t[0]))[:6] + sorted(bad, key=lambda t: -sum(t[0]))[:2]:
    print('   ', lam, 'size', sum(lam), 'SG (two implementations) =', v1, v2, '; longest play 16')
