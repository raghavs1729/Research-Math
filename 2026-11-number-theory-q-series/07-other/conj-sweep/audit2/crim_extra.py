"""Referee: extra checks for R07.
(a) move sets: for every partition of n <= 16, {options by literal row/column removal on cell sets} ==
    {partition of canon(w(lam) minus one letter)} as sets -- i.e. every single-letter deletion is a legal move and
    every legal move is a single-letter deletion.
(b) Corollary 2: count of losing words by length via the transfer matrix, L <= 120, and the recurrence argument.
(c) Corollary 6: number of losing classes of 4-part partitions.
(d) source statements: thm:2parts, thm:3parts (literal), prop:hooksquare, lem:hook, thm:thick-hook, appendix list."""
import itertools
from crim_table import load
E, N, ACC, REP = load(); step = {'E': E, 'N': N}
def run(w, q=0):
    for x in w: q = step[x][q]
    return q
def word(lam):
    w = ''; prev = 0
    for part in reversed(lam): w += 'E' * (part - prev) + 'N'; prev = part
    return w
def unword(w):
    lam = []; c = 0
    for x in w:
        if x == 'E': c += 1
        else: lam.append(c)
    return tuple(reversed(lam))
def canon(w): return w.lstrip('N').rstrip('E')
def parts(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n - k, k): yield (k,) + r
def cells(lam): return {(i, j) for i, p in enumerate(lam) for j in range(p)}
def shape(S):
    rows = sorted({i for i, j in S}); cols = sorted({j for i, j in S})
    ri = {r: k for k, r in enumerate(rows)}; ci = {c: k for k, c in enumerate(cols)}
    T = {(ri[i], ci[j]) for i, j in S}
    lam = [0] * len(rows)
    for i, j in T: lam[i] += 1
    assert all(lam[i] >= lam[i + 1] for i in range(len(lam) - 1)) and T == cells(lam)   # result is again a Young diagram
    return tuple(lam)
cnt = 0
for n in range(1, 17):
    for lam in parts(n):
        S = cells(lam); w = word(lam)
        assert unword(w) == lam and w[0] == 'E' and w[-1] == 'N' and w.count('N') == len(lam) and w.count('E') == lam[0]
        lit = {shape({c for c in S if c[0] != r}) for r in range(len(lam))} | {shape({c for c in S if c[1] != cc}) for cc in range(lam[0])}
        wd = {unword(canon(w[:i] + w[i + 1:])) for i in range(len(w))}
        assert lit == wd, lam
        cnt += 1
print(f'(a) move sets (literal row/column removal vs single-letter deletion + canon) agree for all {cnt} partitions of 1 <= n <= 16')
# (b)
def count(L):
    vec = [0] * 17; vec[E[0]] = 1
    for _ in range(L - 2):
        nv = [0] * 17
        for q in range(17):
            if vec[q]: nv[E[q]] += vec[q]; nv[N[q]] += vec[q]
        vec = nv
    return sum(vec[q] for q in range(17) if ACC[N[q]])
cs = {L: count(L) for L in range(2, 121)}
ok = all(cs[L] == (5 * 2 ** (L - 5) if L % 2 == 0 else 0) for L in range(8, 121))
print('(b) losing canonical words of length L: L=2..12:', [cs[L] for L in range(2, 13)], '; = 5*2^(L-5) (L even), 0 (L odd) for all 8 <= L <= 120:', ok)
print('    (sequence satisfies a linear recurrence of order <= 17+2, so agreement on 19 consecutive L proves it for all L >= 8)')
# (c)
cls = lambda k: '0' if k == 0 else '1' if k == 1 else 'e' if k % 2 == 0 else 'o'
tab = {}
for l4 in range(1, 8):
    for x in range(0, 8):
        for y in range(0, 8):
            for z in range(0, 8):
                key = (cls(l4), cls(x), cls(y), cls(z)); v = ACC[run(word((l4 + x + y + z, l4 + x + y, l4 + x, l4)))]
                assert tab.setdefault(key, v) == v
print('(c) four parts:', sum(tab.values()), 'losing classes of', len(tab))
print('    losing classes:', ' '.join(''.join(k) for k in sorted(k for k, v in tab.items() if v)))
# (d)
P = lambda lam: ACC[run(word(lam))]
print('(d) thm:2parts ([l1,l2] in P iff l1 even, l2 > 0), l1 <= 60:', all(P((a, b)) == (a % 2 == 0) for a in range(1, 61) for b in range(1, a + 1)))
lit3 = lambda a, b, c: a % 2 == 1 and b % 2 == 1 and ((b == 1 and c == 1) or c > 1)
ex3 = [(a, b, c) for a in range(1, 40) for b in range(1, a + 1) for c in range(1, b + 1) if P((a, b, c)) != lit3(a, b, c)]
print('    thm:3parts, literal statement, l1 < 40: exceptions', ex3)
print('    prop:hooksquare ([c,2,1^(r-2)] in P iff r, c even), 2 <= r, c <= 40:', all(P((c, 2) + (1,) * (r - 2)) == (r % 2 == 0 and c % 2 == 0) for r in range(2, 41) for c in range(2, 41)))
print('    lem:hook (hook in P iff min(r,c) >= 2 and r = c mod 2), r, c <= 40:', all(P((c,) + (1,) * (r - 1)) == (min(r, c) >= 2 and (r - c) % 2 == 0) for r in range(1, 41) for c in range(1, 41)))
print('    padding examples of the source: [4,1] in P:', P((4, 1)), ' [4,4,4,1] in N:', not P((4, 4, 4, 1)), ' [5^5,3^7,1^7] in N:', not P((5,) * 5 + (3,) * 7 + (1,) * 7), ' [7^5,5^7,3^7] in P:', P((7,) * 5 + (5,) * 7 + (3,) * 7), ' [4,3] in P:', P((4, 3)), ' [4,3,1,1] in N:', not P((4, 3, 1, 1)), ' [3,1,1] in P:', P((3, 1, 1)))
