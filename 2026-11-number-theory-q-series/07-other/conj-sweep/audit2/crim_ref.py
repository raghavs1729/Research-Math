"""Referee re-check of Theorem 1 of R07 (CRIM).  Independent code; automaton read from the table in R07-crim.md.
Part 1: properties (i), (ii), (iv), (v) by an exhaustive search of configurations (bitmask subsets).
Part 2: direct word-level test of (iv) on all canonical words of length <= 16 (no configurations).
Part 3: brute-force game solver on Young diagrams stored as SETS OF CELLS (rows/columns removed literally,
        pieces re-attached), all partitions of n <= NMAX, compared with the automaton on the boundary word.
"""
import sys, itertools
from crim_table import load
E, N, ACC, REP = load()
step = {'E': E, 'N': N}
def run(w, q=0):
    for x in w: q = step[x][q]
    return q
ACCMASK = sum(1 << q for q in range(17) if ACC[q])
# representatives really reach their states
assert all(run(REP[q]) == q for q in range(17)), 'representatives'
print('representative words reach their states: True; all 17 states reachable')
# ---- Part 1
p_i = (N[0] == 0 and ACC[0])
p_ii = all(ACC[E[q]] == ACC[q] for q in range(17))
print('(i) delta(0,N)=0 and 0 accepting:', p_i, '   (ii) acc(delta(q,E)) == acc(q) for all q:', p_ii)
def img(mask, tab):
    out = 0; q = 0
    while mask:
        if mask & 1: out |= 1 << tab[q]
        mask >>= 1; q += 1
    return out
# configuration: (state of w, set of states of the |w| words w-minus-one-letter, parity of |w|, last letter)
start = (E[0], 1 << 0, 1, 'E')           # w = 'E'
seen = {start}; todo = [start]; viol = []; odd_acc = []; nchecked = 0
while todo:
    q, T, par, last = todo.pop()
    if last == 'N':                       # w starts with E and ends with N
        nchecked += 1
        if ACC[q] != ((T & ACCMASK) == 0): viol.append((q, T, par))
        if ACC[q] and par == 1: odd_acc.append((q, T))
    for x in 'EN':
        c = (step[x][q], img(T, step[x]) | (1 << q), par ^ 1, x)
        if c not in seen: seen.add(c); todo.append(c)
proj = {(q, T) for q, T, par, last in seen}
projp = {(q, T, par) for q, T, par, last in seen}
print('configurations (q,T,parity,last letter):', len(seen), '; projected to (q,T,parity):', len(projp), '; to (q,T):', len(proj))
print('(iv) canonical configurations checked:', nchecked, ' violations:', len(viol))
print('(v) accepted canonical configurations of odd length:', len(odd_acc))
# simpler proof of (v): product of B with the parity bit only
st = {(E[0], 1)}; fr = [(E[0], 1)]; bad = []
while fr:
    q, par = fr.pop()
    for x in 'EN':
        c = (step[x][q], par ^ 1)
        if x == 'N' and ACC[c[0]] and c[1] == 1: bad.append(c)
        if c not in st: st.add(c); fr.append(c)
print('(v\') product B x parity: reachable', len(st), 'pairs; accepted odd canonical:', len(bad))
# stronger: is ANY accepted word starting with E of odd length (canonical or not)?
print('     pairs (accepting state, odd length) reachable at all after a first letter E:', sorted(c for c in st if ACC[c[0]] and c[1] == 1))
# ---- Part 2: direct test of (iv) on words
LW = 16; nw = 0; badw = 0
for L in range(1, LW + 1):
    for mid in itertools.product('EN', repeat=max(L - 2, 0)):
        if L == 1: continue
        w = 'E' + ''.join(mid) + 'N'
        nw += 1
        lhs = ACC[run(w)]
        rhs = not any(ACC[run(w[:i] + w[i+1:])] for i in range(L))
        if lhs != rhs: badw += 1
print(f'(iv) tested directly on all {nw} canonical words of length <= {LW}: violations {badw}')
# ---- Part 3: literal game on cell sets
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
def cells(lam): return frozenset((i, j) for i, p in enumerate(lam) for j in range(p))
def normalise(S):
    """re-attach: compress used row indices and used column indices"""
    rows = sorted({i for i, j in S}); cols = sorted({j for i, j in S})
    ri = {r: k for k, r in enumerate(rows)}; ci = {c: k for k, c in enumerate(cols)}
    return frozenset((ri[i], ci[j]) for i, j in S)
memo = {frozenset(): True}
def losing(S):
    v = memo.get(S)
    if v is not None: return v
    rows = {i for i, j in S}; cols = {j for i, j in S}
    res = True
    for r in rows:
        if losing(normalise(frozenset(c for c in S if c[0] != r))): res = False; break
    if res:
        for cc in cols:
            if losing(normalise(frozenset(c for c in S if c[1] != cc))): res = False; break
    memo[S] = res
    return res
def word(lam):
    w = ''; prev = 0
    for part in reversed(lam): w += 'E' * (part - prev) + 'N'; prev = part
    return w
def parts(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n - k, k): yield (k,) + r
sys.setrecursionlimit(10000)
cnt = mism = oddP = nP = 0
for n in range(0, NMAX + 1):
    for lam in parts(n):
        cnt += 1
        g = losing(cells(lam)); a = ACC[run(word(lam))]
        if g != a: mism += 1
        if g:
            nP += 1
            if lam and (len(lam) - lam[0]) % 2: oddP += 1
print(f'cell-set brute force vs automaton on all {cnt} partitions of n <= {NMAX}: mismatches {mism}; losing {nP}; losing of odd rank {oddP}')
