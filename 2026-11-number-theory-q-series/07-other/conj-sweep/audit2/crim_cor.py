"""Referee re-check of Corollaries 3, 4, 5, 7 of R07 (CRIM).
Brute-force parts use ONLY the game table produced by crim_bf.c (pure game recursion on boundary words of
length <= LMAX; file given as argv[1]) -- never the automaton.  Exact parts use the table printed in R07-crim.md.
"""
import sys, itertools
from crim_table import load
E, N, ACC, REP = load()
step = {'E': E, 'N': N}
LMAX = 26
raw = open(sys.argv[1], 'rb').read()
off = [0]
for L in range(LMAX + 1): off.append(off[-1] + (1 << L) // 8 + 1)
assert off[-1] == len(raw)
def game_losing(w):
    """losing status of canon(w) from the brute-force table (no automaton)"""
    L = len(w); assert L <= LMAX
    v = 0
    for i, x in enumerate(w):
        if x == 'N': v |= 1 << i
    return (raw[off[L] + (v >> 3)] >> (v & 7)) & 1 == 1
def run(w, q=0):
    for x in w: q = step[x][q]
    return q
def tf(w): return tuple(run(w, q) for q in range(17))
def then(f, g): return tuple(g[f[q]] for q in range(17))     # f first, then g
def word(lam):
    w = ''; prev = 0
    for part in reversed(lam): w += 'E' * (part - prev) + 'N'; prev = part
    return w
def conj(lam): return tuple(sum(1 for p in lam if p >= j) for j in range(1, (lam[0] if lam else 0) + 1))
# ---------------------------------------------------------------- literal definitions from the source
def meld(mu, lam):          # mu <> lam = [mu_1+lam_1,...,mu_1+lam_m, mu_1,...,mu_k]
    if not mu: return tuple(lam)
    return tuple(mu[0] + x for x in lam) + tuple(mu)
def durfee(lam):
    d = 0
    while d < len(lam) and lam[d] >= d + 1: d += 1
    return d
def stack(lam, mu):         # lam <| mu ; None if undefined (source: lam_{d+i} = d for 1<=i<=k and lam_d >= d + mu_1)
    d = durfee(lam); k = len(mu); r = len(lam)
    if d == 0: return None
    if d + k > r: return None
    if any(lam[d + i - 1] != d for i in range(1, k + 1)): return None
    if lam[d - 1] < d + (mu[0] if mu else 0): return None
    return tuple(lam[:d]) + tuple(lam[d + i] + mu[i] for i in range(k)) + tuple(lam[d + k:])
def hook(r, c): return (c,) + (1,) * (r - 1)
def thick(r, c, a, b):
    assert r > a >= 2 and c > b >= 2
    return (c,) * a + (b,) * (r - a)
# checks of the encodings on small cases
for mu in [(), (1,), (2, 1), (3, 3, 1), (4, 2)]:
    for lam in [(), (1,), (2, 2), (5, 1, 1), (3, 2, 1)]:
        m = meld(mu, lam)
        assert word(m) == word(mu) + word(lam), (mu, lam)
        assert conj(m) == meld(conj(lam), conj(mu))
print('meld = concatenation of boundary words (and (mu<>lam)\' = lam\'<>mu\'): checked on samples')
assert stack(hook(8, 10), stack(hook(5, 7), hook(3, 3))) is not None      # the source's figure example is defined
# ---------------------------------------------------------------- Corollary 3
c3 = all(tf(x * 2) == tf(x * 4) for x in 'EN')
c3b = all(tf(x * 1) != tf(x * 3) for x in 'EN')
print('Cor 3: E^2 = E^4 and N^2 = N^4 as state maps:', c3, ' (so x^k = x^{k+2} for all k >= 2); x^1 != x^3:', c3b)
CL = (1, 2, 3)              # classes of a positive exponent: 1, even >= 2, odd >= 3
# ---------------------------------------------------------------- Corollary 4: thick hooks
# (a) brute force: status of every thick hook with r + c <= LMAX, and all melds of losing thick hooks of total length <= LMAX
nth = bad_th = 0; losing_th = []
for r in range(3, LMAX):
    for c in range(3, LMAX - r + 1):
        for a in range(2, r):
            for b in range(2, c):
                lam = thick(r, c, a, b); w = word(lam)
                assert w == 'E' * b + 'N' * (r - a) + 'E' * (c - b) + 'N' * a
                nth += 1
                g = game_losing(w)
                if g != ((r - c) % 2 == 0): bad_th += 1
                if g: losing_th.append(w)
print(f'Cor 4 (brute force): {nth} thick hooks with r+c <= {LMAX}: "losing iff r = c mod 2" fails {bad_th} times')
by_len = {}
for w in losing_th: by_len.setdefault(len(w), []).append(w)
nm = badm = 0
def rec(prefix, nparts):
    global nm, badm
    if nparts >= 2:
        nm += 1
        if not game_losing(prefix): badm += 1
    for L, ws in by_len.items():
        if len(prefix) + L <= LMAX:
            for w in ws: rec(prefix + w, nparts + 1)
rec('', 0)
print(f'Cor 4 (brute force): {nm} melds of >= 2 losing thick hooks with total rows+columns <= {LMAX}: not losing: {badm}')
# (b) exact: all exponent classes
img0 = set(); img2 = set(); wrong = 0; ncls = 0
for b in (2, 3):
    for s in CL:
        for t in CL:
            for a in (2, 3):
                f = tf('E' * b + 'N' * s + 'E' * t + 'N' * a)
                same = (s + a - b - t) % 2 == 0          # r = s + a, c = b + t
                if ACC[f[0]] != same: wrong += 1
                if same: ncls += 1; img0.add(f[0]); img2.add(f[2])
print(f'Cor 4 (exact): {ncls} classes of losing thick hooks; automaton status != [r = c mod 2] in {wrong} of 36 classes;')
print(f'               images of state 0: {sorted(img0)}, of state 2: {sorted(img2)}; state 2 accepting: {ACC[2]}')
# ---------------------------------------------------------------- Corollary 5: conj:bloby
# (a) the word formula, from the literal definition of the stack
nform = 0
def parts(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n - k, k): yield (k,) + r
small = [lam for n in range(0, 9) for lam in parts(n)]
for lam in small:
    k = len(lam); l1 = lam[0] if lam else 0
    for r in range(1, 9):
        for c in range(1, 9):
            st = stack(hook(r, c), lam)
            defined = (k <= r - 1 and l1 <= c - 1)
            assert (st is not None) == defined, (r, c, lam)
            if st is not None:
                assert word(st) == 'E' + 'N' * (r - 1 - k) + word(lam) + 'E' * (c - 1 - l1) + 'N'; nform += 1
                assert stack(lam, hook(r, c)) is None or not lam     # "if lam<|mu is defined then mu<|lam is not"
print(f'Cor 5: w(Gamma_(r,c) <| lam) = E N^(r-1-k) w(lam) E^(c-1-lam_1) N and definedness k<=r-1, lam_1<=c-1: {nform} cases checked literally')
# (b) brute force with the game table: all losing lam (word length <= LMAX-2), all losing hooks, stack defined
nb = badb = 0; first_bad = None
for L in range(0, LMAX - 1):
    for mid in itertools.product('EN', repeat=max(L - 2, 0)):
        if L == 1: continue
        w = '' if L == 0 else 'E' + ''.join(mid) + 'N'
        if L > 12 and False: pass
        if not game_losing(w): continue
        k = w.count('N'); l1 = w.count('E')
        for a in range(0, LMAX - 2 - L + 1):
            for b in range(0, LMAX - 2 - L - a + 1):
                r = a + 1 + k; c = b + 1 + l1
                hw = 'E' + 'N' * (r - 1) + 'E' * (c - 1) + 'N'
                if not (r >= 2 and c >= 2 and (r - c) % 2 == 0): continue
                assert game_losing(hw) if len(hw) <= LMAX else True
                nb += 1
                if not game_losing('E' + 'N' * a + w + 'E' * b + 'N'):
                    badb += 1
                    if first_bad is None: first_bad = (r, c, w)
    if L >= 18: break       # keep the enumeration of lam bounded (2^16 words of length 18)
print(f'Cor 5 (brute force): {nb} triples (losing lam with rows+cols <= 18, losing hook, stack defined, total <= {LMAX}): not losing: {badb} {first_bad or ""}')
# hooks: losing iff r, c >= 2 and r = c (mod 2)   (source, Proposition lem:hook)
hk = all(game_losing(word(hook(r, c))) == (r >= 2 and c >= 2 and (r - c) % 2 == 0) for r in range(1, 14) for c in range(1, 14) if r + c <= LMAX)
print('       hooks: losing iff r,c >= 2 and r = c (mod 2), r,c <= 13 (brute force):', hk)
# (c) exact: for each (a, b) in {0..5}^2 with a = b (mod 2): all NONEMPTY canonical w with B(w) => B(E N^a w E^b N)
viol = []; npairs = 0
for a in range(0, 6):
    s = run('E' + 'N' * a)
    init = (E[0], E[s], 'E'); seen = {init}; todo = [init]
    while todo:
        p, t, last = todo.pop()
        if last == 'N' and ACC[p]:
            npairs += 1
            for b in range(a % 2, 6, 2):
                if not ACC[run('E' * b + 'N', t)]: viol.append((a, b, p, t))
        for x in 'EN':
            c2 = (step[x][p], step[x][t], x)
            if c2 not in seen: seen.add(c2); todo.append(c2)
print(f'Cor 5 (exact): product search, a, b in 0..5, a = b (mod 2), nonempty losing w: {npairs} accepting pair-configurations, violations {len(viol)}')
# empty lam: lam* is the hook itself (losing by hypothesis); note the implication is FALSE for w = empty, a = 0 or b = 0:
print('       remark: w = empty, a = 0, b = 2 gives', 'E' + 'EE' + 'N', 'accepted:', ACC[run('EEEN')], '(row [3]; excluded since r = 1)')
# is the parity hypothesis needed?  (a != b mod 2)
viol_odd = 0
for a in range(0, 6):
    s = run('E' + 'N' * a)
    init = (E[0], E[s], 'E'); seen = {init}; todo = [init]
    while todo:
        p, t, last = todo.pop()
        if last == 'N' and ACC[p]:
            for b in range((a + 1) % 2, 6, 2):
                if ACC[run('E' * b + 'N', t)]: viol_odd += 1
        for x in 'EN':
            c2 = (step[x][p], step[x][t], x)
            if c2 not in seen: seen.add(c2); todo.append(c2)
print('       (for a != b mod 2 the stack is never losing -- accepted cases:', viol_odd, ', consistent with conj:main)')
# ---------------------------------------------------------------- Corollary 7: conj:pself, "if" half
# (a) literal generation of all stacks of odd self-conjugate hooks and [1] with rows+cols <= LMAX
stacks = {(1,)}; frontier = [(1,)]
while frontier:
    new = []
    for mu in frontier:
        for r in range(3, LMAX // 2 + 1, 2):
            st = stack(hook(r, r), mu)
            if st is not None and len(st) + st[0] <= LMAX and st not in stacks:
                stacks.add(st); new.append(st)
    frontier = new
rule = {'EN'}; fr = ['EN']
while fr:
    nw = []
    for w in fr:
        for a in range(1, LMAX, 2):
            v = 'E' + 'N' * a + w + 'E' * a + 'N'
            if len(v) <= LMAX and v not in rule: rule.add(v); nw.append(v)
    fr = nw
print(f'Cor 7: {len(stacks)} stacks of odd self-conjugate hooks and [1] with rows+cols <= {LMAX} (literal definition);',
      'same set as the rule w -> E N^a w E^a N (a odd):', {word(s) for s in stacks} == rule)
def psc(lam):
    cj = conj(lam)
    return len(cj) == len(lam) and all((x - y) % 2 == 0 for x, y in zip(lam, cj))
print('       all of them parity self-conjugate (indeed self-conjugate):', all(psc(s) and conj(s) == s for s in stacks),
      '; winning by brute force:', sum(1 for s in stacks if not game_losing(word(s))), 'of', len(stacks))
# left-nested stacks give the same partitions (associativity, as asserted in the source)
assoc_bad = 0; assoc_n = 0
for r1 in range(3, 12, 2):
    for r2 in range(3, r1, 2):
        for r3 in (1, 3, 5):
            if r3 >= r2: continue
            inner = stack(hook(r2, r2), hook(r3, r3) if r3 > 1 else (1,))
            right = stack(hook(r1, r1), inner) if inner else None
            l12 = stack(hook(r1, r1), hook(r2, r2))
            left = stack(l12, hook(r3, r3) if r3 > 1 else (1,)) if l12 else None
            assoc_n += 1
            if left != right: assoc_bad += 1
print(f'       associativity (G1<|G2)<|G3 = G1<|(G2<|G3) on {assoc_n} triples of odd square hooks: failures {assoc_bad}')
# (b) exact
F = {tf('EN')}; fr = [tf('EN')]
while fr:
    nw = []
    for f in fr:
        for a in (1, 3):
            g = then(then(tf('E' + 'N' * a), f), tf('E' * a + 'N'))
            if g not in F: F.add(g); nw.append(g)
    fr = nw
print(f'Cor 7 (exact): closure has {len(F)} state maps; accepted from the start state: {[ACC[f[0]] for f in F]}')
# (c) full conj:pself by brute force (both directions), all partitions with rows+cols <= LMAX and n <= 44
cnt = bad = 0; examples = []
for n in range(1, 45):
    for lam in parts(n):
        if len(lam) + lam[0] > LMAX or not psc(lam): continue
        cnt += 1
        if (not game_losing(word(lam))) != (lam in stacks):
            bad += 1
            if len(examples) < 3: examples.append(lam)
print(f'       full conj:pself, brute force, {cnt} parity self-conjugate partitions (n <= 44, rows+cols <= {LMAX}): counterexamples {bad} {examples}')
