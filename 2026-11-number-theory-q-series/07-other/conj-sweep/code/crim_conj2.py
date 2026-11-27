"""CRIM: (3) one direction of conj:pself; (4) explicit list of losing partitions with four parts."""
import json, itertools, sys
D = json.load(open('../cert/crim_dfa.json')); n = D['n']; acc = D['acc']
d = lambda q, x: D['delta'][f'{q}{x}']
def tf(w):                       # state transformer of a word
    out = []
    for q in range(n):
        for x in w: q = d(q, x)
        out.append(q)
    return tuple(out)
comp = lambda f, g: tuple(g[f[q]] for q in range(n))      # first f then g
# (3) stacks of odd self-conjugate hooks and [1]:  W_0 = EN,  W -> E N^a W E^a N with a = c-1-k, c odd, k = #rows(W) odd => a odd.
#     N^a, E^a for odd a >= 3 act as N^3, E^3 (orbits have preperiod <= 2 and period 2), so a in {1, 3} represents all.
for x in 'EN':
    assert all(tf(x*k) == tf(x*(k+2)) for k in range(2, 8)), 'orbit claim'
F = {tf('EN')}; frontier = [tf('EN')]
while frontier:
    new = []
    for f in frontier:
        for a in (1, 3):
            g = comp(comp(comp(comp(tf('E'), tf('N'*a)), f), tf('E'*a)), tf('N'))
            if g not in F: F.add(g); new.append(g)
    frontier = new
print('(3) transformers of stacks of odd self-conjugate hooks and [1]:', len(F), '; accepted (losing) from the start state:', [acc[f[0]] for f in F])
# numeric test of the full conj:pself for n <= NMAX
def word(lam):
    w = ''; prev = 0
    for part in reversed(lam): w += 'E'*(part-prev) + 'N'; prev = part
    return w
def isP(lam):
    q = 0
    for x in word(lam): q = d(q, x)
    return acc[q]
def conj(lam):
    return tuple(sum(1 for p in lam if p >= j) for j in range(1, lam[0]+1)) if lam else ()
def is_stack(lam):               # stack of odd self-conjugate hooks and [1]
    lam = list(lam)
    while True:
        if lam == [1]: return True
        if not lam: return False
        c = lam[0]; r = len(lam)
        if c != r or c % 2 == 0: return False
        inner = [p-1 for p in lam[1:] if p > 1]
        if not inner: return False      # a bare hook is not a stack ending in [1]
        if inner[0] > c-1 - 0 and False: return False
        lam = inner
def partitions(m, mx=None):
    mx = mx or m
    if m == 0: yield (); return
    for k in range(min(m, mx), 0, -1):
        for rest in partitions(m-k, k): yield (k,) + rest
NMAX = int(sys.argv[1]); bad = []; cnt = 0
for m in range(1, NMAX+1):
    for lam in partitions(m):
        cj = conj(lam)
        if len(cj) != len(lam) or any((a-b) % 2 for a, b in zip(lam, cj)): continue
        cnt += 1
        if (not isP(lam)) != is_stack(lam): bad.append(lam)
print(f'    conj:pself tested on all {cnt} parity self-conjugate partitions of n <= {NMAX}: counterexamples {bad[:5]}')
# (4) four parts: lambda = [l1,l2,l3,l4], word E^{l4} N E^{l3-l4} N E^{l2-l3} N E^{l1-l2} N.
cls = lambda k: '0' if k == 0 else ('1' if k == 1 else ('e' if k % 2 == 0 else 'o'))      # 0, 1, even >= 2, odd >= 3
table = {}
for l4 in range(1, 9):
    for x in range(0, 9):
        for y in range(0, 9):
            for z in range(0, 9):
                lam = (l4+x+y+z, l4+x+y, l4+x, l4)
                key = (cls(l4), cls(x), cls(y), cls(z)); v = isP(lam)
                assert table.setdefault(key, v) == v, 'not determined by classes'
los = sorted(k for k, v in table.items() if v)
print('(4) four parts: losing iff (class(l4), class(l3-l4), class(l2-l3), class(l1-l2)) is one of the following',
      len(los), 'of', len(table), 'classes (0; 1; e = even >= 2; o = odd >= 3):')
print('    ', ' '.join(''.join(k) for k in los))
