"""CRIM under misere play (the player who cannot move wins): DFA for the misere losing positions, with exact proof.
   f(w) = [canon(w) is a misere P-position];  the empty position is a misere N-position.
   Closure property proved: (i) delta(q0,N) = q0; (ii) acc(delta(q,E)) = acc(q); (iii) q0 not accepting;
   (iv) for every word w starting with E and ending with N:  acc(w) <=> no one-letter deletion of w is accepted.
"""
import itertools, json
from collections import deque
from crim_sg import G, canon
SUF, PRE = 8, 10
f = lambda w: G(canon(w), True) == 0
sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
sig = lambda u: tuple(f(u+v) for v in sufs)
rep = {}; order = []; s0 = sig(''); rep[s0] = ''; order.append(s0); trans = {}
i = 0
while i < len(order):
    s = order[i]; u = rep[s]; i += 1
    assert len(u) <= PRE
    for x in 'EN':
        t = sig(u + x)
        if t not in rep: rep[t] = u + x; order.append(t)
        trans[(s, x)] = t
idx = {s: k for k, s in enumerate(order)}; n = len(order)
dE = [idx[trans[(s, 'E')]] for s in order]; dN = [idx[trans[(s, 'N')]] for s in order]; acc = [s[0] for s in order]
print('states:', n)
for k in range(n): print(f'  state {k}: representative {rep[order[k]]!r:10} accepting={acc[k]}  E->{dE[k]}  N->{dN[k]}')
json.dump({'n': n, 'dE': dE, 'dN': dN, 'acc': acc, 'reps': [rep[s] for s in order]}, open('../cert/crim_misere_dfa.json', 'w'))
d = lambda q, x: dE[q] if x == 'E' else dN[q]
c1 = dN[0] == 0; c2 = all(acc[dE[q]] == acc[q] for q in range(n)); c3 = not acc[0]
init = (dE[0], frozenset([0])); seen = {init}; dq = deque([init]); viol = 0
while dq:
    q, T = dq.popleft()
    for x in 'EN':
        q2 = d(q, x); T2 = frozenset([d(t, x) for t in T] + [q])
        if x == 'N' and acc[q2] == any(acc[t] for t in T2): viol += 1
        c = (q2, T2)
        if c not in seen: seen.add(c); dq.append(c)
print(f'(i) {c1} (ii) {c2} (iii) empty word rejected: {c3}; (iv) configurations {len(seen)}, violations {viol}')
print('VERDICT:', 'misere automaton PROVED correct' if c1 and c2 and c3 and viol == 0 else 'FAILED')
# joint outcome classes with the normal-play automaton of R07
B = json.load(open('../cert/crim_dfa.json'))
bd = lambda q, x: B['delta'][f'{q}{x}']
# count canonical words of each length by (normal P?, misere P?)
from collections import Counter
cur = Counter({(bd(0, 'E'), dE[0]): 1}); rows = []
for L in range(2, 21):
    nxt = Counter(); fin = Counter()
    for (p, q), m in cur.items():
        fin[(B['acc'][bd(p, 'N')], acc[dN[q]])] += m
        for x in 'EN': nxt[(bd(p, x), d(q, x))] += m
    rows.append((L, fin[(True, True)], fin[(True, False)], fin[(False, True)], fin[(False, False)])); cur = nxt
print('L: #(normal P & misere P), #(normal P only), #(misere P only), #(neither), over the 2^(L-2) partitions with rows+columns = L')
for r in rows: print('  ', r)
