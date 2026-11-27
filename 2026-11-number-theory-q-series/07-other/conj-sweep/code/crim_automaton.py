"""CRIM (arXiv:2606.16828): learn a DFA B for  f(w) := [canon(w) is a P-position]  on all words w in {E,N}^*,
then PROVE that it is correct by an exact closure check (see R07).

canon(w): strip leading N's and trailing E's (the boundary word of the partition; '' = empty partition).
A move deletes one letter (a row = N, a column = E) and re-canonicalises.
"""
import itertools, sys, json
from crim_struct import isP, canon
SUF, PRE = 8, 8
def f(w): return isP(canon(w))
sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
def sig(u): return tuple(f(u+v) for v in sufs)
# --- learn: Nerode classes of prefixes (breadth-first, closed under appending a letter)
rep = {}; order = []
start = sig('')
rep[start] = ''; order.append(start)
trans = {}
i = 0
while i < len(order):
    s = order[i]; u = rep[s]; i += 1
    assert len(u) <= PRE, 'did not close'
    for x in 'EN':
        t = sig(u + x)
        if t not in rep: rep[t] = u + x; order.append(t)
        trans[(s, x)] = t
idx = {s: k for k, s in enumerate(order)}
n = len(order)
delta = {(idx[s], x): idx[t] for (s, x), t in trans.items()}
acc = [s[0] for s in order]          # sufs[0] == '' : f(u)
reps = [rep[s] for s in order]
print('states:', n)
for k in range(n):
    print(f'  state {k}: representative {reps[k]!r:12} accepting={acc[k]}  E->{delta[(k,"E")]}  N->{delta[(k,"N")]}')
json.dump({'n': n, 'delta': {f'{k}{x}': delta[(k, x)] for k in range(n) for x in 'EN'}, 'acc': acc, 'reps': reps},
          open('../cert/crim_dfa.json', 'w'))
# --- exact verification (no reference to game values below this line)
def run(w):
    q = 0
    for x in w: q = delta[(q, x)]
    return q
ok = True
# (i) leading N's are ignored; (iii) empty word accepted
ok &= (delta[(0, 'N')] == 0) and acc[0]
# (ii) trailing E's are ignored: acc(delta(q,E)) == acc(q) for all states
ok &= all(acc[delta[(q, 'E')]] == acc[q] for q in range(n))
print('(i) delta(q0,N)=q0:', delta[(0,'N')] == 0, ' (iii) acc(q0):', acc[0], ' (ii) acc(qE)=acc(q) for all q:', all(acc[delta[(q,'E')]] == acc[q] for q in range(n)))
# (iv) closure: for every word w that starts with E and ends with N:
#      acc(run(w))  <=>  for no position i is  w-with-letter-i-deleted  accepted.
# BFS over configurations (q, T, parity) with q = run(w), T = { run(w minus one letter) }, for w starting with E.
from collections import deque
init = (delta[(0, 'E')], frozenset([0]), 1)       # w = 'E': deleting it gives '' -> state 0 ; length parity 1
seen = {init}; dq = deque([init]); viol = []; oddP = []
while dq:
    q, T, par = dq.popleft()
    for x in 'EN':
        q2 = delta[(q, x)]; T2 = frozenset([delta[(t, x)] for t in T] + [q]); par2 = par ^ 1
        if x == 'N':                                # the word read so far starts with E and ends with N: canonical
            someP = any(acc[t] for t in T2)
            if acc[q2] == someP: viol.append((q, sorted(T), x))
            if acc[q2] and par2 == 1: oddP.append((q, sorted(T), x))
        c = (q2, T2, par2)
        if c not in seen: seen.add(c); dq.append(c)
print(f'(iv) configurations reached: {len(seen)}; violations of the closure property: {len(viol)}')
print(f'(v) accepted canonical words of odd length (configurations): {len(oddP)}')
print('VERDICT:', 'automaton PROVED correct; every P-position has even rank' if ok and not viol and not oddP else 'FAILED')
