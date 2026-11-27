"""CRIM (arXiv:2606.16828): Moore machines for the three-valued functions
     v+(w) = min(G+(canon w), 2)   (normal-play Grundy value, capped at 2)
     v-(w) = min(G-(canon w), 2)   (misere Grundy value: G-(empty) = 1, otherwise mex of the options)
learned from game values of short words and then PROVED correct by an exact closure check:
  (i)   delta(q0, N) = q0;  (ii) val(delta(q, E)) = val(q) for all q;  (iii) val(q0) = value of the empty position;
  (iv)  for every word w starting with E and ending with N, with V = { val of (w minus one letter) }:
        val(w) = 0 if 0 not in V;  = 1 if 0 in V and 1 not in V;  = 2 if {0,1} subset of V.
By induction on the length this gives val(w) = min(mex, 2) of the true Grundy value for all canonical w.
"""
import itertools, sys, json
from collections import deque
from crim_sg import G, canon
def learn(misere, SUF, PRE):
    f = lambda w: min(G(canon(w), misere), 2)
    sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
    sig = lambda u: tuple(f(u+v) for v in sufs)
    rep = {}; order = []; s0 = sig(''); rep[s0] = ''; order.append(s0); trans = {}
    i = 0
    while i < len(order):
        s = order[i]; u = rep[s]; i += 1
        assert len(u) <= PRE, 'did not close'
        for x in 'EN':
            t = sig(u + x)
            if t not in rep: rep[t] = u + x; order.append(t)
            trans[(s, x)] = t
    idx = {s: k for k, s in enumerate(order)}
    return {'n': len(order), 'dE': [idx[trans[(s, 'E')]] for s in order], 'dN': [idx[trans[(s, 'N')]] for s in order],
            'val': [s[0] for s in order], 'reps': [rep[s] for s in order]}
def verify(M, empty_value):
    n, dE, dN, val = M['n'], M['dE'], M['dN'], M['val']
    d = lambda q, x: dE[q] if x == 'E' else dN[q]
    c1 = dN[0] == 0; c2 = all(val[dE[q]] == val[q] for q in range(n)); c3 = val[0] == empty_value
    init = (dE[0], frozenset([0])); seen = {init}; dq = deque([init]); viol = 0
    while dq:
        q, T = dq.popleft()
        for x in 'EN':
            q2 = d(q, x); T2 = frozenset([d(t, x) for t in T] + [q])
            if x == 'N':
                V = {val[t] for t in T2}
                exp = 0 if 0 not in V else (1 if 1 not in V else 2)
                if val[q2] != exp: viol += 1
            c = (q2, T2)
            if c not in seen: seen.add(c); dq.append(c)
    return c1, c2, c3, len(seen), viol
if __name__ == '__main__':
    out = {}
    jobs = [('misere', True, 1, 8, 13)] + ([('normal', False, 0, 8, int(sys.argv[2]))] if len(sys.argv) > 2 else [])
    for name, mis, emp, SUF, PRE in jobs:
        try: M = learn(mis, SUF, PRE)
        except AssertionError: print(name, ': did not close with prefixes of length <=', PRE); continue
        c1, c2, c3, nconf, viol = verify(M, emp)
        print(f'{name}: {M["n"]} states (longest representative {max(len(r) for r in M["reps"])}); (i) {c1} (ii) {c2} (iii) {c3}; '
              f'(iv) {nconf} configurations, {viol} violations ->', 'PROVED' if c1 and c2 and c3 and viol == 0 else 'FAILED', flush=True)
        out[name] = M
    json.dump(out, open(sys.argv[1], 'w'))
