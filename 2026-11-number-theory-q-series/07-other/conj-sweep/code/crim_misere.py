# CRIM misere play: is the 3-valued function  m(w) = min(G^-(canon(w)), 2)  regular?  (G^- = misere Grundy value, G^-(empty) = 1)
import itertools, sys
from crim_sg import G, canon
SUF, PRE = int(sys.argv[1]), int(sys.argv[2])
f = lambda w: min(G(canon(w), True), 2)
sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
sig = lambda u: tuple(f(u+v) for v in sufs)
rep = {}; order = []; s0 = sig(''); rep[s0] = ''; order.append(s0); trans = {}
i = 0
while i < len(order):
    s = order[i]; u = rep[s]; i += 1
    if len(u) > PRE: print('did not close; classes so far', len(order)); sys.exit()
    for x in 'EN':
        t = sig(u + x)
        if t not in rep: rep[t] = u + x; order.append(t)
        trans[(s, x)] = t
print('closed with', len(order), 'states')
import json
idx = {s: k for k, s in enumerate(order)}
json.dump({'n': len(order), 'delta': {f'{idx[s]}{x}': idx[t] for (s, x), t in trans.items()}, 'val': [s[0] for s in order], 'reps': [rep[s] for s in order]}, open('../cert/crim_misere_dfa.json', 'w'))
