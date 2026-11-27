# CRIM misere play: Nerode classes of (a) the misere P-set {G^- = 0}; (b) the pair (G^+ = 0, G^- = 1) indicator.
import itertools, sys
from crim_sg import G, canon
SUF, PRE, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
if mode == 'P': f = lambda w: G(canon(w), True) == 0
if mode == 'C01': f = lambda w: (G(canon(w), False) == 0 and G(canon(w), True) == 1)
if mode == 'V': f = lambda w: min(G(canon(w), True), 2)
sufs = [''.join(s) for l in range(0, SUF+1) for s in itertools.product('EN', repeat=l)]
sig = lambda u: tuple(f(u+v) for v in sufs)
rep = {}; order = []; s0 = sig(''); rep[s0] = ''; order.append(s0)
i = 0; maxlen = 0
while i < len(order):
    s = order[i]; u = rep[s]; i += 1
    if len(u) > PRE: print(mode, 'did not close; classes so far', len(order), 'max rep length', len(u)); sys.exit()
    for x in 'EN':
        t = sig(u + x)
        if t not in rep: rep[t] = u + x; order.append(t)
print(mode, 'closed with', len(order), 'states; longest representative', max(len(rep[s]) for s in order))
