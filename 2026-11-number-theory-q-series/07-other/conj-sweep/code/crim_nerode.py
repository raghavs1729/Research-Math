# Is the set of P-words of CRIM a regular language?  Estimate the number of Nerode classes of prefixes.
from crim_struct import *
import itertools, sys
SUF = int(sys.argv[1]); PREMAX = int(sys.argv[2])
sufs = []
for l in range(0, SUF+1):
    for s in itertools.product('EN', repeat=l):
        s = ''.join(s)
        if s.endswith('N'): sufs.append(s)
def sig(u):
    return tuple(isP(u+v) for v in sufs)
classes = {}
for l in range(1, PREMAX+1):
    new = 0
    for mid in itertools.product('EN', repeat=l-1):
        u = 'E' + ''.join(mid)
        s = sig(u)
        if s not in classes: classes[s] = u; new += 1
    print(f'prefix length {l}: new classes {new}, total {len(classes)}', flush=True)
