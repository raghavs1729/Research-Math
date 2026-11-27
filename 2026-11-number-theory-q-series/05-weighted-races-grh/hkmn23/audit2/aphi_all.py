"""Second audit: Phi(q) vs (2/3)N(q,1) for every q <= Q (mpmath, own blocks), and the
comparison of every block c_f with the author's raw/t2_closed_400_*.jsonl."""
import sys, json, glob
import mpmath as mp
from sympy import divisors
from aphi import blocks, N1
Q = int(sys.argv[1])
auth = {}
for fn in glob.glob('../raw/t2_closed_400_*.jsonl'):
    for line in open(fn):
        r = json.loads(line); auth[r['f']] = mp.mpf(r['c_f'][0])
maxdiff = 0
for f in range(3, Q + 1):
    b = blocks(f)
    if f in auth:
        maxdiff = max(maxdiff, abs(b - auth[f]))
    elif b != 0:
        print('block missing in author data', f)
print('max |c_f(mine) - c_f(author)| over f <=', Q, '=', mp.nstr(maxdiff, 3))
cz = mp.log(2 * mp.pi) + mp.log(2) - mp.mpf(1) / 2 - mp.zeta(-1, 1, 1) / mp.zeta(-1)
fails = []; minrel = (10, 0)
for q in range(1, Q + 1):
    P = cz + sum(blocks(f) for f in divisors(q) if f > 2)
    thr = mp.mpf(2) / 3 * N1(q)
    if P <= thr:
        fails.append(q)
    else:
        rel = (P - thr) / thr
        if rel < minrel[0]:
            minrel = (rel, q)
print('q with Phi(q) <= (2/3)N(q,1):', fails)
print('smallest relative excess among the others:', mp.nstr(minrel[0], 5), 'at q =', minrel[1])
