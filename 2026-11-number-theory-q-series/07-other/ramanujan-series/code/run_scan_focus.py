"""focused deeper scan: python3 run_scan_focus.py verify.py N seed1,seed2,... id1,id2,..."""
import sys, itertools, json, time
from fractions import Fraction as Fr
from multiprocessing import Pool
from scan import *
from seeds_extra import SEEDS, load_ebisu
import targets as tg
path, N, seeds, ids = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(','), set(sys.argv[4].split(','))
ALL = dict(SEEDS)
TSIG = {}
ns = tg.load_exact(path)
for e in (ns.get('ENTRIES') or ns.get('E')):
    if e['id'] not in ids: continue
    sig = tg.entry_signature(ns, e)
    if isinstance(sig, tuple) and len(sig) == 4:
        TSIG.setdefault((sig[0], sig[1]), []).append(e['id'])
TZ = {}
for (z, cl), i in TSIG.items(): TZ.setdefault(z, []).append(cl)
def den_lcm(cls):
    from math import lcm
    L = 1
    for cl in cls:
        for a, m in cl: L = lcm(L, a.denominator)
    return L
def work(arg):
    name, sh = arg
    params, T = ALL[name]; m = len(params)
    shifts, K = sh[:m], sh[m]
    try: z = z_of_shifts(T, params, shifts, K)
    except ZeroDivisionError: return []
    if z is None or z not in TZ: return []
    cls = TZ[z]; L = den_lcm(cls); want = set(cls); hits = []
    pts = [Fr(i, L) for i in range(L)]
    for bp in itertools.product(pts, repeat=m + 1):
        s = sig_at(T, params, shifts, K, bp[:m], bp[m])
        if s in want:
            hits.append((name, shifts, K, [str(b) for b in bp[:m]], str(bp[m]), str(z), TSIG[(z, s)]))
    return hits
if __name__ == '__main__':
    tasks = []
    for name in seeds:
        m = len(ALL[name][0])
        for sh in itertools.product(range(-N, N + 1), repeat=m + 1):
            if any(sh): tasks.append((name, sh))
    print(len(tasks), 'tasks', file=sys.stderr)
    out = open('scan_focus_hits.jsonl', 'a'); t0 = time.time(); nh = 0
    with Pool(4) as pool:
        for hits in pool.imap_unordered(work, tasks, chunksize=50):
            for h in hits: out.write(json.dumps(h) + '\n'); nh += 1
            out.flush()
    print('done', nh, 'hits', time.time() - t0, file=sys.stderr)
