"""Scan seeds x shifts x base points for signature matches with the survey targets.
usage: python3 run_scan.py <box N for multi-param seeds> <box for 1-2 param seeds> verify_X.py ..."""
import sys, itertools, json, time, os
from fractions import Fraction as Fr
from multiprocessing import Pool
from scan import *
from seeds_extra import SEEDS, load_ebisu
import targets as tg

NBIG, NSMALL = int(sys.argv[1]), int(sys.argv[2])
ALL = dict(SEEDS)
ALL.update(load_ebisu())

# ---- targets
TSIG = {}   # (z, classes) -> [ids]
for path in sys.argv[3:]:
    ns = tg.load_exact(path)
    for e in (ns.get('ENTRIES') or ns.get('E')):
        sig = tg.entry_signature(ns, e)
        if isinstance(sig, tuple) and len(sig) == 4 and isinstance(sig[1], frozenset):
            TSIG.setdefault((sig[0], sig[1]), []).append(e['id'])
        else:
            print('no signature for', e['id'], sig, file=sys.stderr)
TZ = {}
for (z, cl), ids in TSIG.items():
    TZ.setdefault(z, []).append(cl)
print(len(TSIG), 'distinct target signatures;', len(TZ), 'distinct z', file=sys.stderr)


def den_lcm(cls):
    from math import lcm
    L = 1
    for cl in cls:
        for a, m in cl:
            L = lcm(L, a.denominator)
    return L


def work(name):
    params, T = ALL[name]
    m = len(params)
    N = NBIG if m >= 4 else (NSMALL if m <= 2 else 2)
    hits = []
    rng = range(-N, N + 1)
    for sh in itertools.product(rng, repeat=m + 1):
        shifts, K = sh[:m], sh[m]
        if all(s == 0 for s in sh):
            continue
        try:
            z = z_of_shifts(T, params, shifts, K)
        except ZeroDivisionError:
            continue
        if z is None or z not in TZ:
            continue
        cls = TZ[z]
        L = den_lcm(cls)
        Ls = sorted({L, 2 * L} if m <= 2 else {L})
        want = set(cls)
        for LL in Ls:
            pts = [Fr(i, LL) for i in range(LL)]
            for bp in itertools.product(pts, repeat=m + 1):
                base, k0 = bp[:m], bp[m]
                s = sig_at(T, params, shifts, K, base, k0)
                if s in want:
                    hits.append((name, shifts, K, [str(b) for b in base], str(k0), str(z), TSIG[(z, s)]))
    return name, hits


if __name__ == '__main__':
    t0 = time.time()
    names = sorted(ALL, key=lambda nm: len(ALL[nm][0]))
    out = open('scan_hits.jsonl', 'a')
    with Pool(4) as pool:
        for name, hits in pool.imap_unordered(work, names):
            byid = {}
            for h in hits:
                for i in h[6]:
                    byid.setdefault(i, []).append(h)
            if hits:
                print('%-50s %d hits: %s  [%.0fs]' % (name[:50], len(hits), sorted(byid), time.time() - t0), flush=True)
            for h in hits:
                out.write(json.dumps(h) + '\n')
            out.flush()
    print('done', time.time() - t0)
