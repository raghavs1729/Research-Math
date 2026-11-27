"""Stage 2: for each signature hit, try integer offsets of the base point, compute the WZ mate and compare
the term ratio of G(n,0) with that of the target summand (exact rational functions).
usage: python3 match_hits.py hits.jsonl verify_X.py [ids...]"""
import sys, json, itertools, time
from fractions import Fraction as Fr
from multiprocessing import Pool
from wz import *
from seeds_extra import SEEDS, load_ebisu
import targets as tg

ALL = dict(SEEDS)
ALL.update(load_ebisu())
hits_path, tpath = sys.argv[1], sys.argv[2]
only = set(sys.argv[3:])
ns = tg.load_exact(tpath)
TR = {}
for e in (ns.get('ENTRIES') or ns.get('E')):
    if only and e['id'] not in only:
        continue
    sig = tg.entry_signature(ns, e)
    if isinstance(sig, tuple) and len(sig) == 4:
        TR[e['id']] = pari.subst(sig[3], 'k', pari('n'))


def make_F(name, shifts, K, base, k0):
    params, T = ALL[name]
    m = {}
    for i, p in enumerate(params):
        m[p] = lf_add({'1': Fr(base[i])}, {'n': Fr(shifts[i])})
    m['k'] = lf_add({'k': Fr(1), 'n': Fr(K)}, {'1': Fr(k0)})
    return T.subs(m)


def g_ratio(F, R):
    n = pari('n')
    R0 = pari.subst(R, 'k', 0)
    if R0 == 0:
        return None
    Fr0 = pari.subst(F.ratio('n'), 'k', 0)
    return pari.subst(R0, 'n', n + 1) / R0 * Fr0


def try_hit(h):
    name, shifts, K, base, k0, z, ids = h
    ids = [i for i in ids if i in TR]
    if not ids:
        return []
    base = [Fr(b) for b in base]
    k0 = Fr(k0)
    m = len(base)
    out = []
    rng = (0, 1, -1, 2) if m <= 2 else (0, 1)
    n = pari('n')
    for off in itertools.product(rng, repeat=m):
        for koff in (0, 1, 2):
            try:
                F = make_F(name, shifts, K, [b + o for b, o in zip(base, off)], k0 + koff)
                R = wz_mate(F)
                gr = g_ratio(F, R)
            except (NoMate, ZeroDivisionError, AssertionError):
                continue
            except Exception as ex:
                continue
            if gr is None:
                continue
            for i in ids:
                for j in range(-3, 4):
                    if gr == pari.subst(TR[i], 'n', n + j):
                        out.append(dict(id=i, seed=name, shifts=shifts, K=K, base=[str(b + o) for b, o in zip(base, off)],
                                        k0=str(k0 + koff), j=j))
    return out


if __name__ == '__main__':
    hits = [json.loads(l) for l in open(hits_path)]
    if only:
        hits = [h for h in hits if set(h[6]) & only]
    print(len(hits), 'hits to test', file=sys.stderr)
    t0 = time.time()
    seen = set()
    with Pool(4) as pool:
        for res in pool.imap_unordered(try_hit, hits, chunksize=1):
            for r in res:
                key = json.dumps(r, sort_keys=True)
                if key not in seen:
                    seen.add(key)
                    print(key, flush=True)
    print('done', time.time() - t0, file=sys.stderr)
