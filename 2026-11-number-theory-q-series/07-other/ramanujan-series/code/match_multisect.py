"""Stage 2 for multisect_scan.py: exact comparison of the summand of G(m,0) with the target grouped in blocks of q terms.
usage: python3 match_multisect.py hits.jsonl targets.py [ids...]"""
import sys, json, itertools, time
from fractions import Fraction as Fr
from multiprocessing import Pool
from wz import *
from seeds_extra import SEEDS, load_ebisu
import targets as tg

ALL = dict(SEEDS)
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


def block_ratio(r, q, j0):
    """ratio T(m+1)/T(m), T(m) = sum_{i<q} t(q m + j0 + i), as a rational function of n (= m)"""
    n = pari('n')
    rr = [pari.subst(r, 'n', q * n + j0 + i) for i in range(2 * q)]
    def B(off):
        s, p = pari(1), pari(1)
        for i in range(q - 1):
            p = p * rr[off + i]
            s = s + p
        return s
    prod = pari(1)
    for i in range(q):
        prod = prod * rr[i]
    return prod * B(q) / B(0)


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
    return pari.subst(R0, 'n', n + 1) / R0 * pari.subst(F.ratio('n'), 'k', 0)


def try_hit(h):
    name, shifts, K, base, k0, z, ids, q = h
    ids = [i for i in ids if i in TR]
    if not ids:
        return []
    base = [Fr(b) for b in base]; k0 = Fr(k0)
    out = []
    n = pari('n')
    cands = {}
    for i in ids:
        for j0 in range(q):
            br = block_ratio(TR[i], q, j0)
            for j in range(-2, 3):
                cands[(i, j0, j)] = pari.subst(br, 'n', n + j)
    for off in itertools.product((0, 1), repeat=len(base)):
        for koff in (0, 1):
            try:
                F = make_F(name, shifts, K, [b + o for b, o in zip(base, off)], k0 + koff)
                R = wz_mate(F)
                gr = g_ratio(F, R)
            except Exception:
                continue
            if gr is None:
                continue
            for key, cr in cands.items():
                if gr == cr:
                    out.append(dict(id=key[0], q=q, j0=key[1], j=key[2], seed=name, shifts=shifts, K=K,
                                    base=[str(b + o) for b, o in zip(base, off)], k0=str(k0 + koff)))
    return out


if __name__ == '__main__':
    hits = [json.loads(l) for l in open(hits_path)]
    if only:
        hits = [h for h in hits if set(h[6]) & only]
    # dedupe
    seen_h = set(); H = []
    for h in hits:
        key = json.dumps(h)
        if key not in seen_h:
            seen_h.add(key); H.append(h)
    print(len(H), 'hits to test', file=sys.stderr)
    t0 = time.time(); seen = set()
    with Pool(4) as pool:
        for res in pool.imap_unordered(try_hit, H, chunksize=1):
            for r in res:
                key = json.dumps(r, sort_keys=True)
                if key not in seen:
                    seen.add(key); print(key, flush=True)
    print('done', time.time() - t0, file=sys.stderr)
