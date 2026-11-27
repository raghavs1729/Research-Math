"""Multisection scan: does a target series, grouped in blocks of q = 1, 2, 3 consecutive terms, have the canonical
data (z^q, multisected Pochhammer classes) of  n -> F(n,0)  for  F(n,k) = Seed(p + A n; k + x + K n)  with an integer
pattern (A, K) and a base point (p, x) in (1/D) Z ?   (q > 1 corresponds to fractional-step patterns A/q: this is how
Chu's "bisection / trisection series" arise, e.g. T2 = Chu 2018 Ex. 19, and how Theorem 12 arises.)
usage: multisect_scan.py SEED BOX D targets.py out.jsonl [id1,id2,...]"""
import sys, itertools, random, json
import numpy as np
from fractions import Fraction as Fr
from seeds_extra import SEEDS, load_ebisu
from scan import z_of_shifts, sig_at, frac1
import targets as tg

seed, box, D, tpath, outpath = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
only = set(sys.argv[6].split(',')) if len(sys.argv) > 6 else None
ALL = dict(SEEDS)
params, T = ALL[seed]
r = len(params)


def msect(sig, q):
    z, fs = sig
    out = {}
    for al, m in fs:
        for j in range(q):
            h = frac1((al + j) / q)
            out[h] = out.get(h, 0) + m
    return z ** q, frozenset((al, m) for al, m in out.items() if m)


ns = tg.load_exact(tpath)
TARGETS = {}
for e in (ns.get('ENTRIES') or ns.get('E')):
    if only and e['id'] not in only:
        continue
    sig = tg.entry_signature(ns, e)
    if not (isinstance(sig, tuple) and len(sig) == 4):
        continue
    for q in (1, 2, 3):
        TARGETS[(e['id'], q)] = msect((sig[0], sig[1]), q)
rng = random.Random(1)
NCL = D * 840
W = np.array([rng.getrandbits(62) for _ in range(NCL)], dtype=np.int64)
cls_index = lambda fr: int(frac1(Fr(fr)) * NCL) % NCL


def target_hash(sig):
    h = 0
    for a, m in sig[1]:
        h += m * int(W[cls_index(a)])
    return np.int64(((h + 2**63) % 2**64) - 2**63)


grid = np.stack(np.meshgrid(*[np.arange(D)] * (r + 1), indexing='ij'), -1).reshape(-1, r + 1)
zt = {}
for key, (z, fs) in TARGETS.items():
    zt.setdefault(z, []).append(key)
out = open(outpath, 'a')
nh = 0
for sh in itertools.product(range(-box, box + 1), repeat=r + 1):
    shifts, K = sh[:r], sh[r]
    if all(s == 0 for s in sh):
        continue
    try:
        z = z_of_shifts(T, params, shifts, K)
    except ZeroDivisionError:
        continue
    if z is None or z not in zt:
        continue
    H = np.zeros(D ** (r + 1), dtype=np.int64)
    ok = True
    for L, pw in T.gammas:
        cn = sum(L.get(p, Fr(0)) * shifts[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * K
        if cn == 0:
            continue
        cn = int(cn)
        coef = [L.get(p, Fr(0)) for p in params] + [L.get('k', Fr(0))]
        c0 = L.get('1', Fr(0))
        if any((c * D).denominator != 1 for c in coef) or (c0 * D).denominator != 1:
            ok = False
            break
        ci = np.array([int(c) for c in coef], dtype=np.int64)
        l = (grid @ ci + int(c0 * D)) % D
        m = abs(cn)
        tab = np.zeros(D, dtype=np.int64)
        for lv in range(D):
            s = 0
            for j in range(m):
                s += int(W[cls_index((Fr(lv, D) + j) / m)]) if cn > 0 else int(W[cls_index((1 - Fr(lv, D) + j) / m)])
            sgn = pw if cn > 0 else -pw
            tab[lv] = np.int64(((sgn * s + 2**63) % 2**64) - 2**63)
        with np.errstate(over='ignore'):
            H = H + tab[l]
    if not ok:
        continue
    for key in zt[z]:
        idx = np.nonzero(H == target_hash(TARGETS[key]))[0]
        for i in idx[:6]:
            base = [Fr(int(v), D) for v in grid[i]]
            if sig_at(T, params, shifts, K, base[:r], base[r]) == TARGETS[key][1]:
                out.write(json.dumps([seed, list(shifts), K, [str(b) for b in base[:r]], str(base[r]), str(z), [key[0]], key[1]]) + '\n')
                nh += 1
        out.flush()
print(seed, 'box', box, 'hits', nh, flush=True)
