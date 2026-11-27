"""Do T4 or T11 (or their bisections) occur in the Dougall 5F4 / 5H5 / very-well-poised Omega family
(the seed family of Chu-Zhang, Math. Comp. 83 (2014) and Chu, Ramanujan J. 55 (2021)) or in the
Dougall 7F6 family, for ANY integer shift pattern in a box and ANY base point with denominators | 24 ?

Family:  F(n,k) = Seed(p_1 + A_1 n, ..., p_r + A_r n ; k + x + K n),   x = continuous k-offset.
For Dougall5F4 the offset x turns the unilateral 5F4 term into the general bilateral 5H5 / Omega term
   (a+2x ; b+x, c+x, d+x, a+x)  =  (a' ; b', c', d', e'),
so integer patterns (A,B,C,D,K) <-> integer patterns of (a';b',c',d',e'), which are Chu's patterns.
The accelerated series sum_n G(n,0) has summand R(n,0) F(n,0): its canonical data (z and the multiset of
Pochhammer parameters mod 1) are those of n -> F(n,0).  We enumerate the patterns with the right z and then
ALL base points on the grid (1/24)Z^r (vectorised with a random hash of the signed multiset).

usage: chu_family_test.py SEED BOX [D] [G]     (G: targets = Guillera's shifted series)
"""
import sys, itertools, random
import numpy as np
from fractions import Fraction as Fr
from seeds_extra import SEEDS
from scan import z_of_shifts, sig_at, hyp_signature, frac1

D = int(sys.argv[3]) if len(sys.argv) > 3 else 24      # base-point denominators
TARGETS = {
    'T4':  hyp_signature(Fr(27, 32), [Fr(1, 6)] * 2 + [Fr(5, 6)] * 2 + [Fr(1, 2)], [Fr(1, 3), Fr(2, 3), 1, 1, 1]),
    'T11': hyp_signature(Fr(1, 64), [Fr(1, 4)] * 3 + [Fr(3, 4)] * 3, [Fr(1, 2), 1, 1, 1, 1, 1]),
    'T1':  hyp_signature(Fr(1, 4), [Fr(1, 3), Fr(2, 3), 1], [Fr(1, 6), Fr(1, 2), Fr(5, 6)]),
    'T3':  hyp_signature(Fr(1, 9), [1, 1], [Fr(1, 3), Fr(2, 3)]),
    'T7':  hyp_signature(Fr(4, 27), [Fr(1, 2), 1], [Fr(1, 6), Fr(5, 6)]),
    'T12': hyp_signature(Fr(-1, 48), [Fr(1, 4), Fr(3, 4), 1], [Fr(1, 2)] * 3),
    'T2_known': hyp_signature(Fr(1, 9), [Fr(1, 4), Fr(3, 4), 1], [Fr(1, 2)] * 3),   # control: Chu 2018, Ex. 19 (trisected)
    # positive controls: series known to come from Dougall's 5F4 (Guillera; Chu-Zhang 2014)
    'ctl_Guillera_-1/4': hyp_signature(Fr(-1, 4), [Fr(1, 2)] * 5, [1] * 5),
    'ctl_Guillera_1/16': hyp_signature(Fr(1, 16), [Fr(1, 2)] * 3 + [Fr(1, 4), Fr(3, 4)], [1] * 5),
    'ctl_Guillera_27/64': hyp_signature(Fr(27, 64), [Fr(1, 2)] * 3 + [Fr(1, 3), Fr(2, 3)], [1] * 5),
    # actual examples of Chu-Zhang (Math. Comp. 83 (2014)), as quoted by Guillera, arXiv:1611.04385:
    'ctl_ChuZhang_Ex42': hyp_signature(Fr(-16, 27), [1] * 5, [Fr(1, 2), Fr(2, 3), Fr(1, 3), Fr(3, 4), Fr(1, 4)]),
    'ctl_ChuZhang_Ex45': hyp_signature(Fr(-16, 27), [Fr(1, 2)] + [Fr(1, 4)] * 3 + [Fr(3, 4)] * 3, [Fr(1, 3), Fr(2, 3)] + [1] * 5),
    'ctl_Ramanujan_1/4': hyp_signature(Fr(1, 4), [Fr(1, 2)] * 3, [1] * 3),
    'ctl_Ramanujan_-1/4': hyp_signature(Fr(-1, 4), [Fr(1, 2)] * 3, [1] * 3),
    'ctl_Ramanujan_1/64': hyp_signature(Fr(1, 64), [Fr(1, 2)] * 3, [1] * 3),
}


def bisect(sig):
    z, fs = sig
    out = {}
    for a, m in fs:
        for half in (a / 2, (a + 1) / 2):
            h = frac1(half)
            out[h] = out.get(h, 0) + m
    return z * z, frozenset((a, m) for a, m in out.items() if m)


if len(sys.argv) > 4 and sys.argv[4] == 'G':
    # Guillera's half-shifted Ramanujan series (arXiv:1802.02023, Tables 1-6); G8 = our Theorem 2
    sys.path.insert(0, '.')
    import targets as tg
    ns = tg.load_exact('../survey/verify_G.py')
    TARGETS = {}
    for e in ns['ENTRIES']:
        sig = tg.entry_signature(ns, e)
        TARGETS[e['id']] = (sig[0], sig[1])
    TARGETS['ctl_Guillera_1/16'] = hyp_signature(Fr(1, 16), [Fr(1, 2)] * 3 + [Fr(1, 4), Fr(3, 4)], [1] * 5)
    # trisections (a pattern with step 3 in n would produce the series grouped in blocks of three terms)
    def trisect(sig):
        z, fs = sig
        out = {}
        for al, m in fs:
            for j in range(3):
                h = frac1((al + j) / 3)
                out[h] = out.get(h, 0) + m
        return z ** 3, frozenset((al, m) for al, m in out.items() if m)
    for nm in ('G3', 'G8', 'G13'):
        TARGETS[nm + '_trisected'] = trisect(TARGETS[nm])
        TARGETS[nm + '_bisected'] = bisect(TARGETS[nm])
def trisect3(sig):
    z, fs = sig
    out = {}
    for al, m in fs:
        for j in range(3):
            h = frac1((al + j) / 3)
            out[h] = out.get(h, 0) + m
    return z ** 3, frozenset((al, m) for al, m in out.items() if m)


for nm in ('T4', 'T11', 'T1', 'T3', 'T7', 'T12', 'T2_known'):
    if nm not in TARGETS:
        continue
    TARGETS[nm + '_trisected'] = trisect3(TARGETS[nm])
    TARGETS[nm + '_bisected'] = bisect(TARGETS[nm])

seed, box = sys.argv[1], int(sys.argv[2])
params, T = SEEDS[seed]
r = len(params)
rng = random.Random(1)
M = 2 * 3 * 4 * 5 * 6 * 7 * 8          # classes live in (1/(D*M)) Z / Z, enough for |m| <= 8
W = np.array([rng.getrandbits(62) for _ in range(D * 840)], dtype=np.int64)
NCL = D * 840


def cls_index(fr):
    f = frac1(Fr(fr))
    return int(f * NCL) % NCL


def target_hash(sig):
    h = 0
    for a, m in sig[1]:
        h += m * int(W[cls_index(a)])
    return np.int64(((h + 2**63) % 2**64) - 2**63)


# grid of base points: params and the k-offset x, each in (1/D)Z/Z
grid = np.stack(np.meshgrid(*[np.arange(D)] * (r + 1), indexing='ij'), -1).reshape(-1, r + 1)
zt = {}
for nm, (z, fs) in TARGETS.items():
    zt.setdefault(z, []).append(nm)

count = {nm: 0 for nm in TARGETS}
hits = []
npat = 0
for sh in itertools.product(range(-box, box + 1), repeat=r + 1):
    shifts, K = sh[:r], sh[r]
    if K < 0 or all(s == 0 for s in sh):
        continue
    z = z_of_shifts(T, params, shifts, K)
    if z is None or z not in zt:
        continue
    npat += 1
    # factor list: (cn, linear form in (base, x) with integer-ised coefficients (units 1/D), power)
    fac = []
    for L, pw in T.gammas:
        cn = sum(L.get(p, Fr(0)) * shifts[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * K
        if cn == 0:
            continue
        coef = [L.get(p, Fr(0)) for p in params] + [L.get('k', Fr(0))]
        c0 = L.get('1', Fr(0))
        fac.append((int(cn), coef, c0, pw))
    # enumerate base points: only the distinct values of the forms matter; do it on the grid in chunks
    H = np.zeros(D ** (r + 1), dtype=np.int64)
    ok = True
    for cn, coef, c0, pw in fac:
        if any((c * D).denominator != 1 for c in coef) or (c0 * D).denominator != 1:
            ok = False
            break
        ci = np.array([int(c) for c in coef], dtype=np.int64)
        l = (grid @ ci + int(c0 * D)) % D        # l/D = L mod 1
        m = abs(cn)
        tab = np.zeros(D, dtype=np.int64)
        for lv in range(D):
            s = 0
            for j in range(m):
                if cn > 0:
                    s += int(W[cls_index((Fr(lv, D) + j) / m)])
                else:
                    s += int(W[cls_index((1 - Fr(lv, D) + j) / m)])
            sgn = pw if cn > 0 else -pw
            tab[lv] = np.int64(((sgn * s + 2**63) % 2**64) - 2**63)
        with np.errstate(over='ignore'):
            H = H + tab[l]
    if not ok:
        print('skip pattern (non-grid coefficients)', sh)
        continue
    for nm in zt[z]:
        idx = np.nonzero(H == target_hash(TARGETS[nm]))[0]
        for i in idx[:3]:
            base = [Fr(int(v), D) for v in grid[i]]
            sg = sig_at(T, params, shifts, K, base[:r], base[r])
            if sg == TARGETS[nm][1]:
                count[nm] += len(idx)
                hits.append((nm, sh, base))
                # Chu's iteration pattern [da; db, ..., de] of the bilateral / Omega parameters (a+2x; b+x, ..., a+x)
                delta = [shifts[0] + 2 * K] + [shifts[i] + K for i in range(1, r)] + [shifts[0] + K]
                print('HIT', nm, 'shifts', sh, 'base', [str(b) for b in base], 'n_base_points', len(idx),
                      'chu_pattern', delta, 'nonneg' if min(delta) >= 0 else 'has_negative_entry', flush=True)
                break
print(seed, 'box', box, 'patterns with a target z:', npat, 'hits per target:', count)
