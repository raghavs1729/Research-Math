"""Numerical exploration of the four terms of Au's Prop. 2.1 for a matched configuration."""
import sys, json, mpmath as mp
from fractions import Fraction as Fr
from wz import *
from seeds_extra import SEEDS, load_ebisu
from numerics import *
ALL = dict(SEEDS); ALL.update(load_ebisu())

def make_F(name, shifts, K, base, k0):
    params, T = ALL[name]
    m = {}
    for i, p in enumerate(params):
        m[p] = lf_add(lf_parse(str(base[i])), {'n': Fr(shifts[i])})
    m['k'] = lf_add({'k': Fr(1), 'n': Fr(K)}, lf_parse(str(k0)))
    return T.subs(m)

def explore(name, shifts, K, base, k0, dps=40, NN=(20, 40, 80), KK=(2000, 8000), verbose=True):
    mp.mp.dps = dps
    full = [n for n in ALL if n == name or n.startswith(name)][0]
    F = make_F(full, shifts, K, base, k0)
    R = wz_mate(F); assert check_wz(F, R)
    W = WZNum(F, R)
    eps = mp.mpf(10) ** (-dps // 2)
    out = {}
    def safe(f, *a):
        try:
            return f(*a)
        except Exception as ex:
            return mp.nan
    # sum_k F(0,k)
    out['SF0'] = safe(lambda: mp.nsum(lambda k: W.Fv(eps, k), [0, mp.inf]))
    out['SG0'] = safe(lambda: mp.nsum(lambda n: W.Gv(n + eps, 0), [0, mp.inf]))
    out['G(0..3,0)'] = [safe(W.Gv, n + eps, 0) for n in range(4)]
    out['g(0..3)'] = [[safe(W.Gv, n + eps, Kb) for Kb in KK] for n in range(3)]
    out['SFN'] = [safe(lambda N=N: mp.nsum(lambda k: W.Fv(N + eps, k), [0, mp.inf])) for N in NN]
    if verbose:
        print(full, shifts, K, base, k0)
        print('  R =', R if len(str(R)) < 300 else str(R)[:300] + '...')
        for key, v in out.items():
            print('  ', key, '=', v if not isinstance(v, list) else [mp.nstr(x, 15) if not isinstance(x, list) else [mp.nstr(y, 12) for y in x] for x in v])
    return F, R, W, out

if __name__ == '__main__':
    seen = set()
    ids = sys.argv[2:]
    for l in open(sys.argv[1]):
        if not l.startswith('{'): continue
        r = json.loads(l)
        if r['id'] not in ids: continue
        key = (r['seed'], tuple(r['shifts']), r['K'], tuple(r['base']), r['k0'])
        if key in seen: continue
        seen.add(key)
        print('==== target', r['id'], 'j =', r['j'])
        try:
            explore(r['seed'], r['shifts'], r['K'], r['base'], r['k0'])
        except Exception as ex:
            print('   failed', repr(ex)[:200])
