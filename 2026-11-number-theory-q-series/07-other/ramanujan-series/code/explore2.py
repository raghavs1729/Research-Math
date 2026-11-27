import sys, json, mpmath as mp
from explore import *
from numerics import hsum_k, hsum_n
def expl(name, shifts, K, base, k0, dps=40, epss=('1e-8',)):
    mp.mp.dps = dps
    full = [n for n in ALL if n == name or n.startswith(name)][0]
    F = make_F(full, shifts, K, base, k0)
    R = wz_mate(F); assert check_wz(F, R)
    W = WZNum(F, R)
    print(full, shifts, K, base, k0); print('  F =', F); print('  R =', str(R)[:400])
    for e in epss:
        e = mp.mpf(e)
        def safe(f, *a):
            try: return f(*a)
            except Exception as ex: return mp.nan
        SF = safe(hsum_k, W, e)
        SG = safe(hsum_n, W, e)
        g = [safe(W.Gv, n + e, 20000) for n in range(3)]
        SFN = [safe(hsum_k, W, N + e) for N in (5, 20, 60)]
        print('   eps', mp.nstr(e, 3), 'SF', mp.nstr(SF, 20), 'SG', mp.nstr(SG, 20), 'SG/eps', mp.nstr(SG / e, 20))
        print('      g', [mp.nstr(x, 8) for x in g], 'SFN', [mp.nstr(x, 20) for x in SFN], 'SF-SG-SFN', mp.nstr(SF - SG - SFN[-1], 10))
    return F, R, W
if __name__ == '__main__':
    seen=set()
    for l in open(sys.argv[1]):
        if not l.startswith('{'): continue
        r = json.loads(l)
        if r['id'] not in sys.argv[2:]: continue
        key = (r['seed'], tuple(r['shifts']), r['K'], tuple(r['base']), r['k0'])
        if key in seen: continue
        seen.add(key)
        print('==== target', r['id'], 'j =', r['j'])
        try: expl(r['seed'], r['shifts'], r['K'], r['base'], r['k0'])
        except Exception as ex: print('   failed', repr(ex)[:300])
