"""Explore all exact matches of one Guillera target: python3 g_explore.py G3 [seedprefix]"""
import sys, json, mpmath as mp
from explore import *
mp.mp.dps = 30
tid = sys.argv[1]; pref = sys.argv[2] if len(sys.argv) > 2 else ''
recs = []
for l in open('gscan/match_G.out.txt'):
    if l.startswith('{'):
        r = json.loads(l)
        if r['id'] == tid and r['seed'].startswith(pref): recs.append(r)
print(len(recs), 'records')
eps = mp.mpf(10) ** -15
for r in recs:
    try:
        F = make_F(r['seed'], r['shifts'], r['K'], r['base'], r['k0'])
        R = wz_mate(F); assert check_wz(F, R)
        W = WZNum(F, R)
        def safe(f):
            try: return f()
            except Exception as ex: return mp.nan
        SF0 = safe(lambda: hsum_k(W, eps))
        SG0 = safe(lambda: hsum_n(W, eps, 0))
        F00 = safe(lambda: W.Fv(eps, 0))
        SFN = safe(lambda: hsum_k(W, 12 + eps))
        gK = [safe(lambda n=n: W.Gv(n + eps, 4000)) for n in (0, 1)]
        print(r['seed'][:12], r['shifts'], r['K'], r['base'], r['k0'], 'j', r['j'], '| F(0,0)', mp.nstr(F00, 8), 'SF0', mp.nstr(SF0, 12),
              'SG0', mp.nstr(SG0, 12), 'SF12', mp.nstr(SFN, 5), 'G(0,4000),G(1,4000)', [mp.nstr(x, 5) for x in gK], 'degR', pari.poldegree(pari.numerator(R), 'k'))
    except Exception as ex:
        print(r['seed'][:12], r['shifts'], r['K'], r['base'], r['k0'], 'ERR', repr(ex)[:80])
